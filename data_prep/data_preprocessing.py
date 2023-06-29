import os
from pathlib import Path
import tarfile
import gzip
from glob import glob
import shutil
import json
import requests
from requests.exceptions import Timeout
from pandas import merge, read_csv, DataFrame


def get_series_probesets(accessions, metadata_dir):
    """
    Retrieves a dictionary mapping GEO series accessions to their respective probesets.

    Args:
        accessions (list): A list of GEO series accessions.
        metadata_dir (Path): The directory containing the metadata.

    Returns:
        dict: A dictionary mapping GEO series accessions to probesets.
    """

    with open(metadata_dir / "series_platforms.json", encoding="utf-8") as series_platforms_json:
        series_platforms = json.load(series_platforms_json)

    with open(metadata_dir / "platform_probesets.json", encoding="utf-8") as platform_probesets_json:
        platform_probesets = json.load(platform_probesets_json)

    supported_probesets = set()
    for probemap_fname in [Path(f) for f in os.listdir(metadata_dir / "probe_maps")]:
        fname_base = probemap_fname.stem
        probeset = fname_base.split("_", 1)[1]
        supported_probesets.add(probeset)

    series_probesets = {}

    for accession in accessions:
        accession = accession.lower()
        if accession in series_platforms.keys():
            platform = series_platforms[accession]
        else:
            continue
        if (platform in platform_probesets.keys() and
                platform_probesets[platform] in supported_probesets):
            series_probesets[accession] = platform_probesets[platform]
        for probeset in platform_probesets.values():
            if platform in probeset and probeset in supported_probesets:
                series_probesets[accession] = probeset

    return series_probesets


def map_probes(counts_path, probe_map_dir, metadata_dir):
    """
    Maps probe IDs in an expression matrix to gene symbols.

    Args:
        counts_path (str): Path to the expression matrix file.
        probe_map_dir (str): Directory containing probe maps.
        metadata_dir (Path): Directory containing metadata files.

    Returns:
        DataFrame: A Pandas DataFrame with the mapped expression matrix.
    """

    counts_fname = Path(counts_path).name
    accession_id = counts_fname.split("_")[0].lower()
    counts_df = read_csv(counts_path, sep="\t")

    counts_df.columns = counts_df.columns.str.replace('.CEL', '', regex=False)

    probeset = get_series_probesets([accession_id], metadata_dir)[accession_id]
    if not probeset:
        return ""
    counts_df.rename(columns={"probe": probeset}, inplace=True)
    probeset_map_path = Path(probe_map_dir) / f"hsapiens_{probeset}.tsv"
    if not os.path.exists(probeset_map_path):
        print(f"No probe map found: {probeset_map_path}")
    else:
        print(f"Probe map found: {probeset_map_path}")
    probe_map = read_csv(probeset_map_path, sep="\t",
                            dtype={"hgnc_symbol": object, probeset: object})

    counts_df[probeset] = counts_df[probeset].apply(lambda x: str(x).rstrip("_at"))
    probe_map[probeset] = probe_map[probeset].apply(lambda x: str(x).rstrip("_at"))

    counts_df = merge(counts_df, probe_map, on=probeset, how="outer")

    cols = list(counts_df.columns)
    cols = [cols[-1]] + cols[1:-1]
    counts_df = counts_df[cols]
    counts_df.dropna(inplace=True)

    return counts_df


def get_unmapped_counts_paths(data_dir):
    """
    Retrieves paths to unmapped expression matrices.

    Args:
        data_dir (str): Path to the directory containing data files.

    Returns:
        list: A list of paths to unmapped expression matrices.
    """

    unmapped_counts_path_pattern = os.path.join(data_dir, "*_counts_unmapped.tsv")

    return glob(unmapped_counts_path_pattern)


def prep_geo_counts(data_dir: Path, probe_maps_dir: Path, metadata_dir: Path):
    """
    Prepares unmapped expression matrices downloaded from GEO.

    Args:
        data_dir (Path): Path to the directory containing data files.
        probe_maps_dir (Path): Path to the directory containing probe maps.
        metadata_dir (Path): Path to the directory containing metadata files.
    """

    for counts_path in get_unmapped_counts_paths(data_dir):
        counts_df = map_probes(counts_path, probe_maps_dir, metadata_dir)

        if not isinstance(counts_df, DataFrame) or len(counts_df) == 0:
            print(f"Could not map probes for {counts_path}")
            continue

        new_counts_path = counts_path.replace("_unmapped.tsv", ".tsv")
        counts_df.to_csv(new_counts_path, sep="\t", index=False)

        # remove old file
        os.remove(counts_path)


def download_file(url, target_path: Path):
    """
    Downloads a file from a given URL.

    Args:
        url (str): The URL of the file to download.
        target_path (Path): The path where the file should be saved.
    """

    try:
        response = requests.get(url, stream=True, timeout=30)
    except Timeout as exc:
        raise Exception(f"Request timed out: {url}") from exc

    with open(target_path, "wb") as handle:
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)


def untar_and_unzip(tar_path: Path, output_dir: Path, delete_archive=False):
    """
    Extracts files from a tar archive and decompresses gz files contained within it.

    Args:
        tar_path (Path): The path to the tar file.
        output_dir (Path): The directory where files should be extracted.
        delete_archive (bool): If True, the original tar file will be deleted after extraction.
    """

    if tar_path.suffix == (".tar"):
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(path=output_dir)

            # Go through the extracted files and extract gz files
            for member in tar.getmembers():
                if str(member.name).endswith(".gz"):
                    gz_path = output_dir / member.name
                    with gzip.open(gz_path, "rb") as f_in:
                        file_name = gz_path.with_suffix("").name
                        # Remove parts after underscore
                        if "_" in file_name:
                            file_name = file_name.split("_")[0] + ".CEL"
                        with open(output_dir / file_name, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(gz_path)

    if delete_archive:
        # delete tar file
        os.remove(tar_path)

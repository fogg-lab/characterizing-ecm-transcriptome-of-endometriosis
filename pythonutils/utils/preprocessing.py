import pandas as pd
import re
from typing import Dict, List


FIGO_MAP_DF = pd.DataFrame({
    "roman_num": ["I", "II", "III", "IV"],
    "figo_chr": ["figo_stage_1", "figo_stage_2", "figo_stage_3", "figo_stage_4"],
    "figo_num": [1, 2, 3, 4]
})


def cols_to_front(df: pd.DataFrame, cols: List) -> List:
    everything_else = df.columns.difference(cols, sort=False).tolist()
    return df.reindex(columns=cols + everything_else)


def load_matrisome_df(matrisome_list_file: str) -> pd.DataFrame:
    matrisome_df = pd.read_csv(matrisome_list_file, sep='\t')
    matrisome_df.columns = map(lambda s: s.lower().replace(" ", "_"), matrisome_df.columns)
    matrisome_df = matrisome_df.query("division != 'Retired'")
    return matrisome_df


def load_survival_df(survival_data_file: str, event_code: Dict) -> pd.DataFrame:
    survival_df = pd.read_csv(survival_data_file, sep = '\t')
    survival_df = (
        survival_df.assign(vital_status_num = survival_df.vital_status.map(lambda x: event_code[x]))
            .drop(["vital_status"], axis=1)
            .rename({"vital_status_num": "vital_status"}, axis = 1)
            .pipe(cols_to_front, ["sample_name", "survival_time", "vital_status"])
    )
    return survival_df


def transpose_df(df: pd.DataFrame, future_colname_col: str, previous_colname_col: str) -> pd.DataFrame:
    df_t = (
        df.set_index(future_colname_col)                      # set as index so will become column names
            .transpose()
            .rename_axis(None, axis=1)                        # column.name will be set to colname_col, we don't want this
            .reset_index()                                    # former index should now be its own column
            .rename({"index": previous_colname_col}, axis=1)
    )
    return df_t


def decode_figo_stage(df, to="num"):
    if to[0] == "n":
        drop_col = "figo_chr"
    elif to[0] == "c":
        drop_col = "figo_num"
    new_df = (
        df.assign(figo_stage_major_rn = lambda x: x.figo_stage.apply(lambda s: re.findall(r"IV|III|II|I", s)[0]))
            .merge(FIGO_MAP_DF, left_on="figo_stage_major_rn", right_on="roman_num", how = "inner")
            .drop(["figo_stage", "figo_stage_major_rn", "roman_num", drop_col], axis=1)
    )
    return new_df

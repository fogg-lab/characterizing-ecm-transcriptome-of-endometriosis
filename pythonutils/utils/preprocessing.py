import pandas as pd
from typing import Dict, List


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

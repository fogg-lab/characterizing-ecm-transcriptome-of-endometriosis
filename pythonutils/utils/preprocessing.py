import pandas as pd
import re
from typing import Dict, List, Tuple, Optional
from numpy.random import RandomState


FIGO_MAP_DF = pd.DataFrame({
    "roman_num": ["I", "II", "III", "IV"],
    "figo_chr": ["figo_stage_1", "figo_stage_2", "figo_stage_3", "figo_stage_4"],
    "figo_num": [1, 2, 3, 4]
})


def cols_to_front(df: pd.DataFrame, cols: List) -> pd.DataFrame:
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
        survival_df.pipe(lambda df: df[df.vital_status.isin(event_code.keys())])    # Vital status may not be available
            # .assign(vital_status_num = survival_df.vital_status.map(lambda x: event_code[x]))
            .pipe(lambda df: df.assign(vital_status_num = df.vital_status.map(lambda x: event_code[x])))
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


def decode_figo_stage(df: pd.DataFrame, to="num") -> pd.DataFrame:
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


def shuffle_data(df: pd.DataFrame, rand: RandomState) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Assumes index 0 is the label index
    shuffled_df = df.sample(frac=1, random_state=rand)
    x_df = shuffled_df.iloc[:, 1:]
    y_df = shuffled_df.iloc[:, [0]]
    return x_df, y_df


def filter_outliers_IQR(df: pd.DataFrame, filter_col: str, coef: float = 1.5):
    q_1 = df[filter_col].quantile(0.25)
    q_3 = df[filter_col].quantile(0.75)
    iqr = q_3 - q_1
    outlier_thresh = (q_1 - coef * iqr, q_3 + coef * iqr)
    filtered_df = df.copy()
    
    def detect_outlier(val, thresh):
        return (val < thresh[0]) | (thresh[1] < val)
    
    filtered_df["outlier_status"] = (
        filtered_df[filter_col].apply(lambda x: detect_outlier(x, outlier_thresh))
    )
    filtered_df = (
        filtered_df.query("outlier_status != True")
            .drop("outlier_status", axis=1)
    )
    return filtered_df


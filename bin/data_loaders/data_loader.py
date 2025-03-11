import sys
import re
import logging
from typing import List, Union, Dict
from pathlib import Path

import pandas as pd

from ..utils import get_tmt_channels
from ..preprocess_tools import merge_ms1_columns, median_centering

logger = logging.getLogger(__name__)


class DataLoader:
    def load_data(self, use_cols: List[str]):
        pass

    def median_centering_within_batch(self, df_raw: pd.DataFrame) -> pd.DataFrame:
        dfs = []
        tmt_channels = get_tmt_channels(df_raw).columns

        correction_factors_all = {}
        for _, df in df_raw.groupby('Batch'):
            batch_name = df['Batch'].iloc[0]
            df.loc[:, tmt_channels], correction_factors = median_centering(df.loc[:, tmt_channels])
            correction_factors_all[batch_name] = correction_factors
            dfs.append(df)

        df_normalized = pd.concat(dfs, axis=0)

        correction_factors_df = pd.DataFrame.from_dict(correction_factors_all)
        correction_factors_df = correction_factors_df.reset_index()
        correction_factors_df = pd.melt(correction_factors_df, id_vars="index", var_name="Variable", value_name="Value")

        return df_normalized, correction_factors_df
    
    def median_centering_ms1(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalizes samples by multiplying each batch MS1 with its own correction factor.
        Only uses peptides detected in >70% of the batches for computing median to prevent
        low-abundant peptides from dragging down the median in runs with deeper coverage
        :param all_batches: list of evidence dataframes
        :return:
        """

        merged_ms1_df = merge_ms1_columns(df)

        # TODO: Find better solution for selecting peptides for median centering. 
        #       Currently, Sarcoma_Batch20-22 have a z-score bias because they have fewer of the "common" peptides than other batches.
        medians = merged_ms1_df[merged_ms1_df.count(axis=1) > 0.7 * len(merged_ms1_df.columns)].median(axis=0).to_dict()
        mean_median = pd.Series(medians.values()).mean()

        dfs, correction_factors = [], {}
        for _, df in df.groupby('Batch'):
            batch_name = df['Batch'].iloc[0]
            correction_factor = (mean_median / medians[batch_name])

            df['MS1'] = df['Intensity'] * correction_factor
            dfs.append(df)
            correction_factors[batch_name] = correction_factor
        
        df_normalized = pd.concat(dfs, axis=0)

        correction_factors_df = pd.DataFrame(correction_factors, index=[0]).T
        return df_normalized, correction_factors_df

    def scale_ms1_to_reference(self, df: pd.DataFrame, ref_channel_df: pd.DataFrame) -> pd.DataFrame:
        return df
    
    def impute_ms1_intensity(self, df: pd.DataFrame, ref_channel_df: pd.DataFrame) -> pd.DataFrame:
        return df
    
    def redistribute_ms1_intensity(self, df: pd.DataFrame) -> pd.DataFrame:
        return df



def extract_cohort_name(evidence_file_path: Union[str, Path]) -> str:
    """Extract batch name including cohort from a file path, e.g. 
    '/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt' => Sarcoma_Batch1
    """
    # match = re.search(r'Searches/([^/]*)/', evidence_file_path)
    # return match.group(1).replace('/', '_')
    match = re.search(r'([^/]*)/Batch', evidence_file_path)
    return match.group(1).replace('/', '_')


# def extract_cohort_name(evidence_file_path: Union[str, Path]) -> str:
#     # match = re.search(r'[^/]*/', evidence_file_path)
#     match = re.search(r'([^/]*)/Batch', evidence_file_path)
#     return match.group(0).split('/')[0]


def extract_batch_name(evidence_file_path: Union[str, Path]) -> str:
    match = re.search(r'[^/]*/Batch[^_]+', evidence_file_path)
    return match.group(0).replace('/', '_')
    # match = re.search(r'(Batch[A-Za-z]*\d*)', evidence_file_path)
    # return match.group(0)


def extract_experiment_name(evidence_file_path: Union[str, Path]) -> str:
    # match = re.search(r'/(Batch\d+_[^/]*)/', evidence_file_path)
    # match = re.search(r'CPTAC_Dou_[\w]*_(Batch[\w\d]*)', evidence_file_path)
    # we want UCEC
    # match = re.search(r'CPTAC_Dou_[\w]*_(Batch[\w\d]*)', evidence_file_path)
    # match1 = re.search(r'/(Batch[\w\d]*_[^/]*)/', evidence_file_path)
    # match2 = re.search(r'/(CPTAC[\w\d]*_[^/]*)/', evidence_file_path
    match = re.search(r'/(Batch[\w\d]*_[^/]*)/', evidence_file_path)
    # if match1:
    #     return match1.group(0).replace('/', '')

    # if match2:
        # return match2.group(0).replace('/', '')

    # this is last part before combined/txt  --> batch?
    # match = re.search(r'/(Batch[\w\d]*_[^/]*)/', evidence_file_path)
    return match.group(1)


def test_batch_names_equals(experiment_to_batch_name_dict: Dict, df: pd.DataFrame):
    # explain
    evidence_files_batch_names = set(experiment_to_batch_name_dict.keys())
    simsi_output_batch_names = set(df['Experiment'].unique().tolist())

    if evidence_files_batch_names != simsi_output_batch_names:
        missing_in_simsi = evidence_files_batch_names - simsi_output_batch_names
        missing_in_evidence_files = simsi_output_batch_names - evidence_files_batch_names
        raise ValueError(f'Batches in list of evidence files and SIMSI evidence do not match:\nMissing in SIMSI: {missing_in_simsi}\nMissing in evidence_files: {missing_in_evidence_files}')

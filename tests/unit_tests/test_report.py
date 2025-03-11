import numpy as np
import pandas as pd
import pytest

from bin.metrics import Metrics
import bin.report_creation as rc

class TestMetrics:

    def test_rank(self, test_df):
        test_df = Metrics.get_rank(test_df)
        expect_result = pd.DataFrame([(1.0, 3.0, 2.0, 4.0, np.nan, 4),
                                      (5.0, 2.0, 1.0, 4.0, 3.0, 5)],
                                     columns=['rank_0', 'rank_1', 'rank_2',
                                              'rank_3', 'rank_4', 'rank_max'])
        pd.testing.assert_frame_equal(test_df, expect_result)

    def test_fold_change(self, test_df):
        test_df = Metrics.get_fold_change(test_df)
        expect_result = pd.DataFrame([(1e99, 0.1, 10.0, 0.01, np.nan),
                                      (0.0, 122.63, 12262.74, 0.01, 0.02)],
                                     columns=['fc_0', 'fc_1', 'fc_2',
                                              'fc_3', 'fc_4'])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)

    def test_zscore(self, test_df):
        test_df = Metrics.get_zscore(test_df)
        expect_result = pd.DataFrame([(104.162, -0.017, 0.017, -0.033, np.nan),
                                      (-0.812, 0.977, 3.533, -0.574, -0.422)],
                                     columns=['zscore_0', 'zscore_1', 'zscore_2',
                                              'zscore_3', 'zscore_4'])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.001)

    def test_pvalue(self):
        test_df = pd.DataFrame([(104.162, -0.017, 0.017, -0.033, np.nan),
                                      (-0.812, 0.977, 3.533, -0.574, -0.422)],
                                     columns=['zscore_0', 'zscore_1', 'zscore_2',
                                              'zscore_3', 'zscore_4'])
        test_df = Metrics.get_pvalues(test_df)
        expect_result = pd.DataFrame([(0.000, 0.986, 0.986, 0.973, np.nan),
                                      (0.416, 0.328, 0.000, 0.565, 0.673)],
                                     columns=['pvalue_zscore_0', 'pvalue_zscore_1', 'pvalue_zscore_2',
                                              'pvalue_zscore_3', 'pvalue_zscore_4'])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.001)


@pytest.fixture
def test_df():
    # log10 intensity values of 2 genes across 5 samples
    return pd.DataFrame([(100, 1, 2, 0.1, np.nan), (0.4, 3, 5, 0.8, 1)])

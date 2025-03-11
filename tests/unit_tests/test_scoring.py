import numpy as np
import pandas as pd

import bin.z_scoring as scoring

def test_calculate_bounded_zscore():
    df = pd.DataFrame({'Z-score': [1.234, -5.234, 6.123, 2.345]})
    df_with_scores = scoring.calculate_bounded_zscore(df)
    assert list(df_with_scores['score']) == [1.234, -4, 4, 2.345]

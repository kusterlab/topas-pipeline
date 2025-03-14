import math

import numpy as np


def calculate_integer_scores(df):
    score = df.copy(deep=False)

    top20 = math.ceil(df.shape[0] * 0.2)
    last_score = df.shape[0] - 4 * top20
    scores = [top20 * [5], top20 * [4], top20 * [3], top20 * [2], last_score * [1]]
    scores = [score for sublist in scores for score in sublist]
    score["rank_score"] = scores

    score["zscore_score"] = 0
    score.loc[score["Z-score"] > 10, "zscore_score"] = 7
    score.loc[(score["Z-score"] > 5) & (score["Z-score"] < 10), "zscore_score"] = 5
    score.loc[(score["Z-score"] > 3) & (score["Z-score"] < 5), "zscore_score"] = 3
    score.loc[(score["Z-score"] > 2) & (score["Z-score"] < 3), "zscore_score"] = 2
    score.loc[(score["Z-score"] > 1) & (score["Z-score"] < 2), "zscore_score"] = 1

    # check if it works
    score["score"] = score["zscore_score"] + score["rank_score"]
    score = score.drop(["zscore_score", "rank_score"], axis=1)
    return score


def calculate_bounded_zscore(df):
    score = df.copy(deep=False)
    score["score"] = np.maximum(-4, np.minimum(4, score["Z-score"]))
    return score


def calculate_bounded_zscore_all(df):
    score = df.copy(deep=False)
    patient_columns = [c for c in score.columns if c.startswith("zscore_")]
    score[patient_columns] = np.maximum(-4, np.minimum(4, score[patient_columns]))
    score = score.rename(columns=lambda x: x.replace("zscore_", "score_"))
    return score

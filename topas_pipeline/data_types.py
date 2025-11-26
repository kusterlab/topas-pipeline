from enum import Enum


class DataType(str, Enum):
    FULL_PROTEOME = "protein"
    PHOSPHO_PROTEOME = "psite"
    PHOSPHO_SCORE = "phospho_score"

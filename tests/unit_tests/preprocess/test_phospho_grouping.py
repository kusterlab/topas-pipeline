import pandas as pd
import numpy as np
import pytest

from topas_pipeline.preprocess import phospho_grouping as pg


class TestPhosphoGrouping:

    @pytest.mark.parametrize(
        "annotations, expected", 
        [
            #singular case
            ({"imputed"}, "imputed;"),
            ({"quan_OOR"}, "quan_OOR;"),
            ({""}, ""),
            
            #combination case
            ({"imputed", "quan_OOR"}, "imputed;quan_OOR;"),
            ({"imputed", ""}, "partially imputed;"),
            ({"quan_OOR", ""}, "quan_OOR;"),
        ],
    )


    def test_summarize_annotations(self, annotations, expected):
        result= pg.summarize_annotations(annotations)
        assert result == expected


    def test_summarize_annotations_order_independance(self):
        a = {"imputed", "quan_OOR"}
        b = {"quan_OOR", "imputed"}
        assert pg.summarize_annotations(a) == pg.summarize_annotations(b)


    def test_summarize_annotations_raises_invalid_combination_error(self):
        with pytest.raises(ValueError):
            pg.summarize_annotations({"unexpected string"})


    def test_summarize_annotations_raises_invalid_imputation_combination_error(self):
        with pytest.raises(ValueError):
            pg.summarize_annotations({"imputed", "partially imputed"})
            
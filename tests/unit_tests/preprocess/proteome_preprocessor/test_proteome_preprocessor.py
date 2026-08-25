import os
from unittest import mock

import pandas as pd
import pytest

from topas_pipeline.config import Preprocessing, Simsi
from topas_pipeline.preprocess.proteome_preprocessor.proteome_preprocessor import (
    ProteomePreprocessor,
)

MODULE = "topas_pipeline.preprocess.proteome_preprocessor.proteome_preprocessor"
# TODO: Verify the tests written

def make_preprocessing_config(**overrides):
    defaults = dict(
        raw_data_location="/raw",
        fasta_file="/fasta.fasta",
        picked_fdr=0.01,
        fdr_num_threads=4,
        imputation=False,
        normalize_to_reference=False,
        debug=False,
        run_lfq=False,
    )
    defaults.update(overrides)
    return Preprocessing(**defaults)


def make_simsi_config(**overrides):
    defaults = dict(simsi_folder="")
    defaults.update(overrides)
    return Simsi(**defaults)


def make_processor(
    data_types=None, quant_file_formats=None, simsi_overrides=None, **config_overrides
):
    if data_types is None:
        data_types = ["fp", "pp"]
    if quant_file_formats is None:
        quant_file_formats = {"fp": "diann", "pp": "ionquant"}
    if simsi_overrides is None:
        simsi_overrides = {}
    return ProteomePreprocessor(
        results_folder="/results",
        sample_annotation_file="/sample_annotation.csv",
        metadata_annotation_file="/metadata_annotation.csv",
        data_types=data_types,
        quant_strategy="LFQ",
        quant_file_formats=quant_file_formats,
        preprocessing_config=make_preprocessing_config(**config_overrides),
        simsi_config=make_simsi_config(**simsi_overrides),
    )


class TestInit:
    def test_sets_attributes(self):
        processor = make_processor()

        assert processor.results_folder == "/results"
        assert processor.sample_annotation_file == "/sample_annotation.csv"
        assert processor.metadata_annotation_file == "/metadata_annotation.csv"
        assert processor.data_types == ["fp", "pp"]
        assert processor.quant_strategy == "LFQ"
        assert processor.quant_file_formats == {"fp": "diann", "pp": "ionquant"}
        assert processor.preprocessing_config["raw_data_location"] == "/raw"
        assert processor.preprocessing_config["fasta_file"] == "/fasta.fasta"
        assert processor.simsi_config["simsi_folder"] == ""

    def test_preprocessing_func_map(self):
        processor = make_processor()

        assert processor.preprocessing_func["fp"] == processor.preprocess_fp
        assert processor.preprocessing_func["pp"] == processor.preprocess_pp


class TestPreprocess:
    @mock.patch(f"{MODULE}.prep.check_annot")
    @mock.patch(f"{MODULE}.sample_annotation.copy_sample_annotation_file")
    @mock.patch(f"{MODULE}.sample_metadata.copy_metadata_file")
    @mock.patch(f"{MODULE}.filter_sample_annotation")
    @mock.patch(f"{MODULE}.load_sample_annotation")
    def test_loads_and_filters_sample_annotation(
        self,
        mock_load,
        mock_filter,
        mock_copy_metadata,
        mock_copy_annot,
        mock_check_annot,
        tmp_path,
    ):
        raw_df = pd.DataFrame({"Sample name": ["a"]})
        filtered_df = pd.DataFrame({"Sample name": ["a"]})
        mock_load.return_value = raw_df
        mock_filter.return_value = filtered_df

        processor = make_processor(data_types=[])
        processor.preprocess(overwrite=True)

        mock_load.assert_called_once_with("/sample_annotation.csv")
        mock_filter.assert_called_once_with(
            raw_df, remove_qc_failed=True, remove_replicates=True
        )
        pd.testing.assert_frame_equal(
            processor.sample_annotation_df, filtered_df.reset_index()
        )

    @mock.patch(f"{MODULE}.prep.check_annot")
    @mock.patch(f"{MODULE}.sample_annotation.copy_sample_annotation_file")
    @mock.patch(f"{MODULE}.sample_metadata.copy_metadata_file")
    @mock.patch(f"{MODULE}.filter_sample_annotation")
    @mock.patch(f"{MODULE}.load_sample_annotation")
    def test_calls_preprocess_proteome_for_each_data_type(
        self,
        mock_load,
        mock_filter,
        mock_copy_metadata,
        mock_copy_annot,
        mock_check_annot,
    ):
        mock_load.return_value = pd.DataFrame({"a": [1]})
        mock_filter.return_value = pd.DataFrame({"a": [1]})

        processor = make_processor(
            data_types=["fp", "pp"],
            quant_file_formats={"fp": "diann", "pp": "ionquant"},
        )
        with mock.patch.object(
            processor, "get_results_file_list", return_value=[]
        ), mock.patch.object(processor, "preprocess_proteome") as mock_preprocess:
            processor.preprocess(overwrite=True)

        assert mock_preprocess.call_count == 2
        called_data_types = [c.args[0] for c in mock_preprocess.call_args_list]
        assert called_data_types == ["fp", "pp"]

    @mock.patch(f"{MODULE}.prep.check_annot")
    @mock.patch(f"{MODULE}.sample_annotation.copy_sample_annotation_file")
    @mock.patch(f"{MODULE}.sample_metadata.copy_metadata_file")
    @mock.patch(f"{MODULE}.filter_sample_annotation")
    @mock.patch(f"{MODULE}.load_sample_annotation")
    def test_skips_when_preprocessed_file_exists_and_not_overwrite(
        self,
        mock_load,
        mock_filter,
        mock_copy_metadata,
        mock_copy_annot,
        mock_check_annot,
        tmp_path,
    ):
        mock_load.return_value = pd.DataFrame({"a": [1]})
        mock_filter.return_value = pd.DataFrame({"a": [1]})

        results_folder = str(tmp_path)
        (tmp_path / "preprocessed_fp.csv").write_text("dummy")

        processor = make_processor(data_types=["fp"])
        processor.results_folder = results_folder

        with mock.patch.object(
            processor, "get_results_file_list", return_value=[]
        ), mock.patch.object(processor, "preprocess_proteome") as mock_preprocess:
            processor.preprocess(overwrite=False)

        mock_preprocess.assert_not_called()

    @mock.patch(f"{MODULE}.prep.check_annot")
    @mock.patch(f"{MODULE}.sample_annotation.copy_sample_annotation_file")
    @mock.patch(f"{MODULE}.sample_metadata.copy_metadata_file")
    @mock.patch(f"{MODULE}.filter_sample_annotation")
    @mock.patch(f"{MODULE}.load_sample_annotation")
    def test_reprocesses_when_overwrite_true_even_if_file_exists(
        self,
        mock_load,
        mock_filter,
        mock_copy_metadata,
        mock_copy_annot,
        mock_check_annot,
        tmp_path,
    ):
        mock_load.return_value = pd.DataFrame({"a": [1]})
        mock_filter.return_value = pd.DataFrame({"a": [1]})

        results_folder = str(tmp_path)
        (tmp_path / "preprocessed_fp.csv").write_text("dummy")

        processor = make_processor(data_types=["fp"])
        processor.results_folder = results_folder

        with mock.patch.object(
            processor, "get_results_file_list", return_value=[]
        ), mock.patch.object(processor, "preprocess_proteome") as mock_preprocess:
            processor.preprocess(overwrite=True)

        mock_preprocess.assert_called_once()


class TestPreprocessProteome:
    @mock.patch(f"{MODULE}.QuantDataLoaderFactory")
    def test_reuses_existing_preprocessed2_file(
        self, mock_factory, tmp_path
    ):
        results_folder = str(tmp_path)
        preprocessed2_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        preprocessed2_df.to_csv(
            os.path.join(results_folder, "preprocessed_fp2.csv"), index=False
        )

        processor = make_processor(data_types=["fp"])
        processor.results_folder = results_folder

        with mock.patch.object(
            processor, "get_results_file_list"
        ) as mock_get_files, mock.patch.object(
            processor, "preprocess_fp"
        ) as mock_preprocess_fp, mock.patch(
            f"{MODULE}.get_channel_to_sample_id_dict"
        ) as mock_channel_dict, mock.patch(
            f"{MODULE}.sample_mapping.rename_columns_with_sample_ids"
        ) as mock_rename:
            mock_channel_dict.return_value = {}
            mock_rename.return_value = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})

            processor.preprocess_proteome(
                "fp", "diann", pd.DataFrame({"Batch": ["b1"], "Batch FP": ["b1"]})
            )

        mock_get_files.assert_called_once_with("fp", "diann")
        mock_factory.get_loader.assert_not_called()
        mock_preprocess_fp.assert_not_called()
        mock_rename.assert_called_once()

    @mock.patch(f"{MODULE}.get_channel_to_sample_id_dict")
    @mock.patch(f"{MODULE}.sample_mapping.rename_columns_with_sample_ids")
    @mock.patch(f"{MODULE}.QuantDataLoaderFactory")
    def test_loads_and_preprocesses_fp_when_no_preprocessed2(
        self, mock_factory, mock_rename, mock_channel_dict, tmp_path
    ):
        results_folder = str(tmp_path)
        processor = make_processor(data_types=["fp"])
        processor.results_folder = results_folder

        loaded_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        preprocessed_df = pd.DataFrame({"Gene names": ["G1"], "value": [2.0]})
        renamed_df = pd.DataFrame({"Gene names": ["G1"], "value": [2.0]})

        mock_loader_instance = mock.Mock()
        mock_loader_instance.load_and_normalize.return_value = loaded_df
        mock_loader_class = mock.Mock(return_value=mock_loader_instance)
        mock_factory.get_loader.return_value = mock_loader_class
        mock_channel_dict.return_value = {"channel1": "sample1"}
        mock_rename.return_value = renamed_df

        sample_annotation_df = pd.DataFrame(
            {"Batch": ["b1"], "Batch FP": ["b1"]}
        )

        mock_preprocess_fp = mock.Mock(return_value=preprocessed_df)
        processor.preprocessing_func["fp"] = mock_preprocess_fp

        with mock.patch.object(
            processor, "get_results_file_list", return_value=["file1"]
        ) as mock_get_files:
            processor.preprocess_proteome("fp", "diann", sample_annotation_df)

        mock_get_files.assert_called_once_with("fp", "diann")
        mock_factory.get_loader.assert_called_once_with("LFQ")
        mock_loader_class.assert_called_once_with(
            ["file1"],
            "diann",
            sample_annotation_df,
            False,
            False,
            {
                "results_folder": results_folder,
                "data_type": "fp",
                "quant_strategy": "LFQ",
                "simsi_folder": "",
            },
        )
        mock_preprocess_fp.assert_called_once_with(loaded_df)
        mock_channel_dict.assert_called_once()
        mock_rename.assert_called_once_with(
            preprocessed_df, {"channel1": "sample1"}, index_cols=["Gene names"]
        )

        preprocessed2_path = os.path.join(results_folder, "preprocessed_fp2.csv")
        preprocessed_path = os.path.join(results_folder, "preprocessed_fp.csv")
        assert os.path.exists(preprocessed2_path)
        assert os.path.exists(preprocessed_path)

    @mock.patch(f"{MODULE}.get_channel_to_sample_id_dict")
    def test_adds_batch_column_from_batch_fp_when_missing(
        self, mock_channel_dict, tmp_path
    ):
        results_folder = str(tmp_path)
        processor = make_processor(data_types=["fp"])
        processor.results_folder = results_folder

        preprocessed2_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        preprocessed2_df.to_csv(
            os.path.join(results_folder, "preprocessed_fp2.csv"), index=False
        )

        mock_channel_dict.return_value = {}
        sample_annotation_df = pd.DataFrame({"Batch FP": ["b1"]})

        with mock.patch.object(
            processor, "get_results_file_list", return_value=["file1"]
        ), mock.patch(
            f"{MODULE}.sample_mapping.rename_columns_with_sample_ids",
            return_value=pd.DataFrame({"Gene names": ["G1"], "value": [1.0]}),
        ):
            processor.preprocess_proteome("fp", "diann", sample_annotation_df)

        assert "Batch" in sample_annotation_df.columns
        assert sample_annotation_df["Batch"].tolist() == ["b1"]

    @mock.patch(f"{MODULE}.phosphopeptides.group_phosphopeptides_and_normalize")
    @mock.patch(f"{MODULE}.QuantDataLoaderFactory")
    def test_pp_calls_group_phosphopeptides_and_normalize(
        self, mock_factory, mock_group_pp, tmp_path
    ):
        results_folder = str(tmp_path)
        processor = make_processor(
            data_types=["pp"], quant_file_formats={"pp": "ionquant"}
        )
        processor.results_folder = results_folder

        loaded_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        preprocessed_df = pd.DataFrame({"Gene names": ["G1"], "value": [2.0]})
        grouped_df = pd.DataFrame({"Gene names": ["G1"], "value": [3.0]})

        mock_loader_instance = mock.Mock()
        mock_loader_instance.load_and_normalize.return_value = loaded_df
        mock_loader_class = mock.Mock(return_value=mock_loader_instance)
        mock_factory.get_loader.return_value = mock_loader_class
        mock_group_pp.return_value = grouped_df
        processor.preprocessing_func["pp"] = mock.Mock(return_value=preprocessed_df)

        with mock.patch.object(
            processor, "get_results_file_list", return_value=["file1"]
        ):
            processor.preprocess_proteome(
                "pp", "ionquant", pd.DataFrame({"Batch": ["b1"]})
            )

        mock_group_pp.assert_called_once_with(
            results_folder=results_folder,
            sample_annotation_file="/sample_annotation.csv",
            run_lfq=False,
        )

        preprocessed_path = os.path.join(results_folder, "preprocessed_pp.csv")
        assert os.path.exists(preprocessed_path)
        result_df = pd.read_csv(preprocessed_path)
        assert result_df["value"].tolist() == [3.0]


class TestPreprocessFp:
    @mock.patch(f"{MODULE}.id_meta")
    @mock.patch(f"{MODULE}.prep")
    @mock.patch(f"{MODULE}.picked_group")
    def test_calls_pipeline_steps_in_order_and_returns_df(
        self, mock_picked_group, mock_prep, mock_id_meta
    ):
        processor = make_processor()
        input_df = pd.DataFrame({"Gene names": ["G1"]})

        step_dfs = [
            pd.DataFrame({"step": [i]}) for i in range(7)
        ]
        mock_picked_group.picked_protein_grouping.return_value = step_dfs[0]
        mock_prep.filter_data.return_value = step_dfs[1]
        mock_id_meta.create_metadata_columns.return_value = step_dfs[2]
        mock_id_meta.mark_num_peptides.return_value = step_dfs[3]
        mock_id_meta.mark_peptide_id_out_of_range.return_value = step_dfs[4]
        mock_id_meta.mark_detected_in_batch.return_value = step_dfs[5]
        mock_prep.log_transform_intensities.return_value = step_dfs[6]
        final_df = pd.DataFrame({"final": [1]})
        mock_id_meta.mark_quant_out_of_range.return_value = final_df

        result = processor.preprocess_fp(input_df)

        mock_picked_group.picked_protein_grouping.assert_called_once_with(
            input_df, "/results", 0.01, "/fasta.fasta", 4
        )
        mock_prep.filter_data.assert_called_once_with(step_dfs[0], data_type="fp")
        mock_id_meta.create_metadata_columns.assert_called_once_with(step_dfs[1])
        mock_id_meta.mark_num_peptides.assert_called_once_with(step_dfs[2])
        mock_id_meta.mark_peptide_id_out_of_range.assert_called_once_with(step_dfs[3])
        mock_id_meta.mark_detected_in_batch.assert_called_once_with(step_dfs[4])
        mock_prep.log_transform_intensities.assert_called_once_with(step_dfs[5])
        mock_id_meta.mark_quant_out_of_range.assert_called_once_with(step_dfs[6])

        pd.testing.assert_frame_equal(result, final_df)


class TestPreprocessPp:
    def _patch_common(self, mock_prep, mock_id_meta, mock_picked_group):
        save_debug_df = mock.Mock()
        mock_prep.get_save_debug_df_function.return_value = save_debug_df
        mock_picked_group.remap_gene_names.side_effect = lambda df, fasta: df
        mock_id_meta.create_metadata_columns.side_effect = lambda df: df
        mock_prep.filter_data.side_effect = lambda df, data_type: df
        mock_prep.sum_peptide_intensities.side_effect = lambda df, debug, run_lfq: df
        mock_prep.log_transform_intensities.side_effect = lambda df: df
        mock_id_meta.mark_quant_out_of_range.side_effect = lambda df: df
        return save_debug_df

    @mock.patch(f"{MODULE}.id_meta")
    @mock.patch(f"{MODULE}.prep")
    @mock.patch(f"{MODULE}.picked_group")
    def test_skips_imputation_when_disabled(
        self, mock_picked_group, mock_prep, mock_id_meta
    ):
        self._patch_common(mock_prep, mock_id_meta, mock_picked_group)
        processor = make_processor(imputation=False)

        input_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        wide_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        mock_prep.convert_long_to_wide_format.return_value = wide_df

        result = processor.preprocess_pp(input_df)

        mock_prep.impute_data.assert_not_called()
        pd.testing.assert_frame_equal(result, wide_df)

    @mock.patch(f"{MODULE}.id_meta")
    @mock.patch(f"{MODULE}.prep")
    @mock.patch(f"{MODULE}.picked_group")
    def test_runs_imputation_when_enabled(
        self, mock_picked_group, mock_prep, mock_id_meta, tmp_path
    ):
        save_debug_df = self._patch_common(mock_prep, mock_id_meta, mock_picked_group)
        processor = make_processor(imputation=True)
        processor.results_folder = str(tmp_path)

        input_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        imputed_df = pd.DataFrame({"Gene names": ["G1"], "value": [2.0]})
        mock_prep.impute_data.return_value = imputed_df
        wide_df = pd.DataFrame({"Gene names": ["G1"], "value": [2.0]})
        mock_prep.convert_long_to_wide_format.return_value = wide_df

        result = processor.preprocess_pp(input_df)

        mock_prep.impute_data.assert_called_once_with(input_df)
        assert os.path.exists(
            os.path.join(str(tmp_path), "preprocessed_pp_before_imputation.csv")
        )
        save_debug_df.assert_any_call(imputed_df, "_after_imputation")
        pd.testing.assert_frame_equal(result, wide_df)

    @mock.patch(f"{MODULE}.id_meta")
    @mock.patch(f"{MODULE}.prep")
    @mock.patch(f"{MODULE}.picked_group")
    def test_extracts_transferred_spectra_count_columns(
        self, mock_picked_group, mock_prep, mock_id_meta, tmp_path
    ):
        self._patch_common(mock_prep, mock_id_meta, mock_picked_group)
        processor = make_processor(imputation=False)
        processor.results_folder = str(tmp_path)

        input_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        wide_df = pd.DataFrame(
            {
                "value": [1.0],
                "Transferred spectra count sample1": [5],
            }
        )
        mock_prep.convert_long_to_wide_format.return_value = wide_df

        result = processor.preprocess_pp(input_df)

        transfer_path = os.path.join(str(tmp_path), "Transfer_metadata.csv")
        assert os.path.exists(transfer_path)
        assert "Transferred spectra count sample1" not in result.columns
        assert "value" in result.columns

    @mock.patch(f"{MODULE}.id_meta")
    @mock.patch(f"{MODULE}.prep")
    @mock.patch(f"{MODULE}.picked_group")
    def test_no_transfer_metadata_written_when_column_absent(
        self, mock_picked_group, mock_prep, mock_id_meta, tmp_path
    ):
        self._patch_common(mock_prep, mock_id_meta, mock_picked_group)
        processor = make_processor(imputation=False)
        processor.results_folder = str(tmp_path)

        input_df = pd.DataFrame({"Gene names": ["G1"], "value": [1.0]})
        wide_df = pd.DataFrame({"value": [1.0]})
        mock_prep.convert_long_to_wide_format.return_value = wide_df

        result = processor.preprocess_pp(input_df)

        transfer_path = os.path.join(str(tmp_path), "Transfer_metadata.csv")
        assert not os.path.exists(transfer_path)
        pd.testing.assert_frame_equal(result, wide_df)

    @mock.patch(f"{MODULE}.id_meta")
    @mock.patch(f"{MODULE}.prep")
    @mock.patch(f"{MODULE}.picked_group")
    def test_empty_df_skips_transfer_metadata_extraction(
        self, mock_picked_group, mock_prep, mock_id_meta, tmp_path
    ):
        self._patch_common(mock_prep, mock_id_meta, mock_picked_group)
        processor = make_processor(imputation=False)
        processor.results_folder = str(tmp_path)

        input_df = pd.DataFrame({"Gene names": [], "value": []})
        wide_df = pd.DataFrame({"Transferred spectra count sample1": []})
        mock_prep.convert_long_to_wide_format.return_value = wide_df

        result = processor.preprocess_pp(input_df)

        transfer_path = os.path.join(str(tmp_path), "Transfer_metadata.csv")
        assert not os.path.exists(transfer_path)
        pd.testing.assert_frame_equal(result, wide_df)


class TestGetResultsFileList:
    def test_dispatches_to_evidence_loader(self):
        mock_get_files = mock.Mock(return_value=["f1", "f2", "f3"])
        processor = make_processor()
        processor.sample_annotation_df = pd.DataFrame({"a": [1]})

        with mock.patch.dict(
            f"{MODULE}.RESULT_FILE_GETTER_REGISTRY", {"evidence": mock_get_files}
        ):
            result = processor.get_results_file_list("fp", "evidence")

        mock_get_files.assert_called_once_with(
            "fp", "/raw", processor.sample_annotation_df
        )
        assert result == ["f1", "f2", "f3"]

    def test_dispatches_to_diann_loader(self):
        mock_get_files = mock.Mock(return_value=["f1", "f2", "f3"])
        processor = make_processor()
        processor.sample_annotation_df = pd.DataFrame({"a": [1]})

        with mock.patch.dict(
            f"{MODULE}.RESULT_FILE_GETTER_REGISTRY", {"diann": mock_get_files}
        ):
            result = processor.get_results_file_list("fp", "diann")

        mock_get_files.assert_called_once_with(
            "fp", "/raw", processor.sample_annotation_df
        )
        assert result == ["f1", "f2", "f3"]

    def test_dispatches_to_ionquant_loader(self):
        mock_get_files = mock.Mock(return_value=["f1", "f2", "f3"])
        processor = make_processor()
        processor.sample_annotation_df = pd.DataFrame({"a": [1]})

        with mock.patch.dict(
            f"{MODULE}.RESULT_FILE_GETTER_REGISTRY", {"ionquant": mock_get_files}
        ):
            result = processor.get_results_file_list("pp", "ionquant")

        mock_get_files.assert_called_once_with(
            "pp", "/raw", processor.sample_annotation_df
        )
        assert result == ["f1", "f2", "f3"]

    def test_unknown_format_raises_value_error(self):
        processor = make_processor()
        processor.sample_annotation_df = pd.DataFrame({"a": [1]})

        with pytest.raises(ValueError):
            processor.get_results_file_list("fp", "unknown_format")

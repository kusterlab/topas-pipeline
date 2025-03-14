import os
from pathlib import Path
import argparse
import time
import subprocess

import pandas as pd
import pytest

import topas_pipeline.simsi as simsi
import topas_pipeline.config as config
import topas_pipeline.utils

from job_pool import JobPool


class TestMain:
    # Parses command line arguments correctly when provided
    def test_parses_command_line_arguments_correctly(self, mocker):
        mocker.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=argparse.Namespace(config="path/to/config"),
        )
        mocker.patch(
            "topas_pipeline.config.load",
            return_value={
                "preprocessing": {
                    "run_simsi": True,
                    "raw_data_location": "path/to/data",
                },
                "results_folder": "results",
                "sample_annotation": "sample",
                "simsi": {},
                "data_types": [],
                "slack_webhook_url": "dummy_url",
            },
        )
        mocker.patch("os.makedirs")
        mocker.patch("json.dumps")
        mocker.patch("builtins.open", mocker.mock_open())
        mocker.patch(
            "topas_pipeline.simsi.load_sample_annotation",
            return_value=mocker.Mock(to_csv=mocker.Mock()),
        )
        mocker.patch("topas_pipeline.simsi.run_simsi")

        simsi.main(["-c", "path/to/config"])

        config.load.assert_called_once_with("path/to/config")
        simsi.run_simsi.assert_called_once()

    def test_can_skip_simsi(self, mocker):
        mocker.patch(
            "argparse.ArgumentParser.parse_args",
            return_value=argparse.Namespace(config="path/to/config"),
        )
        mocker.patch(
            "topas_pipeline.config.load",
            return_value={
                "preprocessing": {
                    "run_simsi": False,
                    "raw_data_location": "path/to/data",
                },
                "results_folder": "results",
                "sample_annotation": "sample",
                "simsi": {},
                "data_types": [],
            },
        )
        mocker.patch("topas_pipeline.simsi.run_simsi")

        simsi.main(["-c", "path/to/config"])
        simsi.run_simsi.assert_not_called()


class TestRunSimsi:
    # Each data_type in data_types is processed
    def test_each_data_type_processed(self, mocker):
        data_types = ["type1", "type2"]
        kwargs = {"data_types": data_types}

        mocker.patch("topas_pipeline.simsi.JobPool")
        mocker.patch("topas_pipeline.simsi.JobPool.applyAsync")
        mocker.patch("topas_pipeline.simsi.JobPool.checkPool")

        simsi.run_simsi(**kwargs)

        # Assert that applyAsync is called for each data_type
        assert simsi.JobPool.return_value.applyAsync.call_count == len(data_types)
        # Assert that pool is closed at the end
        simsi.JobPool.return_value.checkPool.assert_called_once()


class TestRunSimsiDataType:
    # Function executes without exceptions and sends a success message to Slack
    def test_run_simsi_data_type_success(self, mocker):
        mocker.patch("time.time", side_effect=[0, 10, 20, 30])
        mocker.patch("topas_pipeline.simsi.init_file_logger")
        mocker.patch("topas_pipeline.simsi.send_slack_message")
        mocker.patch(
            "topas_pipeline.preprocess_tools.get_summary_files", return_value=["summary_file"]
        )
        mocker.patch("topas_pipeline.simsi.copy_raw_files")
        mocker.patch(
            "topas_pipeline.meta_input_file.get_meta_input_file_path",
            return_value="meta_input_file",
        )
        mocker.patch("topas_pipeline.simsi.create_meta_input_file")
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"mq_txt_folder": ["folder_Batch1"]}),
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_output_folder", return_value="simsi_output_folder"
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_cache_folder", return_value="simsi_cache_folder"
        )
        mocker.patch("topas_pipeline.simsi.find_matching_summaries_folder", return_value=None)
        mocker.patch("topas_pipeline.simsi.run_simsi_single")
        mocker.patch("topas_pipeline.simsi.store_results_for_reuse")

        simsi.run_simsi_data_type(
            results_folder="results",
            search_result_folder="search_results",
            sample_annotation_file="sample_annotation",
            simsi_folder="simsi",
            tmt_ms_level="2",
            stringencies=1,
            maximum_pep=100,
            num_threads=4,
            slack_webhook_url="dummy_url",
            data_type="typeA",
        )

        simsi.send_slack_message.assert_called_once_with(
            "SIMSI typeA finished", "results", "dummy_url"
        )


class TestGetSimsiFolders:
    def test_returns_correct_path_with_valid_inputs(self):
        simsi_folder = "test_folder"
        data_type = "typeA"
        result_folder_name = "results"
        expected_path = Path("test_folder/simsi_output/TYPEA/results")
        assert (
            simsi.get_simsi_output_folder(simsi_folder, data_type, result_folder_name)
            == expected_path
        )

    def test_valid_simsi_folder_and_data_type(self):
        simsi_folder = "/path/to/simsi"
        data_type = "mzML"
        expected_path = Path("/path/to/simsi/simsi_output/MZML")
        assert simsi.get_simsi_cache_folder(simsi_folder, data_type) == expected_path

    def test_valid_simsi_raw_file_folder(self):
        simsi_folder = "/path/to/simsi"
        data_type = "fp"
        experiment = "Batch1"
        expected_path = "/path/to/simsi/raw_files/FP/Batch1"
        assert (
            simsi.get_simsi_raw_file_folder(simsi_folder, data_type, experiment)
            == expected_path
        )


class TestFindMatchingSummariesFolder:
    def test_finds_correct_summary_folder(self, mocker):
        simsi_cache_folder = mocker.Mock(spec=Path)
        meta_input_file = mocker.Mock(spec=Path)
        meta_input_file.name = "meta_input.json"

        summary_folder = Path("meta_input.json")

        simsi_cache_folder.glob.return_value = [summary_folder]
        mocker.patch("pathlib.Path.is_dir", return_value=True)
        mocker.patch("pathlib.Path.is_file", return_value=True)
        mocker.patch("os.path.getmtime", return_value=1)

        mocker.patch("topas_pipeline.meta_input_file.meta_input_files_equal", return_value=True)

        assert (
            simsi.find_matching_summaries_folder(simsi_cache_folder, meta_input_file)
            == summary_folder
        )


class TestCopyWithSubprocess:
    def test_copy_file_when_source_exists_and_dest_does_not(self, mocker):
        source = "source.txt"
        dest = "dest.txt"

        mocker.patch("os.path.exists", side_effect=lambda x: x == source)
        mocker.patch("os.stat", side_effect=lambda x: mocker.Mock(st_size=100))
        mocker.patch("subprocess.Popen", return_value=mocker.Mock(wait=lambda: 0))

        result = simsi.copy_with_subprocess(source, dest)

        subprocess.Popen.assert_called_once()

        assert result is True

    def test_return_true_when_dest_file_exists_with_same_size(self, mocker):
        source = "source.txt"
        dest = "dest.txt"

        mocker.patch("os.path.exists", return_value=True)
        mocker.patch("os.stat", side_effect=lambda x: mocker.Mock(st_size=100))
        mocker.patch("subprocess.Popen", return_value=mocker.Mock(wait=lambda: 0))

        result = simsi.copy_with_subprocess(source, dest)

        assert not subprocess.Popen.called

        assert result is True

    def test_copy_file_when_source_exists_and_dest_sizes_differ(self, mocker):
        source = "source.txt"
        dest = "dest.txt"

        mocker.patch("os.path.exists", side_effect=lambda x: x in [source, dest])
        mocker.patch(
            "os.stat", side_effect=[mocker.Mock(st_size=100), mocker.Mock(st_size=200)]
        )
        mocker.patch("subprocess.Popen", return_value=mocker.Mock(wait=lambda: 0))

        result = simsi.copy_with_subprocess(source, dest)

        subprocess.Popen.assert_called_once()

        assert result is True

    def test_source_file_does_not_exist(self, mocker):
        source = "non_existent_source.txt"
        dest = "dest.txt"

        mocker.patch("os.path.exists", return_value=False)
        mocker.patch("topas_pipeline.simsi.logger.info")
        mocker.patch("subprocess.Popen", return_value=mocker.Mock(wait=lambda: 0))

        result = simsi.copy_with_subprocess(source, dest)

        assert not subprocess.Popen.called

        assert result is False
        simsi.logger.info.assert_called_once_with(
            f"Could not find source file {source}"
        )


class TestGetCorrectionFactorFiles:
    def test_get_correction_factor_files(
        self, correction_file_mapping_file, queue_file
    ):
        assert simsi.get_correction_factor_files(
            ["Batch2_PP_MASTER"], correction_file_mapping_file, queue_file, Path(".")
        ) == ["UF291262_UE277617.txt"]

    def test_get_correction_factor_files_from_queue_file(
        self, correction_file_mapping_file, queue_file
    ):
        assert simsi.get_correction_factor_files(
            ["Batch3_PP_MASTER"], correction_file_mapping_file, queue_file, Path(".")
        ) == ["UF291262_UE277617.txt"]

    def test_get_correction_factor_files_missing(
        self, correction_file_mapping_file, queue_file
    ):
        with pytest.raises(ValueError, match=r".*Batch34_PP_MASTER.*"):
            simsi.get_correction_factor_files(
                ["Batch34_PP_MASTER"],
                correction_file_mapping_file,
                queue_file,
                Path("."),
            )


class TestFindSimsiEvidenceFile:
    # Correctly identifies and returns the SIMSI evidence file path when all folders and files exist
    def test_correct_evidence_file_path(self, mocker):
        results_folder = "results"
        simsi_folder = "simsi"
        data_type = "type"
        meta_input_file = "meta_input_file_path"
        simsi_cache_folder = "simsi_cache_folder"
        simsi_summaries_folder = "simsi_summaries_folder"
        simsi_evidence_file = f"{simsi_summaries_folder}/p10/p10_evidence.txt"

        mocker.patch(
            "topas_pipeline.meta_input_file.get_meta_input_file_path", return_value=meta_input_file
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_cache_folder", return_value=simsi_cache_folder
        )
        mocker.patch(
            "topas_pipeline.simsi.find_matching_summaries_folder",
            return_value=simsi_summaries_folder,
        )
        mocker.patch("os.path.isfile", return_value=True)

        # Assert
        assert (
            simsi.find_simsi_evidence_file(results_folder, simsi_folder, data_type)
            == simsi_evidence_file
        )


class TestExtractSimsiEvidenceFileFromArchive:
    # Extracts SIMSI evidence file successfully from a valid archive
    def test_extracts_simsi_evidence_file_successfully(self, mocker):
        import io
        import zipfile

        mocker.patch("pathlib.Path.is_file", return_value=True)

        # Mocking the zipfile and its contents
        mock_archive = mocker.Mock()
        mock_archive.namelist.return_value = ["valid_archive/p10/p10_evidence.txt"]
        mock_archive.read.return_value = b"evidence data"

        mock_zipfile = mocker.patch("zipfile.ZipFile")
        mock_zipfile.return_value = mock_archive

        result = simsi.extract_simsi_evidence_file_from_archive("valid_archive/xxx/xxx")

        assert isinstance(result, io.BytesIO)
        assert result.getvalue() == b"evidence data"


class TestStoreResultsForReuse:
    # Successfully renames the summaries directory
    def test_rename_summaries_directory(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tempdir:
            simsi_output_folder = Path(tempdir) / "output"
            simsi_cache_folder = Path(tempdir) / "cache"
            meta_input_file = Path(tempdir) / "meta.txt"
            results_folder_name = "test_results"

            # Create necessary directories and files
            simsi_output_folder.mkdir()
            (simsi_output_folder / "summaries").mkdir()
            (simsi_output_folder / "maracluster_output").mkdir()
            (simsi_output_folder / "SIMSI.log").touch()
            meta_input_file.write_text("meta content")

            # Create cache folder
            simsi_cache_folder.mkdir()

            # Call the function
            simsi.store_results_for_reuse(
                simsi_output_folder,
                simsi_cache_folder,
                meta_input_file,
                results_folder_name,
            )

            # Check if the summaries directory was renamed
            assert (simsi_cache_folder / f"summaries_{results_folder_name}").exists()


class TestRunSimsiSingle:
    # Executes simsi.main with correct configuration parameters
    def test_executes_simsi_main_with_correct_params(self, mocker):
        meta_input_file = Path("/path/to/meta_input_file")
        simsi_output_folder = Path("/path/to/output_folder")
        simsi_cache_folder = Path("/path/to/cache_folder")
        tmt_ms_level = "2"
        stringencies = 3
        tmt_requantify = True
        maximum_pep = 1000
        num_threads = 4

        mock_main = mocker.patch("topas_pipeline.simsi.simsi.main")

        simsi.run_simsi_single(
            meta_input_file,
            simsi_output_folder,
            simsi_cache_folder,
            tmt_ms_level,
            stringencies,
            tmt_requantify,
            maximum_pep,
            num_threads,
        )

        expected_configs = [
            "--meta_input_file",
            str(meta_input_file),
            "--output_folder",
            str(simsi_output_folder),
            "--cache_folder",
            str(simsi_cache_folder),
            "--tmt_ms_level",
            str(tmt_ms_level),
            "--skip_annotated_clusters",
            "--skip_msmsscans",
            "--skip_msms",
            "--stringencies",
            str(stringencies),
            "--maximum_pep",
            str(maximum_pep),
            "--num_threads",
            str(num_threads),
            "--tmt_requantify",
        ]

        mock_main.assert_called_once_with(expected_configs)


class TestCreateMetaInputFile:
    # Correctly reads summary files and extracts experiment names
    def test_reads_summary_files_and_extracts_experiment_names(self, mocker):
        summary_files = ["summary1.txt", "summary2.txt"]
        data_type = "typeA"
        simsi_folder = "simsi_folder"
        meta_input_file = Path("meta_input_file.txt")

        mocker.patch(
            "pandas.read_csv",
            side_effect=[
                pd.DataFrame({"Experiment": ["exp1"]}),
                pd.DataFrame({"Experiment": ["exp2"]}),
            ],
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_raw_file_folder", side_effect=["folder1", "folder2"]
        )
        mocker.patch(
            "topas_pipeline.simsi.get_correction_factor_files", return_value=["file1", "file2"]
        )
        mocker.patch("topas_pipeline.simsi.mi.write_meta_input_file")

        simsi.create_meta_input_file(
            summary_files, data_type, simsi_folder, meta_input_file
        )

        pd.read_csv.assert_any_call("summary1.txt", sep="\t")
        pd.read_csv.assert_any_call("summary2.txt", sep="\t")
        simsi.mi.write_meta_input_file.assert_called_once_with(
            meta_input_file,
            [Path("summary1.txt").parent, Path("summary2.txt").parent],
            ["folder1", "folder2"],
            ["file1", "file2"],
        )


class TestCopyRawFiles:
    # Copies raw files successfully when all files are available
    def test_copies_raw_files_successfully(self, mocker):
        mocker.patch("os.path.exists", return_value=False)
        mocker.patch("os.makedirs")
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame(
                {
                    "Raw file": ["file1", "file2", "Total"],
                    "Experiment": ["exp1", "exp1", "exp1"],
                }
            ),
        )
        mocker.patch("topas_pipeline.simsi.copy_with_subprocess", return_value=True)
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_cache_folder", return_value="/mock/cache/folder"
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_raw_file_folder", return_value="/mock/simsi/folder"
        )

        raw_file_folders = ["/mock/raw/folder1", "/mock/raw/folder2"]
        summary_files = ["/mock/summary/file1.txt"]
        data_type = "mock_data_type"
        simsi_folder = "/mock/simsi/folder"

        simsi.copy_raw_files(raw_file_folders, summary_files, data_type, simsi_folder)

        os.makedirs.assert_called_once_with("/mock/simsi/folder")
        simsi.copy_with_subprocess.assert_any_call(
            "/mock/raw/folder1/file1.raw", "/mock/simsi/folder/file1.raw"
        )
        simsi.copy_with_subprocess.assert_any_call(
            "/mock/raw/folder1/file2.raw", "/mock/simsi/folder/file2.raw"
        )

    def test_skips_copying_existing_mzml_files(self, mocker):
        # Mocking dependencies
        mocker.patch("os.path.exists", return_value=True)
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame(
                {"Raw file": ["file1", "file2"], "Experiment": ["exp1", "exp2"]}
            ),
        )
        mocker.patch("topas_pipeline.simsi.copy_with_subprocess", return_value=True)
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_cache_folder", return_value="/mock/cache/folder"
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_raw_file_folder", return_value="/mock/simsi/folder"
        )

        raw_file_folders = ["/mock/raw/folder1", "/mock/raw/folder2"]
        summary_files = ["/mock/summary/file1.txt"]
        data_type = "mock_data_type"
        simsi_folder = "/mock/simsi/folder"

        # Call the function under test
        simsi.copy_raw_files(raw_file_folders, summary_files, data_type, simsi_folder)

        # Assertions
        simsi.copy_with_subprocess.assert_not_called()

    def test_handles_missing_raw_files(self, mocker):
        # Mocking dependencies
        mocker.patch("os.path.exists", return_value=False)
        mocker.patch("os.makedirs")
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame(
                {
                    "Raw file": ["file1", "file2", "Total"],
                    "Experiment": ["exp1", "exp1", "exp1"],
                }
            ),
        )
        mocker.patch("topas_pipeline.simsi.copy_with_subprocess", return_value=False)
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_cache_folder", return_value="/mock/cache/folder"
        )
        mocker.patch(
            "topas_pipeline.simsi.get_simsi_raw_file_folder", return_value="/mock/simsi/folder"
        )

        raw_file_folders = ["/mock/raw/folder1", "/mock/raw/folder2"]
        summary_files = ["/mock/summary/file1.txt"]
        data_type = "mock_data_type"
        simsi_folder = "/mock/simsi/folder"

        # Call the function under test and assert it raises RuntimeError
        with pytest.raises(
            RuntimeError, match="Failed to find or copy the raw file file1."
        ):
            simsi.copy_raw_files(
                raw_file_folders, summary_files, data_type, simsi_folder
            )


@pytest.fixture
def queue_file(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    queue_file = tmp_path / "queue.csv"
    queue_file_content = """drug, MQ version, pre payload, post payload, threads, experiment, fasta file, raw folder, phospho, protease, mqpar, tmt corr factors, peptide fdr, protein fdr
Batch23_PP_INFORM,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,VA296455_UJ279751.txt,0.01,1
Batch22_PP_INFORM,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,VA296455_UJ279751.txt,0.01,1
Batch1_PP_INFORM_MASTER,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,UF291262_UE277617.txt,0.01,1
Batch2_PP_MASTER,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,UF291262_UE277617.txt,0.01,1
Batch3_PP_MASTER,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,UF291262_UE277617.txt,0.01,1"""

    with open(queue_file, "w") as f:
        f.write(queue_file_content)
    return queue_file


@pytest.fixture
def correction_file_mapping_file(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    correction_file_mapping_file = tmp_path / "correction_factor_files.csv"
    correction_file_mapping_file_content = """experiment	correction_factor_file
Batch23_PP_INFORM	VA296455_UJ279751.txt
Batch22_PP_INFORM	VA296455_UJ279751.txt
Batch1_PP_INFORM_MASTER	UF291262_UE277617.txt
Batch2_PP_MASTER	UF291262_UE277617.txt"""

    with open(correction_file_mapping_file, "w") as f:
        f.write(correction_file_mapping_file_content)
    return correction_file_mapping_file

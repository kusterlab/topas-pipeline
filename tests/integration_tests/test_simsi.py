from pathlib import Path
import topas_pipeline.simsi as simsi
from topas_pipeline.preprocess_tools import MQ_EVIDENCE_COLUMNS
from topas_pipeline.data_loaders.simsi_tmt_loader import SimsiTMTLoader
from topas_pipeline import config

CONFIG_FILE_PATH = './data/test_config.json'


def test_find_simsi_evidence_file():
    configs = config.load(CONFIG_FILE_PATH)
    results_folder = configs["results_folder"]
    simsi_folder = configs["simsi"]["simsi_folder"]
    
    data_type = 'FP'
    simsi_evidence_file = simsi.find_simsi_evidence_file(results_folder, simsi_folder, data_type)
    assert simsi_evidence_file is not None


def test_load_simsi_evidence_file_archived():
    configs = config.load(CONFIG_FILE_PATH)
    results_folder = configs["results_folder"]
    simsi_folder = configs["simsi"]["simsi_folder"]
    raw_data_location = configs["preprocessing"]["raw_data_location"]

    data_type = 'FP'
    evidence_files = [f'{raw_data_location}/Batch01_FP_CPTAC_BRCA/combined/txt/evidence.txt', f'{raw_data_location}/Batch02_FP_CPTAC_BRCA/combined/txt/evidence.txt']
    loader = SimsiTMTLoader(evidence_files, results_folder, simsi_folder, data_type)
    all_batches = loader.load_data(MQ_EVIDENCE_COLUMNS)
    assert len(all_batches) == 2


def test_find_matching_summaries_folder():
    """
    Test that meta_input_file_FP.tsv in a regular summaries folder can be found
    """
    configs = config.load(CONFIG_FILE_PATH)
    results_folder = configs["results_folder"]
    simsi_folder = configs["simsi"]["simsi_folder"]

    simsi_output_folder = Path(f'{simsi_folder}/simsi_output/FP')
    meta_input_file = Path(f'{results_folder}/meta_input_file_FP.tsv')
    summaries_folder = simsi.find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    results_folder_name = results_folder.split('/')[-1]
    assert summaries_folder == Path(f'{simsi_output_folder}/summaries_{results_folder_name}') 


def test_find_matching_summaries_folder_archived():
    """
    Test that meta_input_file_FP.tsv inside a zip archive can be found
    """
    configs = config.load(CONFIG_FILE_PATH)
    results_folder = configs["results_folder"]
    simsi_folder = configs["simsi"]["simsi_folder"]

    simsi_output_folder = Path(f'{simsi_folder}/simsi_output/FP')
    meta_input_file = Path(f'{results_folder}/meta_input_file_FP.tsv')
    summaries_folder = simsi.find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    results_folder_name = results_folder.split('/')[-1]
    assert str(summaries_folder) == str(Path(f'{simsi_output_folder}/summaries_{results_folder_name}'))


def test_find_matching_summaries_folder_not_found():
    """
    Test that trying to find an FP SIMSI search in the PP SIMSI folder will not return any results
    """
    configs = config.load(CONFIG_FILE_PATH)
    results_folder = configs["results_folder"]
    simsi_folder = configs["simsi"]["simsi_folder"]

    simsi_output_folder = Path(f'{simsi_folder}/simsi_output/PP')
    meta_input_file = Path(f'{results_folder}/meta_input_file_FP.tsv')
    summaries_folder = simsi.find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    assert summaries_folder is None
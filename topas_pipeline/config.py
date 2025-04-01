import json
from typing import List

from dataclasses import dataclass, field, asdict


@dataclass
class Simsi:
    run_simsi: bool = True
    simsi_folder: str = ""
    tmt_ms_level: str = "ms2"
    stringencies: int = 10
    tmt_requantify: bool = False
    maximum_pep: int = 1
    num_threads: int = 8
    correction_factor_mapping_file: str = ""


@dataclass
class Preprocessing:
    raw_data_location: str
    fasta_file: str
    picked_fdr: float = 0.01
    fdr_num_threads: int = 8
    program: str = ""
    entity: str = ""
    histologic_subtype: str = ""
    imputation: bool = True
    normalize_to_reference: bool = False
    debug: bool = False
    run_lfq: bool = False
    run_preprocessing: bool = True


@dataclass
class ClinicProc:
    pspFastaFile: str
    pspKinaseSubstrateFile: str
    pspAnnotationFile: str
    pspRegulatoryFile: str
    prot_baskets: str
    extra_kinase_annot: str = ""


@dataclass
class Report:
    samples_for_report: str = "all"
    drug_list_file: str = ""


@dataclass
class Portal:
    update: bool = False
    cohort: str = ""
    url: str = ""
    config: str = ""


@dataclass
class Slack:
    webhook_url: str = ""
    channel: str = ""


@dataclass
class Config:
    results_folder: str
    sample_annotation: str
    metadata_annotation: str
    raw_file_folders: List[str]
    simsi: Simsi
    preprocessing: Preprocessing
    clinic_proc: ClinicProc
    report: Report = Report()
    portal: Portal = Portal()
    data_types: List[str] = field(default_factory=lambda: ["fp", "pp"])
    slack: Slack = Slack()

    def __post_init__(self):
        self.simsi = Simsi(**self.simsi)
        self.preprocessing = Preprocessing(**self.preprocessing)
        self.clinic_proc = ClinicProc(**self.clinic_proc)
        if isinstance(self.report, dict):
            self.report = Report(**self.report)
        if isinstance(self.portal, dict):
            self.portal = Portal(**self.portal)
        if isinstance(self.slack, dict):
            self.slack = Slack(**self.slack)

    def asdict(self):
        return asdict(self)

    def asjson(self):
        return json.dumps(self.asdict(), indent=4)


def load(json_file: str):
    """Reads a JSON file and updates the provided Config dataclass instance."""
    with open(json_file, "r") as f:
        config_data = json.load(f)

    config = Config(**config_data)

    return config


# def merge_dicts(defaults, custom):
#     for key, value in custom.items():
#         if (
#             isinstance(value, dict)
#             and key in defaults
#             and isinstance(defaults[key], dict)
#         ):
#             defaults[key] = merge_dicts(defaults[key], value)  # Recursive merge
#         else:
#             defaults[key] = value  # Overwrite default value
#     return defaults


# def validate_sample_annot(configs):
#     # Check for sample + metadata annotation errors before starting pipeline, currently not used
#     # _ = prep.check_annot(configs.sample_annotation, configs.metadata_annotation, prep.in_metadata)
#     sample_annot_df = load_sample_annotation(configs.sample_annotation)
#     sample_annot_df.to_csv(
#         os.path.join(configs.results_folder, "sample_annotation.csv")
#     )

#     # check_config(configs)


# def check_config(configs):
#
#     # check data types
#     for data_type in configs.data_types:
#         if data_type.upper() not in ['FP', 'PP']:
#             raise ValueError(f'Data type {data_type} not accepted. Accepted data types are `fp` (full proteome) and `pp` (phospho).')
#
#     # check existence of file locations
#     if not os.path.exists(configs.sample_annotation):
#         raise ValueError(f'Sample annotation file {configs.sample_annotation} does not exist in this location.')
#     metadata_file = configs.metadata_annotation
#     if not os.path.isfile(metadata_file):
#             metadata_file = metadata_file.replace('Retrospective_MTBs_Evaluation',
#                                                   'Retrospective_MTBs_Evaluation/Metadata_excel_outdated')
# if not os.path.isfile(metadata_file):
#     raise ValueError(f'Metadata: {metadata_file} not found.')

from topas_pipeline.preprocess.file_processor.diann_file_loader import DiannFileLoader
from topas_pipeline.preprocess.file_processor.evidence_file_loader import (
    EvidenceFileLoader,
)
from topas_pipeline.preprocess.file_processor.ion_quant_file_loader import (
    IonQuantFileLoader,
)
from topas_pipeline.preprocess.file_processor.tmt_evidence_file_loader import (
    TMTEvidenceFileLoader,
)
from topas_pipeline.preprocess.file_processor.simsi_evidence_file_loader import (
    SimsiEvidenceFileLoader,
)


class ResultFileLoaderFactory:
    _file_loaders = {
        "lfq_diann": DiannFileLoader,
        "lfq_evidence": EvidenceFileLoader,
        "lfq_ionquant": IonQuantFileLoader,
        "tmt_evidence": TMTEvidenceFileLoader,
        "simsi_evidence": SimsiEvidenceFileLoader,
    }

    @classmethod
    def get_loader(cls, data_format: str):
        data_format = data_format.lower()
        if data_format not in cls._file_loaders:
            raise ValueError(f"Unsupported quant file format: {data_format}")
        loader = cls._file_loaders[data_format]
        return loader

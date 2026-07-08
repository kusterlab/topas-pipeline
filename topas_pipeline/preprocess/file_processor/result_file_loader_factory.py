from topas_pipeline.preprocess.file_processor.diann_file_loader import DiannFileLoader
from topas_pipeline.preprocess.file_processor.evidence_file_loader import (
    EvidenceFileLoader,
)
from topas_pipeline.preprocess.file_processor.ion_quant_file_loader import (
    IonQuantFileLoader,
)


class ResultFileLoaderFactory:
    _file_loaders = {
        "diann": DiannFileLoader,
        "evidence": EvidenceFileLoader,
        "ionquant": IonQuantFileLoader,
    }

    @classmethod
    def get_loader(cls, quant_file_format: str):
        if quant_file_format not in cls._file_loaders:
            raise ValueError(f"Unsupported quant file format: {quant_file_format}")
        loader = cls._file_loaders[quant_file_format]
        return loader

from topas_pipeline.preprocess.quant_data_loader.lfq_quant_data_loader import (
    LFQQuantDataLoader,
)
from topas_pipeline.preprocess.quant_data_loader.tmt_quant_data_loader import (
    TMTQuantDataLoader,
)
from topas_pipeline.preprocess.quant_data_loader.simsi_tmt_quant_data_loader import (
    SimsiTMTQuantDataLoader,
)

class QuantDataLoaderFactory:
    _loaders = {
        "LFQ": LFQQuantDataLoader,
        "TMT": TMTQuantDataLoader,
        "SIMSI": SimsiTMTQuantDataLoader,
    }

    @classmethod
    def get_loader(cls, quant_strategy: str):
        quant_strategy = quant_strategy.upper()
        if quant_strategy not in cls._loaders:
            raise ValueError(f"Unsupported quant strategy: {quant_strategy}")
        loader_class = cls._loaders[quant_strategy]
        return loader_class

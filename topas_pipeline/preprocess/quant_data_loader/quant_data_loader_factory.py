from topas_pipeline.preprocess.quant_data_loader.lfq_quant_data_loader import (
    LFQQuantDataLoader,
)


class QuantDataLoaderFactory:
    _loaders = {"LFQ": LFQQuantDataLoader}

    @classmethod
    def get_loader(cls, quant_strategy: str):
        if quant_strategy not in cls._loaders:
            raise ValueError(f"Unsupported quant strategy: {quant_strategy}")
        loader_class = cls._loaders[quant_strategy]
        return loader_class

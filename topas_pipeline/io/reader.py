from abc import ABC, abstractmethod
import pandas as pd
import pyarrow.parquet as pq

class Reader(ABC):
    def __init__(self, file_path: str):
        self.file_path = file_path

    @abstractmethod
    def read(self):
        raise NotImplementedError("Subclasses must implement this method")

    @abstractmethod
    def get_header(self):
        raise NotImplementedError("Subclasses must implement this method")


class TxtReader(Reader):
    def __init__(self, file_path: str, **kwargs):
        super().__init__(file_path)
        self.kwargs = kwargs

    def read(self, usecols=None):
        """Read the CSV file and return a pandas DataFrame"""
        self.usecols = usecols if usecols else None
        sep = self.kwargs.get("sep", "\t")
        return pd.read_csv(self.file_path, usecols=usecols, sep=sep)
    
    def get_header(self):
        """Get the header of the CSV file"""
        if self.usecols:
            return self.usecols
        return pd.read_csv(self.file_path, nrows=0, sep=self.kwargs.get("sep", "\t")).columns.tolist()

class ParquetReader(Reader):
    """Reader for parquet files"""
    def __init__(self, file_path: str, **kwargs):
        super().__init__(file_path)
        self.kwargs = kwargs

    def read(self, usecols=None):
        """Read the parquet file and return a pandas DataFrame"""
        self.usecols = usecols if usecols else None
        df = pq.read_table(self.file_path, columns=usecols).to_pandas()
        return df
    
    def get_header(self):
        """Get the header of the parquet file"""
        if self.usecols:
            return self.usecols
        return pq.read_metadata(self.file_path).schema.names

class ReaderFactory:
    _readers = {
        "csv": (TxtReader, {"sep": ","}),
        "tsv": (TxtReader, {"sep": "\t"}),
        "txt": (TxtReader, {"sep": "\t"}),
        "parquet": (ParquetReader, {})
    }

    @classmethod
    def get_reader(cls, file_path: str):
        file_extension = file_path.split(".")[-1].lower()
        if file_extension not in cls._readers:
            raise ValueError(f"Unsupported file type: {file_extension}")
        reader_class, default_kwargs = cls._readers.get(file_extension)
        return reader_class(file_path, **default_kwargs)


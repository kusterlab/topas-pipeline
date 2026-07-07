from topas_pipeline.io.reader import ReaderFactory, TxtReader, ParquetReader
import pytest
import pandas as pd


class TestReaderFactory:
    @pytest.fixture
    def file_paths(self, tmp_path):
        """Fixture to create temporary parquet and tsv files for testing"""
        parquet_file_path = tmp_path / "test.parquet"
        tsv_file_path = tmp_path / "test.tsv"

        data = {"Col1": [1, 4], "Col2": [2, 5], "Col3": [3, 6]}
        df = pd.DataFrame(data)

        df.to_parquet(parquet_file_path)
        df.to_csv(tsv_file_path, sep="\t", index=False)
        return {"parquet": str(parquet_file_path), "tsv": str(tsv_file_path)}

    @pytest.mark.parametrize(
        "file_path, expected_reader_class, expected_kwargs",
        [
            ("data.csv", TxtReader, {"sep": ","}),
            ("data.tsv", TxtReader, {"sep": "\t"}),
            ("data.txt", TxtReader, {"sep": "\t"}),
            ("data.parquet", ParquetReader, {}),
        ],
    )
    def test_factory_returns_correct_reader(
        self, file_path, expected_reader_class, expected_kwargs
    ):
        """Test that the ReaderFactory returns the correct reader class and kwargs based on file extension"""
        reader = ReaderFactory.get_reader(file_path)
        assert isinstance(reader, expected_reader_class)
        assert reader.kwargs == expected_kwargs

    def test_factory_raises_error_for_unsupported_file_type(self):
        with pytest.raises(ValueError) as exc:
            ReaderFactory.get_reader("data.xlsx")
        assert str(exc.value) == "Unsupported file type: xlsx"

    @pytest.mark.parametrize(
        "reader_class, filetype",
        [
            (TxtReader, "tsv"),
            (ParquetReader, "parquet"),
        ],
    )
    def test_read(self, reader_class, filetype, file_paths):
        """Test that the reader reads the file correctly and returns a DataFrame with the expected content"""
        reader = reader_class(file_paths[filetype])
        df = reader.read()
        header = reader.get_header()

        # Check header
        assert header == ["Col1", "Col2", "Col3"]

        # Check Dataframe content
        pd.testing.assert_frame_equal(
            df, pd.DataFrame({"Col1": [1, 4], "Col2": [2, 5], "Col3": [3, 6]})
        )

    @pytest.mark.parametrize(
        "usecols, out_cols",
        [
            (["Col1"], ["Col1"]),
            (["Col1", "Col2"], ["Col1", "Col2"]),
            (None, ["Col1", "Col2", "Col3"]),
        ],
    )
    def test_read_with_use_cols(self, usecols, out_cols, file_paths):
        """Test that the reader reads the file correctly with specified usecols and returns a DataFrame with the expected content"""
        for filetype in file_paths:
            reader = ReaderFactory.get_reader(file_paths[filetype])
            df = reader.read(usecols=usecols)
            assert df.columns.to_list() == out_cols

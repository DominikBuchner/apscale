"""Test oracle: the pandas based fasta_to_parquet that was used before the lightweight worker module."""
import gzip, hashlib
from pathlib import Path

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser


def oracle_fasta_to_parquet(fasta_path: str, output_folder: str, perform_hash: bool) -> None:
    """Function to stream fasta to parquet as tables in the form of
    [header, sequence, count]

    Args:
        fasta_path (str): Path to the fasta file to transform.
        output_folder (str): Output folder to write to
        perform_hash (bool): Wether to hash the sequences prior to parquet generation. Needed if no data aggregation step is used!
    """
    # extract the sample name from the fasta_path
    fasta_name = Path(fasta_path).with_suffix("").stem

    # read the fasta with the simple fasta parser to generate rows
    all_rows = []

    # add a dummy seq to each file to handle empty files greacefully
    all_rows.append([fasta_name, "dummy_hash", "dummy_seq", 0])

    with gzip.open(fasta_path, "rt") as data_stream:
        for header, seq in SimpleFastaParser(data_stream):
            header_data = header.split(";size=")
            seq_hash, seq_count = header_data[0], int(header_data[1].rstrip(";"))

            # normalize first, so that the hash does not depend on case or whitespace
            seq = seq.upper().strip()

            # only compute hash if needed
            if perform_hash:
                seq_hash = hashlib.sha3_256(seq.encode("ascii")).hexdigest()

            all_rows.append([fasta_name, seq_hash, seq, seq_count])

    # put into pandas dataframe
    fasta_data = pd.DataFrame(
        data=all_rows, columns=["sample", "hash", "sequence", "read_count"]
    )

    # add correct datatype annotations
    fasta_data = fasta_data.astype(
        {
            "sample": "string",
            "hash": "string",
            "sequence": "string",
            "read_count": "int64",
        }
    )

    # create an output path, save data to parquet to ingest into duckdb later
    output_path = Path(output_folder).joinpath(f"{fasta_name}.parquet.snappy")
    fasta_data.to_parquet(output_path)

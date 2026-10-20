"""Conversion of one fasta.gz file to a parquet table.

Kept in its own module on purpose: the joblib worker processes import the module of the function they run.
This module only imports what the conversion needs, so every worker stays small (the big libraries of
j_generate_read_table, like pandas and duckdb, are not loaded in the workers).
"""

import gzip, hashlib
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def read_fasta(stream):
    """Yield (header, sequence) for each record of a fasta stream, multi line sequences are joined.

    Args:
        stream: Text stream to read from.
    """
    header, parts = None, []
    for line in stream:
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts)
            header, parts = line[1:], []
        elif header is not None:
            parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def fasta_to_parquet(fasta_path: str, output_folder: str, perform_hash: bool) -> None:
    """Function to stream fasta to parquet as tables in the form of
    [header, sequence, count]

    Args:
        fasta_path (str): Path to the fasta file to transform.
        output_folder (str): Output folder to write to
        perform_hash (bool): Wether to hash the sequences prior to parquet generation. Needed if no data aggregation step is used!
    """
    # extract the sample name from the fasta_path
    fasta_name = Path(fasta_path).with_suffix("").stem

    # add a dummy seq to each file to handle empty files greacefully
    hashes, sequences, read_counts = ["dummy_hash"], ["dummy_seq"], [0]

    with gzip.open(fasta_path, "rt") as data_stream:
        for header, seq in read_fasta(data_stream):
            header_data = header.split(";size=")
            seq_hash, seq_count = header_data[0], int(header_data[1].rstrip(";"))

            # normalize first, so that the hash does not depend on case or whitespace
            seq = seq.upper().strip()

            # only compute hash if needed
            if perform_hash:
                seq_hash = hashlib.sha3_256(seq.encode("ascii")).hexdigest()

            hashes.append(seq_hash)
            sequences.append(seq)
            read_counts.append(seq_count)

    table = pa.table(
        {
            "sample": pa.array([fasta_name] * len(hashes), type=pa.string()),
            "hash": pa.array(hashes, type=pa.string()),
            "sequence": pa.array(sequences, type=pa.string()),
            "read_count": pa.array(read_counts, type=pa.int64()),
        }
    )

    # create an output path, save data to parquet to ingest into duckdb later
    pq.write_table(
        table, Path(output_folder).joinpath(f"{fasta_name}.parquet.snappy"), compression="snappy"
    )

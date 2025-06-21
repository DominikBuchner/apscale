import glob, gzip, pickle, shutil, datetime, os, math, subprocess, duckdb, time
import dask.dataframe as dd
import pyarrow as pa
from zict import File, Buffer, LRU, Func
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Align
import pandas as pd
from apscale.a_create_project import choose_input
from more_itertools import chunked
from joblib import Parallel, delayed
from tqdm import tqdm


def parse_fasta_data(input_path: str):
    """Generator function that takes an input path to a size annotated fasta file
    and yield the samples name, sequence hash, sequence and the sequence count.
    Data is yielded sequences wise

    Args:
        input_path (str): Input path to the fasta file.

    Yields:
        tuple: Yields a tuple with the sample name, sequence has, sequence and the sequence count.
    """
    # extract the sample name
    sample_name = Path(input_path).with_suffix("").with_suffix("").name

    # always yield one empty dummy sequence as first object so empty files are not missing from the read table
    yield sample_name, "empty_seq", "empty_seq", 0

    # extract the data directly from the fasta itself
    with gzip.open(input_path, "rt") as data_stream:
        for header, seq in SimpleFastaParser(data_stream):
            header_data = header.split(";size=")
            seq_hash, seq_count = header_data[0], int(header_data[1])
            seq = seq.upper().strip()
            yield sample_name, seq_hash, seq, seq_count


def generate_buffered_dict(save_path: str, buffer_size: int) -> dict:
    """Function to generate a disk backed dict that can (infinetly) grow.

    Args:
        save_path (str): Path to save the dict to.
        buffer_size (int): How many values to keep in memory before pushing to disk.

    Returns:
        dict: Dict that is backed up if needed
    """
    # define cold and hot storage
    cold = Func(pickle.dumps, pickle.loads, File(Path(save_path)))
    hot = LRU(n=buffer_size, d={})
    buffered_dict = Buffer(hot, cold, n=int(buffer_size * 0.99))

    return buffered_dict


def save_data_lines_to_hdf(save_path: str, data_lines: dict) -> dict:
    """Function to append data lines to the hdf store

    Args:
        save_path (str): Savepath of the output hdf file.
        data_lines (dict): Dict that's holding the data lines

    Returns: Returns an empty dict to free memory space
    """
    # convert the data liens dict to dataframe
    data_lines = pd.DataFrame.from_dict(
        data_lines, orient="index", columns=["hash_idx", "sample_idx", "read_count"]
    )

    # write the data to hdf format
    with pd.HDFStore(
        save_path, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            key="read_count_data", value=data_lines, format="table", data_columns=True
        )

    return {}


def add_dicts_to_hdf(
    save_path: str,
    seq_hash_to_idx: dict,
    seq_hash_to_seq: dict,
    sample_to_idx: dict,
    buffer_size: int,
) -> None:
    """Function to add the idx dictionarys to the hdf store.

    Args:
        save_path (str): Path to save the hdf to.
        seq_hash_to_idx (dict): Dict that maps hash values to idx values
        seq_hash_to_seq (dict): Dict that maps hash values to sequences
        sample_to_idx (dict): Dict that maps sample names to idx values
        buffer_size (int): Buffer size to use to loop over the dictionarys when writing to hdf
    """
    # loop over the chunks of seq hash to idx to push to hdf
    # write the chunk to the hdf store
    with pd.HDFStore(
        save_path, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        for chunk in chunked(seq_hash_to_idx.items(), buffer_size):
            chunk_df = pd.DataFrame(data=chunk, columns=["hash", "hash_idx"]).set_index(
                "hash_idx"
            )
            chunk_df["seq"] = chunk_df["hash"].map(seq_hash_to_seq)
            hdf_output.append(
                key="sequence_data",
                value=chunk_df,
                format="table",
                data_columns=True,
            )

    # write the chunk to the hdf store
    with pd.HDFStore(
        save_path, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        for chunk in chunked(sample_to_idx.items(), buffer_size):
            chunk_df = pd.DataFrame(
                data=chunk, columns=["sample_name", "sample_idx"]
            ).set_index("sample_idx")
            hdf_output.append(
                key="sample_data",
                value=chunk_df,
                format="table",
                data_columns=True,
            )


def index_data_to_hdf(project: str, input_files: str, buffer_size: int) -> str:
    """Function to index the dataset with buffered dicts and saving all data as dataframes to hdf.

    Args:
        project (str): Apscale project to work in.
        input_files (str): List of all input files to process
        buffer_size (int): How many values to keep in memory before flowing data to disk.

    Returns:
        str: Savepath of the hdf storage for easier access later on.
    """
    # generate savepaths for the buffered dicts first
    savepath_seq_hash_to_idx = Path(project).joinpath(
        "11_read_table", "temp", "seq_hash_to_idx"
    )
    savepath_seq_hash_to_seq = Path(project).joinpath(
        "11_read_table", "temp", "seq_hash_to_seq"
    )
    savepath_sample_to_idx = Path(project).joinpath(
        "11_read_table", "temp", "sample_to_idx"
    )

    # generate the buffered dicts
    seq_hash_to_idx = generate_buffered_dict(savepath_seq_hash_to_idx, buffer_size)
    seq_hash_to_seq = generate_buffered_dict(savepath_seq_hash_to_seq, buffer_size)
    sample_to_idx = generate_buffered_dict(savepath_sample_to_idx, buffer_size)

    # define all indeces prior to looping over the files
    hash_idx = 0
    sample_idx = 0
    data_idx = 0

    # spill over to hdf if too full, no lookups needed
    full_data = {}

    # define the hdf savepath, remove data from previous run to not append to it
    hdf_savename = Path(project).joinpath(
        "11_read_table",
        "data",
        "read_data_storage_{}.h5.lz".format(Path(project).stem.replace("_apscale", "")),
    )

    try:
        os.remove(hdf_savename)
    except FileNotFoundError:
        pass

    # loop over all input files
    for file in input_files:
        # loop over the data in the files line by line
        for sample_name, hash, seq, read_count in parse_fasta_data(file):
            # new hashes and samples only have to be added once a new one is found
            if hash not in seq_hash_to_idx:
                seq_hash_to_idx[hash] = hash_idx
                seq_hash_to_seq[hash] = seq
                # increase the hash idx
                hash_idx += 1
            if sample_name not in sample_to_idx:
                sample_to_idx[sample_name] = sample_idx
                sample_idx += 1
            # data has to be added at each iteration
            full_data[data_idx] = (
                seq_hash_to_idx[hash],
                sample_to_idx[sample_name],
                read_count,
            )
            data_idx += 1
            # as soon as buffer size data rows are collected, stream data to hdf
            if data_idx % buffer_size == 0:
                # add data to the hdf store
                save_data_lines_to_hdf(hdf_savename, full_data)
                full_data = {}
    else:
        save_data_lines_to_hdf(hdf_savename, full_data)

    # add the buffered dicts to hdf
    add_dicts_to_hdf(
        hdf_savename, seq_hash_to_idx, seq_hash_to_seq, sample_to_idx, buffer_size
    )

    # remove temporary savefiles
    shutil.rmtree(Path(project).joinpath("11_read_table", "temp"))

    # return to hdf savename for reuse / ordering / reading
    return hdf_savename


def generate_fasta(project: str, hdf_savename: str, chunksize: int) -> None:
    """Function to generate an ordered fasta file that holds all sequences of the dataset.

    Args:
        project (str): Apscale project to work in.
        hdf_savename (str): Path to the hdf file that holds the data.
        chunksize(int): Chunksize to use for saving memory
    """
    read_count_data = dd.read_hdf(
        hdf_savename,
        key="read_count_data",
        columns=["hash_idx", "read_count"],
        chunksize=chunksize,
    )

    # groupby hash_idx to order by abundance
    read_count_data = read_count_data.groupby(by=["hash_idx"]).sum().reset_index()

    # load the sequence data to generate the fasta file
    sequence_data = dd.read_hdf(
        hdf_savename, key="sequence_data", chunksize=chunksize
    ).reset_index()

    # merge on hash idx
    fasta_data = dd.merge(read_count_data, sequence_data, on="hash_idx")

    # sort the values by read count
    fasta_data = fasta_data.sort_values(by=["read_count"], ascending=False)

    # store a fasta order column for later grouping operations
    fasta_data = fasta_data.assign(fasta_order=1)
    fasta_data["fasta_order"] = fasta_data.fasta_order.cumsum() - 1

    # generate a savename for the fasta file
    fasta_savename = Path(project).joinpath(
        "11_read_table",
        "sequences_{}.fasta".format(Path(project).stem.replace("_apscale", "")),
    )

    # open the fasta file to write to
    with open(fasta_savename, "w") as out_stream:
        # write the fasta data
        for part in fasta_data.to_delayed():
            part = part.compute()
            for hash, seq in zip(part["hash"], part["seq"]):
                # skip the empty seq
                if hash != "empty_seq":
                    out_stream.write(f">{hash}\n{seq}\n")

    return fasta_data, fasta_savename


def generate_sequence_groups(
    project: str,
    hdf_savename: str,
    fasta_data: object,
    fasta_savename: str,
    group_threshold: float,
    chunksize: int,
):
    """Function to compute sequence groups with a given threshold (similar to OTUs).

    Args:
        project (str): Apscale project to work in.
        hdf_savename (str): Path to the data storage.
        fasta_data (object): Data from the fasta that was created for the sequences.
        fasta_savename (str): Path to the fasta.
        group_threshold (float): Threshold to group with. Float between 0 and 1.
        chunksize (int): Chunksize to process the resulting sequence groups with.

    Returns:
        _type_: _description_
    """
    # generate the vsearch matchfile for pairwise comparisons
    matchfile_path = Path(project).joinpath(
        "11_read_table",
        "data",
        "matchfile_{}.csv".format(str(group_threshold).replace(".", "")),
    )

    print(
        "{}: Computing matchfile with vsearch.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    with open(matchfile_path, "w") as output:
        # run vsearch and pipe stdout to the file
        subprocess.run(
            [
                "vsearch",
                "--usearch_global",
                fasta_savename,
                "--db",
                fasta_savename,
                "--self",
                "--id",
                str(group_threshold),
                "--userout",
                "-",
                "--quiet",
                "--userfields",
                "query+target+id",
                "--maxaccepts",
                str(0),
            ],
            stdout=output,
        )

    # ingest the data into a dask dataframe
    matchfile = dd.read_csv(
        matchfile_path,
        sep="\t",
        names=["query", "target", "pct_id"],
    )

    # add the hash idx and order to query
    matchfile = (
        dd.merge(
            left=matchfile,
            right=fasta_data[["hash_idx", "hash", "fasta_order"]],
            how="left",
            left_on="query",
            right_on="hash",
        )
        .drop(columns="hash")
        .rename(
            columns={
                "hash_idx": "hash_idx_query",
                "fasta_order": "fasta_order_query",
            }
        )
    )

    # add the hash idx and order to target
    matchfile = (
        dd.merge(
            left=matchfile,
            right=fasta_data[["hash_idx", "hash", "fasta_order"]],
            how="left",
            left_on="target",
            right_on="hash",
        )
        .drop(columns="hash")
        .rename(
            columns={
                "hash_idx": "hash_idx_target",
                "fasta_order": "fasta_order_target",
            }
        )
    )

    print(
        "{}: Creating group lookup table.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # stream data to duck db
    duckdb_savename = Path(project).joinpath(
        "11_read_table",
        "data",
        "group_lookup_table.duckdb",
    )

    # remove old lookup table
    if duckdb_savename.is_file():
        duckdb_savename.unlink()

    # connect to duckdb database
    db_connection = duckdb.connect(duckdb_savename)

    # stream the data to duck db
    for chunk_nr, data_chunk in enumerate(matchfile.to_delayed()):
        # compute the chunk
        data_chunk = data_chunk.compute()
        data_chunk = pa.Table.from_pandas(data_chunk)
        # create the table for the first chunk
        if chunk_nr == 0:
            db_connection.sql("CREATE TABLE matchfile AS SELECT * FROM data_chunk")
        else:
            db_connection.sql("INSERT INTO matchfile SELECT * FROM data_chunk")

    print(
        "{}: Grouping sequences.".format(datetime.datetime.now().strftime("%H:%M:%S"))
    )

    first = True
    sequence_group_data = []
    sequence_group_savename = Path(project).joinpath(
        "11_read_table", "data", "group_mapping.csv"
    )

    # try to remove the old sequence group mapping if it exists
    if sequence_group_savename.is_file():
        sequence_group_savename.unlink()

    # loop over the chunks from fasta_data and sort them into groups
    for fasta_data_chunk in fasta_data.to_delayed():
        # compute the chunk to generate the iterator
        fasta_data_chunk = fasta_data_chunk.compute()

        # loop over the rows of the datachunk
        for hash_idx, _, hash, seq, _ in fasta_data_chunk.itertuples(index=False):
            # look for the matches in the matches table
            all_hash_matches = db_connection.sql(
                f"SELECT * FROM matchfile WHERE hash_idx_query = {hash_idx}"
            ).to_df()
            # if it is the first add it to the groups found table
            if first:
                db_connection.sql(
                    "CREATE TABLE groups_found AS SELECT * from all_hash_matches"
                )
                # continue after the insertion, set first to false
                first = False
                sequence_group_data.append((hash_idx, hash, seq, hash_idx))
                continue

            # now the actual checking happens get all potential matches for the current hash from the groups found
            group_assignment = db_connection.sql(
                f"SELECT * FROM groups_found WHERE hash_idx_target = {hash_idx} ORDER BY fasta_order_query LIMIT 1"
            ).to_df()
            if group_assignment.empty:
                db_connection.sql(
                    "INSERT INTO groups_found SELECT * FROM all_hash_matches"
                )
                sequence_group_data.append((hash_idx, hash, seq, hash_idx))
            else:
                # extract the value from the group assignment
                sequence_group_data.append(
                    (
                        hash_idx,
                        hash,
                        seq,
                        group_assignment["hash_idx_query"].item(),
                    )
                )
        # flush to csv every 10k sequences
        if len(sequence_group_data) >= chunksize:
            # create a dataframe from sequence group data
            sequence_group_data = pd.DataFrame(
                sequence_group_data,
                columns=[
                    "hash_idx",
                    "hash",
                    "seq",
                    "group_idx",
                ],
            )
            sequence_group_data.to_csv(
                sequence_group_savename, sep="\t", index=False, mode="a", header=False
            )
            sequence_group_data = []
    # at the end of the loop flush the sequence group data
    else:
        # create a dataframe from sequence group data
        sequence_group_data = pd.DataFrame(
            sequence_group_data,
            columns=["hash_idx", "hash", "seq", "group_idx"],
        )
        sequence_group_data.to_csv(
            sequence_group_savename, sep="\t", index=False, mode="a", header=False
        )

    # finally close the db connection
    db_connection.close()

    # return sequence group savename to continue
    return sequence_group_savename


def generate_read_table(
    project: str,
    hdf_savename: str,
    to_excel: bool,
    to_parquet: bool,
    fasta_data: object,
) -> None:
    """Function to generate and save read tables to excel and or parquet

    Args:
        project (str): Apscale project to work in.
        hdf_savename (str): Path to the hdf storage that holds all data.
        to_excel (bool): Wether or not to write to excel.
        to_parquet (bool): Wether or not to write to parquet.
        fasta_data (object): fasta data as dask dataframe. holds the sequence data in correct order, so can be used for read_table_generation
    """
    formats = ["excel", "parquet"]

    # collect number of sequences and number of samples
    with pd.HDFStore(hdf_savename, mode="r") as store:
        number_of_sequences = store.get_storer("sequence_data").nrows
        number_of_samples = store.get_storer("sample_data").nrows

    # give user output, ignore the empty seq
    print(
        "{}: Dataset contains {} unique sequences and {} samples.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            number_of_sequences - 1,
            number_of_samples,
        )
    )

    # loop over both possibilities
    for format, setting in zip(formats, (to_excel, to_parquet)):
        # nothing to do if the setting is set to false
        if setting == False:
            continue
        else:
            # compute the batch size for the respective format
            if format == "excel":
                chunksize = 10_000_000 // number_of_sequences
            if format == "parquet":
                chunksize = 100_000_000 // number_of_sequences

        # collect all hashes and sequences
        all_sequence_data = fasta_data.compute()[["hash_idx", "hash", "seq"]]

        ## collect the sample data in chunks
        sample_data = dd.read_hdf(
            hdf_savename, key="sample_data", chunksize=chunksize
        ).reset_index()

        read_count_data = dd.read_hdf(hdf_savename, key="read_count_data")

        # load the data in chunks
        for part_nr, part in enumerate(sample_data.to_delayed()):
            # load the part to dask dataframe and compute it
            part_df = dd.from_delayed(part).compute()

            # create a cross product with current samples and all sequences
            sample_sequence_matrix = pd.merge(
                part_df[["sample_idx", "sample_name"]], all_sequence_data, how="cross"
            )

            # fetch the read counts for the current part
            read_counts_part = read_count_data[
                read_count_data["sample_idx"].isin(part_df["sample_idx"])
            ].compute()

            # merge the read counts into the sample sequence matrix
            sample_sequence_matrix = sample_sequence_matrix.merge(
                read_counts_part, on=["sample_idx", "hash_idx"], how="left"
            )

            # fill the missing values with 0
            sample_sequence_matrix["read_count"] = sample_sequence_matrix[
                "read_count"
            ].fillna(0)

            # pivot the table to create a classic samples x sequences OTU table
            wide_read_table = sample_sequence_matrix.pivot_table(
                index=["hash_idx", "hash", "seq"],
                columns="sample_name",
                values="read_count",
                fill_value=0,
            )

            # flatten column index, reset index for saving
            wide_read_table.columns.name = None
            wide_read_table = wide_read_table.reset_index()

            # sort the dataframe by the order in fasta data
            fasta_order = all_sequence_data[["hash_idx"]].reset_index(drop=True)
            fasta_order["order"] = range(len(fasta_order))

            # merge to the wide read_table for ordering
            wide_read_table = wide_read_table.merge(
                fasta_order, on="hash_idx", how="left"
            )

            # sort and drop the sorting column
            wide_read_table = wide_read_table.sort_values("order").drop(
                columns=["order", "hash_idx"]
            )

            # drop the last row as it contains the empty seq
            wide_read_table = wide_read_table.head(-1)

            # save to excel or parquet
            if format == "excel" and setting:
                savename = Path(project).joinpath(
                    "11_read_table",
                    "{}_read_table_part_{}.xlsx".format(
                        Path(project).stem.replace("_apscale", ""), part_nr
                    ),
                )

                # save to excel
                wide_read_table.to_excel(savename, index=False)

            if format == "parquet" and setting:
                savename = Path(project).joinpath(
                    "11_read_table",
                    "{}_read_table_part_{}.parquet.snappy".format(
                        Path(project).stem.replace("_apscale", ""), part_nr
                    ),
                )

                # save to excel
                wide_read_table.to_parquet(savename, index=False)


def sequence_groups_to_data_storage(
    project: str, sequence_group_savename: str, hdf_savename: str, chunksize: int
):
    # load the group mapping as a dask dataframe, break in 10 MB chunks
    group_mapping = dd.read_csv(
        sequence_group_savename,
        sep="\t",
        names=["hash_idx", "hash", "seq", "group_idx"],
        blocksize=10e6,
    )

    # load the read count data that needs to be updated with the group mappings
    read_count_data = dd.read_hdf(
        hdf_savename, key="read_count_data", chunksize=chunksize
    )

    # add the group idx to the read count data
    sequence_group_read_count_data = dd.merge(
        left=read_count_data,
        right=group_mapping[["hash_idx", "group_idx"]],
        how="left",
        left_on="hash_idx",
        right_on="hash_idx",
    )

    # save the group mapping in human readable format
    group_mapping = group_mapping.rename(
        columns={
            "hash_idx": "Original sequence ID",
            "group_idx": "Maps to group sequence ID",
        }
    )

    # define the output savename, chunk to multiple files if neccessary
    group_mapping_savename = Path(project).joinpath("11_read_table", "data")

    for idx, data_chunk in enumerate(group_mapping.to_delayed()):
        # save the chunk to excel
        data_chunk = data_chunk.compute()
        data_chunk.to_excel(
            group_mapping_savename.joinpath(f"group_mapping_log_part_{idx}.xlsx"),
            index=False,
        )

    # add the sequence groups to the hdf store, replace group idx first, drop duplicates, so only the groups are left
    group_mapping = group_mapping.drop(columns="Original sequence ID").rename(
        columns={"Maps to group sequence ID": "hash_idx"}
    )

    # reorder the columns to make it look like the original
    group_mapping = group_mapping[["hash_idx", "hash", "seq"]]

    # drop the duplicate values
    group_mapping = group_mapping.drop_duplicates(subset="hash_idx", keep="first")

    # give user output
    print(
        "{}: {} sequence groups found.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            group_mapping.shape[0].compute() - 1,
        )
    )

    # convert data types to write to hdf
    group_mapping["hash"], group_mapping["seq"] = group_mapping["hash"].astype(
        object
    ), group_mapping["seq"].astype(object)

    # save to hdf
    group_mapping.to_hdf(
        hdf_savename, key="sequence_group_data", mode="a", compute=True
    )

    # transform sequence group read count data to write the fasta and to push it to hdf
    sequence_group_read_count_data = (
        sequence_group_read_count_data[["group_idx", "sample_idx", "read_count"]]
        .groupby(by=["group_idx", "sample_idx"])
        .sum()
        .reset_index()
        .rename(columns={"group_idx": "hash_idx"})
    )

    # save to hdf
    sequence_group_read_count_data.to_hdf(
        hdf_savename, key="sequence_group_read_count_data", mode="a", compute=True
    )


def main(project=Path.cwd()):
    """Main function to initially create the read table data.

    Args:
        project (str, optional): Path to the apscale project to work in. Defaults to Path.cwd().
    """

    # collect the settings from the excel sheet
    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="11_read_table",
    )
    perform = settings["generate read table"].item()
    to_excel = settings["to excel"].item()
    to_parquet = settings["to parquet"].item()
    group_threshold = settings["sequence group threshold"].item()

    # find out which dataset to load
    prior_step = choose_input(project, "11_generate_read_table")

    # check that at least one hashing step has been performed
    if prior_step == "06_dereplication":
        ## give user output
        print(
            "{}: Please perform at least on data aggregation step before creating the read table.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        return None

    # collect the input for read table generation
    # generate a data folder to save the hdf store
    input_files = glob.glob(
        str(Path(project).joinpath(prior_step, "data", "*.fasta.gz"))
    )

    # generate a temp folder for saving the buffered dicts
    ## create temporal output folder
    for folder in ["temp", "data"]:
        try:
            os.mkdir(Path(project).joinpath("11_read_table", folder))
        except FileExistsError:
            pass

    print(
        "{}: Building data storage.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # index all the data and save as hdf store
    hdf_savename = index_data_to_hdf(project, input_files, 100_000)

    print(
        "{}: Saving sequences to fasta.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # sort the hdf store to generate the fasta file, generate the fasta file
    fasta_data, fasta_savename = generate_fasta(project, hdf_savename, 100_000)

    # create parquet and excel outputs from the hdf
    print(
        "{}: Generating read table(s).".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # generate the read tables if requests
    if perform:
        generate_read_table(project, hdf_savename, to_excel, to_parquet, fasta_data)

    # generate sequence groups
    if group_threshold >= 1:
        # user output
        print(
            "{}: Grouping threshold cannot be larger or equal to one. Step will be skipped.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
    elif group_threshold <= 0 or math.isnan(group_threshold):
        # user output
        print(
            "{}: Grouping threshold cannot 0 or smaller. Step will be skipped.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
    else:
        # user output
        print(
            "{}: Group threshold of {} detected. Sequence groups are calculated.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), group_threshold
            )
        )

        # calculate the sequence groups
        sequence_group_savename = generate_sequence_groups(
            project,
            hdf_savename,
            fasta_data,
            fasta_savename,
            group_threshold,
            chunksize=100_000,
        )

        # update the data storage, generate user log, potentially generate read tables and fasta
        sequence_groups_to_data_storage(
            project, sequence_group_savename, hdf_savename, chunksize=100_000
        )

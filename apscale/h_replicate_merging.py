import datetime, glob, gzip, hashlib
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
from collections import defaultdict
from apscale.a_create_project import choose_input
from Bio.SeqIO.FastaIO import SimpleFastaParser


def merge_replicates(
    input_files: tuple,
    output_name: str,
    minimum_presence: int,
    project=None,
    complevel=None,
):
    """Function to merge n replicate files to one output file.

    Args:
        input_files (tuple): Tuple holding up to n input files.
        output_name (str): Name of the output file to write.
        minimum_presence (int): Sequences have to be present in at least "minimum_presence" replicates.
        project (str, optional): Apscale project to work in. Defaults to None.
        complevel (str, optional): Compression level to use for the output file. Defaults to None.
    """
    # define the initial dict for storing the data
    seq_data = {}

    # count total input sequences and reads here
    total_input_seqs = 0
    total_input_reads = 0

    # hashes as keys, dict as values {seq, size, count}
    # extract the required data from the input file
    for input_file in input_files:
        with gzip.open(input_file, "rt") as in_stream:
            data = SimpleFastaParser(in_stream)
            for header, seq in data:
                # extract the id hash from the header
                header_data = header.split(";size=")
                # account for different size annotations from swarm and vsearch
                seq_id, size = header_data[0], int(header_data[1].split(";")[0])
                # recompute the hash values in case no denoising / clustering is performed
                seq_for_hashing = seq.upper().encode("ascii")
                hash = hashlib.sha3_256(seq_for_hashing).hexdigest()
                # add to the dict
                if hash not in seq_data:
                    seq_data[hash] = {"seq": seq, "size": size, "count": 1}
                    total_input_seqs += 1
                    total_input_reads += size
                else:
                    seq_data[hash]["size"] += size
                    seq_data[hash]["count"] += 1
                    total_input_reads += size

    # transform to dataframe for easier sorting
    output_data = pd.DataFrame.from_dict(data=seq_data, orient="index")

    # transform filtering and sorting if there is anything to sort
    if len(output_data.index) > 0:
        output_data = output_data.loc[output_data["count"] >= minimum_presence]
        output_data = output_data.sort_values(by=["size"], ascending=False)

        # get the data back to dict to easily write the output
        output_data = output_data.to_dict(orient="index")

    # define the output path
    output_path = Path(project).joinpath("09_replicate_merging", "data", output_name)

    # store the total output seqs and reads here
    total_output_seqs = 0
    total_output_reads = 0

    # write the output
    with gzip.open(output_path, "wt") as out_stream:
        for hash in output_data:
            seq, size = output_data[hash]["seq"], output_data[hash]["size"]
            out_stream.write(f">{hash};size={size}\n{seq}\n")
            total_output_seqs += 1
            total_output_reads += size

    # give user output
    # give user output
    print(
        "{}: {}: Combined {} files with {} input sequences and {} input reads resulting in {} output sequences and {} output reads".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            output_name,
            len(input_files),
            total_input_seqs,
            total_input_reads,
            total_output_seqs,
            total_output_reads,
        )
    )

    # return log_data
    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
    log_data = [Path(file).stem for file in input_files] + [
        output_name,
        finished,
        total_input_seqs,
        total_input_reads,
        total_output_seqs,
        total_output_reads,
    ]

    return log_data


def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Will merge replicates and only keep reads that can be found in n replicates.

    Args:
        project (str, optional): Apscale project to work in. Defaults to Path.cwd().
    """
    # collect variables from the settings file
    gen_settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="0_general_settings",
    )
    cores, comp_lvl = (
        gen_settings["cores to use"].item(),
        gen_settings["compression level"].item(),
    )

    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="09_replicate_merging",
    )
    perform = settings["perform replicate merging"].item()
    replicate_del = settings["replicate delimiter"].item()
    minimum_presence = settings["minimum replicate presence"].item()

    # only perform replicate merging if requested
    if perform:
        # select the correct step to pick the input files
        prior_step = choose_input(project, "09_replicate_merging")

        # gather input files
        input = glob.glob(str(Path(project).joinpath(prior_step, "data", "*.fasta.gz")))

        # put a mapping here in case no clustering or denoising is performed --> files have to be renamed accordingly then
        temporary_mapping = [[file_path, file_path] for file_path in input]

        # filter out appendices from previous steps
        appendices = [
            "_PE_trimmed_filtered_dereplicated",
            "_trimmed_filtered_dereplicated",
            "_filtered_dereplicated",
            "_dereplicated",
        ]

        for pair in temporary_mapping:
            for apx in appendices:
                if apx in pair[0]:
                    pair[0] = "{}.gz".format(pair[0].replace(apx, ""))
                    break

        # # generate a defaultdict that maps output_names to input pairs, triplets, etc, depending on the number of replicates
        matches_dict = defaultdict(list)

        # # find matching replicates
        for output_name, file_path in temporary_mapping:
            common_name = "_".join(Path(output_name).name.split("_")[:-1])
            common_name = f"{common_name}.fasta.gz"
            matches_dict[common_name].append(file_path)

        matches_dict = dict(matches_dict)

        # switch keys and values to have direct input for replicate merging
        matches_dict = {tuple(value): key for key, value in matches_dict.items()}

        # perform the replicate merging
        log_data = Parallel(n_jobs=cores)(
            delayed(merge_replicates)(
                matched_files,
                matches_dict[matched_files],
                minimum_presence,
                project,
                comp_lvl,
            )
            for matched_files in matches_dict.keys()
        )

        # extract the number of input files for computing the columns
        nr_of_input_files = len(log_data[0]) - 6
        columns = [f"input file {i}" for i in range(1, nr_of_input_files + 1)] + [
            "output file",
            "finished at",
            "total input sequences",
            "total input reads",
            "total output sequences",
            "total output reads",
        ]
        # put everything in a log dataframe
        log_df = pd.DataFrame(log_data, columns=columns)
        log_df = log_df.sort_values(by=["output file"])

        # write the log file
        log_df.to_excel(
            Path(project).joinpath(
                "09_replicate_merging", "Logfile_09_replicate_merging.xlsx"
            ),
            index=False,
            sheet_name="09_replicate_merging",
        )

        ## add log to the project report
        with pd.ExcelWriter(
            Path(project).joinpath(
                "Project_report_{}.xlsx".format(
                    Path(project).name.replace("_apscale", "")
                )
            ),
            mode="a",
            if_sheet_exists="replace",
            engine="openpyxl",
        ) as writer:
            log_df.to_excel(writer, sheet_name="09_replicate_merging", index=False)

    else:
        ## give user output
        print(
            "{}: Replicate merging is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

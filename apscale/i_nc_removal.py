import datetime, glob, gzip, hashlib
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
from collections import defaultdict
from apscale.a_create_project import choose_input
from Bio.SeqIO.FastaIO import SimpleFastaParser


def max_reads_ncs(file_paths: list) -> dict:
    """Funtion to collect the maximum of reads for each individual sequence from a list of negative control fasta files

    Args:
        file_paths (list): List of all negative controls

    Returns:
        dict: Dict with hash: count for all sequences found in the negative controls
    """
    # define the dict to hold the data
    ncs_dict = {}

    for file in file_paths:
        with gzip.open(file, "rt") as in_stream:
            for header, seq in SimpleFastaParser(in_stream):
                # extract the size from the header
                # extract the id hash from the header
                header_data = header.split(";size=")
                # account for different size annotations from swarm and vsearch
                seq_id, size = header_data[0], int(header_data[1].rstrip(";"))
                # recompute the hash values in case no denoising / clustering is performed
                seq_for_hashing = seq.upper().encode("ascii")
                hash = hashlib.sha3_256(seq_for_hashing).hexdigest()

                # check if the header is in the dict already, then compare, else add
                if hash in ncs_dict.keys():
                    if size > ncs_dict[hash]:
                        ncs_dict[hash] = size
                else:
                    ncs_dict[hash] = size

    return ncs_dict


def substract_nc_reads(
    input_file: str, output_name: str, nc_reads: dict, project=None, comp_lvl=None
):
    """Function to substract the reads found in the negative controls from all remaining input files.

    Args:
        input_file (str): Input file to remove reads from.
        output_name (str): Output name to write to.
        nc_reads (dict): Dict holding the reads found in the negative controls.
        project (_type_, optional): Apscale project to work in. Defaults to None.
        complevel (_type_, optional): Compression level of the output files. Defaults to None.
    """
    # create the output path
    output_path = Path(project).joinpath("10_nc_removal", "data", output_name)
    substracted_sequences = 0
    removed_reads = 0

    # open the input file
    with gzip.open(input_file, "rt") as in_stream:
        with gzip.open(output_path, "wt", compresslevel=comp_lvl) as out_stream:
            for header, seq in SimpleFastaParser(in_stream):
                # extract the id hash from the header
                header_data = header.split(";size=")
                # account for different size annotations from swarm and vsearch
                seq_id, size = header_data[0], int(header_data[1].split(";")[0])
                # recompute the hash values in case no denoising / clustering is performed
                seq_for_hashing = seq.upper().encode("ascii")
                hash = hashlib.sha3_256(seq_for_hashing).hexdigest()
                if hash in nc_reads.keys():
                    # extract the size for substraction
                    new_size = size - nc_reads[hash]
                    substracted_sequences += 1
                    if new_size > 0:
                        # write output
                        removed_reads += nc_reads[hash]
                        out_stream.write(f">{hash};size={new_size}\n{seq}\n")
                    else:
                        removed_reads += size
                else:
                    # just write
                    out_stream.write(f">{hash};size={size}\n{seq}\n")

        ## give user output
        print(
            "{}: {}: {} reads substracted from a total of {} sequences.".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                output_name,
                removed_reads,
                substracted_sequences,
            )
        )
        # return log_data
        finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
        return output_name, finished, removed_reads, substracted_sequences


def main(project=Path.cwd()):
    """Main function to run the negative control removal script

    Args:
        project (str, optional): Path to the apscale project to work in. Defaults to Path.cwd().
    """
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
        sheet_name="10_nc_removal",
    )
    perform = settings["perform nc removal"].item()
    nc_prefix = settings["negative control prefix"].item()

    # only perform is set to true
    if perform:
        # select the correct step to pick the input files
        prior_step = choose_input(project, "10_nc_removal")

        # gather input files
        input = glob.glob(str(Path(project).joinpath(prior_step, "data", "*.fasta.gz")))
        # filter out appendices from previous steps
        appendices = [
            "_PE_trimmed_filtered_dereplicated",
            "_trimmed_filtered_dereplicated",
            "_filtered_dereplicated",
            "_dereplicated",
        ]

        input_output_mapping = []

        for file_path in input:
            for apx in appendices:
                if apx in file_path:
                    output_name = "{}".format(file_path.replace(apx, ""))
                    output_name = Path(output_name).name
                    input_output_mapping.append([file_path, output_name])
                    break
                else:
                    output_name = Path(file_path).name
                    input_output_mapping.append([file_path, output_name])
                    break

        # filter out the negative controls from the mapping
        negative_controls = [
            input_output[0]
            for input_output in input_output_mapping
            if input_output[1].startswith(nc_prefix)
        ]

        # collect the maximum of reads for each sequence in the negative controls
        nc_controls_reads = max_reads_ncs(negative_controls)

        # filter the negative controls from the actual input files, they don't have to be written
        input_output_mapping = [
            input_output
            for input_output in input_output_mapping
            if input_output[0] not in negative_controls
        ]

        # substract those sizes from each of the remaining fasta files if there are any reads in the nc controls
        log_data = Parallel(n_jobs=cores)(
            delayed(substract_nc_reads)(
                input, output, nc_controls_reads, project, comp_lvl
            )
            for input, output in input_output_mapping
        )

        # put everything in a log dataframe
        log_df = pd.DataFrame(
            log_data,
            columns=[
                "output file",
                "finished at",
                "total substracted reads",
                "sequences with substraction",
            ],
        )
        log_df = log_df.sort_values(by=["output file"])

        # write the log file
        log_df.to_excel(
            Path(project).joinpath("10_nc_removal", "Logfile_10_nc_removal.xlsx"),
            index=False,
            sheet_name="10_nc_removal",
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
            log_df.to_excel(writer, sheet_name="10_nc_removal", index=False)
    else:
        ## give user output
        print(
            "{}: Negative control removal is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

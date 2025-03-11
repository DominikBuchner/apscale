import datetime, glob
import pandas as pd
from pathlib import Path
from collections import defaultdict
from apscale.a_create_project import choose_input


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
    pass

    # hashes as keys, dict as values {seq, size, count}


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
        for matched_files in matches_dict.keys():
            merge_replicates(
                matched_files,
                matches_dict[matched_files],
                minimum_presence,
                project,
                comp_lvl,
            )
    else:
        ## give user output
        print(
            "{}: Replicate merging is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

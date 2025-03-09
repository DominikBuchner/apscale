import datetime, glob
import pandas as pd
from pathlib import Path
from collections import defaultdict
from apscale.a_create_project import choose_input


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
        # TODO

        # generate a defaultdict that maps output_names to input pairs, triplets, etc, depending on the number of replicates
        matches_dict = defaultdict(list)

        # find matching replicates
        for file_path in input:
            common_name = "_".join(Path(file_path).name.split("_")[:-1])
            common_name = f"{common_name}.fasta.gz"
            matches_dict[common_name].append(file_path)

        matches_dict = dict(matches_dict)
        print(matches_dict)
    else:
        ## give user output
        print(
            "{}: Replicate merging is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

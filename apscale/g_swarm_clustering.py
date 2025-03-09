import datetime, os, glob
from pathlib import Path
import pandas as pd
from apscale.a_create_project import choose_input


# main function of the swarm clustering script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Will perform swarm clustering on the individual files. If run without denoising will also perform chimera removal and
    exchange the fasta headers with sha 256 hashes for easier processing.

    Args:
        project (str, optional): Path to the apscale project. Defaults to Path.cwd().
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
        sheet_name="08_swarm_clustering",
    )
    perform = settings["perform swarm clustering"].item()

    # only perform clustering if needed
    if perform:
        # select the correct step tp pick the input files
        prior_step = choose_input(project, "08_swarm_clustering")

        # create temporal output folder
        try:
            os.mkdir(Path(project).joinpath("08_swarm_clustering", "temp"))
        except FileExistsError:
            pass

        # gather input files
        input = glob.glob(str(Path(project).joinpath(prior_step, "data", "*.fasta.gz")))

        # since swarm does not work with gzipped data, unzip files first

    else:
        ## give user output
        print(
            "{}: Swarm clustering is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

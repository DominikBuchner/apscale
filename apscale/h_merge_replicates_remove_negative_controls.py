import pandas as pd
from pathlib import Path


def main(project=Path.cwd()) -> None:
    """Main function to merge replicates and remove negative controls from the dataset

    Args:
        project (str, optional): Path to the apscale project. Defaults to Path.cwd().
    """
    # collect the settings from the settings file
    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="7_replicate_negative_controls",
    )

    # define the settings
    replicate_merging, replicate_del, minimum_presence, substract_ncs, nc_prefix = (
        bool(settings["merge replicates"].item()),
        settings["replicate delimiter"].item(),
        int(settings["minimum replicate presence"].item()),
        bool(settings["substract negative controls"].item()),
        settings["negative control prefix"].item(),
    )

    # merge replicates first (optional)
    
    # remove negative controls afterward (optional)

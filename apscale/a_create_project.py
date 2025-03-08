import os, psutil, datetime, gzip
import pandas as pd
from pathlib import Path


def create_project(project_name):
    """Create a new metabarcoding pipeline project with all subfolders."""

    ## try to create the project folder
    try:
        os.mkdir("{}_apscale".format(project_name))
    except FileExistsError:
        print("A project with that name already exists. Please try another name.")
        return None

    ## generate the subfolder structure
    subfolders = [
        "01_raw_data/data",
        "02_demultiplexing/data",
        "03_PE_merging/data",
        "04_primer_trimming/data",
        "05_quality_filtering/data",
        "06_dereplication/data",
        "07_denoising/data",
        "08_swarm_clustering/data",
        "09_replicate_merging/data",
        "10_nc_removal/data",
        "11_read_table",
    ]

    subfolders = [
        Path("{}_apscale".format(project_name)).joinpath(subfolder)
        for subfolder in subfolders
    ]

    for folder in subfolders:
        os.makedirs(folder)

    # generate and populate the settings file, add the project name to the settings file name
    with pd.ExcelWriter(
        Path("{}_apscale".format(project_name)).joinpath(
            "Settings_{}.xlsx".format(Path(project_name).name)
        ),
        mode="w",
        engine="openpyxl",
    ) as writer:
        ## write the 0_general_settings sheet
        df_0 = pd.DataFrame(
            [[int(psutil.cpu_count() - 2), 6]],
            columns=["cores to use", "compression level"],
        )

        df_0.to_excel(writer, sheet_name="0_general_settings", index=False)

        ## write the 03_PE_merging sheet
        df_3 = pd.DataFrame(
            [[25, 199, 5]], columns=["maxdiffpct", "maxdiffs", "minovlen"]
        )

        df_3.to_excel(writer, sheet_name="03_PE_merging", index=False)

        ## write the 04_primer_trimming sheet
        df_4 = pd.DataFrame(
            [["", "", "False"]],
            columns=["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", "anchoring"],
        )

        df_4.to_excel(writer, sheet_name="04_primer_trimming", index=False)

        ## write the 05_quality_filtering sheet
        df_5 = pd.DataFrame(
            [[1, "", ""]], columns=["maxEE", "min length", "max length"]
        )

        df_5.to_excel(writer, sheet_name="05_quality_filtering", index=False)

        ## write the 06_dereplication sheet
        df_6 = pd.DataFrame([[3]], columns=["minimum sequence abundance"])

        df_6.to_excel(writer, sheet_name="06_dereplication", index=False)

        ## write the 07_denoising sheet
        df_7 = pd.DataFrame(
            [["True", 2, "absolute", 4]],
            columns=["perform denoising", "alpha", "threshold type", "size threshold"],
        )

        df_7.to_excel(writer, sheet_name="07_denoising", index=False)

        ## write the 08_swarm clustering sheet
        df_8 = pd.DataFrame(
            [["True"]],
            columns=["perform swarm clustering"],
        )

        df_8.to_excel(writer, sheet_name="08_swarm_clustering", index=False)

        ## write the 09_replicate merging sheet
        df_9 = pd.DataFrame(
            [["True", "_", 2]],
            columns=[
                "perform replicate merging",
                "replicate delimiter",
                "minimum replicate presence",
            ],
        )

        df_9.to_excel(writer, sheet_name="09_replicate_merging", index=False)

        ## write the 10_nc removal sheet
        df_10 = pd.DataFrame(
            [["True", "NC_"]],
            columns=[
                "perform nc substration",
                "negative control prefix",
            ],
        )

        df_10.to_excel(writer, sheet_name="10_nc_removal", index=False)

        ## write the 11 read table sheet
        df_11 = pd.DataFrame(
            [["True", "False", "True", 100]],
            columns=[
                "generate read table",
                "to excel",
                "to_parquet",
                "sequence identity for grouping",
            ],
        )

        df_11.to_excel(writer, sheet_name="11_read_table", index=False)

    ## give user output
    print(
        '{}: "{}_apscale" created as a new project folder.'.format(
            datetime.datetime.now().strftime("%H:%M:%S"), project_name
        )
    )


def empty_file(file_path: str) -> bool:
    """Function to check if an input file is empty. Works on gzip compressed files and regular text files.

    Args:
        file_path (str): Path to the file to be checked

    Returns:
        bool: Returns True if the input file is emtpy, else False
    """
    # convert string to file path
    file_path = Path(file_path)

    # try gzip file first
    try:
        with gzip.open(file_path, "rb") as f:
            data = f.read(1)
            if len(data) == 0:
                return True
            else:
                return False
    except (gzip.BadGzipFile, OSError):
        if os.stat(file_path).st_size == 0:
            return True
        else:
            return False

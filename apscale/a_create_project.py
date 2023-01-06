import os, psutil, datetime
import pandas as pd
from pathlib import Path


def create_project(project_name):
    """Create a new metabarcoding pipeline project with all subfolders"""

    ## try to create the project folder
    try:
        os.mkdir("{}_apscale".format(project_name))
    except FileExistsError:
        print("A project with that name already exists. Please try another name.")
        return None

    ## generate the subfolder structure
    subfolders = [
        "1_raw data/data",
        "2_demultiplexing/data",
        "3_PE_merging/data",
        "4_primer_trimming/data",
        "5_quality_filtering/data",
        "6_dereplication_pooling/data/dereplication",
        "6_dereplication_pooling/data/pooling",
        "7_otu_clustering/data",
        "8_denoising/data",
        "9_lulu_filtering/otu_clustering/data",
        "9_lulu_filtering/denoising/data",
    ]

    subfolders = [
        Path("{}_apscale".format(project_name)).joinpath(subfolder)
        for subfolder in subfolders
    ]

    for folder in subfolders:
        os.makedirs(folder)

    # generate and populate the settings file
    with pd.ExcelWriter(
        Path("{}_apscale".format(project_name)).joinpath("Settings.xlsx"),
        mode="w",
        engine="openpyxl",
    ) as writer:
        # ## write the 3_PE_merging sheet
        df_0 = pd.DataFrame(
            [[int(psutil.cpu_count() - 2), 6]],
            columns=["cores to use", "compression level"],
        )

        df_0.to_excel(writer, sheet_name="0_general_settings", index=False)

        ## write the 3_PE_merging sheet
        df_3 = pd.DataFrame(
            [[25, 199, 5]], columns=["maxdiffpct", "maxdiffs", "minovlen"]
        )

        df_3.to_excel(writer, sheet_name="3_PE_merging", index=False)

        ## write the 4_primer_trimming sheet
        df_4 = pd.DataFrame(
            [["", "", "False"]],
            columns=["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", "anchoring"],
        )

        df_4.to_excel(writer, sheet_name="4_primer_trimming", index=False)

        ## write the 5_quality_filtering sheet
        df_5 = pd.DataFrame(
            [[1, "", ""]], columns=["maxEE", "min length", "max length"]
        )

        df_5.to_excel(writer, sheet_name="5_quality_filtering", index=False)

        ## write the 6_dereplication_pooling sheet
        df_6 = pd.DataFrame([[4]], columns=["min size to pool"])

        df_6.to_excel(writer, sheet_name="6_dereplication_pooling", index=False)

        ## write the 7_otu_clustering sheet
        df_7 = pd.DataFrame([[97, "True"]], columns=["pct id", "to excel"])

        df_7.to_excel(writer, sheet_name="7_otu_clustering", index=False)

        ## write the 8_denoising sheet
        df_8 = pd.DataFrame([[2, 8, "True"]], columns=["alpha", "minsize", "to excel"])

        df_8.to_excel(writer, sheet_name="8_denoising", index=False)

        ## write the 8_denoising sheet
        df_9 = pd.DataFrame(
            [[84, 95, 1, "True"]],
            columns=[
                "minimum similarity",
                "minimum relative cooccurence",
                "minimum ratio",
                "to excel",
            ],
        )

        df_9.to_excel(writer, sheet_name="9_lulu_filtering", index=False)

    ## give user output
    print(
        '{}: "{}" created as a new project.'.format(
            datetime.datetime.now().strftime("%H:%M:%S"), project_name
        )
    )

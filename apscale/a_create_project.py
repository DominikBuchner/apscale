import os, psutil, datetime
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
        "1_raw_data/data",
        "2_demultiplexing/data",
        "3_PE_merging/data",
        "4_primer_trimming/data",
        "5_quality_filtering/data",
        "6_dereplication/data",
        "7_denoising/data",
        "8_esv_table/data",
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

        ## write the 7_denoising sheet
        df_6 = pd.DataFrame([[2, 4, "True"]], columns=["alpha", "minsize"])

        df_6.to_excel(writer, sheet_name="6_denoising", index=False)

    ## give user output
    print(
        '{}: "{}_apscale" created as a new project folder.'.format(
            datetime.datetime.now().strftime("%H:%M:%S"), project_name
        )
    )

import os, subprocess, pickle, datetime, re, glob, gzip, openpyxl, shutil
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from io import StringIO
from joblib import Parallel, delayed
from Bio.SeqIO.FastaIO import SimpleFastaParser
from openpyxl.utils.dataframe import dataframe_to_rows
from Bio.SeqIO.FastaIO import SimpleFastaParser
from apscale.a_create_project import empty_file


# function to perform the remapping of individual read files to individual esv outputs
def fasta_to_esv_tab(esv_file, project=None):
    """Function to remap the individual files tothe generated ESVs."""

    # extract the sample name from the file of the esvs since it is already shortened in the previous step
    sample_name_out = Path(esv_file).with_suffix("").with_suffix("").name

    # parse header and abundance directly from the fasta file if it is not empty
    if not empty_file(esv_file):
        fasta_data = [
            (header.split(";")[0], int(header.split("size=")[-1]))
            for header, _ in SimpleFastaParser(gzip.open(esv_file, "rt"))
        ]
    else:
        fasta_data = []

    # generate the esv tab
    esv_tab = pd.DataFrame(data=fasta_data, columns=["ID", sample_name_out])

    # save the time when finished
    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
    esvs = len(esv_tab)
    reads = esv_tab[sample_name_out].sum()

    ## give user output
    print(
        "{}: {}: Found {} ESVs with a sum of {} reads.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            sample_name_out,
            esvs,
            reads,
        )
    )

    ## pickle log data first for log generation
    with open(
        Path(project).joinpath(
            "8_esv_table", "temp", "{}_log.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump(
            [
                sample_name_out,
                finished,
                esvs,
                reads,
            ],
            log,
        )

    ## pickle otu tab dataframes for otu table generation
    with open(
        Path(project).joinpath(
            "8_esv_table", "temp", "{}_esv_tab.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump(esv_tab, log)


# function to generate the final esv table
def generate_esv_table(project=None):
    """Function to generate the final esv table from the output of the previous step."""
    # find all esv tabs
    esv_tabs = glob.glob(
        str(Path(project).joinpath("8_esv_table", "temp", "*esv_tab.pkl"))
    )
    # load the pickled esv tabs
    esv_tabs = [pickle.load(open(tab_file, "rb")) for tab_file in esv_tabs]
    # set the ID as index
    esv_tabs = [tab.set_index("ID") for tab in esv_tabs]

    # generate the final esv table
    with pd.option_context("future.no_silent_downcasting", True):
        esv_table = pd.concat(esv_tabs, axis=1, join="outer").fillna(0)

    # sort the data with the most abundant esv on top
    esv_table = esv_table.copy()
    esv_table["sum"] = esv_table.sum(axis=1)
    esv_table = esv_table.sort_values(by=["sum"], axis=0, ascending=False)
    esv_table = esv_table.drop(labels=["sum"], axis=1)

    # add unique ID columns
    esv_table = esv_table.rename_axis("unique_ID").reset_index()

    # add an temporary id column
    esv_table.insert(
        0,
        "temporary_ID",
        ["ESV_{}".format(i) for i, _ in enumerate(esv_table["unique_ID"], start=1)],
    )

    # add the sequence data from the fasta file and directly generate the fasta file
    fasta_dict = {}

    # collect all fasta files
    fasta_files = glob.glob(
        str(Path(project).joinpath("7_denoising", "data", "*.fasta.gz"))
    )

    # polulate the dict with the sequences from the fasta file
    for fasta_file in fasta_files:
        for header, seq in SimpleFastaParser(gzip.open(fasta_file, "rt")):
            # parse the header without size
            header = header.split(";")[0]
            if not header in fasta_dict.keys():
                fasta_dict[header] = seq

    # add the sequences to the final ESV table
    esv_table["Seq"] = esv_table["unique_ID"].map(fasta_dict)

    # save the esv table in excel and parquet format
    wb = openpyxl.Workbook(write_only=True)
    ws = wb.create_sheet("ESV table")

    ## save the output line by line for optimized memory usage
    for row in tqdm(
        dataframe_to_rows(esv_table, index=False, header=True),
        total=len(esv_table.index),
        desc="{}: Lines written to ESV table".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        ),
        unit=" lines",
    ):
        ws.append(row)

    ## save the output (otu table)
    print(
        "{}: Saving the ESV table to excel. This may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    wb.save(
        Path(project).joinpath(
            "8_esv_table",
            "{}_ESV_table.xlsx".format(Path(project).stem.replace("_apscale", "")),
        )
    )
    wb.close()

    print(
        "{}: ESV table saved to {}.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            Path(project).joinpath(
                "8_esv_table",
                "{}_ESV_table.xlsx".format(Path(project).stem.replace("_apscale", "")),
            ),
        )
    )

    print(
        "{}: Saving the ESV table to parquet. This may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    esv_table.to_parquet(
        Path(project).joinpath(
            "8_esv_table",
            "{}_ESV_table.parquet.snappy".format(
                Path(project).stem.replace("_apscale", "")
            ),
        ),
        index=False,
    )

    print(
        "{}: ESV table saved to {}.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            Path(project).joinpath(
                "8_esv_table",
                "{}_ESV_table.parquet.snappy".format(
                    Path(project).stem.replace("_apscale", "")
                ),
            ),
        )
    )

    # generate the final fasta file for taxonomic identification
    with open(
        Path(project).joinpath(
            "8_esv_table",
            "{}_ESVs.fasta".format(Path(project).stem.replace("_apscale", "")),
        ),
        "w",
    ) as out_stream:
        for unique_id, seq in zip(esv_table["unique_ID"], esv_table["Seq"]):
            out_stream.write(">{}\n{}\n".format(unique_id, seq))


# main function for the esv table generation script
def main(project=Path.cwd()):
    """Main function of the script. No default values. Will generate an ESV table from the output of
    the two previous steps."""

    # create a temporal output folder
    try:
        os.mkdir(Path(project).joinpath("8_esv_table", "temp"))
    except FileExistsError:
        pass

    # collect variables from the settings file
    gen_settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="0_general_settings",
    )
    cores = gen_settings["cores to use"].item()

    # gather input files for remapping and esv files
    esv_files = sorted(
        glob.glob(str(Path(project).joinpath("7_denoising", "data", "*.fasta.gz")))
    )

    # give user output
    print(
        "{}: Converting {} input fasta files to ESV tabs.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(esv_files)
        )
    )

    # run the remapping in parallel
    Parallel(n_jobs=cores)(
        delayed(fasta_to_esv_tab)(esv_file, project=project) for esv_file in esv_files
    )

    # give user output
    print(
        "{}: Constructing the final ESV table.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # construct the final ESV table
    generate_esv_table(project=project)

    # write the log
    summary_logs = glob.glob(
        str(Path(project).joinpath("8_esv_table", "temp", "*_log.pkl"))
    )
    summary = [pickle.load(open(line, "rb")) for line in summary_logs]

    log_df = pd.DataFrame(
        summary,
        columns=[
            "File",
            "finished at",
            "ESVs",
            "reads",
        ],
    )
    log_df = log_df.sort_values(by="File")
    log_df.to_excel(
        Path(project).joinpath("8_esv_table", "Logfile_8_esv_table.xlsx"),
        index=False,
        sheet_name="8_esv_table",
    )

    ## add log to the project report
    with pd.ExcelWriter(
        Path(project).joinpath(
            "Project_report_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        mode="a",
        if_sheet_exists="replace",
        engine="openpyxl",
    ) as writer:
        log_df.to_excel(writer, sheet_name="8_esv_table", index=False)

    # remove temporary files
    shutil.rmtree(Path(project).joinpath("8_esv_table", "temp"))

    # give user output
    print("{}: Analysis finished.".format(datetime.datetime.now().strftime("%H:%M:%S")))


if __name__ == "__main__":
    main()

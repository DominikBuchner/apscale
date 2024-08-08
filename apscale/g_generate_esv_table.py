import os, subprocess, pickle, datetime, re, glob, gzip, openpyxl, shutil
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from io import StringIO
from joblib import Parallel, delayed
from Bio.SeqIO.FastaIO import SimpleFastaParser
from openpyxl.utils.dataframe import dataframe_to_rows


# function to perform the remapping of individual read files to individual esv outputs
def remap_files_to_esvs(dereplicated_file, esv_file, project=None):
    """Function to remap the individual files tothe generated ESVs."""

    # extract the sample name from the file of the esvs since it is already shortened in the previous step
    sample_name_out = Path(esv_file).with_suffix("").with_suffix("").name

    # run vsearch --search_exact to remap the individual files to the generated ESVs per file
    # capture log and directly picke the output as dataframe for read table generation
    f = subprocess.run(
        [
            "vsearch",
            "--search_exact",
            Path(dereplicated_file),
            "--db",
            Path(esv_file),
            "--output_no_hits",
            "--maxhits",
            "1",
            "--otutabout",
            "-",
            "--quiet",
            "--threads",
            "1",
            "--log",
            Path(project).joinpath(
                "8_esv_table", "temp", "{}_mapping_log.txt".format(sample_name_out)
            ),
        ],
        capture_output=True,
    )

    # directly parse the output to a pandas dataframe for ESV table generation
    esv_tab = pd.read_csv(StringIO(f.stdout.decode("ascii", errors="ignore")), sep="\t")

    ## handle empty outputs correctly
    if not esv_tab.empty:
        esv_tab = esv_tab.set_axis(["ID", sample_name_out], axis=1, copy=False)
    else:
        esv_tab = pd.DataFrame(data=[], columns=["ID", sample_name_out])

    ## collect number of esvs from the output, pickle to logs
    with open(
        Path(project).joinpath(
            "8_esv_table", "temp", "{}_mapping_log.txt".format(sample_name_out)
        )
    ) as log_file:
        content = log_file.read()
        try:
            exact_matches = re.findall("Matching query sequences: (\d+)", content)[0]
        except IndexError:
            exact_matches = re.findall(
                "Matching unique query sequences: (\d+)", content
            )[0]
        version = re.findall("vsearch ([\w\.]*)", content)[0]
        finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## give user output
    print(
        "{}: {}: {} exact matches found ({} reads).".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            sample_name_out,
            exact_matches,
            esv_tab[sample_name_out].sum(),
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
                version,
                exact_matches,
                esv_tab[sample_name_out].sum(),
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

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("8_esv_table", "temp"))

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
    dereplicated_files = sorted(
        glob.glob(str(Path(project).joinpath("6_dereplication", "data", "*.fasta.gz")))
    )
    esv_files = sorted(
        glob.glob(str(Path(project).joinpath("7_denoising", "data", "*.fasta.gz")))
    )

    # give user output
    print(
        "{}: Starting to individually remap {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(dereplicated_files)
        )
    )

    # run the remapping in parallel
    Parallel(n_jobs=cores)(
        delayed(remap_files_to_esvs)(dereplicated_file, esv_file, project=project)
        for dereplicated_file, esv_file in zip(dereplicated_files, esv_files)
    )

    # give user output
    print(
        "{}: Constructing the final ESV table.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # construct the final ESV table
    generate_esv_table(project=project)

    # give user output
    print("{}: Analysis finished.".format(datetime.datetime.now().strftime("%H:%M:%S")))


if __name__ == "__main__":
    main()

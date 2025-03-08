import subprocess, gzip, datetime, glob, os, pickle, shutil, re
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed
from apscale.a_create_project import empty_file


## dereplication function to dereplicate a gzipped fasta file
def dereplication(file, project=None, comp_lvl=None, minimum_seq_abundance=None):
    """Function to dereplicate a gzipped fasta file. Abundance annotations will be
    written to the output fasta."""

    ## extract the filename from the sample path / name and convert to output name
    ## create an output path to write to
    sample_name_out = "{}_dereplicated.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )
    output_path = Path(project).joinpath("06_dereplication", "data", sample_name_out)

    ## run vsearch --derep_fulllength to dereplicate the file
    ## use --log because for some reason no info is written to stderr with this command
    ## relabel to handle different sequencing runs in the otu clustering - seems to only happen if data is downloaded from SRA
    ## write stdout to uncompressed output at runtime
    if not empty_file(file):
        with open(output_path.with_suffix(""), "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--fastx_uniques",
                    Path(file),
                    "--minuniquesize",
                    str(minimum_seq_abundance),
                    "--fastaout",
                    "-",
                    "--quiet",
                    "--fasta_width",
                    str(0),
                    "--log",
                    Path(project).joinpath(
                        "06_dereplication", "temp", "{}.txt".format(sample_name_out)
                    ),
                    "--sizeout",
                    "--relabel",
                    "seq:",
                ],
                stdout=output,
            )

        ## compress the output, remove uncompressed output
        with open(output_path.with_suffix(""), "rb") as in_stream, gzip.open(
            output_path, "wb", comp_lvl
        ) as out_stream:
            shutil.copyfileobj(in_stream, out_stream)
        os.remove(output_path.with_suffix(""))

        ## collect processed and passed reads from the log file
        with open(
            Path(project).joinpath(
                "06_dereplication", "temp", "{}.txt".format(sample_name_out)
            )
        ) as log_file:
            content = log_file.read()
            seqs, unique_seqs = (
                int(re.findall(r"(\d+) seqs", content)[0]),
                int(re.findall(r"(\d+) unique sequences", content)[0]),
            )
            version = re.findall("vsearch ([\w\.]*)", content)[0]
    else:
        with gzip.open(output_path, "wb"):
            seqs, unique_seqs, version = 0, 0, "empty input"

    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## give user output
    print(
        "{}: {}: {} sequences dereplicated into {} unique sequences.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            sample_name_out,
            seqs,
            unique_seqs,
        )
    )

    ## temporarily pickle output for the log file
    with open(
        Path(project).joinpath(
            "06_dereplication", "temp", "{}.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump([sample_name_out, finished, version, seqs, unique_seqs], log)


## main function to call the script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Calls the dereplication over all files, merges them into one file and dereplicates
    the merged file again. Project defaults so current working directory."""

    ## collect variables from the settings file
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

    # collect settings specific for that step
    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="06_dereplication",
    )

    min_abundance = settings["minimum sequence abundance"].item()

    ## collect input files from quality filtering step
    input = glob.glob(
        str(Path(project).joinpath("05_quality_filtering", "data", "*.fasta.gz"))
    )

    print(
        "{}: Starting to dereplicate {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(input)
        )
    )

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath("06_dereplication", "temp"))
    except FileExistsError:
        pass

    ## parallelize the dereplication
    Parallel(n_jobs=cores)(
        delayed(dereplication)(
            file,
            project=project,
            comp_lvl=comp_lvl,
            minimum_seq_abundance=min_abundance,
        )
        for file in input
    )

    ## write log for the dereplication from pkl logs
    summary_logs = glob.glob(
        str(Path(project).joinpath("06_dereplication", "temp", "*.pkl"))
    )
    summary = [pickle.load(open(line, "rb")) for line in summary_logs]

    log_df = pd.DataFrame(
        summary,
        columns=[
            "File",
            "finished at",
            "program version",
            "processed sequences",
            "unique sequences",
        ],
    )
    log_df = log_df.sort_values(by="File")
    log_df.to_excel(
        Path(project).joinpath("06_dereplication", "Logfile_06_dereplication.xlsx"),
        index=False,
        sheet_name="06_dereplication",
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
        log_df.to_excel(writer, sheet_name="06_dereplication", index=False)

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("06_dereplication", "temp"))


if __name__ == "__main__":
    main()

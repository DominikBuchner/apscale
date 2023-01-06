import subprocess, gzip, datetime, glob, os, pickle, openpyxl, shutil, re
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed

## dereplication function to dereplicate a gzipped fasta file
def dereplication(file, project=None, comp_lvl=None, min_size=None):
    """Function to dereplicate a gzipped fasta file. Abundance annotations will be
    written to the output fasta."""

    ## extract the filename from the sample path / name and convert to output name
    ## create an output path to write to
    sample_name_out = "{}_dereplicated.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )
    output_path = Path(project).joinpath(
        "6_dereplication_pooling", "data", "dereplication", sample_name_out
    )

    ## run vsearch --derep_fulllength to dereplicate the file
    ## use --log because for some reason no info is written to stderr with this command
    ## relabel to handle different sequencing runs in the otu clustering - seems to only happen if data is downloaded from SRA
    ## write stdout to uncompressed output at runtime
    with open(output_path.with_suffix(""), "w") as output:
        f = subprocess.run(
            [
                "vsearch",
                "--fastx_uniques",
                Path(file),
                "--fastaout",
                "-",
                "--quiet",
                "--fasta_width",
                str(0),
                "--log",
                Path(project).joinpath(
                    "6_dereplication_pooling", "temp", "{}.txt".format(sample_name_out)
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

    ## generate a new output path to stream min_size dereplicated sequences that are used for pooling and clustering
    ## do so only if min_size if bigger than 1, else use dereplicated files from step before
    if min_size > 1:
        output_path = Path(project).joinpath(
            "6_dereplication_pooling", "data", "pooling", sample_name_out
        )

        with open(output_path.with_suffix(""), "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--fastx_uniques",
                    Path(file),
                    "--fastaout",
                    "-",
                    "--quiet",
                    "--fasta_width",
                    str(0),
                    "--log",
                    Path(project).joinpath(
                        "6_dereplication_pooling",
                        "temp",
                        "{}.txt".format(sample_name_out),
                    ),
                    "--sizeout",
                    "--relabel",
                    "seq:",
                    "--minuniquesize",
                    str(min_size),
                ],
                stdout=output,
            )

    ## collect processed and passed reads from the log file
    with open(
        Path(project).joinpath(
            "6_dereplication_pooling", "temp", "{}.txt".format(sample_name_out)
        )
    ) as log_file:
        content = log_file.read()
        seqs, unique_seqs = (
            re.findall("(\d+) seqs, ", content)[0],
            re.findall("(\d+) unique sequences, ", content)[0],
        )
        version = re.findall("vsearch ([\w\.]*)", content)[0]
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
            "6_dereplication_pooling", "temp", "{}.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump([sample_name_out, finished, version, seqs, unique_seqs], log)


## function to pool the dereplicated reads, pooled reads are dereplicated again
def pooling(file_list, project=None, comp_lvl=None, min_size=None):
    """Function to pool the dereplicated reads, needs a path to dereplicated reads folder
    and a project to work in. Both are passed by the main function."""

    ## write the output file, depending on the value of min_size, delete files that are only used for pooling
    ## and no longer needed for mapping
    if min_size > 1:
        with gzip.open(
            Path(project).joinpath(
                "6_dereplication_pooling",
                "data",
                "pooling",
                "pooled_sequences.fasta.gz",
            ),
            "wb",
            comp_lvl,
        ) as pool:
            number_of_files = len(file_list)
            for i in range(len(file_list)):
                with open(file_list[i], "rb") as to_pool:
                    shutil.copyfileobj(to_pool, pool)
                os.remove(file_list[i])
                print(
                    "{}: Added {} of {} files to the pooled sequences.".format(
                        datetime.datetime.now().strftime("%H:%M:%S"),
                        i + 1,
                        number_of_files,
                    )
                )

    else:
        with gzip.open(
            Path(project).joinpath(
                "6_dereplication_pooling",
                "data",
                "pooling",
                "pooled_sequences.fasta.gz",
            ),
            "wb",
            comp_lvl,
        ) as pool:
            number_of_files = len(file_list)
            for i in range(len(file_list)):
                with gzip.open(file_list[i], "rb") as to_pool:
                    shutil.copyfileobj(to_pool, pool)
                print(
                    "{}: Added {} of {} files to the pooled sequences.".format(
                        datetime.datetime.now().strftime("%H:%M:%S"),
                        i + 1,
                        number_of_files,
                    )
                )

    ## dereplicate the pool with minuniquesize = 2
    print(
        "{}: Dereplicating the pooled sequences for clustering and denoising.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    ## run vsearch --derep_fulllength to dereplicate the file
    output_path = Path(project).joinpath(
        "6_dereplication_pooling",
        "data",
        "pooling",
        "pooled_sequences_dereplicated.fasta.gz",
    )

    with open(output_path.with_suffix(""), "w") as output:
        f = subprocess.run(
            [
                "vsearch",
                "--fastx_uniques",
                Path(project).joinpath(
                    "6_dereplication_pooling",
                    "data",
                    "pooling",
                    "pooled_sequences.fasta.gz",
                ),
                "--fastaout",
                "-",
                "--quiet",
                "--fasta_width",
                str(0),
                "--sizein",
                "--sizeout",
                "--minuniquesize",
                str(2),
            ],
            stdout=output,
        )

    ## compress the output, remove uncompressed output
    with open(output_path.with_suffix(""), "rb") as in_stream, gzip.open(
        output_path, "wb", comp_lvl
    ) as out_stream:
        shutil.copyfileobj(in_stream, out_stream)
    os.remove(output_path.with_suffix(""))


## main function to call the script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Calls the dereplication over all files, merges them into one file and dereplicates
    the merged file again. Project defaults so current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="0_general_settings"
    )
    cores, comp_lvl = (
        gen_settings["cores to use"].item(),
        gen_settings["compression level"].item(),
    )

    settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="6_dereplication_pooling"
    )
    min_size = settings["min size to pool"].item()

    ## collect input files from quality filtering step
    input = glob.glob(
        str(Path(project).joinpath("5_quality_filtering", "data", "*.fasta.gz"))
    )

    print(
        "{}: Starting to dereplicate {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(input)
        )
    )

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath("6_dereplication_pooling", "temp"))
    except FileExistsError:
        pass

    ## parallelize the dereplication
    Parallel(n_jobs=cores)(
        delayed(dereplication)(
            file, project=project, comp_lvl=comp_lvl, min_size=min_size
        )
        for file in input
    )

    ## write log for the dereplication from pkl logs
    summary_logs = glob.glob(
        str(Path(project).joinpath("6_dereplication_pooling", "temp", "*.pkl"))
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
        Path(project).joinpath(
            "6_dereplication_pooling", "Logfile_6_dereplication.xlsx"
        ),
        index=False,
        sheet_name="6_dereplication",
    )

    ## add log to the project report
    with pd.ExcelWriter(
        Path(project).joinpath("Project_report.xlsx"),
        mode="a",
        if_sheet_exists="replace",
        engine="openpyxl",
    ) as writer:
        log_df.to_excel(writer, sheet_name="6_dereplication", index=False)

    ## pool the dereplicated files, if minsize is set to a value > 1 use those for files for pooling
    if min_size > 1:
        files = glob.glob(
            str(
                Path(project).joinpath(
                    "6_dereplication_pooling", "data", "pooling", "*.fasta"
                )
            )
        )
    else:
        files = glob.glob(
            str(
                Path(project).joinpath(
                    "6_dereplication_pooling", "data", "dereplication", "*.fasta.gz"
                )
            )
        )

    pooling(files, project=project, comp_lvl=comp_lvl, min_size=min_size)

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("6_dereplication_pooling", "temp"))


if __name__ == "__main__":
    main()

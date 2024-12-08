import subprocess, gzip, datetime, pickle, glob, os, shutil, re, sys
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed
from apscale.a_create_project import empty_file


## quality filtering function to quality filter the specified file
def quality_filtering(
    file, project=None, comp_lvl=None, maxee=None, min_length=None, max_length=None
):
    """Function to apply quality filtering to a gzipped fastq file. Outputs a gzipped fasta file
    with quality filtered reads. Filtered reads will be discarded."""

    ## extract the filename from the sample path / name and convert to output name
    ## create an output path to write to
    sample_name_out = "{}_filtered.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )
    output_path = Path(project).joinpath("5_quality_filtering", "data", sample_name_out)

    # run vsearch --fastq_filter to apply the quality filtering
    # use --log because for some reason no info is written to stderr with this command
    # write stdout to uncompressed output at runtime, write stderr to a log file
    # run only if input is not empty
    if not empty_file(file):
        with open(output_path.with_suffix(""), "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--fastq_filter",
                    Path(file),
                    "--fastaout",
                    "-",
                    "--quiet",
                    "--fasta_width",
                    str(0),
                    "--log",
                    Path(project).joinpath(
                        "5_quality_filtering", "temp", "{}.txt".format(sample_name_out)
                    ),
                    "--fastq_maxee",
                    str(maxee),
                    "--fastq_minlen",
                    str(min_length),
                    "--fastq_maxlen",
                    str(max_length),
                    "--fastq_qmax",
                    str(64),
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
                "5_quality_filtering", "temp", "{}.txt".format(sample_name_out)
            )
        ) as log_file:
            content = log_file.read()
            kept, discarded = (
                int(re.findall(r"(\d+) sequences kept", content)[0]),
                int(re.findall(r"(\d+) sequences discarded", content)[0]),
            )
            reads = int(kept) + int(discarded)
            version = re.findall("vsearch ([\w\.]*)", content)[0]
    else:
        # generate data for logging if input is empty
        with gzip.open(output_path, "wb"):
            kept, discarded, reads, version = 0, 0, 0, "empty input"

    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## Give user output, if 0 reads are the output handle Zero division exception
    try:
        print(
            "{}: {}: {} of {} reads passed quality filtering ({:.2f}%)".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name_out,
                kept,
                reads,
                int(kept) / int(reads) * 100,
            )
        )
    except ZeroDivisionError:
        print(
            "{}: {}: {} of {} reads passed quality filtering ({:.2f}%)".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name_out,
                0,
                reads,
                0,
            )
        )

    ## temporarily pickle output for the log file, get vsearch version
    with open(
        Path(project).joinpath(
            "5_quality_filtering", "temp", "{}.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump([sample_name_out, finished, version, reads, kept], log)


## main function to call the script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    If default values are desired no arguments are required. Default working directory
    is the current working directory."""

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

    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="5_quality_filtering",
    )

    # check is the settings are setup correctly
    if pd.isna(settings).values.any():
        print(
            "{}: At least on filtering length is missing from the settings file.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        sys.exit()

    maxee, min_length, max_length = (
        settings["maxEE"].item(),
        settings["min length"].item(),
        settings["max length"].item(),
    )

    ## collect the input files from primer trimming step
    input = glob.glob(
        str(Path(project).joinpath("4_primer_trimming", "data", "*.fastq.gz"))
    )

    print(
        "{}: Starting to quality filter {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(input)
        )
    )

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath("5_quality_filtering", "temp"))
    except FileExistsError:
        pass

    ## parallelize the quality filtering
    Parallel(n_jobs=cores)(
        delayed(quality_filtering)(
            file,
            project=project,
            comp_lvl=comp_lvl,
            maxee=maxee,
            min_length=min_length,
            max_length=max_length,
        )
        for file in input
    )

    ## write logfile from pkl log, remove single logs after
    summary_logs = glob.glob(
        str(Path(project).joinpath("5_quality_filtering", "temp", "*.pkl"))
    )
    summary = [pickle.load(open(line, "rb")) for line in summary_logs]

    log_df = pd.DataFrame(
        summary,
        columns=[
            "File",
            "finished at",
            "program version",
            "processed reads",
            "passed reads",
        ],
    )
    log_df = log_df.sort_values(by="File")
    log_df.to_excel(
        Path(project).joinpath(
            "5_quality_filtering", "Logfile_5_quality_filtering.xlsx"
        ),
        index=False,
        sheet_name="5_quality_filtering",
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
        log_df.to_excel(writer, sheet_name="5_quality_filtering", index=False)

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("5_quality_filtering", "temp"))


if __name__ == "__main__":
    main()

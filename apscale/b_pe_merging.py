import subprocess, gzip, glob, pickle, datetime, os, shutil, re, sys
import pandas as pd
from pathlib import Path
from demultiplexer2 import find_file_pairs
from joblib import Parallel, delayed
from tqdm import tqdm
from apscale.a_create_project import empty_file


## file pair: matching forward and reverse reads, project: folder to write to
def pe_merge(
    file_pair,
    project=None,
    comp_lvl=None,
    maxdiffpct=None,
    maxdiffs=None,
    minovlen=None,
):
    """Function to merge to gzipped fastq.gz files via vsearch. Output will be
    a gzipped file containing the merged reads."""

    ## extract the filename from the sample path / name and convert to output name
    ## create an output path to write to
    sample_name_out = "{}_PE.fastq.gz".format(
        "_".join(Path(file_pair[0]).name.split("_")[:-1])
    )
    output_path = Path(project).joinpath("3_PE_merging", "data", sample_name_out)
    log_path = Path(project).joinpath(
        "3_PE_merging", "temp", "{}_log.txt".format(sample_name_out)
    )

    # perform computation on all non empty files
    if not empty_file(file_pair[0]):
        ## write stdout to uncompressed output at runtime, write stderr to a log file
        with open(output_path.with_suffix(""), "w") as output:
            ## run vsearch --fastq_mergepairs to merge the file pair
            f = subprocess.run(
                [
                    "vsearch",
                    "--fastq_mergepairs",
                    Path(file_pair[0]),
                    "--reverse",
                    Path(file_pair[1]),
                    "--fastqout",
                    "-",
                    "--quiet",
                    "--fastq_maxdiffpct",
                    str(maxdiffpct),
                    "--fastq_maxdiffs",
                    str(maxdiffs),
                    "--fastq_minovlen",
                    str(minovlen),
                    "--fastq_allowmergestagger",
                    "--threads",
                    str(1),
                    "--log",
                    Path(log_path),
                ],
                stdout=output,
            )

        ## compress the output, remove uncompressed output
        with open(output_path.with_suffix(""), "rb") as in_stream, gzip.open(
            output_path, "wb", comp_lvl
        ) as out_stream:
            shutil.copyfileobj(in_stream, out_stream)
        os.remove(output_path.with_suffix(""))

        ## read run info from log file
        with open(
            Path(project).joinpath(
                "3_PE_merging", "temp", "{}_log.txt".format(sample_name_out)
            ),
            "rt",
        ) as log_file:
            content = log_file.read()

            reads, merged = (
                re.findall(r"(\d+)  Pairs", content)[0],
                re.findall(r"(\d+)  Merged", content)[0],
            )
    # create an empty output file, generate data needed for logging
    else:
        with gzip.open(output_path, "wb", comp_lvl):
            reads, merged = 0, 0

    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## Give user output, if 0 reads are the output handle Zero division exception
    try:
        print(
            "{}: {}: {} of {} reads merged ({:.2f}%)".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name_out,
                merged,
                reads,
                int(merged) / int(reads) * 100,
            ),
        )
    except ZeroDivisionError:
        print(
            "{}: {}: {} of {} reads merged ({:.2f}%)".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name_out,
                0,
                reads,
                0,
            )
        )

    ## temporarily pickle output for the log file, get vsearch version
    if not empty_file(file_pair[0]):
        f = subprocess.run(["vsearch", "--version"], capture_output=True)
        version = f.stderr.decode("ascii", errors="ignore")
        version = re.findall("vsearch ([\w\.]*)", version)[0]
    else:
        version = "empty input"

    with open(
        Path(project).joinpath(
            "3_PE_merging", "temp", "{}.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump([sample_name_out, finished, version, reads, merged], log)


def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the Settings file.
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
        sheet_name="3_PE_merging",
    )
    maxdiffpct, maxdiffs, minovlen = (
        settings["maxdiffpct"].item(),
        settings["maxdiffs"].item(),
        settings["minovlen"].item(),
    )

    ## collect all files to merge, find matching file pairs
    input = sorted(
        glob.glob(str(Path(project).joinpath("2_demultiplexing", "data", "*.fastq.gz")))
    )

    print(
        "{}: Finding all file pairs in your input. This may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    pairs, singles = find_file_pairs.main(input)

    # give user output if singles are found
    if singles:
        print(
            "{}: Some input files do not have a matching file pair. Please check your input.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        sys.exit()

    print(
        "{}: Found {} matching file pairs in {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(pairs), len(input)
        )
    )
    print(
        "{}: Starting to merge forward and reverse reads.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    ## create folder for temporal output files
    try:
        os.mkdir(Path(project).joinpath("3_PE_merging", "temp"))
    except FileExistsError:
        pass

    ## parallelize the PE merging, find out how many cores to use
    Parallel(n_jobs=cores)(
        delayed(pe_merge)(
            pair,
            project=project,
            comp_lvl=comp_lvl,
            maxdiffpct=maxdiffpct,
            maxdiffs=maxdiffs,
            minovlen=minovlen,
        )
        for pair in pairs
    )

    ## write the log file from pkl log, remove logs after
    summary_logs = glob.glob(
        str(Path(project).joinpath("3_PE_merging", "temp", "*.pkl"))
    )
    summary = [pickle.load(open(line, "rb")) for line in summary_logs]

    ## generate the output dataframe for PE merging
    log_df = pd.DataFrame(
        summary,
        columns=[
            "File",
            "finished at",
            "program version",
            "processed reads",
            "merged reads",
        ],
    )
    log_df = log_df.sort_values(by="File")
    log_df.to_excel(
        Path(project).joinpath("3_PE_merging", "Logfile_3_PE_merging.xlsx"),
        index=False,
        sheet_name="3_PE_merging",
    )

    ## create the general logfile
    log_df.to_excel(
        Path(project).joinpath(
            "Project_report_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        index=False,
        sheet_name="3_PE merging",
    )

    ## remove temporally saved logs from single files
    shutil.rmtree(Path(project).joinpath("3_PE_merging", "temp"))


if __name__ == "__main__":
    main()

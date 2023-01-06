import subprocess, datetime, pickle, glob, os, shutil, sys
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_letters
from io import StringIO
from joblib import Parallel, delayed

## function to trim primers of reads of the specified file
def primer_trimming(file, project=None, p5_primer=None, p7_primer=None, anchoring=None):
    """Function to remove primer sequences from gzipped file via cutadapt. Outputs
    gzipped file with all reads where primers were removed. Untrimmed reads will be
    discarded."""

    ## extract the filename from the sample path / name and convert to output name
    sample_name_out = "{}_trimmed.fastq.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )

    ## if anchoring is True change the cutadapt call
    if anchoring:
        adapter = "^{}...{}".format(p5_primer, str(Seq(p7_primer).reverse_complement()))
    else:
        adapter = "{}...{}".format(p5_primer, str(Seq(p7_primer).reverse_complement()))

    ## run catadapt
    f = subprocess.run(
        [
            "cutadapt",
            "-a",
            adapter,
            "-o",
            str(Path(project).joinpath("4_primer_trimming", "data", sample_name_out)),
            file,
            "--discard-untrimmed",
            "--cores=1",
            "--report=minimal",
        ],
        capture_output=True,
    )

    ## collect processed reads from stderror for the logfile,
    ## handle exception for empty outputs
    log_df = pd.read_csv(StringIO(f.stdout.decode("ascii", errors="ignore")), sep="\t")
    reads, cut_reads = log_df["in_reads"].item(), log_df["out_reads"].item()
    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## give user output
    try:
        print(
            "{}: {}: primers trimmed for {} of {} reads ({:.2f}%)".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name_out,
                cut_reads,
                reads,
                cut_reads / reads * 100,
            )
        )
    except ZeroDivisionError:
        print(
            "{}: {}: primers trimmed for {} of {} reads ({:.2f}%)".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name_out,
                0,
                reads,
                0,
            )
        )

    ## get remaining log information, pickle temporarly to write the log after successfull finish
    py_v = sys.version_info
    py_v = ".".join((str(py_v.major), str(py_v.minor), str(py_v.micro)))

    cutadapt_v = subprocess.run(["cutadapt", "--version"], capture_output=True)
    cutadapt_v = cutadapt_v.stdout.decode("ascii", errors="ignore").rstrip()

    with open(
        Path(project).joinpath(
            "4_primer_trimming", "temp", "{}.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump(
            [sample_name_out, finished, py_v, cutadapt_v, reads, cut_reads], log
        )


def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the Settings file.
    If default values are desired no arguments are required. Default working directory
    is the current working directory."""

    ## collect variables from the Settings file
    gen_settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="0_general_settings"
    )
    cores = gen_settings["cores to use"].item()

    settings = settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="4_primer_trimming"
    )
    p5_primer, p7_primer, anchoring = (
        settings["P5 Primer (5' - 3')"].item(),
        settings["P7 Primer (5' - 3')"].item(),
        settings["anchoring"].item(),
    )

    # check if any primer sequence is missing
    if pd.isna(settings).values.any():
        print(
            "{}: At least one primer sequence is missing in the settings file.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        sys.exit()
    else:
        for primer_sequence in [p5_primer, p7_primer]:
            if (
                all(char in ambiguous_dna_letters for char in primer_sequence)
                and primer_sequence
            ):
                continue
            else:
                print(
                    "{}: The primer entered contains letters that are not part of the IUPAC nucleotide code.".format(
                        datetime.datetime.now().strftime("%H:%M:%S")
                    )
                )
                sys.exit()

    ## collect the input files from PE merging step
    input = glob.glob(str(Path(project).joinpath("3_PE_merging", "data", "*.fastq.gz")))

    print(
        "{}: Starting to trim primers of {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(input)
        )
    )

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath("4_primer_trimming", "temp"))
    except FileExistsError:
        pass

    Parallel(n_jobs=cores)(
        delayed(primer_trimming)(
            file,
            project=project,
            p5_primer=p5_primer,
            p7_primer=p7_primer,
            anchoring=anchoring,
        )
        for file in input
    )

    ## write logfile from pkl log, remove single logs after
    summary_logs = glob.glob(
        str(Path(project).joinpath("4_primer_trimming", "temp", "*.pkl"))
    )
    summary = [pickle.load(open(line, "rb")) for line in summary_logs]

    ## generate log dataframe for primer trimming
    log_df = pd.DataFrame(
        summary,
        columns=[
            "File",
            "finished at",
            "python version",
            "cutadapt version",
            "processed reads",
            "trimmed reads",
        ],
    )
    log_df = log_df.sort_values(by="File")
    log_df.to_excel(
        Path(project).joinpath("4_primer_trimming", "Logfile_4_primer_trimming.xlsx"),
        index=False,
        sheet_name="4_primer_trimming",
    )

    ## add log to the project report
    with pd.ExcelWriter(
        Path(project).joinpath("Project_report.xlsx"),
        mode="a",
        if_sheet_exists="replace",
        engine="openpyxl",
    ) as writer:
        log_df.to_excel(writer, sheet_name="4_primer_trimming", index=False)

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("4_primer_trimming", "temp"))


if __name__ == "__main__":
    main()

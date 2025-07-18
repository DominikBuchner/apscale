import subprocess, datetime, gzip, os, pickle, glob, shutil, re, hashlib, math, sys, powerlaw, contextlib
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from apscale.a_create_project import empty_file
from collections import Counter


# Custom context manager to suppress stdout
# needed to supress prints from the power law package
@contextlib.contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


# function to count the total abundance of a fasta file
def get_total_abundance(file_path: str) -> tuple:
    """Function to count the total abundance of reads in the file provided via path

    Args:
        file_path (str): Path to the input file.

    Returns:
        tuple: Tuple in the form of (file_path, total abundance).
    """
    # read the input file
    file = SimpleFastaParser(gzip.open(file_path, "rt"))
    # extract the abundance information from the header
    total_abundance = 0
    for header, _ in file:
        size = header.split("size=")[-1]
        total_abundance += int(size)

    return file_path, total_abundance


# function to compute the power law cutoff
def get_power_law_cutoff(file_path: str) -> tuple:
    """Function to compute the min_size parameters via a power law distribution.
    --> Assumption that sequencing data follows a power law relationship: A few dominant signals, many rare signals + noise

    Args:
        file_path (str): Path to the input file.

    Returns:
        tuple: Tuple in the form of (file_path, min_size)
    """
    # read the input file line by line and extract the abundance data
    file = SimpleFastaParser(gzip.open(file_path, "rt"))

    abundances = []

    for header, _ in file:
        size = header.split("size=")[-1]
        abundances.append(size)

    # check if a large proportion > 25% of data points are singletons and doubletons -> only fit power law then, else return "normal" threshold of 4
    seq_counts = Counter(abundances)

    try:
        prop = (seq_counts["1"] + seq_counts["2"]) / len(abundances)
    except ZeroDivisionError:
        return file_path, 4

    abundances = np.array(abundances)

    # fit the power law distribution
    if prop > 0.25:
        # only fit power law if we have enough values
        if len(abundances) > 50 and len(np.unique(abundances)) >= 2:
            with suppress_stdout():
                fit = powerlaw.Fit(abundances)
                # only return power law if it actually models the distribution well
                _, p = fit.distribution_compare("power_law", "exponential")
                if p < 0.05:
                    if int(fit.power_law.xmin) >= 4:
                        return file_path, int(fit.power_law.xmin)
                    else:
                        return file_path, 4
                else:
                    return file_path, 4
        else:
            return file_path, 4


## denoising function to denoise all sequences the input fasta with a given alpha and minsize
def denoise(file, project=None, comp_lvl=None, alpha=None, minsize=None):
    """Function to apply denoising to a given gzipped file. Outputs a fasta file with all
    centroid sequences."""
    ## define the name for the output fasta
    ## create an output path to write to
    sample_name_out_1 = "{}_ESVs_with_chimeras.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )
    output_path_1 = Path(project).joinpath("07_denoising", "data", sample_name_out_1)

    ## run vsearch --cluster_unoise to denoise data to ESV level
    ## use --log because for some reason no info is written to stderr with this command
    ## write stdout to uncompressed output at runtime
    ## only run if input file is not empty
    if not empty_file(file):
        with open(output_path_1.with_suffix(""), "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--cluster_unoise",
                    Path(file),
                    "--unoise_alpha",
                    str(alpha),
                    "--minsize",
                    str(minsize),
                    "--sizein",
                    "--sizeout",
                    "--centroids",
                    "-",
                    "--fasta_width",
                    str(0),
                    "--quiet",
                    "--log",
                    Path(project).joinpath(
                        "07_denoising", "temp", "{}.txt".format(sample_name_out_1)
                    ),
                    "--threads",
                    str(1),
                ],
                stdout=output,
                stderr=subprocess.DEVNULL,
            )

        ## compress the output, remove uncompressed output
        with open(output_path_1.with_suffix(""), "rb") as in_stream, gzip.open(
            output_path_1, "wb", comp_lvl
        ) as out_stream:
            shutil.copyfileobj(in_stream, out_stream)
        os.remove(output_path_1.with_suffix(""))

        ## collect processed and passed reads from the log file
        with open(
            Path(project).joinpath(
                "07_denoising", "temp", "{}.txt".format(sample_name_out_1)
            )
        ) as log_file:
            content = log_file.read()
            # only parse if there is anything in the output
            try:
                seqs = int(re.findall(r"(\d+) seqs, min", content)[0])
            except IndexError:
                seqs = 0
            try:
                discarded = int(re.findall(r"(\d+) sequences discarded.", content)[0])
            except IndexError:
                discarded = 0
            try:
                esvs = int(re.findall(r"Clusters: (\d+) Size min", content)[0])
            except IndexError:
                esvs = 0
            try:
                version = re.findall(r"vsearch ([\w\.]*)", content)[0]
            except IndexError:
                f = subprocess.run(["vsearch", "--version"], capture_output=True)
                version = f.stderr.decode("ascii", errors="ignore")
                version = re.findall(r"vsearch ([\w\.]*)", version)[0]
    else:
        with gzip.open(output_path_1, "wb"):
            seqs, discarded, esvs, version = 0, 0, 0, "empty input"

    # add discarded and used reads for denoising to get a correct number for input sequences
    seqs = int(seqs) + int(discarded)

    print(
        "{}: {}: Denoised {} unique sequences into {} ESVs.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), sample_name_out_1, seqs, esvs
        )
    )

    # exchange the names in the output path
    output_path_2 = Path(project).joinpath(
        "07_denoising",
        "data",
        sample_name_out_1.replace("_with_chimeras", "_no_chimeras"),
    )

    ## run vsearch --uchime_denovo to remove chimeric sequences from the ESVs
    # perform only if the input file is not empty
    if not empty_file(output_path_1):
        with open(output_path_2.with_suffix(""), "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--uchime3_denovo",
                    output_path_1,
                    "--relabel",
                    "ESV_",
                    "--nonchimeras",
                    "-",
                    "-fasta_width",
                    str(0),
                    "--quiet",
                    "--sizein",
                    "--sizeout",
                ],
                stdout=output,
                stderr=subprocess.DEVNULL,
            )

        ## compress the output, remove uncompressed output
        with open(output_path_2.with_suffix(""), "rb") as in_stream, gzip.open(
            output_path_2, "wb", comp_lvl
        ) as out_stream:
            shutil.copyfileobj(in_stream, out_stream)
        os.remove(output_path_2.with_suffix(""))

        # remove the files including chimeras afterwards
        os.remove(Path(project).joinpath("07_denoising", "data", sample_name_out_1))
    else:
        with gzip.open(output_path_2, "wb"):
            os.remove(Path(project).joinpath("07_denoising", "data", sample_name_out_1))

    ## collect processed and passed reads from the output fasta, since it is not reported in the log
    f = list(SimpleFastaParser(gzip.open(output_path_2, "rt")))

    # give user output
    print(
        "{}: {}: {} chimeras removed from {} ESV sequences.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            output_path_2.with_suffix("").with_suffix("").name,
            int(esvs) - len(f),
            esvs,
        )
    )

    # calculate final esv count for logging
    final_esvs = int(esvs) - (int(esvs) - len(f))
    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    # generate the final output path only for logging, needed in the next step
    output_path_3 = Path(
        str(output_path_2).replace(
            "_PE_trimmed_filtered_dereplicated_ESVs_no_chimeras", ""
        )
    )

    # generate output sample name
    log_sample_name = Path(project).joinpath(
        "07_denoising",
        "temp",
        "{}_log.pkl".format(output_path_3.with_suffix("").with_suffix("").name),
    )

    with open(
        log_sample_name,
        "wb",
    ) as log:
        pickle.dump(
            [
                output_path_3.with_suffix("").with_suffix("").name,
                finished,
                version if version else "",
                seqs,
                final_esvs,
            ],
            log,
        )


def calculate_hash_headers(file, project=None, comp_lvl=None):
    """Function to transform sequence headers to sha3 256 hashes to easily compare ESV sequences."""
    ## define the name for the output fasta
    ## create an output path to write to
    sample_name_out_1 = "{}_ESVs_with_chimeras.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )

    # exchange the names in the output path
    output_path_2 = Path(project).joinpath(
        "07_denoising",
        "data",
        sample_name_out_1.replace("_with_chimeras", "_no_chimeras"),
    )

    # open the fasta file to a dict
    fasta_data = SimpleFastaParser(gzip.open(output_path_2, "rt"))

    # define all possible appendices to remove for partial pipeline use
    appendices = [
        "_PE_trimmed_filtered_dereplicated_ESVs_no_chimeras",
        "_trimmed_filtered_dereplicated_ESVs_no_chimeras",
        "_filtered_dereplicated_ESVs_no_chimeras",
        "_dereplicated_ESVs_no_chimeras",
    ]

    # generate the final output path
    for apx in appendices:
        if apx in str(output_path_2):
            output_path_3 = Path(str(output_path_2).replace(apx, ""))
            break

    # write the output
    with gzip.open(output_path_3, "wt+") as out_stream:
        for header, seq in fasta_data:
            seq = seq.upper()
            # parse size annotation from fasta header
            size = header.split(";")[-1]
            header = seq.encode("ascii")
            header = hashlib.sha3_256(header).hexdigest()
            out_stream.write(">{};{}\n{}\n".format(header, size, seq))

    os.remove(output_path_2)

    # give user output
    print(
        "{}: {}: Successfully denoised and hashed.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            output_path_3.with_suffix("").with_suffix("").name,
        )
    )


## main function for the denoising script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Will denoise the individual files, perform chimera removal on them and exchange the fasta headers with sha-256 hashes for easier processing.
    """
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
        sheet_name="07_denoising",
    )
    perform, alpha, threshold_type, threshold = (
        settings["perform denoising"].item(),
        settings["alpha"].item(),
        settings["threshold type"].item(),
        settings["size threshold [absolute nr / %]"].item(),
    )

    # check that the data for poisson or cumulative threshold computation is valid
    if threshold_type == "power law":
        derep_settings = pd.read_excel(
            Path(project).joinpath(
                "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
            ),
            sheet_name="06_dereplication",
        )
        derep_min_size = derep_settings["minimum sequence abundance"].item()
        if derep_min_size > 1 and perform == True:
            print(
                "{}: Power law driven minsize is only available if dereplication min size is 1. Please repeat dereplication with a valid value!".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )

            sys.exit()

    # only perform denoising is set to True:
    if perform:
        ## create temporal output folder
        try:
            os.mkdir(Path(project).joinpath("07_denoising", "temp"))
        except FileExistsError:
            pass

        # remove old analysis results to make sure all data is replaced
        for file in Path(project).joinpath("07_denoising", "data").glob("*"):
            file.unlink()

        ## give user output
        print(
            "{}: Starting individual denoising and chimera removal.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        ## collect input files from quality filtering step
        input = glob.glob(
            str(Path(project).joinpath("06_dereplication", "data", "*.fasta.gz"))
        )

        # if threshold type is set to relative, compute the minsize parameter per file
        if threshold_type == "relative":
            print(
                "{}: Computing total read abundance in input files.".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )

            # compute the total abundance
            input = Parallel(n_jobs=cores)(
                delayed(get_total_abundance)(file) for file in input
            )

            # set the individual thresholds accordingly
            input = [(file, math.ceil(size * threshold * 0.01)) for file, size in input]
            # make all integers positive since minsize is a positive int
            input = [(file, size) if size > 1 else (file, 1) for file, size in input]
        if threshold_type == "absolute":
            input = [(file, threshold) for file in input]
        if threshold_type == "power law":
            print(
                "{}: Fitting power law distribution to the each file and computing signal to noise cutoff.".format(
                    datetime.datetime.now().strftime("%H:%M:%S")
                )
            )
            # compute the total abundance
            input = Parallel(n_jobs=cores)(
                delayed(get_power_law_cutoff)(file) for file in input
            )

        ## parallelize the denoising
        Parallel(n_jobs=cores)(
            delayed(denoise)(
                file, project=project, comp_lvl=comp_lvl, alpha=alpha, minsize=size
            )
            for file, size in input
        )

        ## give user output
        print(
            "{}: Constructing final ESV fasta files.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        # parallelize the hash value calculation
        Parallel(n_jobs=cores)(
            delayed(calculate_hash_headers)(file, project=project, comp_lvl=comp_lvl)
            for file, _ in input
        )

        ## write log for the denoising from pkl logs
        summary_logs = glob.glob(
            str(Path(project).joinpath("07_denoising", "temp", "*_log.pkl"))
        )
        summary = [pickle.load(open(line, "rb")) for line in summary_logs]

        log_df = pd.DataFrame(
            summary,
            columns=[
                "File",
                "finished at",
                "program version",
                "unique sequence input",
                "ESVs",
            ],
        )
        log_df["minsize"] = [pair[1] for pair in input]
        log_df = log_df.sort_values(by="File")
        log_df.to_excel(
            Path(project).joinpath("07_denoising", "Logfile_07_denoising.xlsx"),
            index=False,
            sheet_name="07_denoising",
        )

        ## add log to the project report
        with pd.ExcelWriter(
            Path(project).joinpath(
                "Project_report_{}.xlsx".format(
                    Path(project).name.replace("_apscale", "")
                )
            ),
            mode="a",
            if_sheet_exists="replace",
            engine="openpyxl",
        ) as writer:
            log_df.to_excel(writer, sheet_name="07_denoising", index=False)

        ## remove temporary files
        shutil.rmtree(Path(project).joinpath("07_denoising", "temp"))
    else:
        ## give user output
        print(
            "{}: Denoising is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )


if __name__ == "__main__":
    main()

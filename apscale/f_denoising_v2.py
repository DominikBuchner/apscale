import subprocess, datetime, gzip, os, pickle, glob, openpyxl, shutil, psutil, re, sys, hashlib
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from Bio import SeqIO
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from io import StringIO
from tqdm import tqdm
from openpyxl.utils.dataframe import dataframe_to_rows


## denoising function to denoise all sequences the input fasta with a given alpha and minsize
def denoise(file, project=None, comp_lvl=None, alpha=None, minsize=None):
    """Function to apply denoisind to a given gzipped file. Outputs a fasta file with all
    centroid sequences."""
    ## define the name for the output fasta
    ## create an output path to write to
    sample_name_out_1 = "{}_ESVs_with_chimeras.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )
    output_path_1 = Path(project).joinpath("7_denoising", "data", sample_name_out_1)

    ## run vsearch --cluster_unoise to denoise data to ESV level
    ## use --log because for some reason no info is written to stderr with this command
    ## write stdout to uncompressed output at runtime
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
                    "7_denoising", "temp", "{}.txt".format(sample_name_out_1)
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
            "7_denoising", "temp", "{}.txt".format(sample_name_out_1)
        )
    ) as log_file:
        content = log_file.read()
        try:
            seqs, esvs = (
                re.findall("(\d+) seqs, min ", content)[0],
                re.findall("Clusters: (\d+) Size min", content)[0],
            )
        except IndexError:
            seqs, esvs = 0, 0

    print(
        "{}: {}: Denoised {} unique sequences into {} ESVs.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), sample_name_out_1, seqs, esvs
        )
    )

    # exchange the names in the output path
    output_path_2 = Path(project).joinpath(
        "7_denoising",
        "data",
        sample_name_out_1.replace("_with_chimeras", "_no_chimeras"),
    )

    ## run vsearch --uchime_denovo to remove chimeric sequences from the ESVs
    with open(output_path_2.with_suffix(""), "w") as output:
        f = subprocess.run(
            [
                "vsearch",
                "--uchime_denovo",
                output_path_1,
                "--relabel",
                "ESV_",
                "--nonchimeras",
                "-",
                "-fasta_width",
                str(0),
                "--quiet",
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
    os.remove(Path(project).joinpath("7_denoising", "data", sample_name_out_1))

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


def calculate_hash_headers(file, project=None, comp_lvl=None):
    """Function to transform sequence headers to sha3 256 hashes to easily compare ESV sequences."""
    ## define the name for the output fasta
    ## create an output path to write to
    sample_name_out_1 = "{}_ESVs_with_chimeras.fasta.gz".format(
        Path(file).with_suffix("").with_suffix("").name
    )

    # exchange the names in the output path
    output_path_2 = Path(project).joinpath(
        "7_denoising",
        "data",
        sample_name_out_1.replace("_with_chimeras", "_no_chimeras"),
    )

    # open the fasta file to a dict
    fasta_data = SimpleFastaParser(gzip.open(output_path_2, "rt"))

    # generate the final output path
    output_path_3 = Path(
        str(output_path_2).replace(
            "_PE_trimmed_filtered_dereplicated_ESVs_no_chimeras", ""
        )
    )

    # write the output
    with gzip.open(output_path_3, "wt+") as out_stream:
        for _, seq in fasta_data:
            seq = seq.upper()
            header = seq.encode("ascii")
            header = hashlib.sha3_256(header).hexdigest()
            out_stream.write(">{}\n{}\n".format(header, seq))

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

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath("7_denoising", "temp"))
    except FileExistsError:
        pass

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
        sheet_name="6_denoising",
    )
    alpha, minsize = (
        settings["alpha"].item(),
        settings["minsize"].item(),
    )

    ## collect input files from quality filtering step
    input = glob.glob(
        str(Path(project).joinpath("6_dereplication", "data", "*.fasta.gz"))
    )

    ## give user output
    print(
        "{}: Starting individual denoising and chimera removal.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    ## parallelize the denoising
    Parallel(n_jobs=cores)(
        delayed(denoise)(
            file, project=project, comp_lvl=comp_lvl, alpha=alpha, minsize=minsize
        )
        for file in input
    )

    ## give user output
    print(
        "{}: Constructing final ESV fasta files.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    # parallelize the has value calculation
    Parallel(n_jobs=cores)(
        delayed(calculate_hash_headers)(file, project=project, comp_lvl=comp_lvl)
        for file in input
    )


if __name__ == "__main__":
    main("D:\\denoising_test_apscale")

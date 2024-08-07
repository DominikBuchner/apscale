import subprocess, datetime, gzip, os, pickle, glob, openpyxl, shutil, psutil, re, sys
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
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
    output_path = Path(project).joinpath("7_denoising", "data", sample_name_out_1)

    ## run vsearch --cluster_unoise to denoise data to ESV level
    ## use --log because for some reason no info is written to stderr with this command
    ## write stdout to uncompressed output at runtime
    with open(output_path.with_suffix(""), "w") as output:
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
    with open(output_path.with_suffix(""), "rb") as in_stream, gzip.open(
        output_path, "wb", comp_lvl
    ) as out_stream:
        shutil.copyfileobj(in_stream, out_stream)
    os.remove(output_path.with_suffix(""))

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


if __name__ == "__main__":
    main("D:\\denoising_test_apscale")

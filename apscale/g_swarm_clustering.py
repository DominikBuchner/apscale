import datetime, os, glob, gzip, shutil, subprocess, re, pickle, hashlib
import pandas as pd
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from joblib import Parallel, delayed
from apscale.a_create_project import choose_input, empty_file


def unzip_inputs(file: str, project=None) -> None:
    """_summary_

    Args:
        file (str): Path to the file to unzip
        project (_type_, optional): Project to work in. Defaults to None.
    """
    # define the output name
    unzipped_output = Path(project).joinpath(
        "08_swarm_clustering", "temp", Path(file).stem
    )

    # write the output
    with gzip.open(file, "rt") as in_stream:
        with open(unzipped_output, "wt") as out_stream:
            shutil.copyfileobj(in_stream, out_stream)


def swarm_clustering(file, project=None, comp_lvl=None, prior_step=None):
    """Function to perform swarm clustering on a given input. Outputs a fasta file with all centroid sequences.

    Args:
        file (str): Input file to run swarm on.
        project (str, optional): Project to work in. Defaults to None.
        comp_lvl (int, optional): Compression level for gzip to use. Defaults to None.
        prior_step (str, optional): Which step was performed before. Important for renaming files accordingly. Defaults to None.
    """
    # decide which output name to use
    if prior_step == "07_denoising":
        sample_name_out_1 = "{}.gz".format(Path(file).name)
        output_path_1 = Path(project).joinpath(
            "08_swarm_clustering", "data", sample_name_out_1
        )
    elif prior_step == "06_dereplication":
        sample_name_out_1 = "{}_OTUs_with_chimeras.fasta.gz".format(
            Path(file).with_suffix("").with_suffix("").name
        )
        output_path_1 = Path(project).joinpath(
            "08_swarm_clustering", "temp", sample_name_out_1
        )

    # run swarm if the input is not empty
    if not empty_file(file):
        f = subprocess.run(
            [
                "swarm",
                Path(file),
                "-f",
                "-z",
                "-w",
                output_path_1.with_suffix(""),
                "-l",
                Path(project).joinpath(
                    "08_swarm_clustering",
                    "temp",
                    "{}.txt".format(sample_name_out_1),
                ),
            ],
            stdout=subprocess.DEVNULL,
        )

        # compress the output, remove uncompressed output
        with open(output_path_1.with_suffix(""), "rb") as in_stream:
            with gzip.open(output_path_1, "wb", comp_lvl) as out_stream:
                shutil.copyfileobj(in_stream, out_stream)
        os.remove(output_path_1.with_suffix(""))

        # collect processed and passed reads from the log file
        with open(
            Path(project).joinpath(
                "08_swarm_clustering", "temp", "{}.txt".format(sample_name_out_1)
            )
        ) as log_file:
            content = log_file.read()
            input_seqs = int(re.findall(r"in (\d+) sequences", content)[0])
            swarms = int(re.findall(r"Number of swarms:  (\d+)", content)[0])
            version = re.findall(r"Swarm ([\w.\.]*)", content)[0]
    else:
        with gzip.open(output_path_1, "wb"):
            input_seqs, swarms, version = 0, 0, "empty input"

    # remove swarm output files
    try:
        os.remove(
            Path(project).joinpath(
                "08_swarm_clustering", "temp", "{}.txt".format(sample_name_out_1)
            )
        )
    except FileNotFoundError:
        pass

    # give user output
    print(
        "{}: {}: Clustered {} input amplicons into {} swarm clusters.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            Path(file).name,
            input_seqs,
            swarms,
        )
    )

    # save temporary logs
    temporary_log = Path(project).joinpath(
        "08_swarm_clustering",
        "temp",
        "{}_log.pkl".format(Path(file).with_suffix("").name),
    )

    finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    with open(temporary_log, "wb") as log_file:
        pickle.dump([Path(file).name, finished, version, input_seqs, swarms], log_file)

    # remove the input file, no longer needed
    os.remove(file)


def chimera_removal(file: str, project=None, comp_lvl=None) -> None:
    """Function to remove chimeras from swarm clusters with vsearch.

    Args:
        file (str): Input file to perform chimera filtering.
        comp_lvl (int, optional): Compression level to use. Defaults to None.
        project (str, optional): Apscale project to work in. Defaults to None.
    """
    input_name = Path(file).with_suffix("").name
    output_path = Path(project).joinpath(
        "08_swarm_clustering",
        "temp",
        input_name.replace("_with_chimeras", "_no_chimeras"),
    )
    output_path_gzipped = str(output_path).replace(".fasta", ".fasta.gz")

    ## run vsearch --uchime_denovo to remove chimeric sequences from the ESVs
    # perform only if the input file is not empty
    if not empty_file(file):
        with open(output_path, "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--uchime3_denovo",
                    Path(file),
                    "--relabel",
                    "OTU_",
                    "--nonchimeras",
                    "-",
                    "-fasta_width",
                    str(0),
                    "--quiet",
                    "--sizein",
                    "--sizeout",
                ],
                stdout=output,
            )

        ## compress the output, remove uncompressed output
        with open(output_path, "rb") as in_stream:
            with gzip.open(output_path_gzipped, "wb", comp_lvl) as out_stream:
                shutil.copyfileobj(in_stream, out_stream)

        # remove the files including chimeras afterwards
        os.remove(Path(file))
        os.remove(output_path)
    else:
        with gzip.open(output_path_gzipped, "wb", comp_lvl):
            os.remove(file)

    # read the base statistics from the log files
    log_file_name = output_path_gzipped.replace(
        "_OTUs_no_chimeras.fasta.gz", "_log.pkl"
    )

    with open(log_file_name, "rb") as log_file:
        file, finished_clustering, swarm_version, input_amplicons, swarms = pickle.load(
            log_file
        )

    # write the correct output file name to the log here
    # remove all possible appendices to generate the correct output name
    appendices = [
        "_PE_trimmed_filtered_dereplicated",
        "_trimmed_filtered_dereplicated",
        "_filtered_dereplicated",
        "_dereplicated",
    ]

    for apx in appendices:
        if apx in file:
            file = "{}.gz".format(file.replace(apx, ""))
            break

    # add additional data
    finished_chimera_removal = "{}".format(
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    )
    f = subprocess.run(["vsearch", "--version"], capture_output=True)
    vsearch_version = f.stderr.decode("ascii", errors="ignore")
    vsearch_version = re.findall("vsearch ([\w\.]*)", vsearch_version)[0]

    # collect the remaining sequences from the fasta
    fasta_length = len(list(SimpleFastaParser(gzip.open(output_path_gzipped, "rt"))))

    # give user output
    chimeras_removed = swarms - fasta_length

    # give user output
    print(
        "{}: {}: Removed {} chimeric sequences.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            Path(file).name,
            chimeras_removed,
        )
    )

    # update the log files
    with open(log_file_name, "wb") as log_file:
        pickle.dump(
            [
                file,
                finished_clustering,
                swarm_version,
                input_amplicons,
                swarms,
                finished_chimera_removal,
                vsearch_version,
                chimeras_removed,
            ],
            log_file,
        )


def calculate_hash_headers(file, project=None, comp_lvl=None):
    """Function to transform sequence headers to sha3 256 hashes to easily compare centroid sequences.

    Args:
        file (str): File to hash.
        project (str, optional): Project to work in. Defaults to None.
        comp_lvl (int, optional): Compression level to use. Defaults to None.
    """
    input_file_name = Path(file).with_suffix("").with_suffix("").name

    # remove all possible appendices to generate the correct output name
    appendices = [
        "_PE_trimmed_filtered_dereplicated_OTUs_no_chimeras",
        "_trimmed_filtered_dereplicated_OTUs_no_chimeras",
        "_filtered_dereplicated_OTUs_no_chimeras",
        "_dereplicated_OTUs_no_chimeras",
    ]

    for apx in appendices:
        if apx in input_file_name:
            output_file_name = "{}.fasta.gz".format(input_file_name.replace(apx, ""))
            break

    # generate output path
    output_path = Path(project).joinpath(
        "08_swarm_clustering", "data", output_file_name
    )

    # read in the data from the fasta
    fasta_data = SimpleFastaParser(gzip.open(Path(file), "rt"))

    # write the output
    with gzip.open(output_path, "wt", compresslevel=comp_lvl) as out_stream:
        for header, seq in fasta_data:
            seq = seq.upper()
            # parse size annotation from fasta header
            size = header.split(";")[-1]
            header = seq.encode("ascii")
            header = hashlib.sha3_256(header).hexdigest()
            out_stream.write(">{};{}\n{}\n".format(header, size, seq))

    # remove the temporary input
    os.remove(Path(file))

    # give user output
    print(
        "{}: {}: Successfully clustered and hashed.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            input_file_name,
        )
    )


# main function of the swarm clustering script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Will perform swarm clustering on the individual files. If run without denoising will also perform chimera removal and
    exchange the fasta headers with sha 256 hashes for easier processing.

    Args:
        project (str, optional): Path to the apscale project. Defaults to Path.cwd().
    """
    # collect variables from the settings file
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
        sheet_name="08_swarm_clustering",
    )
    perform = settings["perform swarm clustering"].item()

    # only perform clustering if needed
    if perform:
        # select the correct step to pick the input files
        prior_step = choose_input(project, "08_swarm_clustering")

        # create temporal output folder, override if there is old data
        try:
            os.mkdir(Path(project).joinpath("08_swarm_clustering", "temp"))
        except FileExistsError:
            shutil.rmtree(Path(project).joinpath("08_swarm_clustering", "temp"))
            os.mkdir(Path(project).joinpath("08_swarm_clustering", "temp"))

        # gather input files
        input = glob.glob(str(Path(project).joinpath(prior_step, "data", "*.fasta.gz")))

        # since swarm does not work with gzipped data, unzip files first
        Parallel(n_jobs=cores)(delayed(unzip_inputs)(file, project) for file in input)

        # collect the new input for swarm clustering
        input = glob.glob(
            str(Path(project).joinpath("08_swarm_clustering", "temp", "*.fasta"))
        )
        # decide wether to perform just swarm clustering, or chimera removal has to be performed as well
        # if data has been denoised before, only swarm clustering has to be performed
        if prior_step == "07_denoising":
            # parallelize the clustering
            Parallel(n_jobs=cores)(
                delayed(swarm_clustering)(
                    file, project, comp_lvl=comp_lvl, prior_step=prior_step
                )
                for file in input
            )

            # write the project report, remove temporary files
            summary_logs = glob.glob(
                str(Path(project).joinpath("08_swarm_clustering", "temp", "*_log.pkl"))
            )
            summary = [pickle.load(open(line, "rb")) for line in summary_logs]

            log_df = pd.DataFrame(
                summary,
                columns=[
                    "File",
                    "finished at",
                    "program version",
                    "unique sequence input",
                    "swarms",
                ],
            )
            log_df = log_df.sort_values(by="File")
            log_df.to_excel(
                Path(project).joinpath(
                    "08_swarm_clustering", "Logfile_08_swarm_clustering.xlsx"
                ),
                index=False,
                sheet_name="08_swarm_clustering",
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
                log_df.to_excel(writer, sheet_name="08_swarm_clustering", index=False)

            ## remove temporary files
            shutil.rmtree(Path(project).joinpath("08_swarm_clustering", "temp"))
        # else apscale has to perform swarm clustering, chimera removal and calculate the hash values and name the samples accordingly
        else:
            # parallelize the clustering
            Parallel(n_jobs=cores)(
                delayed(swarm_clustering)(
                    file, project, comp_lvl=comp_lvl, prior_step=prior_step
                )
                for file in input
            )
            # gather files for chimera removal
            input_files = glob.glob(
                str(Path(project).joinpath("08_swarm_clustering", "temp", "*.fasta.gz"))
            )

            # perform chimera removal parallelized
            Parallel(n_jobs=cores)(
                delayed(chimera_removal)(file, project=project, comp_lvl=comp_lvl)
                for file in input_files
            )

            # gather the input files after chimera removal
            input_files = glob.glob(
                str(Path(project).joinpath("08_swarm_clustering", "temp", "*.fasta.gz"))
            )

            # calculate the hash headers
            Parallel(n_jobs=cores)(
                delayed(calculate_hash_headers)(
                    file, project=project, comp_lvl=comp_lvl
                )
                for file in input_files
            )

            # write the project report, remove temporary files
            summary_logs = glob.glob(
                str(Path(project).joinpath("08_swarm_clustering", "temp", "*_log.pkl"))
            )
            summary = [pickle.load(open(line, "rb")) for line in summary_logs]

            log_df = pd.DataFrame(
                summary,
                columns=[
                    "File",
                    "finished at",
                    "swarm version",
                    "unique sequence input",
                    "swarms",
                    "finished chimera removal at",
                    "vsearch version",
                    "chimeras removed",
                ],
            )
            log_df = log_df.sort_values(by="File")
            log_df.to_excel(
                Path(project).joinpath(
                    "08_swarm_clustering", "Logfile_08_swarm_clustering.xlsx"
                ),
                index=False,
                sheet_name="08_swarm_clustering",
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
                log_df.to_excel(writer, sheet_name="08_swarm_clustering", index=False)

            ## remove temporary files
            shutil.rmtree(Path(project).joinpath("08_swarm_clustering", "temp"))

    else:
        ## give user output
        print(
            "{}: Swarm clustering is disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

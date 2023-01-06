import os, subprocess, gzip, shutil, datetime, contextlib, joblib, openpyxl, re
import pandas as pd
import numpy as np
from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser
from openpyxl.utils.dataframe import dataframe_to_rows

## function to create a progressmeter for the parallel filtering execution
@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    def tqdm_print_progress(self):
        if self.n_completed_tasks > tqdm_object.n:
            n_completed = self.n_completed_tasks - tqdm_object.n
            tqdm_object.update(n=n_completed)

    original_print_progress = joblib.parallel.Parallel.print_progress
    joblib.parallel.Parallel.print_progress = tqdm_print_progress

    try:
        yield tqdm_object
    finally:
        joblib.parallel.Parallel.print_progress = original_print_progress
        tqdm_object.close()


## function to generate the matchfile
def generate_matchfile(
    project=None, file=None, cores=None, comp_lvl=None, type=None, min_sim=None
):
    """Function to generate a matchfile for the given fasta file. The matchfile will be
    gzipped to save space."""

    ## generate output sample name and filepath
    if type == "otu":
        sample_name_out = "OTU_matchfile.txt.gz"
        output_path = Path(project).joinpath(
            "9_lulu_filtering", "otu_clustering", "data", sample_name_out
        )
        log_path = Path(project).joinpath(
            "9_lulu_filtering", "otu_clustering", "temp", "matchfile.log"
        )
    elif type == "esv":
        sample_name_out = "ESV_matchfile.txt.gz"
        output_path = Path(project).joinpath(
            "9_lulu_filtering", "denoising", "data", sample_name_out
        )
        log_path = Path(project).joinpath(
            "9_lulu_filtering", "denoising", "temp", "matchfile.log"
        )

    ## generate the files with vsearch
    with open(output_path.with_suffix(""), "w") as output:
        f = subprocess.run(
            [
                "vsearch",
                "--usearch_global",
                file,
                "--db",
                file,
                "--self",
                "--id",
                str(min_sim / 100),
                "--iddef",
                str(1),
                "--userout",
                "-",
                "--userfields",
                "query+target+id",
                "--maxaccepts",
                str(0),
                "--query_cov",
                str(0.9),
                "--quiet",
                "--log",
                log_path,
                "--threads",
                str(cores),
                "--maxhits",
                str(10),
            ],
            stdout=output,
        )

    ## convert to dataframe, save to parquet
    with open(output_path.with_suffix(""), "rb") as in_stream, gzip.open(
        output_path, "wb", comp_lvl
    ) as out_stream:
        shutil.copyfileobj(in_stream, out_stream)
    os.remove(output_path.with_suffix(""))


## function to load the otu / esv table and reorder it accoring to the lulu algorithm
def load_data_table(project=None, type=None):
    ## read the table from parquet, set the index to OTU names for easy searching / filtering
    if type == "otu":
        data_table = pd.read_parquet(
            Path(project).joinpath(
                "7_otu_clustering",
                "{}_OTU_table.parquet.snappy".format(Path(project).stem),
            )
        )
    elif type == "esv":
        data_table = pd.read_parquet(
            Path(project).joinpath(
                "8_denoising", "{}_ESV_table.parquet.snappy".format(Path(project).stem)
            )
        )
    data_table = data_table.set_index("ID")
    ## order by occurence, the read count, calculate statistics for calculations later
    ## store the sequence data to readd it to the otu table in the end
    seq_dict = dict(zip(data_table.index, data_table.pop("Seq")))
    data_table["occurence"] = np.count_nonzero(data_table, axis=1)
    data_table["read_count"] = data_table.sum(axis=1, numeric_only=True)
    data_table = data_table.sort_values(by=["occurence", "read_count"], ascending=False)
    occ_statistics = dict(zip(data_table.index, data_table["occurence"]))
    data_table = data_table.drop(["occurence", "read_count"], axis=1)

    return data_table, occ_statistics, seq_dict


def load_matches(project=None, type=None, occ_statistics=None):
    if type == "otu":
        sample_name_out = "OTU_matchfile.txt.gz"
        matches = pd.read_csv(
            Path(project).joinpath(
                "9_lulu_filtering", "otu_clustering", "data", sample_name_out
            ),
            sep="\t",
            names=["query", "match", "similarity"],
        )
    elif type == "esv":
        sample_name_out = "ESV_matchfile.txt.gz"
        matches = pd.read_csv(
            Path(project).joinpath(
                "9_lulu_filtering", "denoising", "data", sample_name_out
            ),
            sep="\t",
            names=["query", "match", "similarity"],
        )
    ## generate a dict of all otus in the table and the potential matches as a list
    matches = matches.set_index("query")

    ##  sort the matches in the order of occurence (or read table)
    sorter = [key for key in occ_statistics.keys() if key in matches.index]
    matches = matches.loc[sorter]

    ## generate the dict with potential hits for every id
    matches = matches.groupby("query", sort=False)["match"].apply(list).to_dict()

    return matches


## function to check the potential daughter sequences according to the lulu algorithm
def check_daughter(
    pot_daughter=None,
    hit_frame=None,
    sub_statistics=None,
    min_rel_co=None,
    min_ratio=None,
):
    ## extract the name of the potential daughter
    pot_daughter_name = pot_daughter.index.item()
    ## extract all potential parents
    pot_parents = [
        hit
        for hit in hit_frame.index
        if sub_statistics[hit] >= sub_statistics[pot_daughter_name]
    ]

    ## loop through all potential parents and check if they are parent to the daughter
    for parent in pot_parents:
        if parent:
            ## generate a dataframe for the calculations, transpose, only select samples where the daughter has reads
            comp = pd.concat([hit_frame.loc[[parent]], pot_daughter]).T
            comp = comp.loc[comp[pot_daughter_name] > 0]
            ## compute relative co-occurence
            try:
                rel_co = len(comp) / np.count_nonzero(comp[parent])
            except ZeroDivisionError:
                rel_co = 0
            ## if relative occurence is sufficient compute abundance ratio
            if rel_co > min_rel_co:
                comp["ab_ratio"] = comp[parent] / comp[pot_daughter_name]
                ## set the potential parent as valid if abundance ratio is sufficient
                if min(comp["ab_ratio"]) > min_ratio:
                    return pot_daughter_name, parent, rel_co, min(comp["ab_ratio"])
                    break
                else:
                    continue
            else:
                continue


## accepts the relations returned by the check daughter function in the form of a dict with
## daughter: parent mappings and return a dict with parent: daughter mappings, where parents
## of parents are already aggregated so that the curation of the OTU table can be done with
## a single groupby instead of iterative sums
def aggregate_relations(relations):
    ## collect all daughters of a given parent in this dict
    all_parents = {}

    ## aggregation logic
    for daughter in reversed(list(relations.keys())):
        parent = relations[daughter]

        ## if the daughter is already a parent
        if daughter in all_parents:
            ## if the parent of that daughter is already a parent, append this parent with the daughter and all daughters of the daughter and
            if parent in all_parents:
                all_parents[parent] += [daughter] + all_parents[daughter]
                ## remove the daughter from the parent, since it listed under a new parent
                del all_parents[daughter]
            else:
                ## if the parent of this daughter is not a parent make it one and append the daughter and all daughters of that daughter
                all_parents[parent] = [daughter] + all_parents[daughter]
        elif daughter not in all_parents:
            ## if the parent is already in the parents simply append the new daughter
            if parent in all_parents:
                all_parents[parent] += [daughter]
            ## if not create a new parent: daughter mapping
            elif parent not in all_parents:
                all_parents[parent] = [daughter]

    ## generate the output as a mapping of all daughter: parent pairs to map to the original dataframe
    output = {}

    for parent in all_parents:
        for daughter in all_parents[parent]:
            output[daughter] = parent

    return output


## function to filter the corrected OTUs / ESVs from the fasta file
def filter_fasta(project=None, type=None, seqs_to_keep=None):
    if type == "otu":
        in_stream = Path(project).joinpath(
            "7_otu_clustering", "{}_OTUs.fasta".format(Path(project).stem)
        )
        out_stream = Path(project).joinpath(
            "9_lulu_filtering",
            "otu_clustering",
            "{}_OTUs_filtered.fasta".format(Path(project).stem),
        )
    elif type == "esv":
        in_stream = Path(project).joinpath(
            "8_denoising", "{}_ESVs.fasta".format(Path(project).stem)
        )
        out_stream = Path(project).joinpath(
            "9_lulu_filtering",
            "denoising",
            "{}_ESVs_filtered.fasta".format(Path(project).stem),
        )

    with open(in_stream, "r") as in_stream:
        with open(out_stream, "w") as out_stream:
            for (header, seq) in SimpleFastaParser(in_stream):
                if header in seqs_to_keep:
                    out_stream.write(">{}\n{}\n".format(header, seq))


## function to write the specific log
def write_log(project=None, results=None, type=None):
    log_df = pd.DataFrame(
        results, columns=["ID", "Error of", "cooccurence", "abundance ratio"]
    )
    log_df["sort"] = log_df["ID"].str.split("_").str[1].astype(float)
    log_df = log_df.sort_values(by=["sort"])
    log_df = log_df.drop(["sort"], axis=1)
    f = subprocess.run(["vsearch", "--version"], capture_output=True)
    version = f.stderr.decode("ascii", errors="ignore")
    version = re.findall("vsearch ([\w\.]*)", version)[0]
    log_df["program version"] = version

    ## add log to the project report
    with pd.ExcelWriter(
        Path(project).joinpath("Project_report.xlsx"),
        mode="a",
        if_sheet_exists="replace",
        engine="openpyxl",
    ) as writer:
        ## write the output
        if type == "otu":
            log_df.to_excel(
                Path(project).joinpath(
                    "9_lulu_filtering",
                    "otu_clustering",
                    "Logfile_9_lulu_filtering_otu.xlsx",
                ),
                sheet_name="9_lulu_filtering_otu",
                index=False,
            )
            log_df.to_excel(writer, sheet_name="9_lulu_filtering_otu", index=False)
        elif type == "esv":
            log_df.to_excel(
                Path(project).joinpath(
                    "9_lulu_filtering", "denoising", "Logfile_9_lulu_filtering_esv.xlsx"
                ),
                sheet_name="9_lulu_filtering_esv",
                index=False,
            )
            log_df.to_excel(writer, sheet_name="9_lulu_filtering_esv", index=False)


def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the Settings file.
    If default values are desired no arguments are required. Default working directory
    is the current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="0_general_settings"
    )
    cores, comp_lvl = (
        gen_settings["cores to use"].item(),
        gen_settings["compression level"].item(),
    )

    ## extract the settings for this module
    settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="9_lulu_filtering"
    )
    min_sim, min_rel_co, min_ratio, to_excel = (
        settings["minimum similarity"].item(),
        settings["minimum relative cooccurence"].item(),
        settings["minimum ratio"].item(),
        settings["to excel"].item(),
    )

    ## create temporal output folder for the log files
    try:
        os.mkdir(Path(project).joinpath("9_lulu_filtering", "otu_clustering", "temp"))
        os.mkdir(Path(project).joinpath("9_lulu_filtering", "denoising", "temp"))
    except FileExistsError:
        pass

    ## user output
    print(
        "{}: Generating matchfiles for OTUs and ESVs.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    ## start with the matchfile generation of otus and esvs
    otu_fasta = Path(project).joinpath(
        "7_otu_clustering", "{}_OTUs.fasta".format(Path(project).stem)
    )
    esv_fasta = Path(project).joinpath(
        "8_denoising", "{}_ESVs.fasta".format(Path(project).stem)
    )

    ## generate the matchlists for otus and esvs
    for file, type in zip((otu_fasta, esv_fasta), ("otu", "esv")):
        generate_matchfile(project, file, cores, comp_lvl, type, min_sim)

    for type in ("otu", "esv"):
        print(
            "{}: Loading the {} table.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
            )
        )
        data_table, occ_statistics, seq_dict = load_data_table(project, type)

        ## load the match list
        print(
            "{}: Loading the match list table.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        matches = load_matches(project, type, occ_statistics)

        print(
            "{}: Starting to filter the {}s.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
            )
        )
        ## run the filtering algorithm parallelized
        with tqdm_joblib(tqdm(desc="Filtering: ", total=len(matches))) as progress_bar:
            results = Parallel(n_jobs=cores)(
                delayed(check_daughter)(
                    data_table.loc[[seq]],
                    data_table.loc[data_table.index.isin(matches[seq])],
                    {seq: occ_statistics[seq] for seq in matches[seq] + [seq]},
                    min_rel_co=min_rel_co / 100,
                    min_ratio=min_ratio,
                )
                for seq in matches
            )

        ## remove none return from parallel execution, check relations between OTUs / ESVs
        results = [hit for hit in results if hit]
        print(
            "{}: Found {} erroneous {}s.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), len(results), type.upper()
            )
        )
        relations = {hit[0]: hit[1] for hit in results}

        ## find parents of parents in the relations, aggregate the result to individual
        ## dauther: parent mappings
        relations = aggregate_relations(relations)

        print(
            "{}: Correcting the {} table.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
            )
        )
        new_idx = [
            relations[otu] if otu in relations else otu for otu in data_table.index
        ]
        data_table.index = new_idx
        data_table.index.name = "ID"
        data_table = data_table.groupby(data_table.index, sort=False).sum()

        ## resort OTU column
        data_table["sort"] = data_table.index.str.split("_").str[1].astype(float)
        data_table = data_table.sort_values(by=["sort"])
        data_table = data_table.drop(["sort"], axis=1)

        ## set index to column
        data_table = data_table.reset_index()
        data_table["Seq"] = data_table["ID"].map(seq_dict)

        ## write table to parquet
        print(
            "{}: Saving the filtered {} table to parquet.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
            )
        )
        if type == "otu":
            outpath = Path(project).joinpath(
                "9_lulu_filtering",
                "otu_clustering",
                "{}_OTU_table_filtered.parquet.snappy".format(Path(project).stem),
            )
        elif type == "esv":
            outpath = Path(project).joinpath(
                "9_lulu_filtering",
                "denoising",
                "{}_ESV_table_filtered.parquet.snappy".format(Path(project).stem),
            )
        data_table.to_parquet(outpath)
        print(
            "{}: {} table saved to {}.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper(), outpath
            )
        )

        ## write table to excel if option is enabled
        if to_excel:
            wb = openpyxl.Workbook(write_only=True)
            ws = wb.create_sheet("{} table".format(type.upper()))

            ## save the output line by line for optimized memory usage
            for row in tqdm(
                dataframe_to_rows(data_table, index=False, header=True),
                total=len(data_table.index),
                desc="{}: Lines written to {} table".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
                ),
                unit=" lines",
            ):
                ws.append(row)

            ## save the output (otu table)
            print(
                "{}: Saving the {} table to excel. This may take a while.".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
                )
            )
            if type == "otu":
                wb.save(
                    Path(project).joinpath(
                        "9_lulu_filtering",
                        "otu_clustering",
                        "{}_{}_table_filtered.xlsx".format(
                            Path(project).stem, type.upper()
                        ),
                    )
                )
            elif type == "esv":
                wb.save(
                    Path(project).joinpath(
                        "9_lulu_filtering",
                        "denoising",
                        "{}_{}_table_filtered.xlsx".format(
                            Path(project).stem, type.upper()
                        ),
                    )
                )
            wb.close()
            print(
                "{}: ESV table saved to {}.".format(
                    datetime.datetime.now().strftime("%H:%M:%S"),
                    Path(project).joinpath(
                        "8_denoising", "{}_ESV_table.xlsx".format(Path(project).stem)
                    ),
                )
            )

        ## rewriting the fasta file
        print(
            "{}: Removing erroneous {}s from fasta file.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
            )
        )
        filter_fasta(project, type, data_table["ID"].to_list())

        ## write the log
        print(
            "{}: Writing log files.".format(
                datetime.datetime.now().strftime("%H:%M:%S"), type.upper()
            )
        )
        write_log(project, results, type=type)

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("9_lulu_filtering", "otu_clustering", "temp"))
    shutil.rmtree(Path(project).joinpath("9_lulu_filtering", "denoising", "temp"))


if __name__ == "__main__":
    main()

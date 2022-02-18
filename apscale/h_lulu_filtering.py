import os, subprocess, gzip, shutil
import pandas as pd
from pathlib import Path

## function to generate the matchfile
def generate_matchfile(project = None, file = None, cores = None, comp_lvl = None, type = None, min_sim = None):
    """Function to generate a matchfile for the given fasta file. The matchfile will be
    gzipped to save space."""

    ## generate output sample name and filepath
    if type == 'otu':
        sample_name_out = 'OTU_matchfile.txt.gz'
        output_path = Path(project).joinpath('9_lulu_filtering', 'otu_clustering', 'data', sample_name_out)
        log_path =Path(project).joinpath('9_lulu_filtering', 'otu_clustering', 'temp', 'matchfile.log')
    elif type == 'esv':
        sample_name_out = 'ESV_matchfile.txt.gz'
        output_path = Path(project).joinpath('9_lulu_filtering', 'denoising', 'data', sample_name_out)
        log_path =Path(project).joinpath('9_lulu_filtering', 'denoising', 'temp', 'matchfile.log')

    ## generate the files with vsearch
    with open(output_path.with_suffix(''), 'w') as output:
        f = subprocess.run(['vsearch',
                            '--usearch_global', file,
                            '--db', file,
                            '--self',
                            '--id', str(min_sim / 100),
                            '--iddef', str(1),
                            '--userout', '-',
                            '--userfields', 'query+target+id',
                            '--maxaccepts', str(0),
                            '--query_cov', str(0.9),
                            '--quiet',
                            '--log', log_path,
                            '--threads', str(cores)], stdout = output)

    ## compress the output
    with open(output_path.with_suffix(''), 'rb') as in_stream, gzip.open(output_path, 'wb', comp_lvl) as out_stream:
            shutil.copyfileobj(in_stream, out_stream)
    os.remove(output_path.with_suffix(''))

def main(project = Path.cwd()):
    """Main function of the script. Default values can be changed via the Settings file.
    If default values are desired no arguments are required. Default working directory
    is the current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '0_general_settings')
    cores, comp_lvl = gen_settings['cores to use'].item(), gen_settings['compression level'].item()

    ## extract the settings for this module
    settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '9_lulu_filtering')
    min_sim, min_rel_co, min_ratio = settings['minimum similarity'].item(), settings['minimum relative cooccurence'].item(), settings['minimum ratio'].item()

    ## create temporal output folder for the log files
    try:
        os.mkdir(Path(project).joinpath('9_lulu_filtering', 'otu_clustering', 'temp'))
        os.mkdir(Path(project).joinpath('9_lulu_filtering', 'denoising', 'temp'))
    except FileExistsError:
        pass

    ## start with the matchfile generation of otus and esvs
    otu_fasta = Path(project).joinpath('7_otu_clustering', '{}_OTUs.fasta'.format(Path(project).stem))
    esv_fasta = Path(project).joinpath('8_denoising', '{}_ESVs.fasta'.format(Path(project).stem))

    for file, type in zip((otu_fasta, esv_fasta), ('otu', 'esv')):
        generate_matchfile(project, file, cores, comp_lvl, type, min_sim)
    ## generate the matchlist for the otus


if __name__ == "__main__":
    main('C:\\Users\\Dominik\\Desktop\\lulu_apscale')

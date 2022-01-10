import subprocess, gzip, datetime, glob, os, pickle, openpyxl, shutil
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed

## dereplication function to dereplicate a gzipped fasta file
def dereplication(file, project = None, comp_lvl = None):
    """Function to dereplicate a gzipped fasta file. Abundance annotations will be
    written to the output fasta."""

    ## extract the filename from the sample path / name and convert to output name
    sample_name_out = '{}_dereplicated.fasta.gz'.format(Path(file).with_suffix('').with_suffix('').name)

    ## run vsearch --derep_fulllength to dereplicate the file
    ## use --log because for some reason no info is written to stderr with this command
    f = subprocess.run(['vsearch',
                        '--derep_fulllength', Path(file),
                        '--output', '-', '--quiet', '--fasta_width', str(0),
                        '--log', Path(project).joinpath('6_dereplication_pooling', 'temp', '{}.txt'.format(sample_name_out)),
                        '--sizeout'], capture_output = True)

    ## write gzipped output so save space
    with gzip.open(Path(project).joinpath('6_dereplication_pooling', 'data', 'dereplication', sample_name_out), 'wb', comp_lvl) as out:
        out.write(f.stdout)

    ## collect processed and passed reads from the log file
    with open(Path(project).joinpath('6_dereplication_pooling', 'temp', '{}.txt'.format(sample_name_out))) as log_file:
        content = log_file.read().split('\n')
        seqs, unique_seqs = content[3].split(' ')[3], content[4].split(' ')[0]
        version = content[0].split(',')[0]
        finished = '{}'.format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## give user output
    print('{}: {}: {} sequences dereplicated into {} unique sequences.'.format(datetime.datetime.now().strftime("%H:%M:%S"), sample_name_out, seqs, unique_seqs))

    ## temporarily pickle output for the log file
    with open(Path(project).joinpath('6_dereplication_pooling', 'temp', '{}.pkl'.format(sample_name_out)), 'wb') as log:
        pickle.dump([sample_name_out, finished, version, seqs, unique_seqs], log)

## function to pool the dereplicated reads, pooled reads are dereplicated again
def pooling(file_list, project = None, comp_lvl = None):
    """Function to pool the dereplicated reads, needs a path to dereplicated reads folder
    and a project to work in. Both are passed by the main function."""

    ## write the output file
    with gzip.open(Path(project).joinpath('6_dereplication_pooling', 'data', 'pooling', 'pooled_sequences.fasta.gz'), 'wt', comp_lvl) as pool:
        number_of_files = len(file_list)
        for i in range(len(file_list)):
            with gzip.open(file_list[i], 'rt') as to_pool:
                shutil.copyfileobj(to_pool, pool)
            print('{}: Added {} of {} files to the pooled sequences.'.format(datetime.datetime.now().strftime("%H:%M:%S"), i + 1, number_of_files))

    ## dereplicate the pool with minuniquesize = 2
    print('{}: Dereplicating the pooled sequences for clustering and denoising.'.format(datetime.datetime.now().strftime("%H:%M:%S")))

    ## run vsearch --derep_fulllength to dereplicate the file
    f = subprocess.run(['vsearch',
                        '--derep_fulllength', Path(project).joinpath('6_dereplication_pooling', 'data', 'pooling', 'pooled_sequences.fasta.gz'),
                        '--output', '-', '--quiet', '--fasta_width', str(0),
                        '--sizein', '--sizeout',
                        '--minuniquesize', str(2)], capture_output = True)

    ## write gzipped output so save space
    with gzip.open(Path(project).joinpath('6_dereplication_pooling', 'data', 'pooling', 'pooled_sequences_dereplicated.fasta.gz'), 'wb', comp_lvl) as out:
        out.write(f.stdout)

## main function to call the script
def main(project = Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Calls the dereplication over all files, merges them into one file and dereplicates
    the merged file again. Project defaults so current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '0_general_settings')
    cores, comp_lvl = gen_settings['cores to use'].item(), gen_settings['compression level'].item()

    ## collect input files from quality filtering step
    input = glob.glob(str(Path(project).joinpath('5_quality_filtering', 'data', '*.fasta.gz')))

    print('{}: Starting to dereplicate {} input files.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(input)))

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath('6_dereplication_pooling', 'temp'))
    except FileExistsError:
        pass

    ## parallelize the dereplication
    Parallel(n_jobs = cores)(delayed(dereplication)(file, project = project, comp_lvl = comp_lvl) for file in input)

    ## write log for the dereplication from pkl logs
    summary_logs = glob.glob(str(Path(project).joinpath('6_dereplication_pooling', 'temp', '*.pkl')))
    summary = [pickle.load(open(line, 'rb')) for line in summary_logs]

    log_df = pd.DataFrame(summary, columns = ['File', 'finished at', 'program version', 'processed sequences', 'unique sequences'])
    log_df = log_df.sort_values(by = 'File')
    log_df.to_excel(Path(project).joinpath('6_dereplication_pooling', 'Logfile_6_dereplication.xlsx'),
                    index = False,
                    sheet_name = '6_dereplication')

    ## add log to the project report
    wb = openpyxl.load_workbook(Path(project).joinpath('Project_report.xlsx'))
    writer = pd.ExcelWriter(Path(project).joinpath('Project_report.xlsx'), engine = 'openpyxl')
    writer.book = wb

    ## write the output
    log_df.to_excel(writer, sheet_name = '6_dereplication', index = False)
    wb.save(Path(project).joinpath('Project_report.xlsx'))
    writer.close()

    ## pool the dereplicated files
    files = glob.glob(str(Path(project).joinpath('6_dereplication_pooling' ,'data', 'dereplication', '*.fasta.gz')))
    pooling(files, project = project, comp_lvl = comp_lvl)

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath('6_dereplication_pooling', 'temp'))

if __name__ == "__main__":
    main()

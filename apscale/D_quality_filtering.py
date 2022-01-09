import subprocess, gzip, datetime, pickle, glob, os, openpyxl, shutil
import pandas as pd
from pathlib import Path
from joblib import Parallel, delayed

## quality filtering function to quality filter the specified file
def quality_filtering(file, project = None, comp_lvl = None, maxee = None, min_length = None, max_length = None):
    """Function to apply quality filtering to a gzipped fastq file. Outputs a gzipped fasta file
    with quality filtered reads. Filtered reads will be discarded."""

    ## extract the filename from the sample path / name and convert to output name
    sample_name_out = '{}_filtered.fasta.gz'.format(Path(file).with_suffix('').with_suffix('').name)

    # run vsearch --fastq_filter to apply the quality filtering
    # use --log because for some reason no info is written to stderr with this command
    f = subprocess.run(['vsearch',
                        '--fastq_filter', Path(file),
                        '--fastaout', '-', '--quiet', '--fasta_width', str(0),
                        '--log', Path(project).joinpath('5_quality_filtering', 'temp', '{}.txt'.format(sample_name_out)),
                        '--fastq_maxee', str(maxee),
                        '--fastq_minlen', str(min_length),
                        '--fastq_maxlen', str(max_length),
                        '--fastq_qmax', str(64)], capture_output = True)

    ## write gzipped output so save space
    with gzip.open(Path(project).joinpath('5_quality_filtering', 'data', sample_name_out), 'wb', comp_lvl) as out:
        out.write(f.stdout)

    ## collect processed and passed reads from the log file
    with open(Path(project).joinpath('5_quality_filtering', 'temp', '{}.txt'.format(sample_name_out))) as log_file:
        content = log_file.read().split('\n')
        kept, discarded = int(content[3].split(' ')[0]), int(content[3].split(' ')[7])
        reads = kept + discarded
        version = content[0].split(',')[0]
        finished = '{}'.format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## Give user output, if 0 reads are the output handle Zero division exception
    try:
        print('{}: {}: {} of {} reads passed quality filtering ({:.2f}%)'.format(datetime.datetime.now().strftime("%H:%M:%S"), sample_name_out, kept, reads, (kept) / reads * 100))
    except ZeroDivisionError:
        print('{}: {}: {} of {} reads passed quality filtering ({:.2f}%)'.format(datetime.datetime.now().strftime("%H:%M:%S"),sample_name_out, 0, reads, 0))

    ## temporarily pickle output for the log file, get vsearch version
    with open(Path(project).joinpath('5_quality_filtering', 'temp', '{}.pkl'.format(sample_name_out)), 'wb') as log:
        pickle.dump([sample_name_out, finished, version, reads, kept], log)

## main function to call the script
def main(project = Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    If default values are desired no arguments are required. Default working directory
    is the current working directory."""

    ## collect variables from the settings file
    gen_settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '0_general_settings')
    cores, comp_lvl = gen_settings['cores to use'].item(), gen_settings['compression level'].item()

    settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '5_quality_filtering')
    maxee, min_length, max_length = settings['maxEE'].item(), settings['min length'].item(), settings['max length'].item()

    ## collect the input files from primer trimming step
    input = glob.glob(str(Path(project).joinpath('4_primer_trimming', 'data', '*.fastq.gz')))

    print('{}: Starting to quality filter {} input files.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(input)))

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath('5_quality_filtering', 'temp'))
    except FileExistsError:
        pass

    ## parallelize the quality filtering
    Parallel(n_jobs = cores)(delayed(quality_filtering)(file, project = project, comp_lvl = comp_lvl, maxee = maxee, min_length = min_length, max_length = max_length) for file in input)

    ## write logfile from pkl log, remove single logs after
    summary_logs = glob.glob(str(Path(project).joinpath('5_quality_filtering', 'temp', '*.pkl')))
    summary = [pickle.load(open(line, 'rb')) for line in summary_logs]

    log_df = pd.DataFrame(summary, columns = ['File', 'finished at', 'program version', 'processed reads', 'passed reads'])
    log_df = log_df.sort_values(by = 'File')
    log_df.to_excel(Path(project).joinpath('5_quality_filtering', 'Logfile_5_quality_filtering.xlsx'),
                    index = False,
                    sheet_name = '5_quality_filtering')

    ## add log to the project report
    wb = openpyxl.load_workbook(Path(project).joinpath('Project_report.xlsx'))
    writer = pd.ExcelWriter(Path(project).joinpath('Project_report.xlsx'))
    writer.book = wb

    ## write the output
    log_df.to_excel(writer, sheet_name = '5_quality_filtering', index = False)
    wb.save(Path(project).joinpath('Project_report.xlsx'))
    writer.close()

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath('5_quality_filtering', 'temp'))

if __name__ == "__main__":
    main()

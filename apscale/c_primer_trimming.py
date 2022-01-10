import subprocess, datetime, pickle, glob, os, openpyxl, shutil
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from io import StringIO
from joblib import Parallel, delayed

## function to trim primers of reads of the specified file
def primer_trimming(file, project = None, p5_primer = None, p7_primer = None, anchoring = None):
    """Function to remove primer sequences from gzipped file via cutadapt. Outputs
    gzipped file with all reads where primers we're removed. Untrimmed reads will be
    discarded."""

    ## extract the filename from the sample path / name and convert to output name
    sample_name_out = '{}_trimmed.fastq.gz'.format(Path(file).with_suffix('').with_suffix('').name)

    ## if anchoring is True change the cutadapt call
    if anchoring:
        adapter = '^{}...{}'.format(p5_primer, str(Seq(p7_primer).reverse_complement()))
    else:
        adapter = '{}...{}'.format(p5_primer, str(Seq(p7_primer).reverse_complement()))

    ## run catadapt
    f = subprocess.run(['cutadapt',
                        '-a', adapter,
                        '-o', str(Path(project).joinpath('4_primer_trimming', 'data', sample_name_out)), file,
                        '--discard-untrimmed', '--cores=1', '--report=minimal'], capture_output = True)

    ## collect processed reads from stderror for the logfile,
    ## handle exception for empty outputs
    log_df = pd.read_csv(StringIO(f.stdout.decode('ascii', errors = 'ignore')), sep = '\t')
    reads, cut_reads = log_df['in_reads'].item(), log_df['out_reads'].item()
    finished = '{}'.format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## give user output
    try:
        print('{}: {}: primers trimmed for {} of {} reads ({:.2f}%)'.format(datetime.datetime.now().strftime("%H:%M:%S"), sample_name_out, cut_reads, reads, cut_reads / reads * 100))
    except ZeroDivisionError:
        print('{}: {}: primers trimmed for {} of {} reads ({:.2f}%)'.format(datetime.datetime.now().strftime("%H:%M:%S"), sample_name_out, 0, reads, 0))

    ## get remainind log information, pickle temporarly to write the log after successfull finish
    py_v = subprocess.run(['python', '--version'], capture_output = True)
    py_v = py_v.stdout.decode('ascii', errors = 'ignore').rstrip()

    cutadapt_v = subprocess.run(['cutadapt', '--version'], capture_output = True)
    cutadapt_v = cutadapt_v.stdout.decode('ascii', errors = 'ignore').rstrip()

    with open(Path(project).joinpath('4_primer_trimming', 'temp', '{}.pkl'.format(sample_name_out)), 'wb') as log:
        pickle.dump([sample_name_out, finished, py_v, cutadapt_v, reads, cut_reads], log)

def main(project = Path.cwd()):
    """Main function of the script. Default values can be changed via the Settings file.
    If default values are desired no arguments are required. Default working directory
    is the current working directory."""

    ## collect variables from the Settings file
    gen_settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '0_general_settings')
    cores, comp_lvl = gen_settings['cores to use'].item(), gen_settings['compression level'].item()

    settings = settings = pd.read_excel(Path(project).joinpath('Settings.xlsx'), sheet_name = '4_primer_trimming')
    p5_primer, p7_primer, anchoring = settings["P5 Primer (5' - 3')"].item(), settings["P7 Primer (5' - 3')"].item(), settings['anchoring'].item()

    ## collect the input files from PE merging step
    input = glob.glob(str(Path(project).joinpath('3_PE_merging', 'data', '*.fastq.gz')))

    print('{}: Starting to trim primers of {} input files.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(input)))

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath('4_primer_trimming', 'temp'))
    except FileExistsError:
        pass

    Parallel(n_jobs = cores)(delayed(primer_trimming)(file, project = project, p5_primer = p5_primer, p7_primer = p7_primer, anchoring = anchoring) for file in input)

    ## write logfile from pkl log, remove single logs after
    summary_logs = glob.glob(str(Path(project).joinpath('4_primer_trimming', 'temp', '*.pkl')))
    summary = [pickle.load(open(line, 'rb')) for line in summary_logs]

    ## generate log dataframe for primer trimming
    log_df = pd.DataFrame(summary, columns = ['File', 'finished at', 'python version', 'cutadapt version', 'processed reads', 'trimmed reads'])
    log_df = log_df.sort_values(by = 'File')
    log_df.to_excel(Path(project).joinpath('4_primer_trimming', 'Logfile_4_primer_trimming.xlsx'),
                    index = False,
                    sheet_name = '4_primer_trimming')

    ## add log to the project report
    wb = openpyxl.load_workbook(Path(project).joinpath('Project_report.xlsx'))
    writer = pd.ExcelWriter(Path(project).joinpath('Project_report.xlsx'), engine = 'openpyxl')
    writer.book = wb

    ## write the output
    log_df.to_excel(writer, sheet_name = '4_primer_trimming', index = False)
    wb.save(Path(project).joinpath('Project_report.xlsx'))
    writer.close()

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath('4_primer_trimming', 'temp'))

if __name__ == "__main__":
    main()

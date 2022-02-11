import sys, os, shutil, openpyxl, psutil, datetime
import pandas as pd
from pathlib import Path

def create_project(project_name):
    """Create a new metabarcoding pipeline project with all subfolder"""

    ## try to create the project folder
    try:
        os.mkdir('{}_apscale'.format(project_name))
    except FileExistsError:
        print('A project with that name already exists. Try another name.')
        return None

    ## generate the subfolder structure
    subfolders = ['1_raw data/data',
                  '2_demultiplexing/data',
                  '3_PE_merging/data',
                  '4_primer_trimming/data',
                  '5_quality_filtering/data',
                  '6_dereplication_pooling/data/dereplication',
                  '6_dereplication_pooling/data/pooling',
                  '7_otu_clustering/data',
                  '8_denoising/data']

    subfolders = [Path('{}_apscale'.format(project_name)).joinpath(subfolder) for subfolder in subfolders]

    for folder in subfolders:
        os.makedirs(folder)

    ## generate and then populate the settings file
    wb = openpyxl.Workbook()
    ws = wb.active
    wb.save(Path('{}_apscale'.format(project_name)).joinpath('Settings.xlsx'))
    wb = openpyxl.load_workbook(Path('{}_apscale'.format(project_name)).joinpath('Settings.xlsx'))
    writer = pd.ExcelWriter(Path('{}_apscale'.format(project_name)).joinpath('Settings.xlsx'), engine = 'openpyxl')
    writer.book = wb
    del wb['Sheet']

    ## write the 3_PE_merging sheet
    df_0 = pd.DataFrame([[int(psutil.cpu_count() * 0.75), 6]],
                        columns = ['cores to use', 'compression level'])

    df_0.to_excel(writer, sheet_name = '0_general_settings', index = False)

    ## write the 3_PE_merging sheet
    df_3 = pd.DataFrame([[25, 199, 5]],
                        columns = ['maxdiffpct', 'maxdiffs', 'minovlen'])

    df_3.to_excel(writer, sheet_name = '3_PE_merging', index = False)

    ## write the 4_primer_trimming sheet
    df_4 = pd.DataFrame([['', '', 'False']],
                        columns = ["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", 'anchoring'])

    df_4.to_excel(writer, sheet_name = '4_primer_trimming', index = False)

    ## write the 5_quality_filtering sheet
    df_5 = pd.DataFrame([[1, '', '']],
                        columns = ["maxEE", "min length", 'max length'])

    df_5.to_excel(writer, sheet_name = '5_quality_filtering', index = False)

    ## write the 6_dereplication_pooling sheet
    df_6 = pd.DataFrame([[4]],
                        columns = ["min size to pool"])

    df_6.to_excel(writer, sheet_name = '6_dereplication_pooling', index = False)


    ## write the 7_otu_clustering sheet
    df_7 = pd.DataFrame([[97, 'True', 'True']],
                        columns = ['pct id', 'to excel', 'to parquet'])

    df_7.to_excel(writer, sheet_name = '7_otu_clustering', index = False)

    ## write the 8_denoising sheet
    df_8 = pd.DataFrame([[2, 8, 'True', 'True']],
                        columns = ['alpha', 'minsize', 'to excel', 'to parquet'])

    df_8.to_excel(writer, sheet_name = '8_denoising', index = False)

    ## save the Settings file againg
    wb.save(Path('{}_apscale'.format(project_name)).joinpath('Settings.xlsx'))
    writer.close()

    ## give user output
    print('{}: "{}" created as a new project.'.format(datetime.datetime.now().strftime("%H:%M:%S"), project_name))

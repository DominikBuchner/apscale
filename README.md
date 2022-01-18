# apscale
Advanced Piepline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

## Introduction
Apscale is a metabarcoding pipeline that handles the most common tasks in metabarcoding
pipelines like paired-end merging, primer trimming, quality filtering, otu clustering and
denoising. It uses a simple command line interface and is configured via a single configuration file.
It automatically uses the available ressources on the machine it runs on while still providing the option
to use less if desired. All modules can be run on their own or as a comprehensive workflow.

A graphical user interface version for apscale is available [here](https://github.com/TillMacher/apscale_gui).

Programs used:
* vsearch (PE merging, quality filtering, otu clustering, denoising)
* cutadapt (primer trimming)

Input:
* demultiplexed gzipped reads

Output:
* log files, project report, OTU/ESV table

## Installation

Apscale can be installed on all common operating systems (Windows, Linux, MacOS).
Apscale requires Python 3.7 or higher and can be easily installed via pip in any command line:

`pip install apscale`

To update apscale run:

`pip install --upgrade apscale`

### Further dependencies - vsearch

Apscale calls vsearch for multiple modules. It should be installed and be in PATH to be executed
from anywhere on the system.

Check the vsearch Github page for further info:

https://github.com/torognes/vsearch

Support for compressed files with zlib is necessary. For Unix based systems this is shipped with
vsearch, for Windows the zlib.dll can be downloaded via:

[zlib for Windows](https://sourceforge.net/projects/mingw-w64/files/External%20binary%20packages%20%28Win64%20hosted%29/Binaries%20%2864-bit%29/zlib-1.2.5-bin-x64.zip/download)

The dll has to be in the same folder as the vsearch executable. If you need help with adding a folder to PATH in windows
please take a look at the first answer on this stackoverflow issue:

[How to add a folder to PATH Windows](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho)

To check if everything is correctly set up please type this into your command line:

`vsearch --version`

It should return a message similar to this:

```
vsearch v2.19.0_win_x86_64, 31.9GB RAM, 24 cores
https://github.com/torognes/vsearch

Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)
VSEARCH: a versatile open source tool for metagenomics
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584

Compiled with support for gzip-compressed files, and the library is loaded.
zlib version 1.2.5, compile flags 65
Compiled with support for bzip2-compressed files, but the library was not found.
```

### Further dependencies - cutadapt

Apscale also calls cutadapt with some modules. Cutadapt should be downloaded and installed
automatically with the Apscale installation. To check this, type:

`cutadapt --version`

and it should return the version number, for example:

`3.5`

## How to use

### Create a new apscale project

Apscale is organized in projects with the following structure:

<pre>
C:\USERS\DOMINIK\DESKTOP\EXAMPLE_PROJECT
├───1_raw data
│   └───data
├───2_demultiplexing
│   └───data
├───3_PE_merging
│   └───data
├───4_primer_trimming
│   └───data
├───5_quality_filtering
│   └───data
├───6_dereplication_pooling
│   └───data
│       ├───dereplication
│       └───pooling
├───7_otu_clustering
│   └───data
└───8_denoising
    └───data
</pre>

A new project can be initialized with the command:

`apscale --create_project NAME`

If you prefer to have your data all in one place you can copy the raw data into 1_raw_data/data.
Demultiplexing won't be handled by Apscale because there are to many different tagging systems out there to implement in a single pipeline.
If you are using inline barcodes you can take a look at https://github.com/DominikBuchner/demultiplexer.
If you are already starting with **demultiplexed** data please copy them into 2_demultiplexing/data.

### Configuring the general settings

Associated with every newly created project, Apscale will generate an Excel sheet in the project folder called "Settings.xlsx".
It is divided into seperate sheets for every module and a 0_general_settings tab.
By default Apscale will set 'cores to use' to all available cores on the system the project is created - 2, but this can be lowered
if the capacity is needed for other processes on your computer.
Apscale only works with compressed data, so it takes compressed data as input and has compressed data as output.
The compression level can be set in the general settings as well. Its default value is 6 since this is the default value of gzip.
The higher the compression level, the longer Apscale will take to process the data and vice versa, so if runtime is an issue, you can
lower the compression level and if disk space is a concern, the compression level can be set to 9 as maximum value.

### Configuring the specific settings

Apscale gives default values for most of its settings. They can be changed if desired, please refer to the manual of vsearch and cutadapt
for further information. The only value Apscale needs from the user is the primers used and the expected length of the fragment excluding which is used
for quality filtering. After these are set, Apscale is ready to run!

### Running Apscale

Navigate to the project folder you would like to process. Apscale can be run from anywhere on the system, but then needs a PATH to the project folder.

`apscale -h`

Will give help on the different functions of apscale.
To run an all-in-one analysis on your dataset run:

`apscale --run_apscale [PATH]`

The PATH argument is optional and only needs to be used if you are not located in the project folder.
This will automatically do PE merging, primer trimming, quality filter, OTU clustering and denoising of the data.
The individual modules can also be run separately (see `apscale -h` for respective commands). A project report will be saved in the project folder as well as an individual
report for the individual steps of the pipeline. Information about the versions of the programs used as well as how many reads where used and passed the module as well as a timestamp when the file finished.

The main output of Apscale will be an OTU table and an ESV table, as well as two .fasta files, which can be used for taxnomic assignment. For example, for COI sequences,
BOLDigger (https://github.com/DominikBuchner/BOLDigger) can be used directly with the output of Apscale to assign taxomoy to the OTUs / ESVs using the Barcode of Life Data system (BOLD) database. Furthermore, the ESV and OTU tables are compatible with TaxonTableTools (https://github.com/TillMacher/TaxonTableTools), which can be used for DNA metabarcoding specific analyses.

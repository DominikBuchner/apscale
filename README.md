# apscale
Advanced Piepline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

## Introduction
Apscale is a metabarcoding pipeline that handles the most common tasks in in metabarcoding
pipelines like paired-end merging, primer trimming, quality filtering, otu clustering and
denoising. It uses a simple command line interface and is configured via a single configuration file.
It automatically uses the available ressources on the machine it runs on while still providing the option
to use less if desired.

## Installation

Apscale can be installed on all common operating systems (Windows, Linux, MacOS).
Apscale requires Python 3.7 or higher and can be easily installed via pip in any command line:

`pip install apscale`

To update apscale run:

`pip install --upgrade apscale`

### Further dependencies

Apscale calls vsearch for multiple modules. It should be installed and be in PATH to be executed
from anywhere on the system.

Check the vsearch Github page for further info:

https://github.com/torognes/vsearch

Support for compressed files with zlib is necessary. For Unix based systems this is shipped with
vsearch, for Windows the zlib.dll can be downloaded via:

[zlib for Windows](https://sourceforge.net/projects/mingw-w64/files/External%20binary%20packages%20%28Win64%20hosted%29/Binaries%20%2864-bit%29/zlib-1.2.5-bin-x64.zip/download)

## How to use

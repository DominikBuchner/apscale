import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale",
    version="4.1.5",
    author="Dominik Buchner",
    author_email="dominik.buchner524@googlemail.com",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DominikBuchner/apscale",
    packages=setuptools.find_packages(),
    license="MIT",
    install_requires=[
        "openpyxl >= 3.0.10",
        "psutil >= 5.8.0",
        "duckdb>=1.3.1",
        "more_itertools >= 10.5.0",
        "demultiplexer2 >= 1.1.6",
        "joblib >= 1.0.0",
        "biopython >= 1.85",
        "cutadapt >= 5.0",
        "tqdm >= 4.56.0",
        "fastparquet >= 0.8.0",
        "numpy >= 2.0.0",
        "pandas >= 2.3.0",
        "powerlaw >= 1.5",
        "pyarrow >= 17.0.0",
        "pygbif >= 0.6.4",
        "pyproj >= 3.4.1",
        "Shapely >= 2.1.1",
        "streamlit >= 1.45.1",
        "tables >= 3.9.0",
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.11",
    entry_points={
        "console_scripts": [
            "apscale = apscale.__main__:main",
        ]
    },
)

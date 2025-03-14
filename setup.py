import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale",
    version="3.0.2",
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
        "pandas >= 1.5.0",
        "demultiplexer2 >= 1.1.1",
        "joblib >= 1.0.0",
        "biopython >= 1.78",
        "cutadapt >= 3.5",
        "tqdm >= 4.56.0",
        "fastparquet >= 0.8.0",
        "pyarrow >= 7.0.0",
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "apscale = apscale.__main__:main",
        ]
    },
)

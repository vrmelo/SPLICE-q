# SPLICE-q

A Python tool for genome-wide **SPLI**cing **E**fficiency **q**uantification from RNA-seq data.

[User Manual](https://github.com/vrmelo/SPLICE-q/wiki)

## Features
- Quantification of individual intron splicing efficiencies from strand-specific RNA-seq data.
- Sensitive to the overlap of genomic elements.
- Fast and user-friendly.

# Installation

SPLICE-q can be installed from pip and from source.

## pip
Using pip is the easiest way to install SPLICE-q.
```bash
 $ pip install SPLICE-q
```

## Development/install from source

```bash
 $ git clone https://github.com/vrmelo/SPLICE-q
 $ cd SPLICE-q
 $ pip install -e .
```

Requirements

- Python 3.6+
- PySam
- InterLap
- NumPy
- Rich

Operating Systems
- Linux, macOS, and Windows 10 Subsystem for Linux.
## Usage

To run SPLICE-q with default parameters, it requires a BAM file and a genome annotation file provided by GENCODE or Ensembl (GTF):
```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf
```
To specify an output file name and location: 

```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf -o outfile.tsv
```

Need help?
```bash
$ SPLICE-q.py -h
```
or check our [User Manual](https://github.com/vrmelo/SPLICE-q/wiki).

## Citation

de Melo Costa, Verônica R., et al. "SPLICE-q: a Python tool for genome-wide quantification of splicing efficiency." BMC bioinformatics 22.1 (2021): 1-14.

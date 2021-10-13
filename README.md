# dcc2
[![PyPI version](https://badge.fury.io/py/dcc2.svg)](https://pypi.org/project/dcc2/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/BIONF/dcc2.svg?branch=master)](https://travis-ci.com/BIONF/dcc2)
![Github Build](https://github.com/BIONF/dcc2/workflows/github%20build/badge.svg)

dcc2 is a tool for compiling core set data for [h1s](httsp://github.com/BIONF/HaMStR) using predicted orthologs from OMA, both OMA-browser and OMA-standalone. Outputs of this tool are 3 (optional 4) folders required for a HaMStR run, including (1) core_orthologs (comprises of OMA orthologous group - OG, or OMA pairs - OP. Each OG/OP has its own directory, where a multiple fasta file and a corresponding profile HMM can be found), (2) genome_dir (contains gene sets of taxa, from which the orthologs are originated), (3) blast_dir (holds the blast databases of those gene sets within genome_dir), and an optional (4) weight_dir (contains feature architecure annotations of all gene sets).

# Table of Contents
* [How to install](#how-to-install)
     * [Install the dcc2 package](#install-the-dcc2-package)
     * [Setup dcc2](#setup-dcc2)
* [Usage](#usage)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*dcc2* is distributed as a python package called *dcc2*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

## Install the dcc2 package
You can install *dcc2* using `pip`:
```
python3 -m pip install dcc2
```

or, in case you do not have admin rights, and don't use package systems like Anaconda to manage environments you need to use the `--user` option:
```
python3 -m pip install --user dcc2
```

and then add the following line to the end of your **~/.bashrc** or **~/.bash_profile** file, restart the current terminal to apply the change (or type `source ~/.bashrc`):

```
export PATH=$HOME/.local/bin:$PATH
```

## Setup dcc2

After installing *dcc2*, you need to run the prepare script to download and parse required OMA browser data.

You can do it by running this command
```
dcc2.prepare -o /output/path/for/oma/data
```

# Usage

For parsing OMA orthologs by using an OMA group ID:
```
dcc2.parseOmaById -g 1 -n HUMAN,THEAM,DESM0 -o /output/path/ -j jobName --cpus 8
```

Or using list of OMA taxa:

```
dcc2.parseOmaBySpec -n HUMAN,ECOLI,YEAST -o /output/path/ -j jobName --annoFas --cpus 8
```

If only 2 OMA taxa are given, you can choose to use OMA pairs instead of OMA groups:
```
dcc2.parseOmaBySpec -n HUMAN,ECOLI -t pair -o /output/path/ -j jobName --annoFAS --cpus 8
```

For parsing an output from OMA-standalone, dcc2 requires:

- the output orthoXML file from OMA,
- a taxon mapping file in tab-delimited format containing 3 columns `<NCBI taxon ID>   <Original taxon name>   <Abbr. taxon name>`, where *original taxon name* is the name that is written in the orthoXML input file and *abbr. taxon name* is its abbreviation species code (for example: HOMSA for *Homo sapiens*)
- protein set of included taxa. This can be either given as a folder, or automatically downloaded from OMA database

```
dcc2.parseOrthoxml -i input.orthoxml -m mapping_file.txt -g /path/to/gene/set -o /output/path/ -j jobName --annoFas --cpus 8
```

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/dcc2/issues/new) or be in touch via email.

# Contributors
- [Vinh Tran](https://github.com/trvinh)
- [Hannah Mülbaier](https://github.com/mueli94)

# Contact
For further support or bug reports please contact: tran@bio.uni-frankfurt.de

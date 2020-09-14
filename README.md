# dcc2

DCCv2 is a tool for compiling core set data for [h1s](httsp://github.com/BIONF/HaMStR) using predicted orthologs from OMA, both OMA-browser and OMA-standalone. Outputs of this tool are 3 (optional 4) folders required for a HaMStR run, including (1) core_orthologs (comprises of OMA orthologous group - OG, or OMA pairs - OP. Each OG/OP has its own directory, where a multiple fasta file and a corresponding profile HMM can be found), (2) genome_dir (contains gene sets of taxa, from which the orthologs are originated), (3) blast_dir (holds the blast databases of those gene sets within genome_dir), and an optional (4) weight_dir (contains feature architecure annotations of all gene sets).

# Table of Contents
* [How to install](#how-to-install)
     * [Install the dcc2 package](#install-the-dcc2-package)
     * [Setup dcc2](#setup-dcc2)
* [Usage](#usage)
* [Bugs](#bugs)
* [How to cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*dcc2* is distributed as a python package called *dcc2*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

## Install the h1s package
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

## Setup HaMStR-oneSeq

After installing *dcc2*, you need to run the prepare script to download and parse required OMA browser data.

You can do it by running this command
```
prepareDcc -o /output/path/for/oma/data
```

# Usage

For parsing OMA orthologs by using an OMA group ID:
```
parseOmaById -g 1 -n HUMAN,THEAM,DESM0 -o /output/path/ -j jobName
```

Or using list of OMA taxa:

```
parseOmaBySpec -n HUMAN,ECOLI,YEAST -o /output/path/ -j jobName --annoFas --cpus 8
```

If only 2 OMA taxa are given, you can choose to use OMA pairs instead of OMA groups:
```
parseOmaBySpec -n HUMAN,ECOLI -o /output/path/ -j jobName -t pair
```

For parsing an output from OMA-standalone, dcc2 requires:
    - the output orthoXML file from OMA,
    - a taxon mapping file in tab-delimited format containing 3 columns (blabla) (blabla) (blabla),
    - protein set of included taxa. This can be either given as a folder, or automatically downloaded from OMA database

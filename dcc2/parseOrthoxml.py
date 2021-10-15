# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This file is part of dcc2.
#
#  dcc2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  dcc2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with phylosophy.  If not, see <http://www.gnu.org/licenses/>.
#
#  Contact: tran@bio.uni-frankfurt.de
#
#######################################################################

import os
import sys
import argparse
import time
from bs4 import BeautifulSoup
from pathlib import Path
import subprocess
from Bio import SeqIO
import multiprocessing as mp
import dcc2.dccFn as dccFn
from datetime import datetime
from tqdm import tqdm

def readFileToDict(file):
    name2id = {}
    name2abbr = {}
    with open(file, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t')
            name2id[fields[1]] = fields[0]
            name2abbr[fields[1]] = fields[2]
    return(name2id, name2abbr)

def concatFasta(fileIn, fileOut):
     cmd = "awk \'/^>/ { print (NR==1 ? \"\" : RS) $0; next } { printf \"%s\", $0 } END { printf RS }\' " + fileIn + " > " + fileOut
     subprocess.call([cmd], shell = True)

def main():
    version = "0.2.3"
    parser = argparse.ArgumentParser(description="You are running dcc2 version " + str(version))
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('additional arguments')
    required.add_argument('-i', '--inFile', help='Input sequence in orthoXML format', action='store', default='', required=True)
    required.add_argument('-o', '--outPath', help='Output directory', action='store', default='', required=True)
    required.add_argument('-g', '--geneSet', help='Path to gene set folder', action='store', default='', required=True)
    required.add_argument('-m', '--mappingFile', help='NCBI taxon ID mapping file', action='store', default='', required=True)
    required.add_argument('-j', '--jobName', help='Job name', action='store', default='', required=True)
    optional.add_argument('-a', '--alignTool', help='Alignment tool (mafft|muscle). Default: mafft', action='store', default='mafft')
    optional.add_argument('-f', '--annoFas', help='Perform FAS annotation', action='store_true')
    optional.add_argument('-l', '--maxGroups', help='Maximum ortholog groups taken into account.', type=int, action='store', default=999999999)
    optional.add_argument('-c', '--cpus', help='Number of CPUs. Default: 4', action='store', default=4, type=int)
    args = parser.parse_args()

    inFile = args.inFile
    dccFn.checkFileExist(inFile)
    dataPath = args.geneSet
    dccFn.checkFileExist(dataPath)
    mappingFile = args.mappingFile
    dccFn.checkFileExist(mappingFile)
    outPath = args.outPath
    try:
        my_abs_path = Path(outPath).resolve(strict=True)
    except FileNotFoundError:
        Path(outPath).mkdir(parents = True, exist_ok = True)
    aligTool = args.alignTool.lower()
    if not (aligTool == "mafft" or aligTool == "muscle"):
        sys.exit("alignment tool must be either mafft or muscle")
    limit = args.maxGroups
    doAnno = args.annoFas
    jobName = args.jobName
    cpus = args.cpus

    start = time.time()
    # create job pool
    pool = mp.Pool(cpus)
    if cpus > (mp.cpu_count()):
        print('Reduce the given number of CPUs to %s' % (mp.cpu_count()))
        pool = mp.Pool(mp.cpu_count())

    ##### read mapping file
    (name2id, name2abbr) = readFileToDict(mappingFile)

    ##### read input file into beatifulsoup object
    print("Reading input XML file...")
    xmlIn = BeautifulSoup(open(inFile),"xml")

    ##### PARSING XML FILE
    ### get list of species together with NCBI taxon IDs and their corresponding genes
    taxonName = {}
    protID = {}
    fasta = {}

    ### create output folders
    Path(outPath + "/genome_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/blast_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/core_orthologs/" + jobName).mkdir(parents = True, exist_ok = True)
    Path(outPath + "/weight_dir").mkdir(parents = True, exist_ok = True)

    ### copy species to genome_dir, blast_dir, weight_dir
    print("Getting gene sets...")
    blastJobs = []
    annoJobs = []
    for spec in xmlIn.findAll("species"):
        specNameOri = spec.get("name")
        if not specNameOri in name2abbr:
            sys.exit("%s not found in %s" % (specNameOri, mappingFile))
        specName = "%s@%s@1" % (name2abbr[specNameOri], name2id[specNameOri])
        Path(outPath + "/genome_dir/" + specName).mkdir(parents = True, exist_ok = True)
        Path(outPath + "/blast_dir/" + specName).mkdir(parents = True, exist_ok = True)
        # get gene set file
        lsCmd = 'ls %s/%s.*' % (dataPath, specNameOri)
        specFile = subprocess.check_output([lsCmd], shell = True).decode(sys.stdout.encoding).strip()
        fileExt = specFile.split(".")[-1]
        # read fasta file to dictionary
        fasta[specName] = SeqIO.to_dict(SeqIO.parse(open(specFile),'fasta'))

        # copy to genome_dir/specName/specName.fa and make smybolic link to blast_dir/specName
        fileInGenome = "%s/genome_dir/%s/%s.fa" % (outPath, specName, specName)
        if not Path(fileInGenome).exists():
            concatFasta(specFile, fileInGenome)
        fileInBlast = "%s/blast_dir/%s/%s.fa" % (outPath, specName, specName)
        if not Path(fileInBlast).exists():
            fileInGenomeMod = "../../genome_dir/%s/%s.fa" % (specName, specName)
            os.symlink(fileInGenomeMod, fileInBlast)
        # write .checked file
        checkedFile = dccFn.openFileToWrite("%s/genome_dir/%s/%s.fa.checked" % (outPath, specName, specName))
        checkedFile.write(str(datetime.now()))
        checkedFile.close()
        # get info for blast
        blastDbFile = "%s/blast_dir/%s/%s.phr" % (outPath, specName, specName)
        if not Path(blastDbFile).exists():
            blastJobs.append([specName, fileInBlast, outPath])
        # get info for FAS annotation
        annoFile = "%s/weight_dir/%s.json" % (outPath, specName)
        if not Path(annoFile).exists():
            annoJobs.append(specFile)

        # save OG members and their spec name to dict
        for gene in spec.findAll("gene"):
            groupID = gene.get("id")
            orthoID = gene.get("protId")
            taxonName[orthoID] = specName
            protID[groupID] = orthoID

    # make blastDB
    print("Creating BLAST databases...")
    if dccFn.is_tool('makeblastdb'):
        msa = pool.map(dccFn.runBlast, blastJobs)
    else:
        print("makeblastdb not found!")

    ### parse ortholog groups
    print("Parsing ortholog groups...")
    alignJobs = []
    hmmJobs = []
    n = 0
    for orthogroup in xmlIn.findAll("orthologGroup"):
        groupID = orthogroup.get("id")
        if groupID:
            n = n + 1
            if (n > limit):
                break
            if groupID.isdigit():
                groupID = "OG_"+str(groupID)
            Path(outPath + "/core_orthologs/" + jobName + groupID).mkdir(parents = True, exist_ok = True)

            # get fasta sequences
            with open(outPath + "/core_orthologs/" + jobName + "/" + groupID + "/" + groupID + ".fa", "w") as myfile:
                for ortho in orthogroup.findAll("geneRef"):
                    orthoID = protID[ortho.get("id")]
                    spec = taxonName[orthoID]
                    orthoSeq = str(fasta[spec][orthoID].seq)
                    myfile.write(">" + groupID + "|" + spec + "|" + orthoID + "\n" + orthoSeq + "\n")

            # get info for MSA
            ogFasta = outPath + "/core_orthologs/" + jobName + "/" + groupID + "/" + groupID
            alignJobs.append([ogFasta, aligTool, groupID])

            # get info for pHMM
            Path(outPath + "/core_orthologs/" + jobName + "/" + groupID + "/hmm_dir").mkdir(parents = True, exist_ok = True)
            hmmFile = "%s/core_orthologs/%s/%s/hmm_dir/%s.hmm" % (outPath, jobName, groupID, groupID)
            flag = 0
            try:
                if os.path.getsize(hmmFile) == 0:
                    flag = 1
            except OSError as e:
                    flag = 1
            if flag == 1:
                hmmJobs.append([hmmFile, ogFasta, groupID])

    ### create MSAs and pHMMs
    print("Calculating MSAs and pHMMs for %s OGs..." % (len(alignJobs)))
    # if dccFn.is_tool(aligTool + " -h"):
    msaOut = []
    for _ in tqdm(pool.imap_unordered(dccFn.runMsa, alignJobs), total=len(alignJobs)):
        msaOut.append(_)
    if dccFn.is_tool('hmmbuild'):
        hmmOut = []
        for _ in tqdm(pool.imap_unordered(dccFn.runHmm, hmmJobs), total=len(hmmJobs)):
            hmmOut.append(_)
    else:
        print("hmmbuild not found!")

    # do FAS annotation
    if doAnno:
        print("Doing FAS annotation...")
        if dccFn.is_tool('fas.doAnno'):
            for specFile in annoJobs:
                dccFn.calcAnnoFas(specFile, outPath, cpus)

    pool.close()
    end = time.time()
    print("Finished in " + '{:5.3f}s'.format(end-start))
    print("Output can be found in %s" % outPath)

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Hannah Mülbaier & Vinh Tran
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
#  Contact: hannah.muelbaier@gmail.com or tran@bio.uni-frankfurt.de
#
#######################################################################

import sys
import os
import argparse
from pathlib import Path
import time
from Bio import SeqIO
import multiprocessing as mp
import dcc2.dccFn as dccFn

def getOGprot(dataPath, omaGroupId, speciesList):
    fileGroups = dccFn.openFileToRead(dataPath + "/oma-groups.txt")
    allGroups = fileGroups.readlines()
    fileGroups.close()
    groupLine = allGroups[int(omaGroupId) + 2].strip().split("\t")
    proteinIds = []
    for prot in groupLine[2:]:
        if prot[0:5] in speciesList:
            proteinIds.append(prot)
    return(proteinIds)

def main():
    version = "1.0.0"
    parser = argparse.ArgumentParser(description="You are running omaParser by OG ID for OMA Browser version " + str(version) + ".")
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('additional arguments')
    required.add_argument('-g', '--OG', help='Input OMA group ID', action='store', default='', required=True)
    required.add_argument('-n', '--name', help='List of OMA species abbr. names', action='store', default='', required=True)
    required.add_argument('-i', '--id', help='List of corresponding taxonomy IDs to OMA species', action='store', default='', required=True)
    required.add_argument('-d', '--dataPath', help='Path to OMA Browser data', action='store', default='', required=True)
    required.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='', required=True)
    required.add_argument('-j', '--jobName', help='Job name', action='store', default='', required=True)
    optional.add_argument('-a', '--alignTool', help='Alignment tool (mafft|muscle). Default: mafft', action='store', default='mafft')
    optional.add_argument('-f', '--annoFas', help='Perform FAS annotation', action='store_true')
    args = parser.parse_args()

    dccFn.checkFileExist(args.dataPath)
    dccFn.checkFileExist(args.dataPath+"/oma-seqs-dic.fa")
    dataPath = str(Path(args.dataPath).resolve())
    omaGroupId = args.OG
    speciesList = str(args.name).split(",")
    speciesTaxId = str(args.id).split(",")
    outPath = str(Path(args.outPath).resolve())
    aligTool = args.alignTool.lower()
    if not aligTool == "mafft" or aligTool == "muscle":
        sys.exit("alignment tool must be either mafft or muscle")
    doAnno = args.annoFas
    jobName = args.jobName

    start = time.time()
    pool = mp.Pool(mp.cpu_count()-2)

    ### create output folders
    print("Creating output folders...")
    Path(outPath + "/genome_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/blast_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/core_orthologs/" + jobName +  "/" + omaGroupId + "/hmm_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/weight_dir").mkdir(parents = True, exist_ok = True)

    ### create spec IDs dict
    specName2id = dict(zip(speciesList, speciesTaxId))

    ### Get genesets
    print("Getting %s gene sets..." % (len(speciesList)))
    dccFn.getGeneset(dataPath, speciesList, speciesTaxId, outPath)

    # read fasta file to dictionary
    fasta = {}
    blastJobs = []
    annoJobs = []
    for i in range(0,len(speciesList)):
        fileName = dccFn.makeOneSeqSpeciesName(speciesList[i], speciesTaxId[i])
        specFile = outPath+"/genome_dir/"+fileName+"/"+fileName+".fa"
        fasta[speciesList[i]] = SeqIO.to_dict(SeqIO.parse(open(specFile),'fasta'))
        # get info for BLAST
        blastDbFile = "%s/blast_dir/%s/%s.phr" % (outPath, fileName, fileName)
        if not Path(blastDbFile).exists():
            blastJobs.append([fileName, specFile, outPath])
        # get info for FAS annotation
        annoFile = "%s/weight_dir/%s.json" % (outPath, fileName)
        if not Path(annoFile).exists():
            annoJobs.append([specFile, outPath])

    ### create blastDBs
    print("Creating BLAST databases for %s taxa..." % len(blastJobs))
    if dccFn.is_tool('makeblastdb'):
        msa = pool.map(dccFn.runBlast, blastJobs)
    else:
        print("makeblastdb not found!")

    ### get OG fasta
    print("Getting protein sequences for OG id %s..." % omaGroupId)
    proteinIds = getOGprot(dataPath, omaGroupId, speciesList)
    dccFn.getOGseq([proteinIds, omaGroupId, outPath, fasta, jobName])

    ### calculate MSAs and pHMMs
    ogFasta = outPath + "/core_orthologs/" + jobName +  "/" + omaGroupId + "/" + omaGroupId
    # do MSAs
    try:
        msaFile = "%s/core_orthologs/%s/%s/%s.aln" % (outPath, jobName, omaGroupId, omaGroupId)
        flag = dccFn.checkFileEmpty(msaFile)
        if flag == True:
            dccFn.runMsa([ogFasta, aligTool, omaGroupId])
    except:
        sys.exit("%s not found or %s not works correctly!" % (ogFasta+".fa", aligTool))
    # do pHMMs
    if dccFn.is_tool('hmmbuild'):
        try:
            hmmFile = "%s/core_orthologs/%s/%s/hmm_dir/%s.hmm" % (outPath, jobName, omaGroupId, omaGroupId)
            flag = dccFn.checkFileEmpty(hmmFile)
            if flag == True:
                dccFn.runHmm([hmmFile, ogFasta, omaGroupId])
        except:
            sys.exit("hmmbuild not works correctly for %s!" % (ogFasta+".fa"))
    else:
        print("hmmbuild not found!")

    ### do FAS annotation
    if doAnno:
        print("Doing FAS annotation...")
        if dccFn.is_tool('annoFAS'):
            anno = pool.map(dccFn.calcAnnoFas, annoJobs)

    pool.close()
    end = time.time()
    print("Finished in " + '{:5.3f}s'.format(end-start))
    print("Output can be found in %s" % outPath)

if __name__ == '__main__':
    main()

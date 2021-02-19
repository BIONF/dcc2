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
from tqdm import tqdm
import dcc2.dccFn as dccFn

def getOriOG(dataPath, omaGroupId):
    fileGroups = dccFn.openFileToRead(dataPath + "/oma-groups.txt")
    allGroups = fileGroups.readlines()
    fileGroups.close()
    groupLine = allGroups[int(omaGroupId) + 2].strip().split("\t")
    return(groupLine)

def getValidSpec(omaGroup, speciesList):
    specList = []
    for spec in speciesList:
        if spec in " ".join(omaGroup):
            specList.append(spec)
        else:
            print('\033[91m%s not found in given OMA group\033[0m' % (spec))
    return(specList)

def getOGprot(omaGroup, speciesList):
    proteinIds = []
    for prot in omaGroup[2:]:
        if prot[0:5] in speciesList:
            proteinIds.append(prot)
    return(proteinIds)

def parseOma(args):
    (omaGroup, omaGroupId, speciesList, speciesTaxId, outPath, jobName, aligTool, doAnno, cpus, dataPath) = args
    # create job pool
    pool = mp.Pool(cpus)
    if cpus > (mp.cpu_count()):
        print('Reduce the given number of CPUs to %s' % (mp.cpu_count()))
        pool = mp.Pool(mp.cpu_count())

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
            annoJobs.append(specFile)

    ### create blastDBs
    print("Creating BLAST databases for %s taxa..." % len(blastJobs))
    if dccFn.is_tool('makeblastdb'):
        # msa = pool.map(dccFn.runBlast, blastJobs)
        blastOut = []
        for _ in tqdm(pool.imap_unordered(dccFn.runBlast, blastJobs), total=len(blastJobs)):
            blastOut.append(_)
    else:
        print("makeblastdb not found!")

    ### get OG fasta
    print("Getting OMA OG id %s..." % omaGroupId)
    proteinIds = getOGprot(omaGroup, speciesList)
    dccFn.getOGseq([proteinIds, omaGroupId, outPath, fasta, specName2id, jobName])

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
            for specFile in annoJobs:
                dccFn.calcAnnoFas(specFile, outPath, cpus)

    pool.close()

def main():
    version = "0.1.0"
    parser = argparse.ArgumentParser(description="You are running dcc2 version " + str(version))
    required = parser.add_argument_group('required arguments')
    required.add_argument('-g', '--OG', help='Input one OMA group ID', action='store', default='', required=True)
    required.add_argument('-n', '--name', help='List of OMA species abbr. names', action='store', default='', required=True)
    required.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='', required=True)
    required.add_argument('-j', '--jobName', help='Job name', action='store', default='', required=True)
    optional = parser.add_argument_group('additional arguments')
    optional.add_argument('-d', '--dataPath', help='Path to OMA Browser data', action='store', default='')
    optional.add_argument('-a', '--alignTool', help='Alignment tool (mafft|muscle). Default: mafft', choices=['mafft', 'muscle'], action='store', default='mafft')
    optional.add_argument('-f', '--annoFas', help='Perform FAS annotation', action='store_true')
    optional.add_argument('-c', '--cpus', help='Number of CPUs. Default: 4', action='store', default=4, type=int)
    args = parser.parse_args()

    omaGroupId = args.OG
    speciesList = str(args.name).split(",")
    outPath = str(Path(args.outPath).resolve())
    aligTool = args.alignTool.lower()
    doAnno = args.annoFas
    jobName = args.jobName
    cpus = args.cpus

    start = time.time()
    # get OMA browser data path
    dataPath = args.dataPath
    if dataPath == '':
        pathconfigFile = os.path.realpath(__file__).replace('parseOmaById.py', 'pathconfig.txt')
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run prepareDcc first!')
        with open(pathconfigFile) as f:
            dataPath = f.readline().strip()
    else:
        dataPath = os.path.abspath(args.dataPath)

    # get valid species list
    print('Checking valid taxa...')
    omaGroup = getOriOG(dataPath, omaGroupId)
    validSpecList = getValidSpec(omaGroup, speciesList)
    if not validSpecList:
        sys.exit('None of the given species could be found in OMA group %s' % omaGroupId)

    # get NCBI Taxon IDs for input species
    speciesTaxId = []
    for spec in validSpecList:
        taxId = dccFn.getTaxonId(dataPath, spec)
        if not taxId:
            sys.exit('%s is not a valid OMA species code!' % spec)
        else:
            speciesTaxId.append(str(taxId))

    # create output folders
    print("Creating output folders...")
    Path(outPath + "/genome_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/blast_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/core_orthologs/" + jobName +  "/" + omaGroupId + "/hmm_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/weight_dir").mkdir(parents = True, exist_ok = True)

    # do parsing

    parseOma([omaGroup, omaGroupId, validSpecList, speciesTaxId, outPath, jobName, aligTool, doAnno, cpus, dataPath])
    end = time.time()
    print("Finished in " + '{:5.3f}s'.format(end-start))
    print("Output can be found at %s" % outPath)

if __name__ == '__main__':
    main()

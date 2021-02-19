# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran & Hannah Mülbaier
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
#  Contact: tran@bio.uni-frankfurt.de or hannah.muelbaier@gmail.com
#
#######################################################################

import sys
import os
import argparse
from pathlib import Path
import time
from Bio import SeqIO
import multiprocessing as mp
from omadb import Client
import dcc2.dccFn as dccFn
from tqdm import tqdm

def gettingOmaGroups(dataPath, speciesSet, nr):
    omaGroupFull = dccFn.openFileToRead(dataPath + "/oma-groups.txt")
    omaGroups = {}
    for i in tqdm(omaGroupFull, desc = 'OMA groups'):
        row = i.split("\t")
        if len(row) != 1:
            ProteinIds = []
            counter = 0
            for j in range(2, len(row)):
                species = row[j][0:5]
                if species in speciesSet:
                    ProteinIds.append(row[j])
                    counter += 1
            if counter >= (len(speciesSet) - int(nr)):
                omaGroups[row[0]] = ProteinIds
    return(omaGroups)

def gettingOmaPairs(speciesList):
    c = Client()
    omaPairs = list(c.pairwise(speciesList[0], speciesList[1], progress = True))
    out = {}
    nr = 1
    for pair in tqdm(omaPairs, desc = 'OMA Pairs'):
        id = 'OP' + pair.rel_type.replace(':','') + "_" + str(nr)
        out[id] = [pair.entry_1.omaid, pair.entry_2.omaid]
        nr = nr + 1
        # if not str(pair.entry_1.oma_group) in out:
        #     out[str(pair.entry_1.oma_group)] = [pair.entry_1.omaid, pair.entry_2.omaid]
        # else:
        #     for id in [pair.entry_1.omaid, pair.entry_2.omaid]:
        #         if id not in out[str(pair.entry_1.oma_group)]:
        #             out[str(pair.entry_1.oma_group)].append(id)
    return(out)

def parseOma(args):
    (speciesList, speciesTaxId, outPath, jobName, nrMissingTaxa, omaType, aligTool, doAnno, cpus, dataPath) = args
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
    print("Loading fasta files...")
    for i in tqdm(range(0,len(speciesList))):
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
        blastOut = []
        for _ in tqdm(pool.imap_unordered(dccFn.runBlast, blastJobs), total=len(blastJobs)):
            blastOut.append(_)
    else:
        print("makeblastdb not found!")

    ### get OGs and their fasta
    if omaType == 'group':
        omaGroups = gettingOmaGroups(dataPath, set(speciesList), nrMissingTaxa)
    elif omaType == 'pair':
        omaGroups = gettingOmaPairs(speciesList)
    print("Getting protein sequences for %s OGs..." % len(omaGroups))
    msaJobs = []
    hmmJobs = []
    for omaGroupId in tqdm(omaGroups):
        # Path(outPath + "/core_orthologs/" + omaGroupId).mkdir(parents = True, exist_ok = True)
        Path(outPath + "/core_orthologs/" + jobName + "/" + omaGroupId + "/hmm_dir").mkdir(parents = True, exist_ok = True)
        # getSeqJobs.append([omaGroups[omaGroupId], omaGroupId, outPath, fasta])  # slower than run sequentially
        dccFn.getOGseq([omaGroups[omaGroupId], omaGroupId, outPath, fasta, specName2id, jobName])

        ogFasta = outPath + "/core_orthologs/" + jobName + "/" + omaGroupId + "/" + omaGroupId
        ### get MSA jobs
        try:
            msaFile = "%s/core_orthologs/%s/%s/%s.aln" % (outPath, jobName, omaGroupId, omaGroupId)
            flag = dccFn.checkFileEmpty(msaFile)
            if flag == True:
                msaJobs.append([ogFasta, aligTool, omaGroupId])
        except:
            sys.exit("%s not found or %s not works correctly!" % (ogFasta+".fa", aligTool))
        ### get pHMM jobs
        try:
            hmmFile = "%s/core_orthologs/%s/%s/hmm_dir/%s.hmm" % (outPath, jobName, omaGroupId, omaGroupId)
            flag = dccFn.checkFileEmpty(hmmFile)
            if flag == True:
                hmmJobs.append([hmmFile, ogFasta, omaGroupId])
        except:
            sys.exit("hmmbuild not works correctly for %s!" % (ogFasta+".fa"))

    ### calculate MSAs and pHMMs
    print("Calculating MSAs and pHMMs for %s OGs..." % len(omaGroups))
    if len(msaJobs) > 0:
        msaOut = []
        for _ in tqdm(pool.imap_unordered(dccFn.runMsa, msaJobs), total=len(msaJobs)):
            msaOut.append(_)
    if len(hmmJobs) > 0:
        if dccFn.is_tool('hmmbuild'):
            hmmOut = []
            for _ in tqdm(pool.imap_unordered(dccFn.runHmm, hmmJobs), total=len(hmmJobs)):
                hmmOut.append(_)
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
    required.add_argument('-n', '--name', help='List of OMA species abbr. names', action='store', default='', required=True)
    required.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='', required=True)
    required.add_argument('-j', '--jobName', help='Job name', action='store', default='', required=True)
    optional = parser.add_argument_group('additional arguments')
    optional.add_argument('-d', '--dataPath', help='Path to OMA Browser data', action='store', default='')
    optional.add_argument('-m', '--missingTaxa', help='Number of allowed missing taxa. Default: 0', action='store', default=0, type=int)
    optional.add_argument('-t', '--omaType', help='OMA type (group|pair). Default: "group"', choices=['group', 'pair'], action='store', default='group')
    optional.add_argument('-a', '--alignTool', help='Alignment tool (mafft|muscle). Default: mafft', choices=['mafft', 'muscle'], action='store', default='mafft')
    optional.add_argument('-f', '--annoFas', help='Perform FAS annotation', action='store_true')
    optional.add_argument('-c', '--cpus', help='Number of CPUs. Default: 4', action='store', default=4, type=int)
    args = parser.parse_args()

    speciesList = str(args.name.upper()).split(",")
    outPath = os.path.abspath(args.outPath)
    jobName = args.jobName
    nrMissingTaxa = args.missingTaxa
    omaType = args.omaType
    aligTool = args.alignTool.lower()
    doAnno = args.annoFas
    cpus = args.cpus

    start = time.time()
    # check oma pair condition
    if omaType == 'pair':
        if len(speciesList) > 2:
            sys.exit('OMA Pair works only with 2 input species!')

    # get OMA browser data path
    dataPath = args.dataPath
    if dataPath == '':
        pathconfigFile = os.path.realpath(__file__).replace('parseOmaBySpec.py', 'pathconfig.txt')
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run prepareDcc first!')
        with open(pathconfigFile) as f:
            dataPath = f.readline().strip()
    else:
        dataPath = os.path.abspath(args.dataPath)

    # get NCBI Taxon IDs for input species
    speciesTaxId = []
    for spec in speciesList:
        taxId = dccFn.getTaxonId(dataPath, spec)
        if not taxId:
            sys.exit('%s is not a valid OMA species code!' % spec)
        else:
            speciesTaxId.append(str(taxId))

    # create output folders
    print("Creating output folders...")
    Path(outPath + "/genome_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/blast_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/core_orthologs/" + jobName).mkdir(parents = True, exist_ok = True)
    Path(outPath + "/weight_dir").mkdir(parents = True, exist_ok = True)

    # do parsing
    parseOma([speciesList, speciesTaxId, outPath, jobName, nrMissingTaxa, omaType, aligTool, doAnno, cpus, dataPath])
    end = time.time()
    print("Finished in " + '{:5.3f}s'.format(end-start))
    print("Output can be found at %s" % outPath)

if __name__ == '__main__':
    main()

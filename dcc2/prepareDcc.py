# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran & Hannah Mülbaier
#
#  This file is part of dcc2.
#
#  phylosophy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  phylosophy is distributed in the hope that it will be useful,
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

from Bio import SeqIO
import time
import json
import sys
import os
import subprocess
from pathlib import Path
import shutil
import gzip
import argparse
import ssl
import urllib.request
import time
from tqdm import tqdm
import dcc2.dccFn as dccFn

def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)
    while True:
        choice = sys.stdin.readline().rstrip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with "yes" or "no" '
                             '(or "y" or "n").\n')

def download_progress(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    if percent > 100:
        percent = 100
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

def download_file(url, file, outDir):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    download_file = urllib.request.URLopener(context=ctx)
    print('Downloading %s' % (url + '/' + file))
    download_file.retrieve(url + '/' + file, outDir + '/' + file, download_progress)
    print(' ... done!')

def downloadFiles(dataPath, keep):
    mainUrl = "https://omabrowser.org/All"
    zipFiles = ['oma-groups.txt.gz', 'oma-seqs.fa.gz']
    specinfo = "oma-species.txt"
    # oma-groups & oma-seqs
    for f in zipFiles:
        try:
            download_file(mainUrl, f, dataPath)
            try:
                print('Extracting...')
                with gzip.open(dataPath + '/' + f, 'rb') as f_in:
                    with open(dataPath + '/' + f.replace('.gz',''), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                if not keep:
                    os.remove(dataPath + '/' + f)
            except:
                sys.exit("Error occurs while extrating %s" % (f))
        except:
            sys.exit("Error occurs while downloading %s/%s" % (mainUrl, f))
    # oma-species
    try:
        download_file(mainUrl, specinfo, dataPath)
    except:
        sys.exit("Error occurs while downloading %s/%s" % (mainUrl, specinfo))

def createDicSpecies(proteins, file):
    sequenceDic = {}
    code = str(proteins.readline()[2:7])
    startline = 0
    lineNr = 0

    for i in tqdm(proteins):
        lineNr += 1
        if code != i[2:7] and i[0] == ">":
            endline = lineNr - 1
            sequenceDic[code] = (startline,endline)
            code = i[2:7]
            startline = lineNr

    sequenceDic[code] = (startline,lineNr)
    json.dump(sequenceDic, file)

def createDicOmaGroup(omaGroups, file):
    groupDic = {}
    for i in tqdm(omaGroups):
        if i[0] != "#":
            line = i.split("\t")
            speciesSet = set()
            groupId = line[0]
            for j in range(2, len(line)):
                species = str(line[j])[0:5]
                speciesSet.add(species)
            groupDic[groupId] = tuple(speciesSet)

    for key in groupDic:
        speciesStr = str(groupDic[key]).replace("(", "")
        speciesStr = speciesStr.replace(")", "")
        speciesStr = speciesStr.replace("'", "")
        speciesStr = speciesStr.replace(" ", "")
        file.write(key + "\t" + speciesStr + "\n")

def prepareDcc(dataPath, keep):
    start = time.time()
    downloadFiles(dataPath, keep)

    allProteins = open(dataPath + "/oma-seqs.fa", "r")
    newFileSpecies = open(dataPath + "/oma-seqs-dic.fa", "w")
    omaGroups = open(dataPath + "/oma-groups.txt", "r")
    newFileOmaGroup = open(dataPath + "/oma-groups-tmp.txt", "w")

    print("############ Indexing proteins...")
    createDicSpecies(allProteins, newFileSpecies)
    print("############ Indexing OMA groups...")
    createDicOmaGroup(omaGroups, newFileOmaGroup)

    newFileSpecies.close()
    newFileOmaGroup.close()
    allProteins.close()
    omaGroups.close()
    end = time.time()
    print("Finished in " + '{:5.3f}s'.format(end-start))

def checkOmaData(dataPath):
    files = ['oma-seqs.fa', 'oma-groups.txt', 'oma-seqs-dic.fa', 'oma-groups-tmp.txt', 'oma-species.txt']
    for f in files:
        dccFn.checkFileExist(dataPath + '/' + f)
        flag = dccFn.checkFileEmpty(dataPath + '/' + f)
        if flag:
            sys.exit('%s is empty!' % (dataPath + '/' + f))

def getPath(dccPath):
    if os.path.exists(dccPath + '/pathconfig.txt'):
        with open(dccPath+'/pathconfig.txt', 'r') as file:
            savedPath = file.read().strip()
            return(savedPath)
    else:
        return('Pathconfig file not found')

def main():
    version = "0.1.0"
    parser = argparse.ArgumentParser(description="You are running dcc2 version " + str(version))
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('additional arguments')
    required.add_argument('-o', '--outPath', help='Directory for saving OMA Browser data', action='store', default='', required=True)
    optional.add_argument('-f', '--force', help='Force override old data (used for reinstalling or updating existing OMA data)', action='store_true', default=False)
    optional.add_argument('-s', '--savePath', help='Save existing OMA Browser path to config file for dcc2', action='store_true', default=False)
    optional.add_argument('-g', '--getPath', help='Get existing OMA Browser path if available', action='store_true', default=False)
    optional.add_argument('-k', '--keep', help='Keep downloaded OMA zip files (not recommended)', action='store_true', default=False)
    args = parser.parse_args()

    dataPath = args.outPath
    Path(dataPath).mkdir(parents = True, exist_ok = True)
    dataPath = os.path.abspath(dataPath)
    force = args.force
    save = args.savePath
    get = args.getPath
    keep = args.keep

    dccPath = os.path.realpath(__file__).replace('/prepareDcc.py','')

    # save existing dataPath to pathconfig file
    if save:
        checkOmaData(dataPath)
        with open(dccPath+'/pathconfig.txt','w') as config:
            config.write(dataPath)
            config.close()
        sys.exit('OMA Browser data can be found at %s. dcc2 is ready to run!' % dataPath)

    # get existing dataPath if available
    savedPath = getPath(dccPath)
    if get:
        if not savedPath == 'Pathconfig file not found':
            checkOmaData(savedPath)
            print('OMA Browser data can be found at\n%s' % savedPath)
        else:
            print('ERROR: OMA Browser data may not yet exist (pathconfig file not found)!')
        sys.exit()

    # download and process oma data
    if not savedPath == 'Pathconfig file not found':
        if not force:
            sys.exit('OMA Browser data can be found already at %s. dcc2 is ready to run!' % savedPath)
        else:
            if dataPath == savedPath:
                print('OMA Browser data found in %s will be deleted and re-downloaded! Enter to continue.' % savedPath)
                if query_yes_no(''):
                    try:
                        shutil.rmtree(savedPath)
                        Path(savedPath).mkdir(parents = True, exist_ok = True)
                        prepareDcc(savedPath, keep)
                    except:
                        sys.exit('Cannot remove %s. Please delete it by yourself and run prepareDcc again!' % savedPath)
            else:
                print('OMA Browser data found in different directory %s! Enter to continue re-downloading into %s.' % (savedPath, dataPath))
                if query_yes_no(''):
                    prepareDcc(dataPath, keep)
    else:
        prepareDcc(dataPath, keep)

    # check and save dataPath to pathconfig file
    checkOmaData(dataPath)
    with open(dccPath+'/pathconfig.txt','w') as config:
        config.write(os.path.abspath(dataPath))
        config.close()
    print("OMA Browser data can be found at %s. dcc2 is ready to run" % dataPath)

if __name__ == '__main__':
    main()

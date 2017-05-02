'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse
import re
import sys
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Get the ancestral consensus sequence from a hmm file')
    
    parser.add_argument('-i', dest='infilepath', metavar='<hmm_file_path>', type=argparse.FileType('rt'),
                   help='path to an hmm file')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()


def getBitscores(inputfile):
    bitscores = []
    for i, line in enumerate(inputfile.readlines()):
        if line.startswith('#'):
            continue
        matchlist = re.split("\s+", line.strip())
        if len(matchlist) >= 19:
            bitscores.append(matchlist[-14]) #not simlpy matchlist[5] because of possible spaces in description/name
        else:
            raise Exception("Error while reading bitscore from " + inputfile.name)
    return bitscores

def getSubjectName(inputname):
    if len(inputname.split('.')) == 2:
        return inputname.split('.')[0]

def getCutoff(bitscore_list, limit=0.9):
    arr = np.asarray(bitscore_list)
    min_score = float(min(arr))
    cutoff_val = min_score*limit
    return cutoff_val


if __name__ == '__main__':

    args = parse_args()
    
    try:
        inputfile = open(args.infilepath.name, 'r')

    except:
        print('IOError occured')
    
    seqname = getSubjectName(os.path.basename(args.infilepath.name))
    bitscores = getBitscores(inputfile)
    cutoff  = getCutoff(bitscores)
    sys.stdout.write('{}\t{}\n'.format(seqname,cutoff))
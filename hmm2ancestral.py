'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse
import re
import sys
from Bio.SeqUtils import GC

def parse_args():
    parser = argparse.ArgumentParser(description='Get the ancestral consensus sequence from a hmm file')
    
    parser.add_argument('-i', dest='infilepath', metavar='<hmm_file_path>', type=argparse.FileType('rt'),
                   help='path to an hmm file')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()


def getConsensus(contig):
    consensus_seq = []
    for i, line in enumerate(inputfile.readlines()):
        matchlist = re.split("\s+", line.strip())
        if len(matchlist) == 26:
            consensus_seq.append(matchlist[22].upper())
    return consensus_seq

def getSubjectName(inputname):
    if len(inputname.split('.')) == 2:
        return inputname.split('.')[0]

def list2fasta(seqname, seqlist, linelen = 60):
    numofbr = len(seqlist)//linelen
    seqlist.insert(0,'BLANK')
    for i in range(1, numofbr+1):
        seqlist.insert(i*linelen + i, '\n')
    del seqlist[0]
    fastaseq = ['>', seqname, '\n'] +  seqlist
    return ''.join(fastaseq).strip()


if __name__ == '__main__':

    args = parse_args()
    
    try:
        inputfile = open(args.infilepath.name, 'r')
#         if not os.path.basename(args.outfilepath.name) == "basename":
#             outputfile = open(args.outfilepath.name, 'w')
#         else:
#             outputfile = open(os.path.join(os.path.dirname(args.outfilepath.name),os.path.basename(args.infilepath.name) + '_consensus.faa'), 'w')
    except:
        print('IOError occured')
    
    seqname = getSubjectName(os.path.basename(args.infilepath.name))
    seqlist = getConsensus(inputfile)
    sys.stdout.write(list2fasta(seqname, seqlist, 60))
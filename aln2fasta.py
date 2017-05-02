'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse
import sys
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='build multifasta files from aln files')
    
    parser.add_argument('-i', dest='infilepath', metavar='<aln_file_path>', type=argparse.FileType('rt'),
                   help='path to an aln file')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()


def makegapless(fasta_sequences):
    for fasta in fasta_sequences:
        fasta.seq = fasta.seq.ungap('-')
    return fasta_sequences

def getSubjectName(inputname):
    if len(inputname.split('.')) == 2:
        return inputname.split('.')[0]

if __name__ == '__main__':

    args = parse_args()
    
    try:
        inputfile = open(args.infilepath.name, 'r')
        fasta_sequences = SeqIO.parse(inputfile,'fasta')
    except:
        print('IOError occured')
    
    seqname = getSubjectName(os.path.basename(args.infilepath.name))
#     import itertools
#     it1, it2 = itertools.tee(fasta_sequences, n=2)
    gapless_fasta_seqs = makegapless(list(fasta_sequences))
    SeqIO.write(gapless_fasta_seqs, sys.stdout, 'fasta')
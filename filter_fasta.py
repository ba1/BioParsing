'''
Created on Apr 17, 2016
#!/usr/bin/python3

#Go through a ncbi genome fna file and filter out the plasmids from the chromosome

@author: bardya
'''


import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Delete all unwanted sequences from a fasta')
    
    parser.add_argument('-i', dest='infilepath', metavar='<fasta_file_path>', type=argparse.FileType('rt'),
                   help='path to an fasta file')
    
    parser.add_argument('-o', dest='outfilepath', metavar='<fasta_file_path>', type=argparse.FileType('w'),
                   help='path to desired output fasta file')
    
    parser.add_argument('-include', dest='include', metavar='<str>', type=str, default=False, nargs='+',
                   help='sequence headers containing any of the specified strings will remain in the fasta')
    
    parser.add_argument('-exclude', dest='exclude',  metavar='<str>',type=str, default=False, nargs='+',
                   help='sequence headers containing any of the specified string will get deleted.')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()


def readfilter(seqdbfile, include=False, exclude=False):
    try:
        seqs = SeqIO.parse(seqdbfile, "fasta")
    except:
        print("ERROR not a fasta file")
        
    def include_test(listofstrings, seq):
        for s in listofstrings:
            if s in seq.description:
                return seq
        return False
            
    def exclude_test(listofstrings, seq):
        print(listofstrings, seq.description)
        for s in listofstrings:
            if s in seq.description:
                return False
        return seq
    
    if include:
        if type(include) is str:
            include=[include]
        seqlst = [include_test(include,s) for s in seqs]

    if exclude:    
        if type(exclude) is str:
            exclude=[exclude]
        seqlst = [exclude_test(exclude,s) for s in seqs]
    
    
    return list(filter(bool, seqlst))


def writefasta(outfile, seqlst):
    count = SeqIO.write(seqlst, outfile, "fasta")
    outfile.close()
    return count
    
if __name__ == '__main__':

    args = parse_args()
    
    try:
        inputfile = open(args.infilepath.name, 'r')
        outputfile = open(args.outfilepath.name, 'w')

    except:
        print('IOError occured')
    
    seqlst = readfilter(args.infilepath.name, include=args.include, exclude=args.exclude)
    writefasta(outputfile, seqlst)
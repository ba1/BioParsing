'''
Created on Nov 3, 2015

@author: bardya
'''
import os
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Given an integer "c" split a multifasta file into n/s files with s sequences.')
    
    parser.add_argument('-i', dest='infilepath', metavar='<fasta_file_path>', type=argparse.FileType('rt'),
                   help='path to an fasta file')
    
    parser.add_argument('-o', '--prefix', dest='outprefix', metavar='<fasta_files_outputpath-prefix>', type=str, default=os.getcwd(),
                   help='prefix (+ path) to desired output fasta files')
    
    parser.add_argument('-c', '--chunksize', dest='chunksize', metavar='<integer>', type=int, required=True,
                   help='Given an integer "c" splits the input multifasta file into n/c files, with c sequences')
    
    parser.add_argument('--version', action='version', version='0.21')
    
    return parser.parse_args()


def readfasta(seqdbfile):
    seqs = SeqIO.parse(seqdbfile, "fasta")
    seqlst = []

    for seq in seqs:
        seqlst.append(seq)

    return seqlst

def chunks(l, n):
    """Generator to yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def writefastas(seqlst, chunksize, outprefix, filesuffix):
    chunkgnr = chunks(seqlst,chunksize)
    
    if os.path.isdir(outprefix): #default file-naming if only a directory path as prefix
        outprefix = os.path.join(outprefix, "chunk")
        
    outdir = os.path.dirname(outprefix)
    
    i = 0
    for chunk in chunkgnr:
        filename = '{}_{}.{}'.format(os.path.basename(outprefix), str(i), filesuffix)
        outpath = os.path.join(outdir, filename)
        outputfile = open(outpath, 'w')
        count = SeqIO.write(chunk, outputfile, "fasta")
        print('[{}]\t{}\t{} Sequences'.format(i, outpath, count))
        outputfile.close()
        i+=1

if __name__ == '__main__':

    args = parse_args()

    try:
        inputfile = open(args.infilepath.name, 'r')
        filesuffix = args.infilepath.name.rsplit('.',1)[1]
    except:
        print('IOError occured, unable to read the input file')
    
    seqlst = readfasta(args.infilepath.name)
    if args.chunksize > len(seqlst) or args.chunksize < 1:
        raise("Chunksize c must be 0< c <= #Sequences in the input file.")
    writefastas(seqlst, args.chunksize, args.outprefix, filesuffix)
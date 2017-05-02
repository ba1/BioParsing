'''
Created on Nov 3, 2015

@author: bardya
'''
import os
import argparse
import random
from Bio import SeqIO 

def parse_args():
    parser = argparse.ArgumentParser(description='From a set of multi fasta files in the input-directory extract 1 or more representative sequences randomly for each')
    
    parser.add_argument('-i', dest='inpath', metavar='<path_to_dir_with_fastas>', type=str, required=True,
                   help='path to a directory containing fasta files. Only textfiles which end with .fa, .fasta, .fastq, .fq, .faa, .fna, or .fnn are considered.')
    
    parser.add_argument('-o', '--prefix', dest='outprefix', metavar='<fasta_files_outputpath-prefix>', type=str, required=True,
                   help='prefix (including path) to desired output fasta files')
    
    parser.add_argument('-s', '--samplesize', dest='samplesize', metavar='<integer>', type=int, default=1,
                   help='Number of sequences randomly drawn from each multi-fasta file')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()


def readfastas(inputpath):
    for fn in os.listdir(path=inputpath):
        if os.path.isfile(os.path.join(inputpath, fn)) and \
        fn.endswith(tuple(['.fa', '.faa', '.fna', '.fnn', '.fasta', '.fastq', '.fq', '.fastn'])):
            
            seqlst = list(SeqIO.parse(os.path.join(inputpath,fn), "fasta"))
            fn, suffix = fn.rsplit('.',1)
            yield seqlst, fn, suffix

def chooseRandom(seqlst, s=1):
    """Draw a random sample from a list of sequence records"""
    reprlist = random.sample(seqlst, s)
    return reprlist

def modify_fasta_header(reprlist, s, mode='description'):
    '''Given an extra string s add this depending on the mode specified to the id or description of each fasta header in list of records'''
    for seqrec in reprlist:
        if mode == 'id':
            seqrec.id += '__{}'.format(s)
        else:
            seqrec.description += ' ##{}##'.format(s)
    return reprlist

def writefastas(seqlst, outdir_prefix, fn, suffix= 'fasta', nametag="_sample",):
    #/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGFasta_representatives/'
    with open(outdir_prefix + fn + nametag + '.' + suffix , 'w') as outputh:
        SeqIO.write(seqlst, outputh, "fasta")

if __name__ == '__main__':

    args = parse_args()
    
    assert os.path.isdir(args.inpath)
    if not os.path.isdir(os.path.dirname(args.outprefix)):
        os.makedirs(os.path.dirname(args.outprefix))
    '/share/project/bardya/Enterobacteriaceae/OMA_prot/Results/HOGFasta'
    
    seqlst_iterator = readfastas(args.inpath)
    for seqlst, fn, suffix in seqlst_iterator:
        reprlst = chooseRandom(seqlst, s=1)
        modify_fasta_header(reprlst, fn, mode='description')
        writefastas(reprlst, args.outprefix, fn, suffix=suffix, nametag="_sample")
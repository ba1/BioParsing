'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse
import sys
from Bio import SeqIO
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Get the 95% confidence threshold length values (2stds) for group of ortholog gen sequences')
    
    parser.add_argument('-i', dest='infilepath', metavar='<fa_file_path>', type=argparse.FileType('rt'),
                   help='path to an multi fasta file')

    parser.add_argument('-ia', dest='ancestral_infilepath', metavar='<fa_file_path>', type=argparse.FileType('rt'), required=False,
                   help='''path to an multi fasta file containing corresponding consensus sequences, if not specified the mode of the length distribution
                        will be outputted at the second column instead of the consensus length''')
    
    parser.add_argument('-m','--mode', dest='mode', metavar="<sd|precalc>", type=str, default="precalc", choices=["sd","'sd'","precalc"],
                        help='''If mode is set to 'precalc' (default) then the third column of the output contains the subtracted sd and the fourth
                        column contains the added number. If mode 'sd' is specified, then the third contains the sd and the fourth contains the average.''')
    
    parser.add_argument('--version', action='version', version='0.11')
    
    return parser.parse_args()


def get95confval(len_distr, mode='precalc'):
    """Determine the mean and 2std of the length distribution of a group
    """
    arr = np.asarray(len_distr)
    m = np.mean(arr)
#     med = int(np.median(arr))
    std = np.std(arr)
    from scipy import stats
    mo = stats.mode(arr)[0][0]
    if mode=='precalc':
        return (mo, mo-2*std, mo+2*std)
    else:
        return (mo, 2*std, m)

def getSubjectName(inputname):
    """Determine first column subject name for output file
    """
    if len(inputname.split('.')) == 2:
        return inputname.split('.')[0]

def fasta_properties(bio_fasta_files):
    """From a multi-fasta determine the length distribution as a list
    """
    len_distr = []
    for fasta in bio_fasta_files:
        sequence = str(fasta.seq)
        len_distr.append(len(sequence.replace('-','')))
    return len_distr

def getItsConsSeqLength(seqname, ancestral):
    """From a multi-fasta containing the consensus sequences
       determine the length of the target groups cons. seq.
    """
    cons_fasta_seqs = SeqIO.index(ancestral,'fasta')
    return len(cons_fasta_seqs[seqname].seq)

if __name__ == '__main__':
    
    args = parse_args()
    
    try:
        inputfile = open(args.infilepath.name, 'r')
        fasta_sequences = SeqIO.parse(inputfile,'fasta')
#         if not os.path.basename(args.outfilepath.name) == "basename":
#             outputfile = open(args.outfilepath.name, 'w')
#         else:
#             outputfile = open(os.path.join(os.path.dirname(args.outfilepath.name),os.path.basename(args.infilepath.name) + '_consensus.faa'), 'w')
    except:
        print('IOError occured')
    
    seqname = getSubjectName(os.path.basename(args.infilepath.name))
    len_list = fasta_properties(fasta_sequences)
    res = get95confval(len_list, mode=args.mode)
    
    if args.ancestral_infilepath:
        cons_len = getItsConsSeqLength(seqname, args.ancestral_infilepath.name)
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(seqname, cons_len, res[1], res[2]))
    else:
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(seqname, *res))

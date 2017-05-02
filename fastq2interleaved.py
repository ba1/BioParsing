'''
Created on Dec 11, 2015

@author: bardya
'''
#This Python script requires Biopython 1.51 or later
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
#Setup variables (could parse command line args instead)


def parse_args():
    parser = argparse.ArgumentParser(description='Given two fastq files join them to a new file with the interleaved format.')
    
    parser.add_argument('-i1', '--forward_reads', dest='infilepath1', metavar='<fastq_file_path>', type=argparse.FileType('rt'),
                   help='path to an forward reads fastq file')
    
    parser.add_argument('-i2', '--reverse_reads',  dest='infilepath2', metavar='<fastq_file_path>',type=argparse.FileType('rt'),
                   help='path to an reverse reads fastq file')
    
    parser.add_argument('-o', '--out', dest='outfilepath', metavar='<fastq_file_output>', type=argparse.FileType('wt'),
                   help='prefix (+ path) to desired output fasta files')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()

if __name__ == '__main__':
    
    ##TEST
#     file_f = "SRR001666_1.fastq"
#     file_r = "SRR001666_2.fastq"
#     file_out = "SRR001666_interleaved.fastq"
    ##
    
    args = parse_args()
    
    handle = open(args.outfilepath.name, "w")
    count = 0

    f_iter = FastqGeneralIterator(open(args.infilepath1.name,"rU"))
    r_iter = FastqGeneralIterator(open(args.infilepath2.name,"rU"))
    for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter,r_iter):
        assert f_id == r_id
        count += 2
        #Write out both reads with "/1" and "/2" suffix on ID
        handle.write("@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n" \
                     % (f_id, f_seq, f_q, r_id, r_seq, r_q))
    handle.close()
    print("%i records written to %s" % (count, args.outfilepath.name))
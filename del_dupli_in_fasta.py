'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Delete all duplicate entries (header+sequence) in fasta. If only sequence identical, add "| duplicate" to header.')
    
    parser.add_argument('-i', dest='infilepath', metavar='<fasta_file_path>', type=argparse.FileType('rt'),
                   help='path to an fasta file')
    
    parser.add_argument('-o', dest='outfilepath', metavar='<fasta_file_path>', type=argparse.FileType('w'),
                   help='path to desired output fasta file')
    
    parser.add_argument('-m', dest='mode', metavar='<header|sequence>', type=str, choices=["header", "Header", "sequence", "Sequence"],
                   default="header", help='mode headers checks for "headers and then sequence". Mode sequence searches only for sequence duplicates')
    
    parser.add_argument('-k', dest='keep_flag', action="store_true",
                        help='with this options nothing gets deleted. Headers get count number attached to end of the line to make them unique.')
    
    parser.add_argument('-rn', dest='rename_flag', action="store_true",
                        help='with this options nothing gets deleted. Headers get replaced by an integer reflecting the count')
    
    
    
    parser.add_argument('--version', action='version', version='0.12')
    
    return parser.parse_args()


def readfasta(seqdbfile, keep_flag=False, rename_flag=False):
    from collections import Counter
    try:
        seqs = SeqIO.parse(seqdbfile, "fasta")
    except:
        seqs = SeqIO.parse(seqdbfile, "clustal")
    
    seqlst = []
    dupcount = 0
    modcount = 0
    

    for seq in seqs:
        currIDlst = [e.id for e in seqlst]
        if rename_flag:
            seq.id = ">" + str(Counter(seqlst)[str(seq.id)] + 1)
            modcount += 1
            continue
        if seq.id in currIDlst:
            ind = currIDlst.index(seq.id)
            if keep_flag:
                seq.id = str(seq.id) + "_" + str(Counter(seqlst)[str(seq.id)] + 1)
                modcount += 1
                continue
            if seqlst[ind].seq == seq.seq:
                dupcount += 1
                continue
            else:
                seq.id = str(seq.id) + "_" + str(Counter(seqlst)[str(seq.id)] + 1)
                modcount += 1
        seqlst.append(seq)

    stats_dict = {"delentries":dupcount, "numofseqs":len(seqlst), "modentries":modcount}
    return seqlst, stats_dict


def printStats(stats_dict):
    outp = """
#Entries remaining in output:\t{numofseqs}
#Entries deleted:\t{delentries}
#Headers modified:\t{modentries}
    """.format(**stats_dict)
    print(outp)

def writefasta(outfile, seqlst):
    count = SeqIO.write(seqlst, outfile, "fasta")
    outfile.close()

if __name__ == '__main__':

    args = parse_args()
    
    try:
        inputfile = open(args.infilepath.name, 'r')
        outputfile = open(args.outfilepath.name, 'w')
#         if not os.path.basename(args.outfilepath.name) == "basename":
#             outputfile = open(args.outfilepath.name, 'w')
#         else:
#             outputfile = open(os.path.join(os.path.dirname(args.outfilepath.name),os.path.basename(args.infilepath.name) + '_consensus.faa'), 'w')
    except:
        print('IOError occured')
    
    seqlst, stats_dict = readfasta(args.infilepath.name, keep_flag=args.keep_flag, rename_flag=args.rename_flag)
    printStats(stats_dict)
    writefasta(outputfile, seqlst)
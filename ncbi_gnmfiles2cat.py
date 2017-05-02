'''
Created on Oct 20, 2015
This script takes the genome files in a specific directory and concatenates them if they
have the same file ending. The end file therefore contains all sequence information of 
all plasmids/replicons/contigs in a single file outputted in the mother directory.

@author: bardya
'''
import os
import argparse
import subprocess
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Concatenate genome sequences of different replicons into one for a specified directory')
    
    parser.add_argument('-i', dest='dirpath', metavar='<path_to_directory>',
                   help='path to a directory containing the files to concatenate - dirname will become the output file prefix')
    
    parser.add_argument('-s', '--suffix', dest='suffix', metavar='<suffix>',
                   help='suffix specified without the dot e.g. "faa", or "ffn"')
    
    parser.add_argument('-m', '--m', dest='mode', metavar='<mode>', choices=['single','batch'], default='batch',
                   help='if set the specified dir will be regarded mother dir')
    
    parser.add_argument('-sr', '--subdirregex', dest='subdir_regex', metavar='<subdir_regex>', default='.*',
                   help='specify the regex to match the subdir name, all by default')
    
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args()

def concatenate(suffix, genomesdir, subdir_regex, depth):
    if depth == 1:
        for path, subdirs, files in os.walk(genomesdir):
            for subdir in subdirs:
                    if not re.match(subdir_regex, subdir):
                        continue
                    outfilepath = os.path.join(genomesdir, subdir + "." + suffix)
                    outfile = open(outfilepath, 'w')
                    cmd = ['cat', os.path.join(genomesdir, subdir, '*' + suffix)]
                    cat = subprocess.call(' '.join(cmd), shell=True, stdout=outfile, stderr=subprocess.PIPE)
                    #cat.wait()
                    #cat.communicate()[0]
                    #print(outfile.name, outfilepath)
                    outfile.close()
    else:
        if not re.match(subdir_regex, genomesdir):
            import sys
            print("Specified directory does not match defined regex. Quitting script.", file=sys.stderr) 
        outfilepath = os.path.join(genomesdir, os.path.basename(path) + "." + suffix)
        outfile = open(outfilepath, 'w')
        cat = subprocess.Popen(['cat', os.path.join(genomesdir, genomesdir) + '*'+ suffix], stdout=outfile, stderr=subprocess.PIPE)
        cat.communicate()[0]
        outfile.close()

if __name__ == '__main__':

    args = parse_args()
    
    if args.dirpath.endswith('/'):          #the specified path shouldn't contain a slash at the end
        args.dirpath = args.dirpath[:-1]
    
    if not os.path.isdir(args.dirpath):
        raise Exception("Specified path is not a directory.")
    depth = 1
    if args.mode == 'single':
        depth = 0
    concatenate = concatenate(args.suffix, args.dirpath, args.subdir_regex, depth)
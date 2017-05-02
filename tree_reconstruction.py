'''
Created on Dec 14, 2015

@author: bardya
'''

import os,glob
import argparse

class Range:
    """A simple helper class to allow for a range restriction of floats in argparse
    by extending the 'choices' functionality"""
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return '{0}-{1}'.format(self.start, self.end)

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate a species tree')
    
    parser.add_argument('-o', '--outpath', dest='outpath', metavar='<path/to/output>', type=str, default="",
                   help='path to desired output directory')
    
    parser.add_argument('-a', '--msa-alignment', dest='align', action="store_true",
                        help='calculate ms-alignments for all sequences in subdir /<wdpath>/fa')
    
    parser.add_argument('-cat', '--concat-alignments', dest='cat', action="store_true",
                        help='Generate supermatrix by concatenating with perl script')
    
    parser.add_argument('-p', '--phylip-conversion', dest='phylip', action="store_true",
                        help='convert all alignments in aln subdirectory to phylip format')
    
    parser.add_argument('-t', '--tree-recons', dest='tree', action="store_true",
                        help='construct gene trees with RaXML')
    
    parser.add_argument('-hmm', '--hmm-models', dest='hmm', action="store_true",
                        help='Generate hmm with hmmer3')
    
    parser.add_argument('-hmmprof', '--hmm_profiles', dest='hmmprof', action="store_true",
                        help='Generate hmm profiles with augustus')
    
    parser.add_argument('-wd', '--wdpath', dest='wdpath', metavar='<string>', type=str, required=True,
                   help='Path to the working directory')
    
    parser.add_argument('-dg', '--degap', dest='degap', metavar='<string>', type=float, choices=[Range(0.0, 100.0)],
                   help='Degapping of concatenated alignment')
    
    parser.add_argument('--no_multiprocessing', dest='no_multiprocessing', action="store_true",
                         help='If this set, no multiprocessing will be done ')
    
    parser.add_argument('--version', action='version', version='0.2')
    
    return parser.parse_args()




def buildMSA(inpath, outpath, f=False):
    """build MSA with MAFFT-linsi
    """
    fileset = glob.glob(os.path.join(inpath, '*.fa'))
    if not f:
        parallelizer(fileset, buildMSA, inpath, outpath)
    else:
        if not os.path.isfile(os.path.join(outpath, os.path.basename(f).replace('.fa','.'+align_format))):
            os.system("/share/applications/mafft/mafft-7.213/bin/linsi --thread -1 --quiet {0} > {1}".format(
                f, os.path.join(outpath, os.path.basename(f).replace('.fa','.'+align_format))))


def convert2Phylip(outpath, outpath1, f=False):
    """convert alignment f to phylip f, uses clustalw, 
    therefore expexts aln file format
    """
    fileset =  glob.glob(os.path.join(outpath, '*.aln'))
    if not f:
        parallelizer(fileset, convert2Phylip, inpath, outpath)
    else:
        os.system('/usr/bin/clustalw -convert -output=Phylip -infile={0} -outfile={1}'.format(
                f, os.path.join(outpath1, os.path.basename(f).replace('.aln', '.phy'))))


def concatAlignment(inpath,outfilepath):
    """concatenate fastas in directory to supermatrix
    """
    #With Perl Script
    #os.system("time perl /home/bardya/usr/scripts/perl_scripts/concat_alignments_dmp.pl -in {} -out {}".format(
    #            inpath, outfilepath))

    #With FastConCat
    #os.system("cd {} && perl /home/holger/src/FASconCAT_v1.0.pl -s -i -p".format(inpath))
    #Not in phylip format
    os.system("cd {} && perl /home/bardya/usr/bin/FASconCAT/FASconCAT-G_v1.02.pl -s".format(inpath))
    os.system( "mv {} ../ && mv {} ../".format(os.path.join(inpath,"FcC_supermatrix.fas"),os.path.join(inpath,"FcC_info.xls")));

def degapConcatedAlignment(infilepath, t=0):
    """degap the concatenated alignment in the supermatrix file. Exclude columns with less than fraction t rows
    """
    os.system("time perl /home/bardya/usr/scripts/perl_scripts/degapper.pl -in {} -limit {}".format(
                infilepath, t))
    #os.system("time perl /home/bardya/usr/scripts/perl_scripts/RBF_FaPhylip.pl -fileend .proc -indir {} -outdir {}".format(
    #            os.path.dirname(infilepath)))
    
def prepareTrees(outpath, cwdpath):
    pass
    # prepare RXaml Trees
#     if args.phylip:
#         for f in glob.glob(os.path.join(outpath, '*.phy')):
#             os.system('/share/applications/standard-RAxML/bin/raxmlHPC -f a -# 100 -x 12345 -m PROTGAMMAAUTO -s {0} -p 12345 -n {1} -w {2}'.format(
#             f, os.path.basename(f).replace('.phy',''), cwdpath))
# //    elif args.cat:
#         os.system('/share/applications/standard-RAxML/bin/raxmlHPC -f a -# 100 -x 12345 -m PROTGAMMAAUTO -s {0} -p 12345 -n {1} -w {2}'.format(
#             f, os.path.basename(f).replace('.phy',''), cwdpath))

def generateHmms(outpath, outpath1, f=False):
    """build hmm models with HMMER3
    """
    fileset = glob.glob(os.path.join(outpath, '*.'+align_format))
    if not f:
        parallelizer(fileset, generateHmms, outpath, outpath1)
    else:
        if not os.path.isfile(os.path.join(outpath1,os.path.basename(f).replace('.'+align_format,'.hmm'))):
            os.system('/usr/bin/hmmbuild {1} {0}'.format(
                f, os.path.join(outpath1,os.path.basename(f).replace('.'+align_format,'.hmm'))))


def generateHmmProfiles(outpath, outpath2, f=False):
    """generate profile files with augustus from msa files
    """
    fileset = glob.glob(os.path.join(outpath, '*' + align_format))
    if not f:
        parallelizer(fileset, generateHmmProfiles, outpath, outpath2)
    else:
        if not os.path.isfile(os.path.join(outpath2,os.path.basename(f).replace('.'+align_format, '.prfl'))):
            os.system('perl /home/bardya/usr/bin/augustus-3.2/scripts/msa2prfl.pl {0} > {1}'.format(
                f, os.path.join(outpath2,os.path.basename(f).replace('.'+align_format, '.prfl'))))


def parallelizer(fileset, my_function, path1, path2):
    """Is called by the function actually recalled to do the job (saves one function)"""
    
    import multiprocessing
    #processes = []
    #for f in fileset:
    #    processes.append(multiprocessing.Process(target=my_function, args=(path1,path2), kwargs={"f":f}))
    global multip
    pool = multiprocessing.Pool(processes=multip) #this forces the value of cpu_count() as number of simultaneous workers
    for filepath in fileset:
        pool.apply(my_function, (path1, path2), dict(f=filepath))

                  
if __name__ == '__main__':
    args = parse_args()
    #args = parser.parse_args(['',''])    #for testing purposes
    
    os.chdir(args.wdpath) #has to be specified otherwise error
    cwdpath = os.getcwd()
    
    if not args.outpath:
        if not args.phylip:
            align_format='fas'
        else:
            align_format='aln'
        outpath = os.path.join(cwdpath, align_format)
        outpath1 = os.path.join(cwdpath, 'hmm')
        outpath2 = os.path.join(cwdpath, 'prfl')
        
    if not os.path.exists(outpath):
        os.makedirs(outpath, 0o777)
    if not os.path.exists(outpath1):
        os.makedirs(outpath1, 0o777)
    if not os.path.exists(outpath2):
        os.makedirs(outpath2, 0o777)
    
    inpath = os.path.join(cwdpath,'fa')
    
    if args.no_multiprocessing:
        multip = int(args.no_multiprocessing) #used global
    else:
        multip = None
    
    if args.align:
        buildMSA(inpath, outpath)
    if args.hmm:
        generateHmms(outpath, outpath1)
    if args.hmmprof:
        generateHmmProfiles(outpath, outpath2)
    if args.phylip:
        convert2Phylip(outpath, outpath)
    if args.cat:
        concatAlignment(outpath, os.path.join(cwdpath, "FcC_supermatrix.fas"))
    if args.degap:
        degapConcatedAlignment(os.path.join(cwdpath, "FcC_supermatrix.fas"), t=args.degap)
    if args.tree:
        prepareTrees(outpath, cwdpath)
    

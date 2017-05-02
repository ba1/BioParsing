'''
Created on Dec 14, 2015

@author: bardya
'''
import argparse
def parse_args():
    parser = argparse.ArgumentParser(description='Given a set of filters get ncbi assembly data and ftp addresses')
    
    parser.add_argument('-o', '--outputfile', dest='outputfile', metavar='<path/to/output>', type=str,
                        help='path to desired output file')
    
    parser.add_argument('-d', '--field-delimiter', dest='delimiter', metavar='<character|string>', type=str, default='\t',
                        help='specify the field delimiter in the assembly summary file and the output')
    
    parser.add_argument('--usetable', dest='download_list_given', metavar='</path/to/assembly_list>', type=str,
                        help='you can use a prepared assembly list as well, all parameters for filtering will be ignored')
    
    parser.add_argument('-s', '--subject', dest='subject', metavar='<string>', type=str, required=False,
                        help='Expects a substring for column 8 of the assembly_summary.txt file, namely the organisms name')

    parser.add_argument('-l', '--latest', dest='latest', action="store_true",
                        help='filter for ftp addresses of the latest assembly only')
    
    parser.add_argument('-c', '--complete', dest='complete', action="store_true",
                        help='filter for assemblies with status "complete genome" only')
    
    parser.add_argument('-chr', '--chromosome', dest='chromosome', action="store_true",
                        help='filter for assemblies with status "chromosome" or higher')

    parser.add_argument('-sc', '--scaffold', dest='scaffold', action="store_true",
                        help='filter for assemblies with status "scaffold" or higher')
    
    parser.add_argument('-co', '--contig', dest='contig', action="store_true",
                        help='filter for assemblies with status "contig" or higher')
    
    parser.add_argument('-fm', '--filematch', dest='filematch', type=str, 
                        default='*.gz',
                        help='specifies a UNIX-like regular expression, only matching filenames are downloaded')

    parser.add_argument('-u','--url', dest='url', metavar='<URL>', type=str,
                        default="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", 
                        help='Specifies the download URL of the assembly summary path')
    
    parser.add_argument('-p','--pathsonly', dest='pathsonly', action="store_true",
                        help='if set prints only the ftp path')
    
    parser.add_argument('-v','--verbose', dest='verbose', action="store_true",
                        help='Print information on download process to standard out')
    
    parser.add_argument('-odl', '--download_prefix', dest='destination', metavar='</path/prefix>', type=str,
                        default='',
                        help='path and prefix to desired output directory for ftp downloaded assembly files')
    
    parser.add_argument('--version', action='version', version='0.4')
    
    return parser, parser.parse_args()


#awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths

def get_assembly_summary_file(URL):
    import urllib.request
    filename, headers = urllib.request.urlretrieve(URL)
    assembly_summary_txt = open(filename)
    return assembly_summary_txt

def filterLinesForFields(assembly_summary_txt, delimiter, filterlist=None):
    filtered_lines = []
    headers = []
    
    for __, line in enumerate(assembly_summary_txt.readlines()):
        if line.startswith('#'):
            headers.append(line)
            continue
        fields = line.split(delimiter)
        if filterlist:
            filtered_line = filter(fields, delimiter, filterlist)
        else:
            filtered_line = delimiter.join(fields)
        if filtered_line:
            filtered_lines.append(filtered_line)
    
    return headers, filtered_lines       

def filter(fields, delimiter, filterlist):
    """two lists are given and according to positions they are zipped into tuples, 
    if a value exists in the filter list (can also be string or a list of strings)
    than it has to be in it"""
    comparison_pairs = list(zip(fields, filterlist))
    for cp in comparison_pairs:
        if cp[1]:
            if isinstance(cp[1], str):
                if not cp[1] in cp[0]:          #here insert regular expression pattern comparison
                    return False
            elif isinstance(cp[1], list):
                if not cp[0] in cp[1]:
                    return False
    
    return delimiter.join(fields)
    
def create_filterlist(subject, latest, complete, chromosome, scaffold, contig):
    filterlist = [False] * 19
    if subject:
        filterlist[7] = subject
    if latest:
        filterlist[10] = "latest"
    if complete:
        filterlist[11] = "Complete Genome"
    if chromosome:
        filterlist[11] = ["Chromosome", "Complete Genome"]
    if scaffold:
        filterlist[11] = ["Scaffold", "Chromosome", "Complete Genome"]
    if contig:
        filterlist[11] = ["Scaffold", "Complete Genome", "Chromosome", "Contig"]
    
    return filterlist

def generateOutput(outputfile, pathsonly, delimiter, headers, filtered_lines, download_list_given):
    
    if not headers:
        print("No header line found.")
    elif len(headers)>1:
        header = headers[0]
        print("Several potential header lines found, assuming first line to be correct.")
    else:
        header = headers[0]
    
    filtered_lines.insert(0,header)   
    #if pathsonly:
        #filtered_lines = [line.split(delimiter)[-1] for line in filtered_lines]
    
    if outputfile:
        with open(outputfile, 'w') as outfile:
            outfile.writelines(filtered_lines)
    elif not download_list_given:
        for line in filtered_lines:
            print(line.rstrip())

def ftpDownload(filtered_lines, delimiter, destination, filematch, verbose=False):
    """for all filtered_lines in the assembly list use the URL in field 19 to calculate
    the results"""
    import ftplib
    ftpurlist = [line.rstrip().split(delimiter)[-1] for line in filtered_lines]

    try:
        serveraddress = 'ftp.ncbi.nlm.nih.gov'
        ftp = ftplib.FTP(serveraddress)
        ftp.login()
        
        #filematch = '*.gz' #'*.gz'
        addr_prefix = "ftp://ftp.ncbi.nlm.nih.gov"
        
        for ftpurl in ftpurlist:
            if verbose:
                print("Downloading from {} to {}".format(ftpurl, destination))
                getFiles(ftp, ftpurl, addr_prefix, filematch, destination)
                print("Done")
            else:
                getFiles(ftp, ftpurl, addr_prefix, filematch, destination)
        ftp.quit()
    except:
        print("FTP Download Error. Check addresses and destination paths -->", ftpurl)
        import traceback
        tb = traceback.format_exc()
        print(tb)
    
#wget -r -b -c -nH -N ftp://ftp.patricbrc.org/patric2/patric3/genomes/* -A '.faa' --cut-dirs=3
  
def getFiles(ftp, ftpurl, addr_prefix, filematch, destination):
    import os
    directory = ftpurl.replace(addr_prefix,'',1)
    subjectname = os.path.basename(ftpurl)
    ftp.cwd(directory)
    if not os.path.exists(os.path.join(destination,subjectname)):
        os.makedirs(os.path.join(destination,subjectname))
    else:
        return
    for filename in ftp.nlst(filematch):
        fhandle = open(os.path.join(destination, subjectname, filename), 'wb')
        ftp.retrbinary('RETR ' + filename, fhandle.write)
        fhandle.close()
                     
if __name__ == '__main__':
    parser, args = parse_args()
    #args = parser.parse_args(['-s', 'Acinetobacter', '-c', '-l', 
    #                          '-odl', '/share/project/bardya/genomes/Acinetobacter', 
    #                          '-o', '/share/project/bardya/genomes/Acinetobacter/filtered_assembly_summary.txt'])
    if not args.download_list_given:
        assembly_summary_txt = get_assembly_summary_file(args.url)
        filterlist = create_filterlist(args.subject, args.latest, args.complete, args.chromosome, args.scaffold, args.contig)
        headers, filtered_lines = filterLinesForFields(assembly_summary_txt, args.delimiter, filterlist=filterlist)
    else:
        assembly_summary_txt = open(args.download_list_given, 'r')
        headers, filtered_lines = filterLinesForFields(assembly_summary_txt, args.delimiter, filterlist=None)
    #Toggle here if download or not
    if not args.pathsonly:
        ftpDownload(filtered_lines, args.delimiter, args.destination, args.filematch, verbose = args.verbose)
    generateOutput(args.outputfile, args.pathsonly, args.delimiter, headers, filtered_lines, args.download_list_given)
    

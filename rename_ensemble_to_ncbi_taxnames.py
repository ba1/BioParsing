'''
Created on Jul 17, 2015

@author: bardya
'''
import re


mapdict = {}

with open("/home/bardya/usr/bin/megan/class/resources/files/acinetobacter_species_EnsemblBacteria.txt", 'r') as mapfile:
    mapdict = dict(line.strip().split("\t", 1)[::-1] for line in mapfile) #returns the 2-element iterable to the dict function

with open('/home/bardya/usr/bin/ncbi-blast-2.2.30+/db/acinetobacter_ensemble_prot.faa', 'r') as infile2, open('/home/bardya/usr/bin/ncbi-blast-2.2.30+/db/acinetobacter_ensemble_ncbitax_prot.faa', 'w') as outfile:
    for i,line in enumerate(infile2):
        for word in re.findall(r"\w+", line):
            if word in mapdict:                
                line = line.replace(word, mapdict[word])
        outfile.write(line)
#So we first split our input line into individual words, then we perform a fast lookup of the words 
#within the keys of our dictionary. Thereby, only if this word exists, the replace function for this line is called.       
        
#Some explanations: The re module is optimzed for string searching. re.findall returns all non-overlapping matches of <pattern> in <string>, as a list of strings, scans left-to-right.
#Alternatively: re.sub
#>>> re.sub(r'\sAND\s', ' & ', 'Baked Beans And Spam', flags=re.IGNORECASE)
#'Baked Beans & Spam'

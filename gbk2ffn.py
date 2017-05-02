'''
Created on Dec 15, 2015

@author: bardya
'''

from Bio import SeqIO
gb_file = "/share/project/bardya/Acinetobacter/ncbi_ftp/GCF_001077655.1_ASM107765v1/GCF_001077655.1_ASM107765v1_genomic.gbff"
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
    # now do something with the record
    print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
    print(repr(gb_record.seq))
    
gb_feature = gb_record.features[26]
print(gb_feature)
x = gb_feature.extract(gb_record.seq)
print(x)
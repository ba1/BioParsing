'''
Created on Mar 9, 2017

@author: bardya
'''

from Bio import SeqIO
gbk_filename = "/home/bardya/usr/data/Oleivorans-like/AC-SGAC_2117.gbk_mod"
faa_filename = "/home/bardya/usr/data/Oleivorans-like/AC-SGAC_2117.faa"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank"):
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type == "CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">{} from {}\n{}\n".format(
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
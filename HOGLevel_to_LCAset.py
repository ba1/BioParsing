'''
Created on Dec 20, 2016

@author: bardya
'''
# HOG_herit = [0,1,2]
#  
# hogset = set()
#  
# for i in HOG_herit:
#     fh = open("/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGLevel_Gains/{}".format(i), 'r')
#     hogset |= set([line.strip() for line in fh if line.strip()])
#     fh.close()
#  
#     fh2 = open("/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGLevel_Losses/{}".format(i), 'r')
#     hogset -= set([line.strip() for line in fh2 if line.strip()])
#     fh2.close()
#  
# LCA_set = []
# for i in hogset:
#     i_sample = i.replace('.fa', '_sample.fa')
#     with open('/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGFasta_random_representative/{}'.format(i_sample), 'r') as fh3:
#         for line in fh3:
#             if line.startswith('>'):
#                 protein_id = line.split(' ', 1)[0][1:]
#                 LCA_set.append(protein_id)
#  
# with open('/share/project/bardya/Enterobacteriaceae/OMA_prot/HOG_LCA/{}_LCA_set.txt'.format(HOG_herit[-1]), 'w') as LCA_seth:
#     for protein_id in LCA_set:
#         LCA_seth.write(protein_id + '\n')


HOG_herit = [0]
  
  
hogset = set()
fh = open("/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGLevel_Gains/{}".format(HOG_herit[0]), 'r')
hogset |= set([line.strip() for line in fh if line.strip()])
fh.close()
  
LCA_set = []
  
for i in hogset:
    i_sample = i.replace('.fa', '_sample.fa')
    with open('/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGFasta_random_representative/{}'.format(i_sample), 'r') as fh3:
        for line in fh3:
            if line.startswith('>'):
                protein_id = line.split(' ', 1)[0][1:]
        LCA_set.append(protein_id)
  
with open('/share/project/bardya/Enterobacteriaceae/OMA_prot/{}_Gains.txt'.format(HOG_herit[0]), 'w') as LCA_seth:
    for protein_id in LCA_set:
        LCA_seth.write(protein_id + '\n')



# HOG_herit = [20]
# 
# 
# hogset = set()
# fh = open("/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGLevel_Losses/{}".format(HOG_herit[0]), 'r')
# hogset |= set([line.strip() for line in fh if line.strip()])
# fh.close()
# 
# LCA_set = []
# 
# for i in hogset:
#     i_sample = i.replace('.fa', '_sample.fa')
#     with open('/share/project/bardya/Enterobacteriaceae/OMA_prot/HOGFasta_random_representative/{}'.format(i_sample), 'r') as fh3:
#         for line in fh3:
#             if line.startswith('>'):
#                 protein_id = line.split(' ', 1)[0][1:]
#         LCA_set.append(protein_id)
# 
# with open('/share/project/bardya/Enterobacteriaceae/OMA_prot/{}_Losses.txt'.format(HOG_herit[0]), 'w') as LCA_seth:
#     for protein_id in LCA_set:
#         LCA_seth.write(protein_id + '\n')

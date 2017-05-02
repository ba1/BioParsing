'''
Created on Dec 12, 2016

@author: bardya
'''

MODE="BRITE"

linelst=[]

khan = open("/share/project/bardya/dbs/kegg/ko02044.keg", 'r')

for line in khan:
    if MODE == "PATHWAY":
    
        if line.startswith(('A', 'C', 'D')):
            if line.startswith('A'):
                super = line.split('>')[1].split('<')[0]
                continue
            if line.startswith('C'):
                #print(line)
                meta = line[1:].strip()
                #print('meta', meta)
                continue
            else:
                line = line[1:].strip()
                linelst.append('{}\t{}\t{}\t{}\n'.format(line[0], line[1], meta, super))
                        
    else:
        
        if line.startswith(('A', 'B', 'C')):
            if line.startswith('A'):
                super = line.split('>')[1].split('<')[0]
                continue
            if line.startswith('B'):
                #print(line)
                meta = line[1:].strip()
                #print('meta', meta)
                continue
            else:
                line = line[1:].strip()
                line = line.rsplit("  ",1)
                linelst.append('{}\t{}\t{}\t{}\n'.format(line[0], line[1], meta, super))
        
        
        
with open("/share/project/bardya/dbs/kegg/keggK2brite.map", "w") as hn:
    for line in linelst:
        hn.write(line)

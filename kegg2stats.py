'''
Created on Dec 14, 2016

Purpose: Given a list of keggKO Results from "Detail Page". Create a map which contains 
further information besides protein ID (e.g. HOG membership)

Purpose2: For individual lists containing this secondary information extract all its 
genes by id and extract all its annotated genes and pathways

@author: bardya
'''

import os
import subprocess
import csv
import sqlite3
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='generate fasta files of all groups, naming scheme: OG<#species>_<groupID>.fa')

    parser.add_argument('-omawd', dest='oma_workpath', metavar='<oma_working_directory_path>', required= True,
                   help='path to the OMA working directory')

    parser.add_argument('-omaout', dest='oma_output_dirname', metavar='<oma_output_dirname>', default = "Results",
                   help='base directory of the OMA output')

    parser.add_argument('-m', '--mode', dest='mode', metavar='<hog|og|hogaware>', type=str, 
                        choices=["hog", "hogs", "HOG", "HOGs", "og", "ogs", "OG", "OGS", "HOGAware","HOGAWARE","hogaware","hogAware"],
                   default="ogs", help='based on selected mode parse the OrthologousGroups.orthoxml or the HierarchicalGroups.orthoxml file located in OMA output directory.')

    parser.add_argument('-f', '--force', dest='force_flag', metavar='<force overwrite flag>', action='store_const', const=1, default=0,
                   help='if set, the output in the directory "./Bins" will be overwritten')

    parser.add_argument('--no-stats', dest='stats_flag', metavar='<stats_to_stdout_flag>', action='store_const', const=0, default=1,
                   help='if set, script does not give out a statistical overview to stdout')
    
    parser.add_argument('--no-accessory', dest='accesory_flag', metavar='<produce accessory genomes flag>', action='store_const', const=1, default=0,
                   help='if set, script gives out the accessory genomes into the directory "./Accessory" relativ to omawd')

    parser.add_argument('-t', '--speciestree', dest='nwcktree', metavar='<path/to/tree.file>', type=argparse.FileType('rt'),
                   help='path to a file containing the species tree in string representation')
    

    parser.add_argument('--version', action='version', version='0.1')

    return parser.parse_args()


def clearCheckPath(outdir, force=0):
    if os.path.isdir(outdir):
        if force:
            import shutil
            shutil.rmtree(outdir, ignore_errors=True)
        else:
            raise IOError("Output Directory already exiting. Specify '-f' option to force overwrite")

    os.makedirs(outdir)


def clearCheckFile(filepath, force=0):
    if os.path.isfile(filepath):
        if force:
            import shutil
            shutil.rmtree(filepath, ignore_errors=True)
        else:
            raise IOError("Output File already exiting. Specify '-f' option to force overwrite")

def createSqlTable(con, filepath, tablename, columnnames, sep='\t', mode='fail', primary_key=()):
    '''creates a database table '''
    
    try:
        df = pd.read_csv(filepath, sep=sep, names=columnnames, index_col=False)
        df.to_sql(tablename, con, if_exists=mode, index=False)
    
    except ValueError as e:
        print(tablename +  ' already exists in the database.')
    except Exception as e:
        print('Problem with creation of ' + tablename)
        print(''' Error Message:\n'''.format(e))
    #if primary_key:
    #    con.execute('ALTER TABLE {} ADD PRIMARY KEY {};'.format(tablename, tuple(primary_key)))
    
    con.commit()    
    
#def createProtein2Ktable():
# def concatAllKeggResults(listofcsvs, keggjobid2organism):
#     '''Concatenate textfile outputs of KEGG Blastkoala analysis to a single large file with structure specied-id \t kegg result line
#     '''
#     subprocess.call('''for i in *.csv; do id=$(echo $i | sed 's/\.csv$//'); while read line; do echo "$id    $line" >> protein2keggK.map; done < $i; done''', shell=True)
#     
#     subprocess.call('''while read line; do id="$(echo -e "$line" | cut -f1)"; id="$(echo -e "$id" | sed 's/\./__/')"; name="$(echo -e "$line" | cut -f2,3)"; name="$(echo -e "$name" | sed 's/\t/ /')"; sed -i "s/$id/$name/" protein2keggK.map; done < keggoutput2organism.map''', shell=True)
# 
#     subprocess.call('''cut -f1,2,3,4 protein2keggK.map > protein2keggK_no_seconday.map''', shell=True)
# 
# 
# def addAdditionalInfo(additional_info_map):
#     ''' Adds an additional info column to each line of the protein2keggK_no_seconday.map 
#     '''
#     
#     subprocess.call('''for i in HOG*; do ids=$(grep ">" $i | cut -d' ' -f1 | cut -c2-25); while read id; do res=$(grep -m1 "$id" ../../kegg/protein2keggK_no_seconday.map); K=$(echo "$res" | cut -f3); if [[ ! -z "$K" ]]; then categ=$(grep -m1 $K /share/project/bardya/dbs/kegg/keggK2pathway.map); else categ=$(echo -e '\t\t\t'); fi; echo -e "$i"'\t'"$id"'\t'"$categ"; done < "$ids"; done > results_kegg_all.map''', shell=True)
#
#
# def OrganismProtein2HOGmap():
#     suprocess.call('''for i in *.fa; do headers=$(grep ">" $i); while read -r line; do id=$(echo "$line" | cut -d ' ' -f1 | cut -c2-50 |  tr -d '[:space:]'); organism=$(echo "$line" | rev | cut -d' ' -f1 | rev | tr -d '\[\]' | tr -d '[:space:]'); echo "$organism        $id     $i"; done <<< "$headers"; done > OrganismProtein2HOG.map''')     
#    
    
if __name__ == '__main__':
    
    KEGGKO2PATHWAY_MAP = '/share/project/bardya/dbs/kegg/keggK2pathway.map'
#created with kegg_hierarch2map.py --> Each line contains tab separated
#KO    Gene name/ description/ EC    Pathway id + Pathway name    Category
#K00844  HK; hexokinase [EC:2.7.1.1]     01200 Carbon metabolism [PATH:ko01200]  Metabolism

    PROTEIN2KEGG_MAP = '/share/project/bardya/Enterobacteriaceae/kegg/protein2keggK0.map'
#-->Each line contains tab separated information on
#OrganismID    ProteinID    KEGG Results (tab separated detail-output of BlastKoala)
#GCF_000736695__1_XBKBD_3526_PRJEB4325_protein   WP_002211347.1  K02518  infA; translation initiation factor IF-1        73


    ORGANISMPROTEIN2HOG_MAP = '/share/project/bardya/Enterobacteriaceae/OMA_prot/Results/HOGFasta/OrganismProtein2HOG.map2'
#created with #suprocess.call('''for i in *.fa; do headers=$(grep ">" $i); while read -r line; do id=$(echo "$line" | cut -d ' ' -f1 | cut -c2-50 |  tr -d '[:space:]'); organism=$(echo "$line" | rev | cut -d' ' -f1 | rev | tr -d '\[\]' | tr -d '[:space:]'); echo "$organism        $id     $i"; done <<< "$headers"; done > OrganismProtein2HOG.map''')     
#-->Each line contains tab separated information on
#OrganismID    ProteinID    HOG_ID
#GCF_000252955__1_ASM25295v1_protein    WP_041573901.1    HOG10000.fa
    
    
    LOSSGAINHOG2Level = '/share/project/bardya/Enterobacteriaceae/OMA_prot/level2hoglossgain.map'
#created with cd HOGLevel_Gains && for i in *; do while read line; do echo "$i   $line"; done < $i >> ../LOSSMAP; done && > sed -i 's/^/G\t/' ../GAINMAP
#cd HOGLevel_Losses && for i in *; do while read line; do echo "$i   $line"; done < $i >> ../LOSSMAP; done && > sed -i 's/^/L\t/' ../LOSSMAP
#cd .. && cat GAINMAP LOSSMAP > level2hoglossgain.map

    
    #con = sqlite3.connect(":memory:")
    con = sqlite3.connect("mykeggdb")
    cur = con.cursor()
    
    
    createSqlTable(con, KEGGKO2PATHWAY_MAP, 'keggK2pathway', 
                   ['KO', 'Description', 'Pathway', 'Category'])
    
    
    createSqlTable(con, PROTEIN2KEGG_MAP, 'OrganismProtein2KeggResults', 
                   ['Organism', 'ProteinID', 'KO_primary', 'Description', 'SimScore_primary', 'KO_secondary', 'SimScore_secondary' ], mode='append')
    
    
    createSqlTable(con, ORGANISMPROTEIN2HOG_MAP, 'OrganismProtein2HOG',
                   ['Organism', 'ProteinID', 'HOG' ])
    
    createSqlTable(con, LOSSGAINHOG2Level, 'lossgainhog2level',
        ['Event', 'Level', 'HOG' ])
    
    
    
    con.execute('''CREATE TABLE IF NOT EXISTS OrgProHog2Keggres AS
                   SELECT OrganismProtein2KeggResults.Organism, 
                          OrganismProtein2KeggResults.ProteinID,
                          OrganismProtein2HOG.HOG,
                          OrganismProtein2KeggResults.KO_primary, 
                          OrganismProtein2KeggResults.Description, 
                          OrganismProtein2KeggResults.SimScore_primary, 
                          OrganismProtein2KeggResults.KO_secondary, 
                          OrganismProtein2KeggResults.SimScore_secondary 
                   FROM OrganismProtein2KeggResults LEFT OUTER JOIN OrganismProtein2HOG ON 
                        ((OrganismProtein2KeggResults.Organism = OrganismProtein2HOG.Organism) AND
                         (OrganismProtein2KeggResults.ProteinID = OrganismProtein2HOG.ProteinID) );
                ''')
    con.commit()
    
    sql = '''CREATE TABLE IF NOT EXISTS HogKo2Pathway AS 
                SELECT OrgProHog2Keggres.HOG, keggK2pathway.KO, keggK2pathway.Pathway, keggK2pathway.Category
                                    FROM OrgProHog2Keggres, keggK2Pathway
                                    WHERE 
                                    OrgProHog2Keggres.KO_primary = keggK2pathway.KO;'''
    
    con.execute(sql)
    
    sql = '''CREATE TABLE IF NOT EXISTS gainsAtlevel AS 
                SELECT * FROM HogKo2Pathway LEFT OUTER JOIN lossgainhog2level ON (HogKo2Pathway.HOG = lossgainhog2level.HOG) WHERE lossgainhog2level.Event = 'G';'''
    
    con.execute(sql)

    
    sql = '''CREATE TABLE IF NOT EXISTS lossesAtlevel AS 
                SELECT * FROM HogKo2Pathway LEFT OUTER JOIN lossgainhog2level ON (HogKo2Pathway.HOG = lossgainhog2level.HOG) WHERE lossgainhog2level.Event = 'L';'''
    
    con.execute(sql)
    
    
#     for i in range(77):
#         sql = 'select HOG, KO, Pathway, Category from gainsatlevel where Level={};'.format(i)
#         df = pd.read_sql(sql, con)
#         df.to_csv('/share/project/bardya/Enterobacteriaceae/kegg/only_relevant/kegg_HOGLevel_Gains/{}'.format(i), sep='\t', header=False, index=False)
#          
#     for i in range(77):
#         sql = 'select HOG, KO, Pathway, Category from lossesatlevel where Level={};'.format(i)
#         df = pd.read_sql(sql, con)
#         df.to_csv('/share/project/bardya/Enterobacteriaceae/kegg/only_relevant/kegg_HOGLevel_Losses/{}'.format(i), sep='\t', header=False, index=False)        
#  

   
    hogrange = [0,1,2] 
    losses = ()
    df = pd.DataFrame()
    for i in hogrange:
        sql = 'select HOG from lossesatlevel where Level = {};'.format(i)
        df = pd.read_sql(sql, con)
        losses += tuple(set(df['HOG']))        
      
    r = str(tuple(hogrange))
    sql = 'select HOG, KO, Pathway, Category from gainsatlevel where Level in {} and not HOG in {};'.format(r, losses)
    df = pd.read_sql(sql, con)
    df.to_csv('/share/project/bardya/Enterobacteriaceae/kegg/only_relevant/kegg_HOGLevel_LCA/{}_LCA_set.txt'.format(hogrange[-1]), sep='\t', header=False, index=False)       
      
     

    con.close()


#grep -v "GCF_001281565__1_ASM128156v1_protein" protein2keggK0.map | grep -v "GCF_000009065__1_ASM906v1_protein" | grep -v "GCF_000834335__1_ASM83433v1_protein" | grep -v "GCF_000314875__2_ASM31487v2_protein" | grep -v "GCF_000330865__1_ASM33086v1_protein" | grep -v "GCF_000240185__1_ASM24018v2_protein" | grep -v "GCF_000742135__1_ASM74213v1_protein" | grep -v "GCF_000788015__1_ASM78801v1_protein" | grep -v "GCF_000025565__1_ASM2556v1_protein" | grep -v "GCF_000026345__1_ASM2634v1_protein" | grep -v "GCF_000012005__1_ASM1200v1_protein" | grep -v "GCF_000006945__1_ASM694v1_protein" | grep -v "GCF_000252995__1_ASM25299v1_protein > protein2keggK0.map        
        
        
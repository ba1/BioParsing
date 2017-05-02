'''
Created on Jan 16, 2017

@author: bardya
'''

import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse
import statsmodels.sandbox.stats.multicomp as sm
import seaborn as sns
import os

def parse_args(args):
    parser = argparse.ArgumentParser(description='From a set of multi fasta files in the input-directory extract 1 or more representative sequences randomly for each')
    
    parser.add_argument('-t', '--test-set', dest='testset', metavar='<path_to_csv>', type=str, #required=True,
                   help='path to a directory containing the csv')
    
    parser.add_argument('-r', '--reference-set', dest='refset', metavar='<path_to_csv>', type=str, #required=True,
                   help='path to a directory containing the csv')
    
    parser.add_argument('-a', '-alpha', '-thresh', '--threshold', dest='thresh', metavar='<float [0:1]>', type=float, default=0.05,#required=True,
                   help='a floating point number indicating threshold for the p value and FDR correction')
    
    
    parser.add_argument('-mode', dest='mode', type=int,
                        help='specify the number of the column in both input sets containing the "terms" you want to test for enrichment')
    
    parser.add_argument('--version', action='version', version='0.1')
    
    return parser.parse_args(args)


def generate2x2ContingencyTables(df):
    
    df = df.fillna(0).T

    N, n = df.sum(axis=1).tolist()
    
    for c in df.columns:
        df['Non-'+ str(c)] = pd.Series([N - df.loc['Reference_set'][c], n - df.loc['Test_set'][c]], index=df.index)

        yield df[[c, 'Non-' + c]]



if __name__ == '__main__':
    #args = parse_args()
    
    basedir = '/share/project/bardya/Enterobacteriaceae/kegg/only_relevant/kegg_HOGLevel_LCA/'
    ref = '2_LCA_set.uniq.txt'
    test = '3_Losses.uniq.txt'
    
    outfile = ref.split('.',1)[0] + '_vs_' + test.rsplit('.', 1)[0] + '.csv'
    
    args = parse_args(['-r',os.path.join(basedir,ref), 
                       '-t',os.path.join(basedir,test),
                      '-alpha', '0.05'])    #for testing purposes
    
    keggK2pathway = pd.read_csv('/share/project/bardya/dbs/kegg/keggK2pathway.map', names=['KORTH', 'Description', 'Pathway', 'Category'], sep='\t')
    
    
    testset_df = pd.read_csv(args.testset, sep='\t', header=None, names=["HOG", "KORTH", "KODESC", "CATEG"])
    refset_df = pd.read_csv(args.refset, sep='\t', header=None, names=["HOG", "KORTH", "KODESC", "CATEG"])
    
    term = "KODESC"
    
    s1 = testset_df.groupby(term, as_index=True).size();
    s2 = refset_df.groupby(term, as_index=True).size();
    
    s2 = s2[s2 > 2]
    df = pd.DataFrame(dict(Test_set = s1, Reference_set = s2))

    fisherResults = {}
    id_value = 0
    
    for ct in generate2x2ContingencyTables(df):
        
        oddsratio, pvalue = stats.fisher_exact(ct, alternative="two-sided")        
        #oddsratio, pvalue = stats.fisher_exact(ct)
        
        if (pvalue <= args.thresh):
            if (pvalue <= args.thresh):
                if oddsratio >= 1:
                    sig = 'under'.upper()
                else:
                    sig = 'over'.upper()
        else:
            sig = '-'
        

        if term == "KODESC":
            termlabel = "Pathway"
            termname = ct.columns[0].split(' ',1)[1].split(" [")[0]
            id_label = "k ID"
            id_value = "k" + str(ct.columns[0].split(' ',1)[0])
            id_start = str(ct.columns[0].split(' ',1)[0])
        elif term == "CATEG":
            termlabel = "Category"
            termname = ct.columns[0]
            id_label = "ID"
            id_value += 1
            id_start = ct.columns[0]
            
        genelist_test = ", ".join(testset_df[testset_df[term].str.startswith(id_start)]["HOG"].tolist())
        genelist_ref  = ", ".join(refset_df[refset_df[term].str.startswith(id_start)]["HOG"].tolist())
        
        fisherResults[id_value] = {id_label: id_value,
                            termlabel: termname,
                            'p-Value': pvalue,
                            '#Test': ct.loc['Test_set', ct.columns[0]],
                            '#Ref': ct.loc['Reference_set', ct.columns[0]],
                            '#notAnnotTest': ct.loc['Test_set', ct.columns[1]],
                            '#notAnnotRef': ct.loc['Reference_set', ct.columns[1]],
                            'Over/Under': sig,
                            'TestSeqs': genelist_test,
                            'RefSeqs': genelist_ref }

    
    fr = pd.DataFrame(fisherResults).T

    benjamini = sm.multipletests(fr['p-Value'], method = 'fdr_bh', alpha=args.thresh)

    fr = pd.concat([fr, pd.Series(benjamini[1], name='FDR', index=fr.index)], axis=1) #p-adjusted
    
    fr = pd.concat([fr, pd.Series(benjamini[0], name='FDR_TEST', index=fr.index)], axis=1) #is_rejected
    
#     sns.set(color_codes=True)
#     sns.distplot(fr['p-Value'], kde=False, bins=20)
#     sns.plt.show()
     
    fr_filtered = fr[fr['p-Value'] <= args.thresh]
    
    fr_filtered.to_csv(os.path.join(basedir, outfile + "ix"), 
                       columns=[termlabel,'FDR','p-Value','#Test','#Ref','#notAnnotTest','#notAnnotRef','Over/Under','TestSeqs','RefSeqs'], 
                       header=True, index_label=id_label, sep='\t')
    
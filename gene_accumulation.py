'''
Created on Jan 4, 2017

@author: bardya
'''

import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import random
import re

#Starting from a Matrix with 1/0 indicating presence/absence of a taxon in a orthologous group
#With OMA the thresholds for pairwise orthology detection are: 
#    Alignment Length minimum 61% of the shorter Sequence 
#    Alignment Score minimum 183
#    Stable Pair Tolerance Value 1.53
#    Minimum Sequence Length: 40aa

    
def readCsvToMap(filepath, fieldnames, delimiter='\t'):
    fh = open(filepath, 'r')
    fh_filtered = [row for row in fh if not row.startswith('#')]
    
    assert len(fieldnames) >= 2
    
    reader = csv.DictReader(fh_filtered, fieldnames=fieldnames, delimiter=delimiter)
    
    xtoymap = {}
    
    if len(fieldnames) == 2:  
        for linedict in reader:
            for k in list(linedict.keys()):
                if k == None or k == '':
                    del linedict[k]
                    continue
                linedict[k] = linedict[k].strip()
            
            xtoymap[linedict[fieldnames[0]]] = '\t'.join(linedict[f] for f in fieldnames[1:])
    
    if len(fieldnames) == 3:
        for linedict in reader:
            if not linedict[fieldnames[0]] in xtoymap.keys():
                xtoymap[linedict[fieldnames[0]]] = {linedict[fieldnames[1]]: linedict[fieldnames[2]]} #.split(' ', 1)[0]
            else:
                #add this to existing subdict
                xtoymap[linedict[fieldnames[0]]][linedict[fieldnames[1]]] = linedict[fieldnames[2]]
        
    return xtoymap

def rarefaction_simultaneous(presabs_df, taxonset):
    
    taxonset = list(taxonset) #make sure it is a list
    random.shuffle(taxonset) #works in-place

    presabs_df = presabs_df[taxonset] #reordered and filtered
    presabs_df = presabs_df.dropna(how='all') #delete all rows with no content
    
    size_vector_pan = ['pan']
    size_vector_cor = ['cor']
    size_vector_sin = ['sin']
    
    joined_df_pan = pd.DataFrame()
    joined_df_cor = pd.DataFrame()
    joined_df_sin = pd.DataFrame()
    
    for num, i in enumerate(presabs_df.columns):
        if num == 0:
            joined_df_pan = presabs_df[i].dropna().to_frame()
            joined_df_cor = presabs_df[i].dropna().to_frame()
            joined_df_sin = presabs_df[i].dropna().to_frame()
            
        s = presabs_df[i].dropna()
        joined_df_pan = pd.merge(joined_df_pan, s.to_frame(), left_index=True, right_index=True, how='outer')
        size_vector_pan.append(float(len(joined_df_pan)))
            
        joined_df_cor = pd.merge(s.to_frame(), joined_df_cor, left_index=True, right_index=True, how='inner')
        size_vector_cor.append(float(len(joined_df_cor)))
           
        joined_df_sin = pd.merge(joined_df_sin, s.to_frame(), left_index=True, right_index=True, how='outer', indicator=True)
        singleton_df = joined_df_sin[joined_df_sin._merge == 'right_only']
        
        joined_df_sin.drop('_merge', axis=1, inplace=True)
        if num == 0:
            size_vector_sin.append(float(len(joined_df_sin)))
            
        else:
            #size_vector_sin.append(float(len(singleton_df)))
            size_vector_sin.append(size_vector_pan[-1]-size_vector_pan[-2])
            
    return [size_vector_pan, size_vector_cor, size_vector_sin]

def prepare_all_dfs(presabs_matrix, TAXONSET, REPS):
    
    artifcolnames = ['cond'] + [float(i) for i in range(len(TAXONSET) + 1)[1:]]

    pan_cor_sin_df = pd.DataFrame(columns=artifcolnames)
        
    for i in range(REPS):
        svec_lst = rarefaction_simultaneous(presabs_matrix, TAXONSET)
        svec_df = pd.DataFrame(np.array(svec_lst), columns=artifcolnames)
        pan_cor_sin_df = pan_cor_sin_df.append(svec_df, ignore_index=True)
    
    pan_cor_sin_df = pan_cor_sin_df.set_index('cond')
    pan_cor_sin_series = pan_cor_sin_df.stack() #flattened
    pan_cor_sin_series.index.names = ['cond', 'genomes'] #name indices
     
    pan_cor_sin_df = pan_cor_sin_series.to_frame() #convert to frame

    pan_cor_sin_df = pan_cor_sin_df.reset_index(level=['genomes', 'cond']) #reset indices to columns
        
    pan_cor_sin_df.columns = ["cond","genomes","genes"]
    pan_cor_sin_df[['genomes', 'genes']] = pan_cor_sin_df[['genomes','genes']].apply(pd.to_numeric) #change column types
    pan_cor_sin_df["cond"] = pan_cor_sin_df["cond"].astype('category')

    return pan_cor_sin_df

        
def exp_decay1(x, a, b, c):
    return a * np.exp(b*x) + c

def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

def power_law(x, a, b, c):
    return a * x**b + c

def power_law_decay(x, a, b, c):
    return a / (x*b) + c

def fit_exp_linear(t, y, C=0):
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    return A, K

def fit_exp_nonlinear(t, y, func):
    if func == exp_decay:
        initial_guess = [2000, 0.65, 100]
    elif func == power_law:
        initial_guess = [1, -1, 1]
    
    opt_parms, parm_cov = sp.optimize.curve_fit(func, t, y, maxfev=20000, p0=initial_guess)
        
    A, K, C = opt_parms
    r2 = calc_rsquared(t,y,func, opt_parms)
    
    alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

    n = len(y)    # number of data points
    p = len(opt_parms) # number of parameters
    
    dof = max(0, n - p) # number of degrees of freedom
    
    # student-t value for the dof and confidence level
    tval = sp.stats.distributions.t.ppf(1.0-alpha/2., dof) 
    
    return A, K, C, r2, zip(range(n), opt_parms, np.diag(parm_cov)), tval
    
def calc_rsquared(t, y, func, opt_parms):
    #residual sum of squares (ss_tot) with
    residuals = y - func(t, *opt_parms)
    ss_res = np.sum(residuals**2)
    #You can get the total sum of squares (ss_tot) with
    ss_tot = np.sum((y-np.mean(y))**2)
    #r_squared-value with,
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getSubTaxonlist(id2taxon_map, regex, one_per_species=True, exclude_sp=False):
        
    TAXONSET = []
    SPECIESLIST = []
    if one_per_species:
    
        for assemblyid, taxon in id2taxon_map.items():
            match = re.search(regex, taxon)
            if match:
                spname = ' '.join(taxon.split(' ')[0:2])
                if exclude_sp:
                    if not spname in SPECIESLIST and not spname[1] == 'sp.':
                        TAXONSET.append(assemblyid)
                        SPECIESLIST.append(' '.join(taxon.split(' ')[0:2]))
                else:
                    if not spname in SPECIESLIST or spname[1] == 'sp.':
                        TAXONSET.append(assemblyid)
                        SPECIESLIST.append(' '.join(taxon.split(' ')[0:2]))
        return TAXONSET
    
    else:
        for assemblyid, taxon in id2taxon_map.items():
            match = re.search(regex, taxon)
            if match:
                spname = ' '.join(taxon.split(' ')[0:2])
                if exclude_sp:
                    if not spname[1] == 'sp.':
                        TAXONSET.append(assemblyid)
                else:
                    TAXONSET.append(assemblyid)
        return TAXONSET



def produceOGMatrix(ogmatrix_df, outpath, TAXONSET):
    try:
        ogmatrix_df_taxonset = ogmatrix_df.copy()
        ogmatrix_df_taxonset = ogmatrix_df_taxonset[TAXONSET]
        ogmatrix_df_taxonset[ogmatrix_df_taxonset > 1] = 1
#         ogmatrix_df_taxonset = ogmatrix_df_taxonset.replace(0, np.nan)
#         ogmatrix_df_taxonset.dropna(how='all', inplace=True)
#         ogmatrix_df_taxonset.fillna(value=0, inplace= True)
        ogmatrix_df_taxonset.to_csv(outpath, header=False, sep='\t', index=False)
        
        import subprocess
        excode = subprocess.call(["sed -i 's/\t//g' {}".format(outpath),], shell=True)
    
    except Exception as e:
        print("Failed to create correct file to " + outpath)
        print(''' > Error Message:\n\n'''.format(e))


if __name__ == '__main__':
    
    #### DATA IO PART 1####
    
    #filename = '/share/project/bardya/Enterobacteriaceae/OMA_prot/Results/OrthologousMatrix.txt'
    fieldnames = ['assemblyID', 'Taxon']
    id2taxon_map = readCsvToMap('/share/project/bardya/Enterobacteriaceae/OMA_prot/assembly2strain_protein.map', fieldnames)
    
    fieldnames = ['Taxon', 'segNum', 'geneID']
    taxon2seqNum2geneid_map = readCsvToMap('/share/project/bardya/Enterobacteriaceae/OMA_prot/Results/Map-SeqNum-ID.txt', fieldnames)

    fh = open('/share/project/bardya/Enterobacteriaceae/OMA_prot/Results/OrthologousMatrix.txt', 'r')
    ogmatrix = pd.read_csv(fh, sep='\t', header='infer', comment='#' )
    
    presabs_matrix = ogmatrix.copy()
    presabs_matrix[presabs_matrix > 1] = 1
    presabs_matrix = presabs_matrix.replace(0, np.nan)
    
    
    #### DATA IO PART 2####
    
    
    #### PAN CORE AND SINGLETONS DATA EXTRACTION ####
    
    REPS = 100
    #regex = '.*Photorhabdus.*'
    #TAXONSET = getSubTaxonlist(id2taxon_map, regex, one_per_species=True, exclude_sp=False)
    regex = '.*Xenorhabdus.*'
    TAXONSET = getSubTaxonlist(id2taxon_map, regex, one_per_species=False, exclude_sp=False)
    #TAXONSET.extend(getSubTaxonlist(id2taxon_map, regex, one_per_species=False, exclude_sp=False))
    print('Analysis started for {} Taxa...'.format(len(TAXONSET)))

    #produceOGMatrix(ogmatrix, '/share/project/bardya/Enterobacteriaceae/OMA_prot/OrthologousMatrix_only_1-0.txt', TAXONSET)
    
    pan_core_sin_df = prepare_all_dfs(presabs_matrix, TAXONSET, REPS)

    #x = pan_core_sin_df[pan_core_sin_df["cond"] == "cor"].drop('cond', 1)
    #x.to_csv('/home/bardya/Desktop/testcore.tsv', sep='\t', index=False, header=False)
    
    
    #### NON-LINEAR DATA FITTING ####
    
    fitted_y = {}
    
    t = np.linspace(1, np.ceil(len(TAXONSET)), num=np.ceil(len(TAXONSET))) # Generate x_coords based on these
    
    for cond in ['pan', 'cor', 'sin']:
        
        #running_median = [list(pan_core_sin_df[(pan_core_sin_df['genomes'] == k) & (pan_core_sin_df['cond'] == cond) ].median())[1] for k in range(len(TAXONSET)+1)[1:]]
        #print("running median {}: {}".format(cond, running_median))
        #[ax.fig.annotate(p[0], (p[1]+50, p[1]), color='r') for p in zip(range(len(TAXONSET)+1)[1:], running_median)]
        
        if cond == 'pan':
            func = power_law
        else:
            func = exp_decay
        
        m = pan_core_sin_df[pan_core_sin_df.cond == cond ].drop('cond', 1)
        noisy_y = m["genes"]
        A, K, C, r2, confzip, tval  = fit_exp_nonlinear(m["genomes"], noisy_y, func)
        
        #print(cond.upper(), A, K, C, r2)
        fitted_y[cond] = func(t, A, K, C)
        
        for i, p, var in confzip:
            sigma = var**0.5
            print('{0} p{1}: {2} ±{3}'.format(cond.upper(), i, p, sigma*tval))
        
        print('{0} r-squared: {1}'.format(cond.upper(),r2))
        
        
    #### DATA VISUALIZATION ####
    
    sns.set(color_codes=True)
    palette1={'pan':'b', 'sin':'g', 'cor':'r'}
    palette2={'pan':'darkblue', 'sin':'darkgreen', 'cor':'darkred'}
    
    ax = sns.lmplot(x="genomes", y="genes", 
                    hue="cond", 
                    data=pan_core_sin_df, 
                    #x_estimator=np.median, 
                    fit_reg=False, 
                    #ci=100, 
                    scatter=True, #x_jitter=True
                    palette=palette1, 
                    markers=["o", "s", "v"], 
                    legend=True) #use col="cond" instead to split into three plots with shared y-axis

    ax.set(xlabel='Number of Genomes added', 
           ylabel='Number of Ortholog Clusters', 
           title='Genome Rarefaction Analysis ({} Random Shuffles)'.format(REPS), 
           xlim=1)    
    
    
    for cond, fit_y in fitted_y.items():
        plt.plot(t, fit_y, c=palette2[cond], linewidth=2)
    
    
    sns.plt.xlim(0.75,)
    sns.plt.ylim(0,)
    
    
#     #From PanGP
#     fit_y = exp_decay1(t, 4005.52, -0.81, 1767.86)
#     #plt.plot(t, fit_y, c='m', linewidth=2)
#     
#     #From mycurvefit.com
#     fit_y = exp_decay1(t, 5045.184,-0.8953288, 1845.755)
#     #plt.plot(t, fit_y, c='m', linewidth=2)

#     
# #     # Linear Fit (Note that we have to provide the y-offset ("C") value!! COR
# #     m = pan_core_sin_df[pan_core_sin_df.cond == 'cor'].drop('cond', 1)
# #     noisy_y = m["genes"]
# #     A, K = fit_exp_linear(m["genomes"], noisy_y, 632)
# #     fit_y = exp_decay(t, A, K, 632)
# #     plt.plot(t, fit_y, linewidth=2)
    
   
    #sns.plt.show()
    plt.savefig("/share/project/bardya/Enterobacteriaceae/OMA_prot/pan_cor_sin_Xeno_all_taxa.svg", format='svg', dpi=1200)
    #plt.savefig("/share/project/bardya/Enterobacteriaceae/OMA_prot/pan_cor_sin_Xeno_Photo_single_species.svg", format='svg', dpi=1200)
'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse

#Info:
#The primary entry point for the ElementTree library is the parse() function,
#which can take a filename or a file-like object. This function parses the entire document at once.
#If memory is tight, there are ways to parse an XML document incrementally instead.

NAMESPACE = "{http://orthoXML.org/2011/}"

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

    parser.add_argument('--no-accessory', dest='noaccesory_flag', metavar='<produce accessory genomes flag>', action='store_const', const=1, default=0,
                   help='if set, script does not give out the accessory genomes into the directory "./Accessory" relativ to omawd')

    parser.add_argument('-t', '--speciestree', dest='nwcktree', metavar='<path/to/tree.file>', type=argparse.FileType('rt'),
                   help='path to a file containing the species tree in string representation')
    
    parser.add_argument('-tm', '--oma2taxmap', dest='oma2taxmap', metavar='<path/to/oma2tax.map>', type=argparse.FileType('rt'),
                   help='path to a file containing a tab separated mapping of OMA genome ids to taxon names, each mapping in a new line')
#'/share/project/bardya/Enterobacteriaceae/OMA_prot/assembly2strain_protein.map'

    parser.add_argument('--version', action='version', version='0.1')

    return parser.parse_args()



def parseXMLTree(inpt, fromstring=False):
    try:
        from lxml import etree
    except ImportError:
        import xml.etree.ElementTree as etree

    if not fromstring:
        doc = etree.parse(inpt)
        root = doc.getroot()
    else:
        root = etree.fromstringlist(inpt)

    return root


def parseDoc(root, omadir):
    """
    Devide the OrthoXML file into child nodes of species and their genes and
    Orthologous Groups and their corresponding genes. Then create two dictionaries
    one mapping the name of each species to its genes and the other mapping the
    OG group id to its genes
    :param root: the root of the xml document
    :param omadir: the working directory where OMA generated its output
    """
    species = []
    groups = []
    for child in root:
        if child.tag == NAMESPACE + "species":
            species.append(child)
        if child.tag == NAMESPACE + "groups":
            groups.append(child)

    gene2species_map, specnames = parseSpecies(species, omadir)
    groups_dict = parseGroups(groups)                   #{gr1ID:[gene1Id,gene2Id,..], gr2ID:[..], ..}

    return gene2species_map, groups_dict


def parseDocHOGAware(root, omadir):
    """
    Devide the OrthoXML file into child nodes of 1) species and their genes and 2)
    Orthologous Groups and their corresponding genes. Then create two dictionaries
    the first maps the name of species on its genes and the other maps the OG group
    by its id onto its genes
    :param root: the root of the xml document
    :param omadir: the working directory where OMA generated its output
    """
    species = []
    groups = []
    for child in root:
        if child.tag == NAMESPACE + "species":
            species.append(child)
        if child.tag == NAMESPACE + "groups":
            groups.append(child)
    
    gene2species_map, specnames = parseSpecies(species, omadir)
    groups_dict, HOG_groups_dict = parseGroupsHOGAware(groups, all_alias='LUCA', all_list=specnames)                   #{gr1ID:[gene1Id,gene2Id,..], gr2ID:[..], ..}

    return gene2species_map, groups_dict, HOG_groups_dict

def parseSpecies(species, omadir):
    gene2species_map = {}        #for every gene point at species name, memory efficient
    specnum = 0
    specnames = []
    for child in species:
        specnum += 1
        specname = child.attrib["name"]
        specnames.append(specname)
        genes = child.findall(".//"+ NAMESPACE + "gene")
        dbpath = child.find("./" + NAMESPACE + "database").attrib["version"]
        seq_dict = loadDatabase(omadir, dbpath)
        for gene in genes:
            gene2species_map[gene.attrib["id"]] = {"id":gene.attrib["id"],
                                            "protId":gene.attrib["protId"],
                                            "species":specname,
                                            "seq":seq_dict[gene.attrib["protId"]]}
    return gene2species_map, specnames

def parseGroups(groups):
    groups_dict = {}
    for group in groups:
        orthogroups = group.findall("./" + NAMESPACE + "orthologGroup")
        for og in orthogroups:
            genes = og.findall(".//" + NAMESPACE + "geneRef")
            groups_dict[og.attrib["id"]] = [gene.attrib["id"] for gene in genes] #dict-value = a list of geneids #dict-key: og group id

    return groups_dict


def parseGroupsHOGAware(groups, all_alias='', all_list=''):
    """
    Groups all orthologous groups of the highest level (these are HOGs) into Bins
    with the same level attribute, i.e. same list of species
    """
    from collections import defaultdict
    
    groups_dict = defaultdict(tuple)
    
    HOG_groups_dict = defaultdict(dict)
    
    groups_count_dict = defaultdict(list)
    ak=""
    hep=0
    c = 0

    for group in groups:
#         orthogroups = group.findall(".//" + NAMESPACE + "orthologGroup")

        orthogroups = group.findall("./" + NAMESPACE + "orthologGroup") #only highest level
        
        for og in orthogroups:
            level = og.findall(NAMESPACE + "property")[0].attrib["value"]

            numsubogs = len(og.findall(NAMESPACE + "orthologGroup"))
            numsubpogs = len(og.findall(NAMESPACE + "paralogGroup"))

            genes = og.findall(".//" + NAMESPACE + "geneRef")
            gene_ids = tuple([gene.attrib["id"] for gene in genes])

            #print("LEVEL={}\tID={}\t#GENES={}".format(level, og.attrib['id'], len(gene_ids)))
            
            groups_dict[level] += gene_ids # '+=' operator extends the tuple
            HOG_groups_dict[level][og.attrib['id']] = gene_ids
            #groups_count_dict[level].extend([numsubogs])

#ONLY TESTING START
#             if og.attrib['id'] == "57":
#                 print("IDCOUNT:", c)
#                 ak = level
#
#             if ak:
#                 if hep < len(HOG_groups_dict[ak]["57"]):
#                     print(len(HOG_groups_dict[ak]["57"]), "KA{}AK".format(og.attrib["id"]), gene_ids)
#                     if level == ak:
#                         print("JA")
#                     hep = len(HOG_groups_dict[ak]["57"])
#                     c+=1
#                     if c>5:
#                         exit()
#     for k,v in groups_dict.items():
#         print('\n')
#         print(k,':',len(v))
#
#     print('\n***************************************************************\n')
#
#     for k,v in groups_count_dict.items():
#         print('\n')
#         print(k,':',sum(v))
#ONLY TESTING END

    if all_list and all_alias:
        sall_list = '/'.join(sorted(all_list))
        groups_dict[sall_list] = groups_dict[all_alias]
        HOG_groups_dict[sall_list] = HOG_groups_dict[all_alias]
        del groups_dict[all_alias]
        del HOG_groups_dict[all_alias]

#     import pprint
#     pprint.pprint(groups_count_dict)

    #return groups_dict #dict-key: level, #dict-value = list of gene ids
    return groups_dict, HOG_groups_dict


def produceAccessory(gene2species_map, groups_dict, omadir, force=0):
    lstAllGrpdGenes = []
    for grpID, genes in groups_dict.items():
        lstAllGrpdGenes += genes

    print("Number of Total Genes: {}\nNumber of Genes groups in {}s : {} with {} double occurences.".format(len(gene2species_map),MODE, len(lstAllGrpdGenes), len(lstAllGrpdGenes) - len(set(lstAllGrpdGenes))))

    accessory_geneids = sorted(list(set(gene2species_map.keys()) - set(lstAllGrpdGenes)))
    species_accessory_genomes = {}
    for geneID in accessory_geneids:
        if gene2species_map[geneID]["species"] in species_accessory_genomes:
            species_accessory_genomes[gene2species_map[geneID]["species"]].append(gene2species_map[geneID])
        else:
            species_accessory_genomes[gene2species_map[geneID]["species"]] = [gene2species_map[geneID]]

    outdir = os.path.join(omadir, "Accessory")
    clearCheckPath(outdir, force=force)
    for species, gene_dict_lst in species_accessory_genomes.items():
        filename = species + "_singletons.fa"
        with open(os.path.join(outdir, filename), 'w') as outfile:
            outfile.write(genelist2fasta(gene_dict_lst))

##def checkSpeciesNumber(groups_dict, gene2species_map, specnum):
##    res = dict(zip(list(range(specnum)),{{}}*specnum))
##    res_single_copy = dict(zip(list(range(specnum)),{{}}*specnum))
##
##    resname =
##    resname_single_copy =
##
##
##    for groupID, geneIDS in groups_dict.items():
##        namelist = [gene2species_map[geneID][1] for geneID in geneIDS]                 # 1 for second tuple element, of k=geneID v=(gene, specname, database)
##        unique_names = set(namelist)
##
##        if len(unique_names) == len(namelist):
##            res_single_copy[len(unique_names)][frozenset(namelist)][groupID] = [gene2species_map[geneID] for geneID in geneIDS]}                 #produce e.g. res[15]= [32, 123, 14] , i.e. list of group numbers
##
##        res[len(unique_names)][frozenset(namelist)][groupID] = [gene2species_map[geneID] for geneID in geneIDS]}
##
##    return res

def clearCheckPath(outdir, force=0):
    if os.path.isdir(outdir):
        if force:
            import shutil
            shutil.rmtree(outdir, ignore_errors=True)
        else:
            raise IOError("Output Directory already exiting. Specify '-f' option to force overwrite")

    os.makedirs(outdir)


def produceOGFastas(gene2species_map, groups_dict, omadir, force=0, stats=1):
    """Produce fasta files for each OG or HOG into a subdirectory *_Bins
    """
    outdir = os.path.join(omadir,"Bins" + "_" + MODE)
    clearCheckPath(outdir, force=force)

    singcpy_groups_dict = groups_dict.copy()
    group_by_species_dict = {}

    for groupID, geneIDS in groups_dict.items():
        gene_dict_lst = [gene2species_map[geneID] for geneID in geneIDS]
        namelist = [g["species"] for g in gene_dict_lst]
        unique_names = set(namelist)

        if not str(sorted(namelist)) in group_by_species_dict:
            group_by_species_dict[str(sorted(namelist))] = {groupID:gene_dict_lst}
        else:
            group_by_species_dict[str(sorted(namelist))][groupID] = gene_dict_lst

        if not len(unique_names) == len(namelist):
            del singcpy_groups_dict[groupID]

        filename = "{}{}_{}.fa".format(MODE, len(unique_names), groupID) #also add the species name to distinguish

        with open(os.path.join(outdir, filename), 'w') as outfile:
            outfile.write(genelist2fasta(gene_dict_lst))

    if stats:
        produceStats(groups_dict, categ="ALL")
        produceStats(singcpy_groups_dict, categ="SCO")
        produceStatsB(group_by_species_dict, categ="SPEC")


def produceHOGAwareFastas(gene2species_map, groups_dict, omadir, intnodes, intnodes_order, force=0, stats=1):
    """Produce fasta files for each HOG level into a subdirectory called Bins_HOGAware/NODE_*.fa
    """
    stats = False
    outdir = os.path.join(omadir, "Bins" + "_" + MODE)
    clearCheckPath(outdir, force=force)

    for i in intnodes_order:

        singcpy_groups_dict = groups_dict.copy()
        group_by_species_dict = {}

        for level, geneIDS in groups_dict.items():        #group dict contains a list of genes as its values
            gene_dict_lst = [gene2species_map[geneID] for geneID in geneIDS]
            namelist = [g["species"] for g in gene_dict_lst]
            unique_names = set(namelist)                 # get the corresponding species names for each gene in the group
            #assert sorted(list(unique_names)) == level.split('/') #is this meaningful? It could be that there are no specific genes of a taxa located here

            if not str(sorted(namelist)) in group_by_species_dict:
                group_by_species_dict[str(sorted(namelist))] = {level:gene_dict_lst}
            else:
                group_by_species_dict[str(sorted(namelist))][level] = gene_dict_lst

            if not len(unique_names) == len(namelist):
                del singcpy_groups_dict[level]

            filename = "{}_{}.fa".format('NODE', i)

            with open(os.path.join(outdir, filename), 'w') as outfile:
                outfile.write(genelist2fasta(gene_dict_lst))    #also add the species name to distinguish

    if stats:
        produceStats(groups_dict, categ="ALL")
        produceStats(singcpy_groups_dict, categ="SCO")
        produceStatsB(group_by_species_dict, categ="SPEC")


def produceHOGLevels(HOG_groups_dict, intnodes, gene2species_map, force=0):
    """For each level defined by the id of the internal node, give out a file listing all corresponding
    HOGs linewise with the #genes it containes and the number of unique taxa they represent, #leaves in subtree
    and a phyletic pattern indicating presence/absence of the gene in the alphabetically ordered list of taxa"""

    #How it works:
    #Iterate through all levels and their assigned HOGs. For each HOG resolve
    #the species which belongs to the geneID stored and store in a set 
    #The difference of the set when compared to the HOG-Level set of species
    #indicates the missing species
    #The subfunction ReducedLossCount will determine if any subset of this set
    #of missing species is itself an inner node of the tree


    levelwise_speccount_dict = {i: len(i.split('/')) for i in HOG_groups_dict.keys()}
    levelwise_speccount_dict["LUCA"] = max(levelwise_speccount_dict.values()) #re-introduce the "LUCA" level, the set now contains both 

    ##PHYLETIC PATTERN START###
    from collections import OrderedDict
    maxkey = max(HOG_groups_dict.keys(), key=lambda x: len(x))
    maxkey = maxkey.split('/')

    phyletic_pattern = OrderedDict()
    [phyletic_pattern.update({x:0}) for x in maxkey] #maxkey is alphabetically sorted

    print("## Phyletic Pattern - Order of Taxa ##")
    print(phyletic_pattern)
    print('\n')
    ##PHYLETIC PATTERN END ###


    outdir = os.path.join(omadir, "HOGLevel")
    outdir2 = os.path.join(omadir, "HOGLevel_Losses")

    clearCheckPath(outdir, force=force)
    clearCheckPath(outdir2, force=force)

    inv_intnodes = {v: k for k, v in intnodes.items()} # key:SpeciesLevel, value=NodeID
    
    
    intnodes_species_sets = [set(k.split('/')) for k in inv_intnodes.keys()] #species sets of all internal nodes

    from collections import defaultdict
    global loss_count_dict
    global loss_dict
    loss_dict = defaultdict(list)
    loss_count_dict = defaultdict(lambda: 0)

#     def getReducedLossCount(intnodes_species_sets, missing_spec):
#         from itertools import combinations
#         loss_count=len(missing_species)
#         #python idea search from other perspective
#         def reduceLosses(missing_species, loss_count):
#             for i in reversed(range(2, len(missing_species)+1)):
#                 pwrset = [set(comb) for comb in combinations(missing_species,i)]
#                 for s in pwrset:
#                     if s in intnodes_species_sets:
#                         missing_species = missing_species - s
#                         loss_count -= len(s) - 1
#                         loss_count_dict["/".join(sorted(list(s)))]+=1
#                         return reduceLosses(missing_species, loss_count)
#
#             return missing_species,loss_count
#
#         ms, tot_loss_count = reduceLosses(missing_species, loss_count)
#         for i in ms:
#             loss_count_dict[i]+=1
#
#         return tot_loss_count


    def getReducedLossCount(intnodes_species_sets, missing_species, og_id, debug=False):
        """
        """
        totloss_cnt=0

        if not missing_species:  #empty-set speed-up
            return totloss_cnt

        for intnode_set in sorted(intnodes_species_sets, key=lambda x: len(x), reverse=True): #sort descending
                  
            if all(species in missing_species for species in intnode_set): #if all species of intnodeset also in missing_species_set
                missing_species = missing_species - intnode_set
                totloss_cnt += 1
                intnode_name = "/".join(
                                        sorted(
                                               list(intnode_set)
                                               )
                                        )
                
                loss_count_dict[intnode_name] += 1
                loss_dict[intnode_name].append(og_id)
                #if debug:
                #    print("MISSING", missing_species)
                #    print(intnode_set)
                #    print(totloss_cnt)
            if len(missing_species) < 2 : #speed-up if reduced to one missing species
                break

        for i in missing_species:
            totloss_cnt+=1
            #if debug:
                #print(i, totloss_cnt)
            loss_count_dict[i]+=1
            loss_dict[i].append(og_id)

        return totloss_cnt
    

    for level, ogdict in HOG_groups_dict.items():
        
        hog_to_number_map = dict()
        species_in_level = set(level.split('/'))
        for og_id, gene_id_list in ogdict.items():

            species_in_hog = [gene2species_map[gene_id]['species'] for gene_id in gene_id_list]
            uniq_species_in_hog = set(species_in_hog)

            missing_species = species_in_level - uniq_species_in_hog

            #if og_id=="12991":
            #    #print(intnodes_species_sets, "\n","\n")
            #    #print("ALL MISSING", missing_species)
            #    tot_losses = getReducedLossCount(intnodes_species_sets, missing_species, og_id, debug=True)
            #else:
            tot_losses = getReducedLossCount(intnodes_species_sets, missing_species, og_id)

            og_phyl_pat = phyletic_pattern.copy()
            for speciesname in species_in_hog:
                og_phyl_pat[speciesname]=1   #change to +=1 for non-binary vector with paralog infos

            phyltcfngrpr = "".join([str(v) for v in og_phyl_pat.values()])

            hog_to_number_map[og_id] = (len(uniq_species_in_hog),len(species_in_hog), tot_losses, phyltcfngrpr)
        
        
        #sorted_level = "/".join(sorted(level.split('/'), key=lambda x: x.lower()))
        with open(os.path.join(omadir, 'HOGLevel', str(inv_intnodes[level])), 'w') as hn:
            
            hn.write('{}\n'.format('\n'.join(['HOG{}.fa\t{}\t{}\t{}\t{}\t{}'.format(
                                                k, v[0], v[1], levelwise_speccount_dict[level], v[2], v[3]) for k,v in hog_to_number_map.items()])))

    
    for level, ogdict in HOG_groups_dict.items():
        
        with open(os.path.join(omadir, 'HOGLevel_Losses', str(inv_intnodes[level])), 'w') as ln:
            
            ln.write('{}\n'.format('\n'.join(['HOG{}.fa'.format(x) for x in loss_dict[level]])))
    
    intnode_lossdict = {inv_intnodes[k]:v for k,v in loss_count_dict.items() if len(k.split('/')) > 1}
    
    leafnode_lossdict = {k:v for k,v in loss_count_dict.items() if len(k.split('/')) < 2}

    #print(intnode_lossdict)
    #print("\n\n\n")
    #print(leafnode_lossdict)
    #print("\n\n\n\n")
    #print(len(loss_dict.keys()))
    for level, og_list in loss_dict.items():
        if len(level.split('/')) < 2:
            print('\n'.join(['L\t{}\tHOG{}.fa'.format(level, x) for x in og_list]))


    return intnode_lossdict, leafnode_lossdict

def produceHOGGainLossTree(HOG_groups_dict, intnode_lossdict, leafnode_lossdict, intnodes, tree, leafname_dict, outdir):

    import ete3

    for n in tree.traverse():
        if n.is_leaf() and leafname_dict:
            n.name = leafname_dict[n.name]

    schema_names = ete3.COLOR_SCHEMES.keys()

    def mylayout(n):
        if not n.is_leaf():
            #print("\n\n\nintnode_lossdict", intnode_lossdict.keys())
            shortn = int(n.name.split('_')[1])

            number_of_hogs = len(HOG_groups_dict[intnodes[shortn]].keys())
            number_of_gains = number_of_hogs
            try:
                number_of_losses = int(intnode_lossdict[shortn])
            except:
                number_of_losses = 0

            #thisface = ete3.faces.BarChartFace([int(number_of_losses), int(number_of_gains)],
            #                                   colors=ete3.COLOR_SCHEMES["set1"],
            #                              #["Tomato","DarkCyan"],
            #                              labels=["#Losses","#Gains"],
            #                              width=20, height=40, min_value=0, max_value=4100)
            #thisface.opacity = 0.8


            thisfaceA = ete3.faces.TextFace('L: ' + str(int(number_of_losses)) + '\n' + 'G: ' + str(int(number_of_gains)))
            thisfaceA.opacity = 0.8
            ete3.faces.add_face_to_node(thisfaceA, n, column=0, position="branch-right")
            ete3.faces.add_face_to_node(ete3.AttrFace("name"), n, column=1)
        else:
            pass
            #ete3.faces.add_face_to_node(ete3.AttrFace("name"), n, column=0)

   # for n in tree.traverse():
   #     if not n.is_leaf():

    ts = ete3.TreeStyle()
    ts.layout_fn = mylayout
    #ts.rotation = 90
    ts.allow_face_overlap = True
    ts.title.add_face(ete3.faces.TextFace("HOG Gain Loss Eventcount", fsize=20), column=0)
    tree.render(os.path.join(outdir,"img_faces.png"), tree_style = ts)


def produceStats(groups_dict, categ= ''): #{grid1:genelist1, grid2:genelist2}
    count_dict = {}

    for k,v in groups_dict.items():
        if len(v) in count_dict:
            count_dict[len(v)].append(v)
        else:
            count_dict[len(v)] = [k]
    print("Category {}".format(categ))
    print("\nTotal #Genes")
    print(len(groups_dict.keys()))
    print("\nStatistical Overview: #Groups_per_Group_Size")
    for size in sorted(list(count_dict.keys())):
        print("{}:\t{}".format(size,len(count_dict[size])))
    #print("Largest Group (#Genes)")
    #print("{}:\t{}".format(max(count_dict)),count_dict[max(count_dict)])
    #print("Smallest Group (#Genes)")
    #print("{}:\t{}".format(min(count_dict)),count_dict[min(count_dict)])
    print()
    print("\nAll Groups")
    print("{}:\t{}\t{}".format("Group_Id","#Genes","Genes"))
    for k,v in groups_dict.items():
        print("{}:\t{}\t{}".format(k,len(v), str([i for i in v]) ) )

def produceStatsB(groups_dict, categ= ''):
    print("All groups grouped by species in group")
    print("{}:\t{}\t{}\t{}\t{}".format("Group of Species","#Members","#Species","#Groups","Group_Ids"))
    lis = []
    for k,v in groups_dict.items(): # {"spname1,spname2,spname3,..": {}] }
#         list2d = [t for t in v]
#         import itertools
#         merged = list(itertools.chain.from_iterable(list2d))
        print()
        lis.append((set(k.split(',')),len(k.split(',')),len(set(k.split(','))),len(v),str(list(v.keys()))))

    slis = sorted(lis,key=lambda x: x[1])
    for entry in slis:
        print("{}:\t{}\t{}\t{}\t{}".format(*entry))

#OUTPUTSUUGGEST : {HOG ID, #totalSeq, #totalSpecies, SpeciesA:[#seqs], SpeciesN:[#seqs]}

def getInternalNodes(nwcktreepath):
    '''Parse the tree file (argument) and give each internal node a name
    Save all internal nodes in a dictionary called intnodes with its number
    as the key and a sorted list of all its leaves as the value. Save then
    the value as a single string separating the leaves with forward slashes
    and also give out the labeled tree in newick format in a textfile

    Return intnodes dictionary and the order of the int_nodes as traversed
    in the first step.
    '''

    from ete3 import Tree

    try:
        handle = open(nwcktreepath, 'r')
        nwcktree = handle.read().rstrip()
    except:
        print("IOError importing tree at" + nwcktreepath)

    unrooted_tree = Tree(nwcktree)

    i = 0
    orderings = []
    intnodes = dict()
    for node in unrooted_tree.traverse("preorder"):
        if not node.is_leaf():
            node.name = "NODE_%d" %i
            orderings.append((i, node.dist))
            intnodes[i] = sorted([subnode.name for subnode in node.get_descendants() if subnode.is_leaf()])
            i += 1

    intnodes_order = [e[0] for e in orderings]
    #intnodes_order = [e[0] for e in sorted(orderings, key=lambda x: x[1], reverse=True)]
    #print(intnodes_order    )

    for k,v in intnodes.items():
        intnodes[k] = '/'.join(v)

    #for k,v in intnodes.items():
    #   print("DICT\t{}\t{}".format(v,k))


    unrooted_tree.write(format=1, outfile=nwcktreepath + '_labeledInternalNodes')

    return unrooted_tree, intnodes, intnodes_order


def loadDatabase(omadir, dbpath):
    with open(os.path.join(omadir, dbpath), 'r') as db:
        for i in range(4):
            next(db)
        import itertools
        it = itertools.chain('<root>', db, '</root>')
        root = parseXMLTree(it, fromstring=True)
        genes = root.findall('.//E')
        seq_dict = {}
        for entry in genes:
            seq_dict[entry.find('./ID').text] = entry.find('./SEQ').text

        return seq_dict


def genelist2fasta(gene_dict_lst):
    res_str = ""
    for gene in gene_dict_lst:
        res_str += ">{}\n{}\n".format(gene["protId"],gene["seq"])

    return fasta_linebreak(res_str.rstrip())


def fasta_linebreak(s, seqrowlength=70):
    fastastring_formatted = []
    for line in s.split("\n"):
        if not line.startswith(">"):
            n = 70
            linesplit = [line[i:i+n] for i in range(0, len(line), n)]
            formattedline = "\n".join(linesplit) + "\n"
            fastastring_formatted.append(formattedline)
        else:
            fastastring_formatted.append(line + "\n")

    return "".join(fastastring_formatted)

if __name__ == '__main__':
    #argumentlist = ['--sum', '7', '-1', '42']
    #args = parse_args(argumenlist)

    args = parse_args()

    if args.mode.lower() == "hogs" or args.mode.lower() == "hog":
        filename = "HierarchicalGroups.orthoxml"
        MODE = "HOG"                                         #global Variable!!! be aware
    elif args.mode.lower() == "hogaware" or args.mode.lower() == "hogsaware":
        filename = "HierarchicalGroups.orthoxml"
        MODE = "HOGAware"
    else:
        filename = "OrthologousGroups.orthoxml"
        MODE = "OG"

    import csv
    if args.oma2taxmap and os.path.isfile(args.oma2taxmap.name):
        with open(args.oma2taxmap.name, mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            leafname_dict = {rows[0]:rows[1] for rows in reader if not rows[0].startswith('#')}
    else:
        leafname_dict = None

    omadir = args.oma_workpath
    omaresdir = os.path.join(omadir, args.oma_output_dirname) #Path to 'Results' directory within OMA wd
    omacachedir = os.path.join(omadir,"Cache", "DB")

    if not os.path.isdir(omadir):
        raise IOError("Problem finding specified OMA working dir")
    else:
        if not (os.path.isdir(omacachedir)):
            raise IOError("Problem finding specified OMA Cache dir and/or Cache/DB dir within the OMA working dir")
        else:
            if not (os.path.isdir(omaresdir)):
                raise IOError("Problem finding specified OMA output dir within the OMA working dir")
            else:
                orthoxml_file = open(os.path.join(omaresdir,filename), 'r')

    root = parseXMLTree(orthoxml_file)

    gene2species_map, groups_dict = parseDoc(root, omadir)

    if args.nwcktree and MODE =="HOGAware":
        gene2species_map, groups_dict, HOG_groups_dict = parseDocHOGAware(root, omadir)
        unrooted_tree, internal_nodes, internal_nodes_order = getInternalNodes(args.nwcktree.name)
        intnode_lossdict, leafnode_lossdict = produceHOGLevels(HOG_groups_dict, internal_nodes, gene2species_map, force=args.force_flag)
        produceHOGGainLossTree(HOG_groups_dict, intnode_lossdict, leafnode_lossdict, internal_nodes, unrooted_tree, leafname_dict, omadir)
        #produceHOGAwareFastas(gene2species_map, groups_dict, omadir, internal_nodes, internal_nodes_order, force=args.force_flag, stats=args.stats_flag)
    else:
        gene2species_map, groups_dict, = parseDoc(root, omadir)
        produceOGFastas(gene2species_map, groups_dict, omadir, force=args.force_flag, stats=args.stats_flag)

    if not args.noaccesory_flag:
        produceAccessory(gene2species_map, groups_dict, omadir, force=args.force_flag)

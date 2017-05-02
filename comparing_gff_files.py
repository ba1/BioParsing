'''
Created on Jul 22, 2015

@author: bardya
'''

#Description of the gff file format: http://www.sequenceontology.org/gff3.shtml

class Gene:
    def __init__(self, seqid = '', source = '', genetype = '', start = 0, end = 0, score = '', strand = '', phase = '', attributes = ''):
        self.seqid = seqid
        self.source= source
        self.genetype = genetype
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self.length = self.end-self.start+1 if not (self.start == 0 and self.end ==0) else 0
    
    def __eq__(self, other): 
        return self.__dict__ == other.__dict__
    
class Genome:
    def __init__(self):
        self.chromosomes = dict()
    
    def __len__(self):
        return len(self.chromosomes)
    
    def numberofgenes(self):
        l = 0
        for k,v in self.chromosomes.items():
            l += len(v)
        return l

def parsegff(filename):
    genome = Genome()
    for line in filename.readlines():
        if not line.startswith('#'):
            if line.startswith('>'): #if the fasta sequence is appended to the gff file, like prokka does e.g,
                break
            seqid,source,gentype,start,end,score,strand,phase,attributes = line.strip().split('\t')
            g = Gene(seqid,source,gentype,start,end,score,strand,phase,attributes)
            if g.seqid in genome.chromosomes:
                genome.chromosomes[g.seqid].append(g)
            else:
                genome.chromosomes[g.seqid] = [g]
    return genome

if __name__ == "__main__":
    fileB = open("/home/bardya/usr/data/15_03/a.equi/quast_results/results_2015_07_22_15_51_36/predicted_genes/polished_assembly_A_equi_genemark_genes.gff",'r')
    fileA = open("/home/bardya/usr/data/15_03/a.equi/prokka/A.equi.gff",'r')
    genomeA = parsegff(fileA)
    genomeB = parsegff(fileB)
    identicalgenes = []
    samestartpos = []
    samestartpos_careful = []
    sameendpos = []
    sameendpos_careful = []
    notfound = []
    for chromA,genelistA in genomeA.chromosomes.items():
        for geneA in genelistA: #for every gene of each chromosome of genomeA
            for chromB,genelistB in genomeB.chromosomes.items():
                for geneB in genelistB:
                    if geneA.start == geneB.start and geneA.end == geneB.end:
                        identicalgenes.append(geneA)
                    elif geneA.start == geneB.start:
                        if abs(geneA.length - geneB.length) <= (min(geneB.length,geneA.length)*0.15):
                            samestartpos_careful.append(geneA)
                        else:
                            samestartpos.append(geneA)
                    elif geneA.end == geneB.end:
                        if abs(geneA.length - geneB.length) <= (min(geneB.length,geneA.length)*0.15):
                            sameendpos_careful.append(geneA)
                        else:
                            sameendpos.append(geneA)

    for chromA,genelistA in genomeA.chromosomes.items():
        for geneA in genelistA: #for every gene of each chromosome of genomeA
            if not geneA in identicalgenes:
                notfound.append(geneA)
    
    print("#Total Genes in Genome A: {}".format(genomeA.numberofgenes()))
    print("""
  #Identical genes: {}
  #Genes with same start pos only: {}
  #Genes with same end pos only: {}
  #Genes not identical: {}
  #Genes with same start pos and total length deviation less than 15% of shorter gene's length: {}
  #Genes with same end pos and total length deviation less than 15% of shorter gene's length: {}  
  """.format(len(identicalgenes),len(samestartpos),len(sameendpos),len(notfound), len(samestartpos_careful),len(sameendpos_careful)))
    print("\n".join([x.attributes for x in samestartpos]))
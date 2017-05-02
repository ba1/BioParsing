'''
Created on Jul 6, 2015

@author: bardya
'''

import re
import string


def clearup(s, chars):
    return re.sub('[%s]' % chars, '', s).upper()

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
     
if __name__ == "__main__":
    s = '''        1 ctcaaacaat tataaatttc tttttattgt tgagtattta ggtcttgcac cgcgataatt
       61 ttgaactgca taggtagaag ctcaatttta attgctcata tcataagcaa agatcgaaat
      121 aatatgcagt atttttaatg aaaaccactc taatggttat ttccaagcat aaaaaaacac
      181 ccccaacatc agttgagggt gtttagaata atgagctggc gatgacttac tctcacatgg
      241 gaaaccccac actaccatca gcgctaagag gtttcacttc tgagttcggg aagggatcag
      301 gtggttcact cttgctatgg tcgccagcac aactggtatg gatacttgca ttggtcttat
      361 tggttgccga tgcgttttcc aaatctttac agatgggctg attgagtctt actttaattt
      421 tgttcatttt agctaaacat ataactaaat caagttgctt tgcatattaa tgaatcgatt
      481 gatgccatat atacaactgc ttgggtgttg tatagtcaag cctcacgagc aattagtatt
      541 ggtcagcttc acatatcact atgcttccac atccaaccta tcaacgtcct agtctcgaac
      601 ggctctttag aggacataaa gtcctaggga aatcttatct tgaggtaggc ttcccgctta
      661 gatgctttca gcggttatcc cttccgaaca tagctacccg gcgatgcgac tggcgtcaca
      721 accggtacac cagaggttcg tccactctgg tcctctcgta ctaggagcag atcctctcaa
      781 atttccagcg cccacggtag atagggaccg aactgtctca cgacgttcta aacccagctc
      841 gcgtacctct ttaaatggcg aacagccata cccttgggac ctgcttcagc cccaggatga
      901 gatgagccga catcgaggtg ccaaacaccg ccgtcgatat gaactcttgg gcggtatcag
      961 cctgttatcc ccagagtacc ttttatccgt tgagcgatgg cccttccata cagaaccacc
     1021 ggatcactaa gacctacttt cgtacctgct cgacttgtgg gtctcgcagt taagcgcgct
     1081 tttgccttta tactctacgc gtgatttccg accacgctga gcgcaccttc gtactcctcc
     1141 gttactcttt aggaggagac cgccccagtc aaactaccca ccagacatgg tcctcgcccc
     1201 ggattacggg gcagagttag aacctcaaca ttaccagggt ggtatttcaa ggacggctcc
     1261 attggaacta gcgttccaac ttcaaagcct cccacctatc ctacacaagt aaggtcaaag
     1321 ttcaatgtca agctgcagta aaggttcacg gggtctttcc gtctagccgc gggtacactg
     1381 catcttcaca gcgatttcga tttcactgag cctctgctgg agacagcgcc gccatcatta
     1441 tgccattcgt gcaggtcgga acttacccga caaggaattt cgctacctta ggaccgttat
     1501 agttacggcc gccgtttact ggggcttcga tcaagagctt cgcttacgct aaccccatca
     1561 attaaccttc cagcaccggg caggcatcac accctatacg tccactttcg tgtttgcaga
     1621 gtgctatgtt tttaataaac agttgcagcg gcctggtttc tgtggctgtc atcagctcag
     1681 gaagcaagtt ccatcaccaa caacagcgta ccttctcccg aagttacggt accattttgc
     1741 ctagttcctt cagcagagtt ctctcaagcg ccttggtcta ctcgacctga ccacctgtgt
     1801 cggtttaggg tacgattcct gtgtaactga agcttagaga cttttcctgg aagcatggta
     1861 tcagccactt cgctgtacaa gtacagcttg ctatcagatc tcagcataga gcaccccgga
     1921 tttgcctaag atgcatgcct acttcctttc acctggacaa ccaacgccag gctgacttaa
     1981 ccttctccgt cctctcatcg cattacacag aagtattgga atattaacca atttcccatc
     2041 gactacgcct ttcggcctcg ccttaggggt cgactcaccc agccccgatt aacgttggac
     2101 tggaaccctt ggtctttcgg cgaacgggtt tttcacccgt tttgtcgtta ctcacgtcag
     2161 cattcgcact tctgatacct ccagcatact tctcaataca ccttcatcgg cttacagaac
     2221 gctcccctac cacttgactt aaaagtcaaa tccgcagctt cggcacatag ttttagcccc
     2281 gttacatctt ccgcgcaggc cgactcgact agtgagctat tacgctttct ttaaagggtg
     2341 gctgcttcta agccaacctc ctagctgtct atgccttccc acatcgtttc ccacttaact
     2401 atgattttgg ggccttagct ggcggtctgg attgttttcc tcttgactac ggacgttagc
     2461 acccgcagtc tgtctcccgg atagtactca taggtattcg gagtttgcat cggtttggta
     2521 agtcgggatg accccctagc cgaaacagtg ctctaccccc tatggtattc gtccgaggcg
     2581 ctacctaaat agctttcggg gagaaccagc tatcaccagg cttgattagc ctttcacccc
     2641 tatccacaag tcatcccctg gcttttcaac gacagtgggt tcggtcctcc agttagtgtt
     2701 acccaacctt caacctgctc atggatagat cgcctggttt cgggtctata cccagcaact
     2761 aaacgcccta ttaagactcg atttctctac ggctccccta tacggttaac ctcgctactg
     2821 aatataagtc gctgacccat tatacaaaag gtacgcagtc accgaacaag tcggctccca
     2881 ctgcttgtat gcatgcggtt tcaggatcta tttcactccc ctcacagggg ttcttttcgc
     2941 ctttccctca cggtactggt tcactatcgg tcagtcagga gtatttagcc ttggaggatg
     3001 gtccccccat attcagacaa ggtttcacgt gcctcgccct actcgtcatc attatgtgtg
     3061 ccctttcgtg tacgggaata tcaccctcta cgttcgcact tcccagagcg ttccactaaa
     3121 acacacataa cttaatgggc tgatccccgt tcgctcgccg ctactaaggg aatctcaatt
     3181 gatttctttt cctaagggta ctgagatgtt tcacttcccc tcgttcgcct tgcaacacta
     3241 tgtattcatg ttgcaatacc taccttaaag taggtgggtt cccccattca gaaatctccg
     3301 gatcaaagga tatttgccgc ctccccggag cttatcgcag gctattacgt ctttcatcgc
     3361 ctctgactgc caaggcatcc accacatgca cttaattact tgactataca accccaaaca
     3421 gtcgttatta cctacaagga gtaataagac atgatcaatg attcctcatc aatctcttac
     3481 agtttgaagt actgtgtatt taaacactgt acagcttcaa tctaattcat ataccaaaac
     3541 gcttgattca gttaatttgc tagttctcaa tgaattccaa gattagaatt acttcttctc
     3601 ttttcattat tgagtgaaca atttatttca gactcaattt tgccaatctg ttaatgaata
     3661 aacatgcctt cgtcaggtca tgcttaatac cgtgatactt aaatcacaga agttaataaa
     3721 ccaagatcta aatctttatt tactaatttc tgtaatccga acttctctta agttctggtg
     3781 gagactagga gagtcgaact cctgacctcc tgcgtgcaaa gcaggcgctc taccaactaa
     3841 gctaagtccc cagcttacat catcagttat gtatctttct gtctataacc ttaatcagtc
     3901 aaagtcatgg tgggtctgac aagacttgaa cttgtgaccc cacgcttatc aagcgtgtgc
     3961 tctaaccaac tgagctacag accctcagat acatctatga agaacaactt gttgtggatt
     4021 cttaccaatc gtcaatcttt cgttaaggag gtgatccagc cgcaggttcc cctacggcta
     4081 ccttgttacg acttcacccc agtcatcggc cacaccgtgg taagcgtcct ccttgcggtt
     4141 agactaccta cttctggtgc aacaaactcc catggtgtga cgggcggtgt gtacaaggcc
     4201 cgggaacgta ttcaccgcgg cattctgatc cgcgattact agcgattccg acttcatgga
     4261 gtcgagttgc agactccaat ccggactacg atcggctttt tgagattagc atcctatcgc
     4321 taggtagcaa ccctttgtac cgaccattgt agcacgtgtg tagccctggc cgtaagggcc
     4381 atgatgactt gacgtcgtcc ccgccttcct ccagtttgtc actggcagta tccttaaagt
     4441 tcccgacatt actcgctggc aaataaggaa aagggttgcg ctcgttgcgg gacttaaccc
     4501 aacatctcac gacacgagct gacgacagcc atgcagcacc tgtatgtaag ttcccgaagg
     4561 caccaatcca tctctggaaa gttcttacta tgtcaaggcc aggtaaggtt cttcgcgttg
     4621 catcgaatta aaccacatgc tccaccgctt gtgcgggccc ccgtcaattc atttgagttt
     4681 tagtcttgcg accgtactcc ccaggcggtc tacttatcgc gttagctgcg ccactaaagc
     4741 ctcaaaggcc ccaacggcta gtagacatcg tttacggcat ggactaccag ggtatctaat
     4801 cctgtttgct ccccatgctt tcgcacctca gcgtcagtgt taggccagat ggctgccttc
     4861 gccatcggta ttcctccaga tctctacgca tttcaccgct acacctggaa ttctaccatc
     4921 ctctcccaca ctctagctaa ccagtatcga atgcaattcc caagttaagc tcggggattt
     4981 cacatttgac ttaattagcc gcctacgcgc gctttacgcc cagtaaatcc gattaacgct
     5041 tgcaccctct gtattaccgc ggctgctggc acagagttag ccggtgctta ttctgcgagt
     5101 aacgtccact atctgaaggt attaacttca gtagcctcct cctcgcttaa agtgctttac
     5161 aaccataagg ccttcttcac acacgcggca tggctggatc aggcttgcgc ccattgtcca
     5221 atattcccca ctgctgcctc ccgtaggagt ctgggccgtg tctcagtccc agtgtggcgg
     5281 atcatcctct cagacccgct acagatcgtc gccttggtag gcctttaccc caccaactag
     5341 ctaatccgac ttaggctcat ctattagcgc aaggtccgaa gatcccctgc tttctcccgt
     5401 aggacgtatg cggtattagc attcctttcg aaatgttgtc ccccactaat aggcagattc
     5461 ctaagcatta ctcacccgtc cgccgctaag tgatagtgca agcaccatca ctccgctcga
     5521 cttgcatgtg ttaagcctgc cgccagcgtt caatctgagc catgatcaaa ctcttcagtt
     5581 aaaaatcatt ttgcacctta ttaaagacaa ggtgccaatt ctggctcatc aattttctga
     5641 cttaaatttc gctcaaataa acttcgagta atttaaacca atcaatcaat gataatattt
     5701 cgatcaatca atcagtaaaa atccacacaa gttgttcttc ataatctctt aatgattatg
     5761 aatgacggtg gttatggcgt aatgcgtggc attcaaaata actactttgg tggtcgtcag
     5821 tatttcaatg aattgcatac accagactac aaacttcttg gcgaatctat gggtgttaag
     5881 agctggaaag tcggtagtgc agatgaattt aaaactgtta ttcaag'''
    print(fasta_linebreak(clearup(s, string.punctuation + string.digits + string.whitespace.replace("\n",""))))
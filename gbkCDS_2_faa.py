
seqs = []


def read_fasta(filename):
    flag = False
    hub = ''
    with open(filename) as file:
        for line in file:
            if line.startswith("     CDS") and not flag:
                hub += line
                flag = True
            elif not line.startswith("     CDS") and flag and line.startswith("                     "):
                hub += line
                flag = True
            elif line.startswith("     CDS") and flag:
                seqs.append(hub)
                hub = line
                flag = True
            elif not line.startswith("     CDS") and not line.startswith("                     ") and flag:
                seqs.append(hub)
                hub = ""
                flag = False
            elif not line.startswith("     CDS") and not line.startswith("                     ") and not flag:
                continue

read_fasta("/home/bardya/usr/data/15_03/goettig_test/Acinetobacter_calcoaceticus_P5077_2117/raw_data/A._calcoaceticus_P5077_2117.gbk")

def cdsentry2fasta(cdsentry):
    translation = "".join(cdsentry.split("/translation")[1].strip().split("/")[0].replace('"','').replace("=","").split())
    product = cdsentry.split("/product")[1].strip().replace('"','').replace("=","").replace("\n                     "," ").replace("/","| ")
    return translation, product

fastastring = ""

for seq in seqs:
    transl, prod = cdsentry2fasta(seq)
    fastastring += ">{}\n{}\n".format(prod,transl)

fastastring_formatted = ""
for line in fastastring.split("\n"):
    if not line.startswith(">"):
        n = 70
        linesplit = [line[i:i+n] for i in range(0, len(line), n)]
        formattedline = "\n".join(linesplit) + "\n"
        fastastring_formatted += formattedline
    else:
        fastastring_formatted += line + "\n"

output = open("/home/bardya/usr/data/15_03/goettig_test/Acinetobacter_calcoaceticus_P5077_2117/raw_data/A._calcoaceticus_P5077_2117.gbk.faa", 'w')
output.write(fastastring_formatted)
output.close()
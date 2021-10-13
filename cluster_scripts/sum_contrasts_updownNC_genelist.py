##script used to combine multiple contrast files into a matrix
import os, sys

inp = open(sys.argv[1], "r")
FCind = int(sys.argv[2]) #indice where Fold change is
pind = int(sys.argv[3]) #indice where pvalue is
output1 = open(str(sys.argv[1])+"_upclassfile.txt","w")
output2 = open(str(sys.argv[1])+"_dwnclassfile.txt","w")
output3 = open(str(sys.argv[1])+"_neg.txt", "w")
output4 = open(str(sys.argv[1])+"_up.txt", "w")
output5 = open(str(sys.argv[1])+"_dn.txt", "w")

def add_data_to_dict(inp,Dpos,Dneg,FCind,pind):
    header= inp.readline()
    for line in inp:
            L = line.strip().split("\t")
            gene = L[0].split(".")[0]
            FC = L[FCind]
            adj_pvalue = L[pind]
            if adj_pvalue == 'NA':
                pass
            else:
                if float(adj_pvalue) <=0.05:
                    if float(FC) >= 1:
                        if gene not in Dpos:
                            Dpos[gene] = 1
                        else:
                            print (gene, "gene listed twice")
                    elif float(FC) <= -1:
                        if gene not in Dneg:
                            Dneg[gene] = 1
                        else:
                            print (gene, "gene listed twice")
                else:
                    if float(FC) <= 0.1:
                        if float(FC) >= -0.1:
                            if gene not in Dpos:
                                Dpos[gene] = 0
                            if gene not in Dneg:
                                Dneg[gene] = 0
                
Dup = {}
Ddn = {}

add_data_to_dict(inp,Dup,Ddn,FCind,pind)
inp.close()
print (Dup,Ddn)
j=0
k=0
for i in Dup.values():
    if i == 0:
        j=j+1
    if i == 1:
        k=k+1
print (j, "number of genes in up-reg, neg class")
print (k, "number of genes in up-reg, pos class")
j=0
k=0
for i in Ddn.values():
    if i == 0:
        j=j+1
    if i == 1:
        k=k+1
print (j, "number of genes in dn-reg, neg class")
print (k, "number of genes in dn-reg, pos class")

#write class data to each gene
output1.write("Gene\tClass\n")
for gene in Dup:
    data= Dup[gene]
    #for data in data_list:
    string= str(data)
    output1.write(gene + "\t" + "%s" % string + "\n")
    if data == 0:
        output3.write(gene + "\n")
    elif data == 1:
        output4.write(gene + "\n")
    else:
        pass
output1.close()
output3.close()
output4.close()
output2.write("Gene\tClass\n")
for gene in Ddn:
    data= Ddn[gene]
    #for data in data_list:
    string= str(data)
    output2.write(gene + "\t" + "%s" % string + "\n")
    if data == 1:
        output5.write(gene + "\n")
    else:
        pass
output2.close()
output5.close()
#file.close()
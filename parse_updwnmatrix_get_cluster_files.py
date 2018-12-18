#parse a matix file (up, down, or NC) to get clusters for each condition and negative cluster
#write each positive cluster and negative cluster to its own file
import sys
inp_matrix= open(sys.argv[1], 'r')
header = inp_matrix.readline()
def all_same(items):
    return all(x == items[0] for x in items)

def get_matrix_dict(inp, D):
    nlist= []
    for line in inp:
        if line.startswith("AT"):
            L= line.strip().split('\t')
            gene= L[0]
            data= L[1:]
            if gene not in D:
                D[gene]= data
            else:
                print (gene, "already in dict")
            if all_same(data) == True:
                nlist.append(gene)
            else:
                pass
    return D, nlist
    
D= {}
gene_dict, neg_list = get_matrix_dict(inp_matrix, D)
print (gene_dict, "number of negative genes: ", len(neg_list))
oup3= open("negative_gene_cluster.txt", "w")
oup3.write("Gene\tClass\n")
for gene in neg_list:
    oup3.write("%s\t0\n" % gene)
oup3.close()
inp_matrix.close()
#get up/down cluster files
header = header.strip().split("\t")[1:]
for i in header:  
    x= header.index(i)
    up_list=[]
    dwn_list=[]
    for gene in gene_dict:
        lista = gene_dict[gene]
        value= lista[x]
        try:
            value= int(value)
            if value == 1:
                up_list.append(gene)
            elif value == -1:
                dwn_list.append(gene)
            else:
                pass
        except ValueError:
            print ("Could not convert data to an integer.")
            if value.lower() == 'up':
                up_list.append(gene)
            elif value.lower() == 'down':
                dwn_list.append(gene)
            elif value.lower() == 'dwn':
                dwn_list.append(gene)
            else:
                pass
                
    oup1= open(str(i)+"_cluster_up.txt", "w")
    oup1.write("Gene\tClass\n")
    oup2= open(str(i)+"_cluster_dwn.txt", "w")
    oup2.write("Gene\tClass\n")
    for gene in up_list:
        oup1.write("%s\t1\n" % gene)
    for gene in dwn_list:
        oup2.write("%s\t1\n" % gene)
    oup1.close()
    oup2.close() 
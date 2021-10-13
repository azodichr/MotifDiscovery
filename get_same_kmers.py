# get overlap of 2 kmer lists

import os, sys

dir1= sys.argv[1]
#kmers= open(sys.argv[2],"r")
out= open("all_clusters_cre_rank.txt","w")

# def get_kmers(kmers, lista):
#     for line in kmers:
#         L=line.strip().split('\t')
#         k=L[0]
#         lista.append(k)
#     return lista

def add_to_dict(inp, D, curr_list, count):
    header= inp.readline()
    for line in inp:
        L=line.strip().split('\t')
        k2=L[0]
        rank=L[1]
        if count != 0:
            if k2 not in curr_list:
                newlist=[]
                curr_list.append(k2)
                for i in range(count):
                    newlist.append("NA")
                if k2 not in D:
                    D[k2]=newlist+[str(rank)]
                else:
                    D[k2].append(str(rank))
            else:
                D[k2].append(str(rank))
        else:
            if k2 not in D:
                D[k2]=[str(rank)]

            
        
    return D, curr_list

def add_to_dict2(inp, D, curr_list, count):
    header= inp.readline()
    for line in inp:
        L=line.strip().split('\t')
        k2=L[0]
        rank=L[6]
        if count != 0:
            if k2 not in curr_list:
                newlist=[]
                curr_list.append(k2)
                for i in range(count):
                    newlist.append("NA")
                if k2 not in D:
                    D[k2]=newlist+[str(rank)]
                else:
                    D[k2].append(str(rank))
            else:
                D[k2].append(str(rank))
        else:
            if k2 not in D:
                D[k2]=[str(rank)]

            
        
    return D, curr_list

kmer_list=[]
#final_list= get_kmers(kmers, kmer_list)
#kmers.close()
kdict={}
title_list=[]
current_k_list=[]
count=0
for file in os.listdir(dir1):
    if file.endswith("_imp_avgrank_RF.txt") or file.endswith("_imp"):
        name = file.strip().split("_imp")[0]
        title_list.append(name)
        inp = open(dir1 + "/" + file)
        kdict, current_k_list= add_to_dict(inp, kdict, current_k_list, count)
        inp.close()
        count= count+1
    if file.endswith("_imp_scaled.txt"):
        name = file.strip().split("_imp")[0]
        title_list.append(name)
        inp = open(dir1 + "/" + file)
        kdict, current_k_list= add_to_dict2(inp, kdict, current_k_list, count)
        inp.close()
        count= count+1

print(kdict)
titlestr= "\t".join(title_list)
out.write("kmer\t%s\n" % (titlestr))
for key in kdict:
    data = kdict[key]
    if(len(set(data))==1):
        pass
    else:
        out.write("%s\t" % (key))
        for d in data:
            out.write("%s\t" % (d))
        out.write("\n")
out.close()
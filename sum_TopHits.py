##script used to combine TopHits
#python sum_contrasts.py <start directory> <output file name>
import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
start_dir = sys.argv[1]
ending= "_TopHits.txt.txt" #string that file ends with
sum_matrix = open(sys.argv[2],"w") #output
TYPE= str(sys.argv[3])

def add_data_to_dict(inp, name_dict, k_dict, count, TYPE, name):
    header= inp.readline()
    if name not in name_dict:
        name_dict[name]={}
    if TYPE== "cisbp":
        for line in inp:
            L=line.strip().split('\t')
            if len(L) > 2:
                k2=L[0]
                rank=L[1]
                cisBP=L[3]
                PCC=L[5]
                fam=L[9]
                if k2 not in k_dict:
                    k_dict[k2]=[str(cisBP), float(PCC), str(fam)]

                if k2 not in name_dict[name]:
                    k1= Seq(k2)
                    if str(k1.reverse_complement()) not in name_dict[name]:
                        name_dict[name][k2]=[float(rank)]
            else:
                k2=L[0]
                rank=L[1]
                if k2 not in name_dict[name]:
                    name_dict[name][k2]=[float(rank)]
                    
        dict1= name_dict[name] #get dictionary associated with name
        dict2= {key: rank for rank, key in enumerate(sorted(dict1, key=dict1.get), 1)} #convert to rank
        name_dict[name]= dict2 #associate name with new ranked dictionary
        
        
                
    elif TYPE == "dap":
        for line in inp:
            L=line.strip().split('\t')
            if len(L) > 2:
                k2=L[0]
                rank=L[1]
                cisBP=L[2]
                PCC=L[4]
                fam=L[5]
                if k2 not in k_dict:
                    k_dict[k2]=[str(cisBP), float(PCC), str(fam)]

                if k2 not in name_dict[name]:
                    k1= Seq(k2)
                    if str(k1.reverse_complement()) not in name_dict[name]:
                        name_dict[name][k2]=[float(rank)]
                        
            else:
                k2=L[0]
                rank=L[1]
                if k2 not in name_dict[name]:
                    name_dict[name][k2]=[float(rank)]
                    
        dict1= name_dict[name] #get dictionary associated with name
        dict2= {key: rank for rank, key in enumerate(sorted(dict1, key=dict1.get), 1)} #convert to rank
        name_dict[name]= dict2 #associate name with new ranked dictionary
    
    elif TYPE == "both":
        for line in inp:
            L=line.strip().split('\t')
            if len(L) > 7:
                k2=L[0]
                rank=L[1]
                cisBP=L[7]
                PCC=L[9]
                fam=L[10]
                cisBP2=L[12]
                PCC2=L[14]
                fam2=L[18]
                if k2 not in k_dict:
                    k_dict[k2]=[str(cisBP), float(PCC), str(fam), str(cisBP2), float(PCC2), str(fam2)]

                if k2 not in name_dict[name]:
                    k1= Seq(k2)
                    if str(k1.reverse_complement()) not in name_dict[name]:
                        name_dict[name][k2]=[float(rank)]
                        
            else:
                k2=L[0]
                rank=L[1]
                if k2 not in name_dict[name]:
                    name_dict[name][k2]=[float(rank)]
                    
        dict1= name_dict[name] #get dictionary associated with name
        dict2= {key: rank for rank, key in enumerate(sorted(dict1, key=dict1.get, reverse=True), 1)} #convert to rank
        name_dict[name]= dict2 #associate name with new ranked dictionary
    
    
    else:
        print("need TYPE, dap or cisbp")
        sys.exit(1)
    
    return (name_dict, k_dict)

#loop through directory for each file to add input and each filename
namedict={}
current_k_dict={}
count=0
title_list = []
titlelist1= []
dir2 = start_dir + "/"

for file in os.listdir(dir2):
    if file.endswith(str(ending)):
        name = file.strip().split(ending)[0]
        print (name)
        titlelist1.append(name)
        inp = open(dir2 + "/" + file)
        namedict, current_k_dict = add_data_to_dict(inp, namedict, current_k_dict, count, TYPE, name)
        inp.close()
        count= count+1
        inp.close()
###################################################################
print(current_k_dict)
print(namedict)

k_list= list(current_k_dict.keys())
print(len(k_list))
for k1 in k_list:
    k1= Seq(k1)
    #find and remove reverse compliment
    if str(k1.reverse_complement()) in k_list:
        if k1 != k1.reverse_complement():
            k_list.remove(str(k1.reverse_complement()))
print(len(k_list))

print(titlelist1)
title_str= "\t".join(titlelist1)

#write heading for gene and each filename
if TYPE == "both":
    sum_matrix.write("Motif\tDAP\tPCC\tDAPfam\tcisBP\tPCC\tcisBPfam\t%s\n" % title_str)
    na_str= "NA\t"
    #write data to each gene
    for k in k_list:
        k1= Seq(k)
        sum_matrix.write(k + "\t")
        kdata=current_k_dict[k]
        for x in kdata:
            data2= str(x)
            sum_matrix.write(data2 + "\t")
        for name in titlelist1:
            if name in namedict.keys():
                dict1= namedict[name]
                if k in dict1.keys():
                    data=dict1[k]
                    data2= str(data)
                    sum_matrix.write(data2 + "\t")
                elif str(k1.reverse_complement()) in dict1.keys():
                    data=dict1[str(k1.reverse_complement())]
                    data2= str(data)
                    sum_matrix.write(data2 + "\t")
                else:
                    sum_matrix.write(na_str)
        sum_matrix.write("\n")

else:
    sum_matrix.write("Motif\tcisBP\tPCC\tfam\t%s\n" % title_str)
    na_str= "NA\t"
    #write data to each gene
    for k in k_list:
        k1= Seq(k)
        sum_matrix.write(k + "\t")
        kdata=current_k_dict[k]
        for x in kdata:
            data2= str(x)
            sum_matrix.write(data2 + "\t")
        for name in titlelist1:
            if name in namedict.keys():
                dict1= namedict[name]
                if k in dict1.keys():
                    data=dict1[k]
                    data2= str(data)
                    sum_matrix.write(data2 + "\t")
                elif str(k1.reverse_complement()) in dict1.keys():
                    data=dict1[str(k1.reverse_complement())]
                    data2= str(data)
                    sum_matrix.write(data2 + "\t")
                else:
                    sum_matrix.write(na_str)
        sum_matrix.write("\n")

sum_matrix.close()

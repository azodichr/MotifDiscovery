# parse kmers10.gff with list of kmers to get kmer-specific gff for overlap

import os, sys

gff_file= open(sys.argv[1])
kmers= open(sys.argv[2])
out= open(str(sys.argv[2])+"_gff.txt","w")

kmer_list=[]
for line in kmers:
    L=line.strip().split('\t')
    k=L[0]
    kmer_list.append(k)
    
for line in gff_file:
    L2=line.strip().split('\t')
    #print(L2)
    if len(L2) > 3:
        k2= L2[2]
        if k2 in kmer_list:
            out.write(line)
        else:
            pass
    else:
        pass
        
gff_file.close()
kmers.close()
out.close()
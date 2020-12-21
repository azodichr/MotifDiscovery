import pandas as pd
import numpy as np
import sys, os
import timeit
from math import sqrt
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import fisher_exact
start = timeit.default_timer()

start_dir= sys.argv[1]
output= open(sys.argv[2], "w")

## function to get rank of importance score
def get_rank(inp):
	df = pd.read_csv(inp, sep="\t") # read in file as data frame
	df['pct_rank'] = df['mean_imp'].rank(pct=True) #get percentile rank of importance
	df= df.rename(columns={"Unnamed: 0": "kmer"})
	return df


## function to get kmers and importance. 
## Make sure for each cluster you are only getting 1 reverse complement
def get_kmer(df, name, D, klist1):
	klist=[]
	for i in range(len(df)) : 
		#print(df.loc[i, "kmer"], df.loc[i, "pct_rank"])
		km= df.loc[i, "kmer"]
		k1 = Seq(km)
		imp= df.loc[i, "pct_rank"]
		if str(k1) in klist1.keys():
			count= klist1[str(k1)]
			count= count+1
			klist1[str(k1)]=count
			if str(k1) in klist or str(k1.reverse_complement()) in klist:
				pass
			else:
				klist.append(km)
				if name not in D.keys():
					D[name]=[(km,imp)]
				else:
					D[name].append((km, imp))
		elif str(k1.reverse_complement()) in klist1.keys():
			k2= k1.reverse_complement()
			count= klist1[str(k2)]
			count= count+1
			klist1[str(k2)]=count
			if str(k2) in klist or str(k2.reverse_complement()) in klist:
				pass
			else:
				klist.append(str(k2))
				if name not in D.keys():
					D[name]=[(str(k2),imp)]
				else:
					D[name].append((str(k2), imp))
		else:
			klist1[km]=1
			if str(k1) in klist or str(k1.reverse_complement()) in klist:
				pass
			else:
				klist.append(km)
				if name not in D.keys():
					D[name]=[(km,imp)]
				else:
					D[name].append((km, imp))


	return D,klist1

## loop through files as input to get_kmer
kdict={}
klist1={}	
for file in os.listdir(start_dir):
	if file.endswith("_imp"):
		name=str(file)
		print(name)
		# get rank here
		df= get_rank(file)
		#print(df.head(4))
		# get kmers in a dict here
		#inp= open(file,"r")
		kdict,klist1= get_kmer(df, name, kdict, klist1)
		#inp.close()
		
print(kdict)
print("length of kmer list",len(klist1.keys()))

# sort dictionary based on imp values

kdict_keylist= sorted(kdict,key=lambda k: kdict[k][1], reverse=True)
print(kdict_keylist)



# sort kmer count dictionary based on count
klist2= sorted(klist1.items(), key=lambda x: x[1], reverse=True)
print(klist2)
kliststr= "\t".join([i[0] for i in klist2])
print(kliststr)
output.write("file\t%s\n" % kliststr)

# loop through dict and write output matrix

for key in kdict_keylist:
		output.write('%s\t' % key)
		tuplist = kdict[key]
		firstelements= [str(a_tuple[0]) for a_tuple in tuplist]
		secelements= [a_tuple[1] for a_tuple in tuplist]
		for k in klist2:
			km= k[0]
			if km in firstelements:
				ind= firstelements.index(km)
				val= secelements[ind]
				output.write('%.2f\t' % val)
			else:
				output.write('NA\t')
		
		output.write('\n')

# write count
output.write("count\t")
for k2 in klist2:
	countfin= k2[1]
	output.write('%s\t' % countfin)
output.write('\n')

output.close()
	
		
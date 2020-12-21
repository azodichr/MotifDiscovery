"""
Sums all result files from tomtom output
output is kmers in the same consensus sequence

Required input:
  -dir : path to tomtom_out folders
  
OPTIONAL:
  -pval : float point number that gives pval cutoff. Default is 0.05
  -eval : float point number that gives eval cutoff. Default is none
  
OUTPUT:
  fasta file with kmers that map to the same motif <motif_name.fa>
"""
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq

#DEFAULT
PVAL=0.05
EVAL="none"
for i in range(1,len(sys.argv),2):
	if sys.argv[i] == "-dir":
		startdir= sys.argv[i+1]
	if sys.argv[i] == "-pval":
		PVAL= sys.argv[i+1]
	if sys.argv[i] == "-eval":
		EVAL= sys.argv[i+1]
		
def get_kmer(inp, D, PVAL, EVAL):
	header= inp.readline()
	for line in inp:
		L=line.strip().split("\t")
		if len(L)==10:
			q= L[0]
			t= L[1]
			pv=L[3]
			ev=L[4]
			o=L[9]
			if o == "-":
				k1=Seq(t)
				t=str(k1.reverse_complement())
			else:
				pass
			if EVAL=="none":
				if float(pv) < PVAL:
					if t not in D.keys():
						D[t]=[q]
					else:
						if q not in D[t]:
							D[t].append(q)
				else:
					pass
			else:
				if float(ev) < EVAL:
					if t not in D.keys():
						D[t]=[q]
					else:
						if q not in D[t]:
							D[t].append(q)
				else:
					pass
	return D

motif_D={}
for dir in os.listdir(startdir):
	if dir.endswith("_tomtom_out"):
		startdir1=startdir + "/" +dir + "/"
		inp= open(startdir1+"tomtom.tsv","r")
		motif_D= get_kmer(inp, motif_D, PVAL, EVAL)
		inp.close()

#write fastas
print(motif_D)
for motif in motif_D:
	kmers= motif_D[motif]
	output= open(str(motif)+".fa","w")
	count=1
	for k in kmers:
		output.write(">%s_%d\n%s\n" %(motif,count,k))
		count=count+1
	output.close()

	
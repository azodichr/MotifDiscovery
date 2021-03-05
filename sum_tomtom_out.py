"""
Sums all result files from tomtom output
output is kmers in the same consensus sequence

Required input:
  -dir : path to tomtom_out folders
  
OPTIONAL:
  -pval : float point number that gives pval cutoff. Default is 0.05
  -eval : float point number that gives eval cutoff. Default is none
  -out: name of summary file if desired
  
OUTPUT:
  fasta file with kmers that map to the same motif <motif_name.fa>
  summary file with motif and all associated CREs
"""
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq

#DEFAULT
PVAL=0.05
EVAL="none"
OUT="NA"
for i in range(1,len(sys.argv),2):
	if sys.argv[i] == "-dir":
		startdir= sys.argv[i+1]
	if sys.argv[i] == "-pval":
		PVAL= sys.argv[i+1]
	if sys.argv[i] == "-eval":
		EVAL= sys.argv[i+1]
	if sys.argv[i] == "-out":
		OUT= sys.argv[i+1]
		
def get_kmer(inp, D, PVAL, EVAL, name, D2):
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
					if name not in D2.keys():
						D2[name]=[t]
					else:
						if t not in D2[name]:
							D2[name].append(t)
				else:
					pass
			else:
				if float(ev) < EVAL:
					if t not in D.keys():
						D[t]=[q]
					else:
						if q not in D[t]:
							D[t].append(q)
					if name not in D2.keys():
						D2[name]=[t]
					else:
						if t not in D2[name]:
							D2[name].append(t)
				else:
					pass
	return D,D2

motif_D={}
name_D={}
for dir in os.listdir(startdir):
	if dir.endswith("_tomtom_out"):
		name= str(dir).strip().split('_tomtom')[0]
		startdir1=startdir + "/" +dir + "/"
		inp= open(startdir1+"tomtom.tsv","r")
		motif_D,name_D= get_kmer(inp, motif_D, PVAL, EVAL, name, name_D)
		inp.close()

#write fastas
print(motif_D,name_D)
if OUT != "NA":
	outputM= open(OUT+"_motif-kmerlist.txt",'w')
	for motif in motif_D:
		kmers= motif_D[motif]
		kmerstr= "\t".join(kmers)
		outputM.write('%s\t%s\n' % (motif,kmerstr))
		output= open(str(motif)+".fa","w")
		count=1
		for k in kmers:
			output.write(">%s_%d\n%s\n" %(motif,count,k))
			count=count+1
		output.close()
	outputM.close()
	output2= open(OUT+"_cluster-motif_matrix.txt","w")
	motiflist= motif_D.keys()
	motifstr="\t".join(motiflist)
	output2.write('cluster\t%s\n' % motifstr)
	for n in sorted(name_D.keys()):
		output2.write('%s\t' % n)
		motifs=name_D[n]
		for m in motiflist:
			if m in motifs:
				output2.write('%d\t' % 1)
			else:
				output2.write('%d\t' % 0)
		output2.write('\n')
	output2.close()
				
else:
	for motif in motif_D:
		kmers= motif_D[motif]
		kmerstr= "\t".join(kmers)
		output= open(str(motif)+".fa","w")
		count=1
		for k in kmers:
			output.write(">%s_%d\n%s\n" %(motif,count,k))
			count=count+1
		output.close()

	
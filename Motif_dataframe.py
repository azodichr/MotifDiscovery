"""
PURPOSE:
Given a list of k-mers and fasta file of genes to consider, make a dataframe of presence or absense of each motif in the promoter region.
Before running, set path to Miniconda:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH


INPUT:
  -fasta        File with two columns (sep by '\t') Col 1 = Class Name; Col 2 = path to FASTA file
  -k            List of kmers to start with (/mnt/home/azodichr/ML_Python/6mers.txt or 5mers.txt)

OUTPUT:
  -fasta_k_df   Dataframe of presence (1) or absense (0) of k in fasta promoters.
"""

import pandas as pd
import numpy as np
import sys


for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-fasta":
        FASTA = sys.argv[i+1]
      if sys.argv[i] == '-k':
        km = sys.argv[i+1]


def make_DF(FASTA, km):
  from Bio import SeqIO
  from Bio.Seq import Seq
  
  numpy_header = ['Class']
  
  kmers = []
  for l in open(km, 'r'):
    kmers.append(l.strip("\n").strip("\r").split('\t')[0])
  try:
    kmers.remove('feature')
  except:
    print("feature is not in k-mer list - proceed")

  for i in kmers:
    numpy_header.append(i)  
  



  
  for line in open(FASTA,'r'):
    dataframe = np.zeros([1,len(kmers)+1])       # This fits the np df into the pd df - the plus 1 is for the Class!
    genes = ['Skip_this_line']  #index for pandas df
    info = line.strip().split("\t")
    name = info[0]
    file_path = info[1]
    fasta = open(file_path, 'r')
    
    for seq_record in SeqIO.parse(fasta, 'fasta'):
      header = seq_record.id
      genes.append(header)
      seq = str(seq_record.seq)
      gene_array =np.array([name])       # Array of P/A (1/0) for each gene - starts with '1' For Positive Class
      for ki in kmers:
        kmer = Seq(ki)
        if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
        else:
          gene_array = np.append(gene_array, 0)

      dataframe = np.vstack((dataframe,gene_array))
    

    pd_DF= pd.DataFrame(dataframe, index=genes, columns=numpy_header)  # Turn numpy into pandas DF
    pd_DF= pd_DF.drop("Skip_this_line",0)

    SAVE = file_path.split("/")[-1]+"_"+ km.split("/")[-1][:-21] + "_df.txt" 
    pd_DF.to_csv(SAVE, sep="\t")



make_DF(FASTA, km)


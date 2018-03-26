"""
PURPOSE:
Find k+mers enriched in your positive dataset and run RandomForest Classifier to determine how well those k+mers predict your classes

Before running, set path to Miniconda:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH

*When submitting jobs ask for 8 nodes! 

INPUT:
  -pos_file FASTA file with positive examples
  -neg_file FASTA file with negative examples - the script will run ML on 50 random samples to get balanced test
  -pos      String for what codes for the positive example (Default = 1)
  -neg      String for what codes for the negative example (Default = 0)
  -k        List of kmers to start with (/mnt/home/azodichr/ML_Python/6mers.txt or 5mers.txt)
  -pval     P-value cut off for Fisher's exact test (Default = 0.01)
  -n        Number of draws on negative data set you would like to run
  -save     Save name (will overwrite results files if the same as other names in the directory youre exporting to)
  -score    Default: F-measure. Can change to AUC-ROC using "-score roc_auc"
  -feat     Default: all (i.e. everything in the dataframe given). Can import txt file with list of features to keep.


OUTPUT:
  -SAVE_df_pPVAL.txt       Dataframe that goes into SK-learn for ML.
  -SAVE_FETresults.txt    Results for features enriched with fisher's exact test: Feature, # Pos Examples with Feature, # Neg examples with feature, Pvalue
  -SAVE_RF_results.txt    Results from RF runs
  -RESULTS.txt            Final results get added to this file: Run Name, # Features, # Reps (different Neg Datasets), CV, F_measure, StDev, SE
"""
import pandas as pd
import numpy as np
import sys, os



def Make_DF(POS, NEG, K, SAVE):
  from Bio import SeqIO
  from Bio.Seq import Seq
  from scipy.stats import fisher_exact

  #Put all kmers/kmer pairs into list
  km = []
  for l in open(K, 'r'):
    km.append(l.strip("\n"))

  numpy_header = ['Class']
  for i in km:
    numpy_header.append(i)  
  
  dataframe = np.zeros([1,len(km)+1])       # This fits the np df into the pd df - the plus 1 is for the Class!
  positive_present = {}.fromkeys(km, 0)     # Count occurence of each feature in positive examples
  negative_present = {}.fromkeys(km, 0)     # Count occurence of each feature in negative examples

  #Open positive and negative fasta files
  p = open(POS, 'r')
  n = open(NEG, 'r')

  num_pos = 0
  num_neg = 0
  genes = ['Skip_this_line']  #index for pandas df

  for seq_record in SeqIO.parse(p, 'fasta'):
    num_pos += 1
    header = seq_record.id
    genes.append(header)
    seq = str(seq_record.seq)
    gene_array =np.array([1])       # Array of P/A (1/0) for each gene - starts with '1' For Positive Class
    for ki in km:
      if " " in ki:                         #Checks to see if motif is a pair - pairs are separated by a space
        k1 = Seq(ki.split(" ")[0])
        k2 = Seq(ki.split(" ")[1])
        if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq: 
          gene_array = np.append(gene_array, 1)
          positive_present[ki] = positive_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)

      else:                                #If no separation by a space, assumes you're looking at singletons.
        kmer = Seq(ki)
        if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
          positive_present[ki] = positive_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)
    dataframe = np.vstack((dataframe,gene_array))
  
  for seq_record in SeqIO.parse(n, 'fasta'):
    num_neg += 1 
    header = seq_record.id
    genes.append(header)
    seq = str(seq_record.seq)         
    gene_array =np.array([0])               # Array of P/A (1/0) for each gene - starts with '0' For Negative Class
    for ki in km:
      if " " in ki:                         #Checks to see if motif is a pair - pairs are separated by a space
        k1 = Seq(ki.split(" ")[0])
        k2 = Seq(ki.split(" ")[1])
        if str(k1) in seq or str(k1.reverse_complement()) in seq and str(k2) in seq or str(k2.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
          negative_present[ki] = negative_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)

      else:                                 #If no separation by a space, assumes you're looking at singletons.
        kmer = Seq(ki)
        if str(kmer) in seq or str(kmer.reverse_complement()) in seq:
          gene_array = np.append(gene_array, 1)
          negative_present[ki] = negative_present[ki]+1
        else:
          gene_array = np.append(gene_array, 0)
    dataframe = np.vstack((dataframe,gene_array))

  DF= pd.DataFrame(dataframe, index=genes, columns=numpy_header, dtype=int)  # , dtype=int # Turn numpy into pandas DF
  DF= DF.drop("Skip_this_line",0)
  print(DF)
  out_name = SAVE + "_df.txt"
  DF.to_csv(out_name, sep = "\t")
  


if __name__ == '__main__':

  neg = '0'
  pos = '1'

  for i in range (1,len(sys.argv),2):

    if sys.argv[i] == '-neg':              #String for negative class : Default = 0
      NEG = sys.argv[i+1]
    if sys.argv[i] == "-pos":              #String for positive class : Default = 1
      POS = sys.argv[i+1]
    if sys.argv[i] == '-save':
      SAVE = sys.argv[i+1]
    if sys.argv[i] == '-k':
      K = sys.argv[i+1]
  
  Make_DF(POS, NEG, K, SAVE)


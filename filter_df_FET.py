"""
PURPOSE:
Take a machine learning dataframe and filter out features that aren't enriched in the positive class
Path to Miniconda:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH


INPUT:
  -df       Unfiltered df
  -pos      FASTA file with positive examples
  -neg      FASTA file with negative examples - the script will run ML on 50 random samples to get balanced test
  -pval     P-value cut off for Fisher's exact test (Default = 0.01)
  -FDR      Default: N. Designate (Y/N) if you want to run FDR correction during enrichment test


OUTPUT:
  -df_FET       Dataframe that goes into SK-learn for ML.

"""

import pandas as pd
import numpy as np
import sys, os
from math import sqrt
from scipy.stats import fisher_exact

PVAL = 0.01
neg = '0'
pos = '1'
cla = 'Class'


for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-pos':         #Fasta file for positive examples
    POS = sys.argv[i+1]
  if sys.argv[i] == '-neg':         #Fasta file for negative examples
    NEG = sys.argv[i+1]
  if sys.argv[i] == '-pval':             #Default is 0.01
    PVAL = float(sys.argv[i+1])
  if sys.argv[i] == '-df':
    DF = sys.argv[i+1]
  if sys.argv[i] == "-fdr":
    FDR = sys.argv[i+1]


df = pd.read_csv(DF, sep='\t', index_col = 0, header=0)
kmers = list(df)
kmers.remove(cla)

enriched = []

for k in kmers:
  #if k == 'ABF2':
    temp = df.groupby([cla, k]).size().reset_index(name="Count")
    try:
      TP = temp.loc[(temp[cla] == 1) & (temp[k] == 1), 'Count'].iloc[0]
    except:
      TP = 0
    try:
      TN = temp.loc[(temp[cla] == 0) & (temp[k] == 0), 'Count'].iloc[0]
    except:
      TN = 0
    try:
      FP = temp.loc[(temp[cla] == 0) & (temp[k] == 1), 'Count'].iloc[0]
    except:
      FP = 0
    try:
      FN = temp.loc[(temp[cla] == 1) & (temp[k] == 0), 'Count'].iloc[0]
    except:
      FN = 0

    #print(TP, FP, TN, FN)
    oddsratio,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater')
      
    if pvalue <= PVAL:
      enriched.append(k)
print(enriched)

exit()
  
if FDR == "N":
  for k in positive_present:
    try:
      count += 1
      TP = positive_present[k]            #Positive Examples with kmer present
      FP = negative_present[k]            #Negative Examples with kmer present
      TN = num_neg-negative_present[k]    #Negative Examples without kmer
      FN = num_pos-positive_present[k]    #Positive Examples without kmer

      oddsratio,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater')
      outFISH.write('\n%s\t%d\t%d\t%.7f' % (k, (positive_present[k]),(negative_present[k]),pvalue))
      if pvalue <= PVAL:          # Remove unenriched features from dataframe
        enriched_kmers[k] = pvalue
      if pvalue > PVAL:
        DF = DF.drop(k, 1)
      if count%10000==0:
        print("Completed " + str(count) + " features")

    except ValueError:
      count += 1 
      outFISH.write('\n%s\t%d\t%d\t1.0' % (k, (positive_present[k]),(negative_present[k])))

elif FDR == "Y":
  fdr_dict = {}
  for k in positive_present:
    try:
      count += 1
      TP = positive_present[k]            #Positive Examples with kmer present
      FP = negative_present[k]            #Negative Examples with kmer present
      TN = num_neg-negative_present[k]    #Negative Examples without kmer
      FN = num_pos-positive_present[k]    #Positive Examples without kmer

      oddsratio,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater')
      outFISH.write('\n%s\t%d\t%d\t%.7f' % (k, (positive_present[k]),(negative_present[k]),pvalue))
      pvalue_i = str(pvalue)
      fdr_dict[k] = pvalue_i

    except ValueError:
      count += 1 
      outFISH.write('\n%s\t%d\t%d\t1.0' % (k, (positive_present[k]),(negative_present[k])))
  
  outFISH.close()
  f_path = os.getcwd()+ "/" + SAVE + "_FETresults.txt"
  
  R=('R --vanilla --slave --args '+ os.getcwd() + " " + SAVE + " " + f_path+'< /mnt/home/azodichr/GitHub/MotifDiscovery/FDR.R')
  os.system(R)
  
  fdr_file = os.getcwd() +"/" + SAVE + "_FETresults_FDR.csv"

  for l in open(fdr_file,'r'):
    if "feature" in l:
          pass
    else:
      kmer, poscount, negcount, pval, adjp = l.strip().split(",")
      if float(adjp) <= PVAL:          #Remove unenriched features from dataframe
        enriched_kmers[kmer] = adjp
      elif float(adjp) > PVAL:
        DF = DF.drop(kmer, 1)
      if count%10000==0:
        print("Completed " + str(count) + " features")

else:
  print("Please include ''-FDR Y/N'' in command to designate if you want FDR correction")




if __name__ == '__main__':
  
  Make_DF(K, PVAL, SAVE)




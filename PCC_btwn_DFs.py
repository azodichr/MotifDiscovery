"""
PURPOSE:
Given two ML dataframes, find the PCC for the most similar feature from df2 
for every feature in df1


INPUT:
  -df1      
  -df2        
  -save  

OUTPUT:
  -save.txt

AUTHOR: Christina Azodi

REVISIONS:   Submitted 8/24/2017

"""
import pandas as pd
import numpy as np
import sys, os

DF2_KEEP = 'na'

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-df1":
    DF1 = sys.argv[i+1]
  if sys.argv[i] == "-df2":
    DF2 = sys.argv[i+1]
  if sys.argv[i] == "-df2_keep":
    DF2_KEEP = sys.argv[i+1]
  if sys.argv[i] == '-save':
    SAVE = sys.argv[i+1]

if len(sys.argv) <= 1:
  print(__doc__)
  exit()



# Get list of positive genes
df1 = pd.read_csv(DF1, sep='\t', index_col = 0)
df1_feats = list(df1)
df1_feats.remove('Class')
if DF2_KEEP == 'na':
  df2 = pd.read_csv(DF2, sep='\t', index_col = 0)
  df2_feats = list(df2)
  df2_feats.remove('Class')
else:
  keep_df2 = open(DF2_KEEP).readlines()
  keep_df2 = map(str.strip, keep_df2)
  print(keep_df2)
  df2 = pd.read_csv(DF2, sep='\t', index_col = 0, usecols=keep_df2)
  df2_feats = list(df2) 
print(df1_feats)
print(df2_feats)


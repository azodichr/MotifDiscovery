import sys
import numpy as np
import pandas as pd


for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-df":                # Df of pCRE criteria
    DF = sys.argv[i+1]
  if sys.argv[i] == "-key":               # Key with column name \t If best value is high or low
    KEY = sys.argv[i+1]


if len(sys.argv) <= 1:
  print(__doc__)
  exit()


# Load Criteria Dataframe

df = pd.read_csv(DF, sep='\t', header =0, index_col = 0)

# Get all the column names you want to rank
key = {}
with open(KEY, 'r') as f:
  for l in f:
    column, up_down = l.strip().split('\t')
    key[column] = up_down

# Adds any GO term columns since those will have different names depending on the cluster
for c in list(df):
  if "GO_" in c and "_%" in c:
    key[c] = 'high'

ranks = pd.DataFrame(index = df.index.values, columns = key.keys())
for col in key:
  if key[col] == 'low':
    r_low = df[col].rank(ascending = True)
    ranks[col] = r_low
  elif key[col] == 'high':
    r_high = df[col].rank(ascending = False)
    ranks[col] = r_high
  elif key[col] == 'binary':
    print('hello')
    r_high = df[col].rank(ascending = False)
    ranks[col] = r_high
  else:
    print("error with: " + col)


#print(df.head())
print(ranks.head())
name = DF + '_Rank'
ranks.to_csv(name, sep = '\t', na_rep = 'na')
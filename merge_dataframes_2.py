import sys
import pandas as pd


HOW = 'inner'
SAVE = 'test'

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-df":
    DF = sys.argv[i+1]
  if sys.argv[i] == "-key":
    KEY = sys.argv[i+1]
  if sys.argv[i] == "-how":
    HOW = sys.argv[i+1]
  if sys.argv[i] == "-save":
    SAVE = sys.argv[i+1]


key = {}
with open(KEY, 'r') as f:
  for l in f:
    prefix, file_path = l.strip().split('\t')
    key[prefix] = file_path


df = pd.read_csv(DF, sep='\t', header=0, index_col = 0)

for d in key:
  df_temp = pd.read_csv(key[d], sep='\t', header=0, index_col = 0)
  df_temp.drop('Class', 1, inplace=True) 
  df_temp.columns = [d + '_' + str(col) for col in df_temp.columns]
  df = pd.merge(df, df_temp, how = HOW, left_index = True, right_index=True)

#print(df.shape)
#print(df.head(20))

save = SAVE + "_df.txt"
df.to_csv(save, sep='\t')

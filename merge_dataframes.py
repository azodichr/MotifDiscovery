"""
PURPOSE:
Merge Machine Learning Dataframes together.
Adds a pre-fix to features added to the base data frame.

INPUTS:
  REQUIRED:
  -df     Base dataframe to use (no prefix added to these features)
  -key    Tab-delimited file with col1 = prefix and col2 = path_to_dataframe
  -save   Save name
  
  OPTIONAL:
  -how    How to merge (inner or outer merge) Default = inner
 
OUTPUT:
  -SAVE_df.txt

"""
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

print('Merging %i dataframes' % (len(key)+1))

df = pd.read_csv(DF, sep='\t', header=0, index_col = 0)
print('Start shape: %s' % str(df.shape))
for d in key:
  print('Merging %s' % d)
  df_temp = pd.read_csv(key[d], sep='\t', header=0, index_col = 0)
  df_temp.drop('Class', 1, inplace=True) 
  df_temp.columns = [d + '_' + str(col) for col in df_temp.columns]
  df = pd.merge(df, df_temp, how = HOW, left_index = True, right_index=True)

df = df.groupby(df.index).first()
#df['index'] = df.index
#df.drop_duplicates(subset=['index','Class'],inplace=True)
#df.set_index = df['index']
#df = df.drop('index', axis=1)

print('Final shape: %s' % str(df.shape))
#print(df.shape)
#print(df.head(20))

save = SAVE + "_df.txt"
df.to_csv(save, sep='\t')

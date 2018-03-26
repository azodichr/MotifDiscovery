""" 
Compare pCREs
"""

import sys, os
import numpy as np
import pandas as pd
from collections import defaultdict

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-path":
    PATH = sys.argv[i+1]
  if sys.argv[i] == "-string":
    STRING = sys.argv[i+1]

# Read all pCREs into a dictionary of lists
pCREs = defaultdict(list)

for f in os.listdir(PATH):
  if STRING in f and '_imp.csv' in f and 'RF' in f:
    name = f.strip('_imp.csv')
    with open(PATH + f) as file_:
      for l in file_:
        pCREs[name].append(l.strip().split(',')[0])

pCREs1 = pCREs.copy()

for data in pCREs:
  for data2 in pCREs1:
    if data != data2:
      overlap = list(set(pCREs1[data]) & set(pCREs1[data2]))
      length = len(overlap)
      print('%s:%s\t%i\t%s' % (data, data2, length, overlap))
  del pCREs1[data]





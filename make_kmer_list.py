import sys, os
import itertools
from Bio.Seq import Seq

#Makes list of all possible kmers
bases = ['A','T','G','C']
num = sys.argv[1]
km = [''.join(p) for p in itertools.product(bases, repeat=int(num))]

name = num + "mers.txt"
out = open(name, 'w')
l = []
final_list = []

for i in km:
  forw = str(Seq(i))
  rev = str(Seq(i).reverse_complement())
  l.append(forw)
  l.append(rev)
  
  if forw not in final_list and rev not in final_list:
    final_list.append(forw)
    out.write(forw + "\n")

print("Number of %s-mers with reverse complements removed: %i" % (num, len(final_list)))

pairs = ['_'.join(k) for k in itertools.combinations(final_list, 2)]

print("Number of %s-mer pair combinations with reverse complements removed: %i" % (num, len(pairs)))

pair_name = num + "mer_pairs.txt"
out2 = open(pair_name,'w')
for p in pairs:
  out2.write(p + "\n")
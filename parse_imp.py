"""
Sorts importance files output by RandomForest_v2.0 and related SciKit-learn ML scripts and allows for other selection.

Required input:
  -f : path to file or path to directory with multiple imp.txt files

Other options:
  -n : Gives top n most important features
  -p : Gives top percent p most important features
  -value : default = True, if False then don't print pvalue in output

"""
import os, sys
import operator
n = "n"
cutoff = "n"
p = "n"
value = "True"
f = "help"

for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-f':    #Path to imp.txt file or to directory with files
    f = sys.argv[i+1]
  if sys.argv[i] == '-n':    #Return the top n 
    n = int(sys.argv[i+1])
  if sys.argv[i] == '-cutoff':    #Return all features with imp over cutoff
    cutoff = sys.argv[i+1]
  if sys.argv[i] == '-p':    #Return the top p percent
    p = sys.argv[i+1]
  if sys.argv[i] == '-value':    #Return the top p percent
    value = sys.argv[i+1]


def sort(f):
  dic = {}
  for l in open(f, 'r'):
    kmer, val = l.strip().split("\t")
    dic[kmer] = float(val)
  sorted_dic = sorted(dic.items(), key=operator.itemgetter(1), reverse = True)

  if n == p == "n":
    name = f + "_sort"
    out = open(name, 'w')
    if value == "True":
      for i in sorted_dic:
        out.write("%s\t%s\n" % (i[0], i[1]))
    if value == "False" or value == "false" or value == "f":
      for i in sorted_dic:
        out.write("%s\n" % (i[0]))

  elif n != "n":
    name = f + "_top" + str(n)
    out = open(name, 'w')
    kmer_list = sorted_dic[0:n]
    if value == "True":
      for i in kmer_list:
        out.write("%s\t%s\n" % (i[0], i[1]))
    if value == "False" or value == "false" or value == "f":
      for i in kmer_list:
        out.write("%s\n" % (i[0]))

  elif p != "n":
    name = f + "_top" + str(p) + "perc"
    out = open(name, 'w')
    top = int(float(len(sorted_dic)) * float(p) * 0.01)
    kmer_list = sorted_dic[0:top]
    if value == "True":
      for i in kmer_list:
        out.write("%s\t%s\n" % (i[0], i[1]))
    if value == "False" or value == "false" or value == "f":
      for i in kmer_list:
        out.write("%s\n" % (i[0]))



if ".txt" in f:
  print("Parsing given file")
  sort(f)

else:
  print("Parsing all .imp files in directory")
  for j in os.listdir(f):
    if j.startswith(".") or not "_imp.txt" in j:
        pass
    else:
      print(j)
      sort(j)


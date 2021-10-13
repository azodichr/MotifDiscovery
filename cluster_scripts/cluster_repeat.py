import sys, os

f = sys.argv[1]
use = 'a'
use = sys.argv[2]

if use == 'a':
  clusters = ['UNN','UNU','NUN','NUU','NNU','DNN','DND','NDN','NDD','NND']
elif use == 'p':
  clusters = ['UNN','UNU','NUN','NNU','DNN','DND','NND']

with open(f, 'r') as file:
  lines = file.readlines()
lines = [item.strip() for item in lines]

out = open(f, 'w')

for c in clusters:
  for l in lines:
    out.write(l.replace('XXX', c))
    out.write('\n')



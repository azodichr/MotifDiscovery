#parse FET results to get kmer lists as input for overlap

import os, sys

start_dir = str(sys.argv[1])

for file in os.listdir(start_dir):
    if file.endswith("_imp_avgrank_RF.txt") or file.endswith("_imp_avgrank_SVM.txt") or file.endswith("_imp"):
        name= str(file).strip()
        save_file= name+"_CREonly.txt"
        inp= open(file, 'r')
        output= open(start_dir+'/'+save_file, 'w')
        header= inp.readline()
        for line in inp:
            L= line.strip().split('\t')
            cre= L[0]
            if cre.startswith("GSM") or cre.startswith("Coords") or cre.startswith("no"):
                pass
            else:
                output.write('%s\n' % (cre))
        output.close()
        inp.close()
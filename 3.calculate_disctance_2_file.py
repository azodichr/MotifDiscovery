

import sys, re, os, math, getopt, random
GLOBALS = {}

# Packages required for pccrange.
from  TAMO            import MotifTools
from  TAMO.util       import Arith
from  TAMO.MotifMetrics import ProbeSet

from TAMO.Clustering.MotifCompare import *
from TAMO.Clustering              import MotifCompare
from TAMO                         import MotifTools
from   TAMO.MotifTools import Motif, print_motifs


class MyWriter:

    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = file(filename, 'a')

    def write(self, text):
        self.stdout.write(text)
        self.logfile.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()
    def flush(self):
        self.stdout.flush()
        
    
def main():
    #os.system("%s-%s.dm" % (getarg('file1'),getarg('file2')))
    motifs1 = getarg('motifs1')
    motifs2 = getarg('motifs2')
    print len(motifs1),len(motifs2)
    n=0
    outname = getarg('file1').split('/')[-1]+'-'+getarg('file2').split('/')[-1]+'.dm'
    oup=open('/%s/%s' %('/'.join(getarg('file1').split('/')[:-1]),outname),"w")
    for i in range(len(motifs1)):      
        Dmat = computeDpairs(motifs2,[motifs1[i]],VERBOSE=0,DFUNC=DFUNC)
        oup.write('%s\n' % ("\t".join(map(lambda x: "%0.2f" % x,Dmat[0].values()))))
        n+=1
    oup.close()

# This function is not in the normal TAMO distribution on HPCC. I found it on
# calculon in /home/zou/lib/python2.5/site-packages/TAMO/Clustering/MotifCompare.py.
# My guess is that it is a custom written script by Cheng, because I cannot 
# find any evidence of the function in any other copy of TAMO or its documentation.

def pccrange(self,other,Srange,Orange):
    from numpy import array
    from scipy.stats import pearsonr
    #import statistics
    '''Utility function: compute diff of '''
    '''self and other from self+Sstart, width of other'''
    POW     = math.pow
    Dtot    = 0
    divisor = float(len(Srange))
    '''Computes distance'''
    pcc=0
    
    for si, oi in zip(Srange,Orange):
        Sl=[]
        Ol=[]
        for L in ACGT:
            Sl.append(POW(2,self.logP[si][L]))
            Ol.append(POW(2,other.logP[oi][L]))
        pcc+= pearsonr(Sl, Ol)[0]
    _dist=1-pcc/divisor
    return _dist

# This is another custom function written by Cheng and found in the same place
# as pccrange. Alex Seddon, 27 Jun 2012

def computeDpairs(motifs1,motifs2,VERBOSE=0,DFUNC=DFUNC):
    N = len(motifs1)
    N2 =len(motifs2)
    dmat = {}

    Nhalf = int(.293*N)
    for i in range(N2):
        dmat[i]={}
        if VERBOSE:
            if i == Nhalf: sys.stdout.write('|'); sys.stdout.flush()
            else:          sys.stdout.write('.'); sys.stdout.flush()
        A = motifs2[i]
        for j in range(N):
            B = motifs1[j]
            # minshortestoverhangdiff finds the allignment of two motifs that
            # minimizes the distance measurement specified. 
            D=minshortestoverhangdiff(A,B,5,DFUNC=DFUNC) 
            # Symmetric accross the diagonal of the matrix.
            dmat[i][j] = D              
    #End of Pretty status bar
    if VERBOSE: print                               
    return dmat

def parse_opts():
    global GLOBALS
    global DFUNC
    short_opts = 'i:j:'
    long_opts = ['dfunc=']
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts,long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()
    print opts
    GLOBALS['args'] = args
    GLOBALS['motifs1'] = []
    GLOBALS['motifs2'] = []

    DFUNCtxt = None
    for opt,value in opts:
        if opt == '-i':                   GLOBALS['motifs1'] = MotifTools.txt2motifs(value)
        if opt == '-i':                   GLOBALS['file1'] = value
        if opt == '-j':                   GLOBALS['motifs2'] = MotifTools.txt2motifs(value)
        if opt == '-j':                   GLOBALS['file2'] = value
        if opt == '--dfunc':              DFUNCtxt = value


     # Deal with DFUNC and DMAX
    if DFUNCtxt == 'pccrange':
        DFUNC = pccrange
    else:
        if DFUNCtxt == 'NCB':
            #_DFUNC = MotifCompare.negcommonbits
            _DFUNC = MotifCompare.negcommonbitsrange   
        else:
            try:
                exec ("_DFUNC = MotifCompare.%s"%DFUNCtxt)
            except:
                usage("No such distance metric: %s"%DFUNCtxt)
        if _DFUNC:  set_dfunc(_DFUNC)

def usage(txt=''):
    '''
    Place information about command-line arguments here
    '''
    if txt: print "Error: %s"%txt
    print 'Usage: %s -m motifs'%(
        re.sub('^.*/','',sys.argv[0]))
    print ""

    print '         [-i]        motif1'
    print '         [-j]        motif2'
    print '         [--dfunc <function>]  Distance metric.  Examples are NCB, and diffrange.  Any'
    print '                               "range" function in MotifCompare.py is acceptible'
    print ''
    print 'Some useful examples:'
    print ''
    print '   UPGMA.py -i <motiffile.tamo>  -j   --dfunc diffrange (default)'

    
    sys.exit(1)

       
def getarg(varname):
    global GLOBALS
    if GLOBALS.has_key(varname):   return GLOBALS[varname]
    else:                          return None


def set_dfunc(_dfunc,):
    '''
    set_dfunc(dfunc, dmax)

    Set the distance/divergence/difference metric and threshold, and synchronize
    between this module and MotifCompare

    Examples:
    set_dfunc(MotifCompare.negcommonbitsrange, -8.0)
    set_dfunc(MotifCompare.diffrange, 0.23)
    '''
    global DFUNC
    DFUNC = _dfunc
    MotifCompare.DFUNC = _dfunc



   
if __name__ == '__main__': 
    parse_opts()
    main()  

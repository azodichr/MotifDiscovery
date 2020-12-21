import os,sys
from TAMO            import MotifTools
from TAMO.MotifTools import Motif, save_motifs
import motility 
import math
from operator import itemgetter
import bisect
import time
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."
    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L

def score_range(start,end,threshold):
	# Start = 0, end = max_score
	global ATbias,GCbias
	
	d={} # Dictionary containing an index for each score in the range specified.
		 # The key for each index is the weight of all the sites above that 
		 # score. (weight = # sites above score/ # total sites)
		 
	for i in frange(start,end,(end-start)/float(10)):
		d[i]=pwm.weight_sites_over(i,ATbias,GCbias)
		
	if i < end:
		#print "i",i,end
		d[end]= pwm.weight_sites_over(end,ATbias,GCbias)

	dsort=sorted(d.items(), key=itemgetter(1)) # Sorted list of scores and
											   # weights as tuples.
	
	klist=map(lambda x: x[0],dsort) # List of scores in score range
	vlist=map(lambda x: x[1],dsort) # List of weights
	#print klist,vlist
	
	bi=bisect.bisect_left(vlist,threshold) # index with weights above the 
										   # threshold
	
	print bi
	for i in vlist:
		print i,threshold,i-threshold
	if bi >=len(vlist) or bi== 0:
		return None
	else:
		return klist[bi],klist[bi-1] # returns the scores that give weights at
					                 # and just below the threshold

def run_search(seq_file,pwm,score):
	
	# Search for the motif in every gene in the FASTA file.
	inp=open(seq_file,"r")
	line=inp.readline()
	while line:
		if line.startswith(">"):
			name=line[1:].strip()
			block=''
			line=inp.readline()
			
			# Block is the promoter sequence for one gene in the FASTA
			while not line.startswith(">") and line:
				block+=line.strip()
				line=inp.readline()
			#print block,score,name
			parse_block(block,pwm,score,name)
 

def parse_block(seq_ori,pwm,score,name):
	global oup
	threshold = score
	
	# Checks that there are no Ns in the sequence file.
	dic={"A":1,"T":1,"G":1,"C":1}
	seq_list=filter(dic.has_key,list(seq_ori))
	seq="".join(seq_list)
	if len(seq)!=len(seq_list):
		oup.write("wrong")
	#print len(pwm)
	
	if len(pwm) < len(seq):
	
		# Searches sequence using the PWM.
		matches = pwm.find(seq,threshold)
		
		'''
		print "\nfound %d match(es) to pwm in '%s' at %f:" % (len(matches), name,
                                                    threshold)
		for (start, end, o, match) in matches:
			print '\t%d to %d in %d orientation; match is %s' % (start, end, o, match)
		'''
		pdic={}
		
		# Start, end, o is for orientation, and match is the sequence
		for (start, end, o, match) in matches:
			#print "Matches",len(matches),threshold,name
			
			score = pwm.calc_score(match) # Get the score for the match
			
			if score<= pwm.max_score(): # If the score is <= max_score
				
				if pdic.has_key(score):
					oup.write('%s\t%d\t%d\t%d\t%s\t%s\t%s\n' % (name,start, end, o, match,score,pdic[score]))
				else:
					pdic[score] = pwm.weight_sites_over(score)
					oup.write('%s\t%d\t%d\t%d\t%s\t%s\t%s\n' % (name,start, end, o, match,score,pdic[score]))
				#print name,score,pdic[score]

			else:
				score= pwm.max_score()
				if pdic.has_key(score):
					oup.write('%s\t%d\t%d\t%d\t%s\t%s\t%s\n' % (name,start, end, o, match,score,pdic[score]))
				else:
					pdic[score] = pwm.weight_sites_over(score)
					oup.write('%s\t%d\t%d\t%d\t%s\t%s\t%s\n' % (name,start, end, o, match,score,pdic[score]))


if __name__=="__main__":
    
    ##
    # Parse the arguments for the mapping
	file=sys.argv[1] # TAMO sequence file
	threshold= math.pow(10,-float(sys.argv[2])) 
	maxthreshold = float(sys.argv[3]) # for strong score, using 0.9*max score
	ATbias=float(sys.argv[4]) # 0.33
	GCbias=float(sys.argv[5]) # 0.17
	seq_file=sys.argv[6] # FASTA file of the sequence
	tar_dir="" # Target directory for the output file
	#
	##
	
	##
	#
	for i in range(1,len(sys.argv)):
		if sys.argv[i]== "-d":
			tar_dir=sys.argv[i+1].rstrip("/")
	print tar_dir
	ml=MotifTools.txt2motifs(file)
	n=0
	new_list=[]
	
	##
	# Looks at each motif from the TAMO file. Uses the find function from 
	# motility to find the sequences with that motif.
	for Ikey in range(len(ml)):
		#print m.ll 
		

		time1=time.time()
		
		m=ml[Ikey] # Pull out the motif from the motif list.
		save_motifs([m],file+'_'+str(Ikey)) # Save the motif as a file.
		
		##
		# Create the PWM matrix from the TAMO motif.
		matrix=[]
		for i in range(len(m.ll)):
			position=[]
	 		for nt in ["A","C","G","T"]:
				position.append(m.ll[i][nt])
			matrix.append(position)
		#print matrix
		pwm = motility.PWM(matrix)
		#
		##
		
		##
        # Print the motif index and the max and min PWM score, and the IUPAC.
		print "START",Ikey
		score_max = pwm.max_score() # Maximum score for the motif.
		score_min = pwm.min_score() # Minimum score for the motif.
		start2=0
		stop2=0
		print score_min,score_max,len(m.oneletter)
		#
		##
		
		##
        # 
		if len(m.oneletter)<18: # Checks for motifs under 18 nt.
		
			tmp1=score_range(0,score_max,threshold) # Gives scores with weights
			                                        # at and just below the
			                                        # threshold.
			                                        
			if tmp1!=None:
			
				start1,stop1=tmp1 
				
				print "run1",start1,stop1
				tmp2=score_range(start1,stop1,threshold)
				if tmp2!=None:
					start2,stop2= tmp2
					print "run2",start2,stop2
			if start2!=0:
				score=(start2+stop2)/2
				filename=file+"_"+str(Ikey)+"_"+str(threshold)+".out.pvalue"
				if tar_dir=="":
					fo=filename
				else:
					fo=tar_dir+"/"+filename[filename.rfind("/")+1:]
				print "OUT1",Ikey,str(Ikey),filename,fo
				oup=open(fo,"w")
				oup.write("%s\t%s\t%s\n" % (Ikey,m.oneletter,len(m.oneletter)))
				run_search(seq_file,pwm,score)
				oup.close()
				
			else:
				score= score_max*maxthreshold
				filename=file+"_"+str(Ikey)+"_"+str(maxthreshold)+"max"+".out.pvalue"
				if tar_dir=="":
					fo=filename
					print "pre1",fo
				else:
					fo=tar_dir+"/"+filename[filename.rfind("/")+1:]
					print "pre",tar_dir,filename[filename.rfind("/")+1:],fo
				print "OUT2",Ikey,str(Ikey),filename,fo
				oup=open(fo,"w")
				oup.write("%s\t%s\t%s\n" % (Ikey,m.oneletter,len(m.oneletter)))
				run_search(seq_file,pwm,score)
				oup.close()
				#os.system("rm %s" % fo)
		else:
			score="too Long"

		print "time",time.time()-time1
		
		#oup.write("%s\n" % score)
		#break
		

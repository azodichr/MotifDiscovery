# MotifDiscovery
Original pipeline used python for processing dataframes and R for running Random Forest. For those instructions see bottom of the document (Python & R Pipeline).

New pipeline uses SciPy to determine enrichment, Numpy and Pandas for dataframe management, and SciKit Learn for RandomForest

## Updates to the Pipeline:
December 2020 : Added scripts to use MEME suite in kmer correlation instead of TAMO

October 15, 2018 : Added scripts to generate a tamo file and map to cis-bp

May 19 2016 : Added option to run multiple test correction during kmer enrichment. Add "-FDR Y" to your command line to do so.

April 27 2016 : Added MapEnrichedKmers.py to directory. This is the starting point for incorporating DAP-Seq data into the pipeline. Script is being actively developed by cbazodi and will likely be renamed!

March 6 2016 : Change RandomForest_v2.0.py to call the random forest script from RF_scikit.py. That way ML runs are consistent between methods.

Febrary 2 2016 : Add RF scikit.py. Which allows you to run RF ML on any dataframe given. Can select scoring type ('f1' = F-measure; 'roc_auc' = AUC-ROC  : default = f1) and give list of features to include (default is to use the whole dataframe)

Nov 25 2015 : Simplify and standardize outputs (gives stdev and SE), also changed the random_state generator in RandomForestClassifier from an int to the default which is np.random.

Nov 24 2015 : Make Fisher's Exact Test 1-Tailed (only looking at enrichment in positive class) & remove Training/Testing data split since using cross-validation

Nov 13 2015 : Alter script so that enriched kmers are lengthened by 1 bp until they are no longer enriched

Oct 26 2015 : Switch from Python+R Pipeline to running everything in Python, this included changing to include reverse complement information, and to run the ML using 20 sets of random negative example genes. 

## Get positive and negative examples from log FC data set

1. Get up-reg gene list (logFC >= 1, adj. p-val < 0.05), dwn-reg gene list (logFC <= -1, adj. p-val < 0.05), and negative gene list (logFC <= 0.1, adj. p-val > 0.05; logFC >= -0.1, adj. p-val > 0.05). 

File needs to have a header, genes in first column and a log FC and p-value. Must indicate what columns are logFC and p-value

    python ~/Github/MotifDiscovery/sum_contrasts_updownNC_genelist.py <logFC data with p-val> <indice with logFC> <indice with p-value>

2. Set Up Your Files:

Inside directory for Pairwise experiment make directories for FASTA files and Motif Lists:
mkdir FastaFiles
mkdir MotifLists

Put cluster file in FastaFiles dir and get promoter sequence:

    cd FastaFiles/
    cp [pos_examples]

    python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta [promoter sequences] -name [pos examples]

For arabidopsis you can use: /mnt/home/azodichr/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa

Put negative example file in FastaFiles dir and get promoter sequence.

    cp [neg_examples]
    python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta [promoter sequences] -name [neg examples]

## Get enriched kmers

3. Get enriched kmer dataframe using Fisher's Exact Test

       python ~/Github/MotifDiscovery/pCRE_Finding_FET.py -pos <pos fasta> -neg <neg fasta> -k ~/1-herb_CRE_project/motifs/6mer.txt -save <name of output files>
        
     optional:
      
        -pos_str  String for what codes for the positive example (Default = 1)
        -neg_str  String for what codes for the negative example (Default = 0)
        -k        List of kmers to start with (/mnt/home/azodichr/ML_Python/6mers.txt or 5mers.txt)
        -pval     P-value cut off for Fisher's exact test (Default = 0.01)
        -FDR      Default: N. Designate (Y/N) if you want to run FDR correction during enrichment test
        
     Output:
     
        -SAVE_df_pPVAL.txt       Dataframe that goes into SK-learn for ML
        
   submit as qsub file
   
        python ~/Github/parse_scripts/qsub_hpc.py -f submit -c FET.runcc -u john3784 -w 239 -m 12 -wd ~/4-Trichome_project/FastaFiles/

## Use ML Pipeline (most recent version) here to get class predictions (see ML-pipeline https://github.com/bmmoore43/ML-Pipeline)

# Sequence similarity between pCREs and CIS-BP/DAP-Seq motifs

## Using MEME-suite

1. Download and install MEME-suite

    1. Download: http://meme-suite.org/doc/download.html
    2. Installation instructions: http://meme-suite.org/doc/install.html
    
2. Installation on mac:

    1. use anaconda to install python 2.7
    
    Create python 2.7 environment
    
        conda create --name py27 python=2.7
        
    Activate environment
    
        source activate py27
        
    2. change to meme directory and check dependencies
    
           cd meme-5.2.0
           cd scripts/
           perl dependencies.pl
        
    3. install any/all depedencies. try
    
           sudo cpan <library>
        
    or 
        
         sudo port install <library>

    4. configure and make (in meme-5.2.0 folder)
    
           ./configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
        
           make test
        
           make install
        
    5. check installation in bash profile
    
           nano ~/.bash_profile
        
3. Convert previous TAMO files to meme files (this is the case for TFBM and DAPseq sites) and kmers to memes

    1. Add source/motif name to tamo file
    
           python ~/Desktop/post_doc/scripts/parse_tamo_get_meme.py <tamo file> <original index file>
           
           python ~/Desktop/post_doc/scripts/parse_tamo_get_meme.py Athaliana_TFBM_v1.01.tm.index.direct.index.tm.txt Athaliana_TFBM_v1.01.tm.index.direct.index
           
    2. convert tamo to meme
    
           tamo2meme <tamo file> > <meme output>
           
           /Users/Beth/Desktop/Github/meme-5.2.0/scripts/tamo2meme Athaliana_TFBM_v1.01.tm.index.direct.index.tm.tamo > Athaliana_TFBM_v1.01.tm.index.direct.index.tm.tamo_meme.txt
           
     3. convert kmers to memes using iupac2meme from meme suite
     
     no arguments, just run in folder where imp files are
     
            python ~/Desktop/post_doc/scripts/kmer_files_2_meme.py
            
     a meme file is created for each imp file
     
4. correlate kmer meme files to TFBM or DAPseq meme files using tomtom

        ~/Desktop/Github/meme-5.2.0/src/tomtom [options] <query motif file> <target motif file>+
        
        ~/Desktop/Github/meme-5.2.0/src/tomtom -o Sp_0.7_2091_mRNA5utr_tomtom_out Sp_brachy_clusters_0.7_2091.txt_v3.2.2.txt_mRNA5utr.fa_random1_df_p0.01.txt_nodups.txt_LogReg_imp_meme.txt Athaliana_TFBM_v1.01.tm.index.direct.index.tm.tamo_meme.txt
        
5. sum tomtom output to make consensus sequence from kmers that map to the same motif

   default uses a p-value of less than 0.01 as cutoff. Otherwise cutoff can be input as -pval or -qval

        python ~/Desktop/Github/MotifDiscovery/sum_tomtom_out.py -dir ./
     
   output is fasta file that can be used as input to ggseqlogo.R to form consensus

## Motif PCC distance using TAMO

1. generate TAMO file based on motif sequence
module load TAMO

    python generate_PWM.py [motif list]
    
    python generate_PWM.py kmers10.txt

2. command will divide the tamo files and also generate the command_line files(named “runcc”). Tamo 2 is the Athaliana file

        python /mnt/home/seddonal/scripts/5_motif_merging/pcc_merge_CC.py create_cc_2 -t [tamo_file_1] -t2 [tamo_file_2]

CISBP: 

    python pcc_merge_CC.py create_cc_2 -t kmers10.txt.tamo -t2 Athaliana_TFBM_v1.01.tm.index.direct.index.tm

File from: 
    
    /mnt/home/mjliu/kmer_5/Athaliana_TFBM_v1.01.tm.index.direct.index
    
DAPSeq: 

    python /mnt/home/seddonal/scripts/5_motif_merging/pcc_merge_CC.py create_cc_2 -t kmers10.txt.tamo -t2 DAP_motifs.txt.tm

File from: 

    /mnt/research/ShiuLab/14_DAPseq/PWM_to_tamo/DAP_motifs.txt.tm
    
 
3. Run the command line from step 2 using qsub_hpcc

        python /mnt/home/shius/codes/qsub_hpc.py -f submit -c runcc -mo TAMO
        
 check for failed jobs:
 
        python pcc_merge_CC.py check_outputs -c runcc
 
 output is a distance matrix, where the lower the number, the closer the distance (more similar the sequence)
 * Make sure all CIS-BP jobs finish running before moving on - the job files will overwrite eachother
** Fix this bug later!

4. Merge the outputs into one matrix

        python pcc_merge_CC.py combine_distance_matrix_2 -t [tamo_file_1] -t2 [tamo_file_2]
        
   
    CISBP: 

        python pcc_merge_CC.py combine_distance_matrix_2 -t kmers10.txt.tamo -t2 Athaliana_TFBM_v1.01.tm.index.direct.index.tm

    DAPSeq: 

        python ~mjliu/script_from_Alex/pcc_merge_CC.py combine_distance_matrix_2 -t kmers10.txt.tamo -t2 DAP_motifs.txt.tm

5. In Excel add the column and row labels. Colums are the motifs from t2 in order of "Athaliana_TFBM_v1.01.tm.index.direct.index" and “DAP_motifs.txt.tm_index” respectively; the row represent the motifs from "t1"; the order is based on “kmers.txt”

final PCC distance files: 

    kmers10.txt.tamo-Athaliana_TFBM_v1.01.tm.index.direct.index.tm.dm_mod
    kmers10.txt.tamo-DAP_motifs.txt.tm_mod

6. get enriched TF families
 
            python TF_with_between_correlation_low_PCC_average.py [TF family file] [merged distance matrix from step 4]
            
            TF family file: Athaliana_TFBM_v1.01.tm.index.direct.index_subset.txt









# Old Python & R Pipeline (useful for doing paired kmer enrichment)

Pairwise_kmers.py: Contains functions to make lists of paired kmers, make data frames of what genes contain those motifs, and run Fisher's Exact test to determine enrichment of those kmers/kmer pairs in the positive genes. 
RandomForest.R: Runs Random Forest on input dataframe. 10 replicates and 10 fold cross validation. 

## What you need:
•	File with all your positive examples (naming scheme will be based off the name of the positive example file, so make sure that makes sense and isn't too long)
•	File with all your negative examples
•	Fasta file with all gene promoter regions. 

If you already have fasta files of positive and neg examples, skip steps 2-3.

*Anytime you log in to HPC and want to use the pipeline you have to first run:

    export   PATH=/mnt/home/azodichr/miniconda3/bin:$PATH

Scripts: /mnt/home/azodichr/GitHub/MotifDiscovery/

## 1. Set Up Your Files:
1. Inside directory for Pairwise experiment make directories for FASTA files and Motif Lists:
  - mkdir FastaFiles
  - mkdir MotifLists

2. Put cluster file in FastaFiles dir and get promoter sequence:
  - cd FastaFiles/
  - cp [pos_examples] .
  - python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta [promoter sequences] -name [pos examples]
  - *For arabidopsis you can use: /mnt/home/azodichr/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa*

3. Put negative example file in FastaFiles dir and get promoter sequence. 
  - cp [neg_examples] .     #For Random Forest you want a 1:1 ratio of positive and negative examples, if 73 genes in cluster, randomly select 73 negative examples.
  - python /mnt/home/shius/codes/FastaManager.py -f getseq2 -fasta [promoter sequences] -name [neg examples]
  - *For arabidopsis you can use: /mnt/home/azodichr/01_DualStress_At/TAIR10_upstream1000_Alex_20101104.mod.fa*

4. Make singleton and paired kmer list. Reverse complement sequences are separated by “.”, pairs separated by a space (k = length of kmer you want):
  - cd ../MotifLists
  - python Pairwise_kmers.py -f make_pairs2 –k 5
  - python Pairwise_kmers.py -f make_pairs2 –k 6


## 2. Make presence/absence dataframes and do enrichment
The work flow here is 1) Make df with all kmers/pairs. 2) Make list of enriched kmers/pairs. 3) Remake df with just those enriched kmers/pairs.

1. Make data frame with presence or absence of all kmer/kmer pair:
  - python Pairwise_kmers.py -f make_df –k [ListOfKmers] -p [positive fasta files] -n [negative fasta files]
  - *If you want to add DNA structure information to your prediction ask me. It didn't add much to my prediction...*

2. If you want all the motifs move on to step 4, otherwise parse your motifs using Fisher's Exact Test:
  - python Pairwise_kmers.py -f parse_df –df [output df from step 5]
  OPTIONAL: -pval <Default is 0.05>

3. Re-make data frame with only enriched motifs:
  - python Pairwise_kmers.py -f make_df –k [output from step 6, ending in: “_sig_0.05.txt”] -p [positive fasta files] -n [negative fasta files]

## 3. Run Random Forest
####Method 1: Using scikit learn in python 
See 'If you want to make a data table and run RandomForest in one step' in the Python Pipeline above

####Method 2: Using R

If randomForest is not in your library yet, see *Getting RandomForest onto HPC.
  - export R_LIBS_USER=~/R/library
  - R --vanilla --slave --args [df*] < RandomForest.R
  - *Can use output df from step 1 or 3*

This will output two files:
  - .imp.txt: Open in exel, sort by "Mean Decrease Accuracy" - Make sure you shift the column headers over by one- they skip the motif name heading...
  - .Results.txt: Output with F-measure, stdev, sterror, and 95% confidence intervals.

*Getting RandomForest onto HPC:
  - Rscript -e "install.packages(‘LIBRARY_NAME',lib='~/R/library',contriburl=contrib.url('http://cran.r-project.org/'))”
  - export R_LIBS_USER=~/R/library      *you will need to run this line every time you run RandomForest.R*
  - R
  - Sys.getenv("R_LIBS_USER")
  - library(“LIBRARY_NAME")
  - q()



"""

Compiles information about pCREs into one database and performs statistical tests (Fisher's Exact & Welch's T-test)
to look for differences between pCRE features in cluster vs. non-cluster gene promoters. 

##### Input #####
Required:
-df                 Dataframe with Col 1 = genes, Col 2 = Class, Col 3-... = pCREs

Optional:
-enrich             FET results from pCRE finding pipeline (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/01_Kmer_Features/NDD_Ras_FET01_FETresults.txt)
-imp_RF             RF Machine learning importance file 
-imp_SVC            SVC Machine learning importance file
-map                GFF file of all pCREs of interest (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/all_kmers.gff)
-DHS                Also required (-map)
-DAP                Also required (-map)
-CNS                Also required (-map)
-WSC                Key file with histone markers, activation/repression status, and path to file. See azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/02_Histone/histone_overlap_files.txt). Also required (-map)
-Histone            Txt file with histone overlap info (Col 1 = histone marker name, Col 2 = activator/repressor, Col 3 = full path to overlap file)
                    See: /mnt/home/azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/02_Histone/histone_overlap_files.txt
-TFBM_CISBP         PCCD between kmers and CIS-BP motifs (mjliu/kmer_5/Athaliana_TFBM_v1.01.tm.index.direct.index.tm) using output from Ming's pipeline (~mjliu/kmer_5/kmer_5.sh)
-TFBM_DAPSeq        PCCD between kmers and DAP-Seq motifs (ShiuLab/14_DAPseq/PWM_to_tamo/DAP_motifs.txt.tm) using output from Ming's pipeline (~mjliu/kmer_5/kmer_5.sh)
-FS                 pCREs selected by decision tree feature selection rather than enrichement.
-GO                 Results from GO analysis output from http://go.princeton.edu/. Uses only sig GO terms (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/06_GOTerms/query_NNU.txt)
-arules             Dataframe with association rules as features (/mnt/home/azodichr/01_CombinedStress/Prasch_HD/07_6clusters/01_Kmer_Features/NNU_Ras_FET01_df_p0.01.txt_onlyARs_sup0.2_conf0.2_lift)
-Dist               Distribution of pCREs compared to shuffled sequences (azodichr/01_CombinedStress/Prasch_HD/07_6clusters/02_kmer_mapping/03_pCRE_Distribution/observe_vs_random_100_pCREs_NNU)

##### Possible Output Columns #####
pvalue              Enrichement Score 
percent_positives   Support score - i.e. Percent of cluster genes that contain the kmer in the promoter
Imp_RF              Importance score in RF predictive models
Imp_SVC             Importance score in SVC predictive models
Mean Copy#_pos_pres Average copy number in promoter region of cluster genes where kmer is present
Mean_Copy#_neg_pres Average copy number in promoter region of non-cluster genes where kmer is present
Copy#_Tstat         T-statistic from Welch's T-test. + = greater in pos, - = greater in negative 
Copy#_pval          Significance (2-tailed)
Mean_WSC_pos        Average within species nucleotide diverity of kmer in positive example genes
Mean_WSC_neg        Average within species nucleotide diverity of kmer in negative example genes
WSC_tstat           T-statistic from Welch's T-test. + = greater in pos, - = greater in negative 
WSC_pval            Significance (2-tailed)
DAP_odds            Fisher's Exact Test comparing pCREs overlapping with DAP-Deq peaks (200bp each) in cluster and non-cluster genes
DAP_pval            Significance score for DAP FET
DHS_Odds            Fisher's Exact Test comparing pCREs overlapping with DHS sites in cluster and non-cluster genes
DHS_pval            Significance score for DHS FET
CNS_Odds            Fisher's Exact Test comparing pCREs overlapping with CNS sites in cluster and non-cluster genes
CNS_pval            Significance score for CNS FET
FS_ovlp             Does the k-mer overlap with one of the 6-mers identified to predict the cluster using feature selection
arule               Does the pCRE belong to a significant association rule (Check out README.txt file to see pipeline for finding sig. arules)
CISBP_TFBM          Best match known TFBM from the CIS-BP Database
CISBP_TFBM_PCCD     Similarity (PCC-Distance - lower = better match) to best match known TFBM from the CIS-BP Database
DAPSeq_TFBM         Best match known TFBM from DAP-Seq Database
DAPSeq_TFBM_PCCD    Similarity (PCC-Distance - lower = better match) to best match known TFBM from DAP-Seq Database
GO_TERM             Created for each significant GO term, number of genes with a kmer present in the GO term
GO_TERM_%           Created for each significant GO term, percent of genes in the GO category that have the k-mer

"""
import sys
import numpy as np
np.set_printoptions(threshold=np.inf)
import pandas as pd
from scipy import stats

DF = ENRICH = IMP_RF = IMP_SVC = MAP = DHS = CNS = GO = WSC = TFBM = DAP = FS = TFBM_CISBP = TFBM_DAPSeq = HIST = ARULES = DIST = False
#act_repress = 'activation'

def ovrp_enriched(e_df, hits):
  """ Function to calculate enrichement of overlap with feature in cluster compared to non-cluster genes 

  Fisher's Exact Test
                                 Overlap with Feature   |  No overlap with feature
    pCRE in Cluster gene      |       t1                            t2
    pCRE in Non-Cluster gene  |       t3                            t4

"""
  subset = e_df[e_df['kmer'].isin(kmers)]
  subset = subset[subset.details.str.contains("overlap=NA") == False]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]

  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['t1'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['t3'])

  hits = pd.merge(hits, grouped_pos, how = 'left', right_index=True, left_index=True)
  hits = pd.merge(hits, grouped_neg, how = 'left', right_index=True, left_index=True)

  hits['t2'] = hits['Hits_pos_pres'] - hits['t1']
  hits['t4'] = hits['Hits_neg_pres'] - hits['t3']
  hits = hits.fillna(0)

  #if act_repress == "repression":
    #hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.t2, r.t1], [r.t4, r.t3]]), axis=1)
  #else:
    #print("pos enriched:")
  hits['temp'] = hits.apply(lambda r: stats.fisher_exact([[r.t1, r.t2], [r.t3, r.t4]]), axis=1)
  hits[['Odds', 'pval']] = hits['temp'].apply(pd.Series)

  return hits



for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-df":
        DF = sys.argv[i+1]
      if sys.argv[i] == '-enrich':
        ENRICH = sys.argv[i+1]
      if sys.argv[i] == '-imp_RF':
        IMP_RF = sys.argv[i+1]
      if sys.argv[i] == "-imp_SVC":
        IMP_SVC = sys.argv[i+1]
      if sys.argv[i] == '-map':
        MAP = sys.argv[i+1]
      if sys.argv[i] == "-DHS":
        DHS = sys.argv[i+1]
      if sys.argv[i] == "-CNS":
        CNS = sys.argv[i+1]
      if sys.argv[i] == "-WSC":
        WSC = sys.argv[i+1]
      if sys.argv[i] == "-TFBM_CISBP":
        TFBM_CISBP = sys.argv[i+1]
      if sys.argv[i] == "-TFBM_DAPSeq":
        TFBM_DAPSeq = sys.argv[i+1]
      if sys.argv[i] == "-DAP":
        DAP = sys.argv[i+1]
      if sys.argv[i] == "-Histone":
        HIST = sys.argv[i+1]
      if sys.argv[i] == "-FS":
        FS = sys.argv[i+1]
      if sys.argv[i] == "-arules":
        ARULES = sys.argv[i+1]
      if sys.argv[i] == "-GO":
        GO = sys.argv[i+1]
      if sys.argv[i] == "-Dist":
        DIST = sys.argv[i+1]

if len(sys.argv) <= 1:
  print(__doc__)
  exit()


####### Load Dataframe  #######

df_used = pd.read_csv(DF, sep='\t', header =0, index_col = 0)

kmers = np.delete(df_used.columns.values, 0)
pos_genes = df_used[df_used.Class==1].index.values
neg_genes = df_used[df_used.Class==0].index.values

df = pd.DataFrame(index = kmers, columns=None)


####### pCRE Criteria #######

if ENRICH:
  enrich_df = pd.read_csv(ENRICH, sep='\t', header =0, index_col = 0, usecols = ['feature', 'percent_postives', 'pvalue'])
  df = pd.merge(df, enrich_df, how = 'left', right_index=True, left_index=True)

if IMP_RF:
  imp_RF = pd.read_csv(IMP_RF, sep='\t', header =0, index_col = 0, names = ['Imp_RF'])
  df = pd.merge(df, imp_RF, how = 'left', right_index=True, left_index=True)

if IMP_SVC:
  imp_SVC = pd.read_csv(IMP_SVC, sep='\t', header =0, index_col = 0, names = ['Imp_SVC'])
  df = pd.merge(df, imp_SVC, how = 'left', right_index=True, left_index=True)

if MAP:
  map_df = pd.read_csv(MAP, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer'])
  subset = map_df[map_df['kmer'].isin(kmers)]
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]

  count_pos = subset_pos.groupby(['kmer','gene']).size()
  grouped_pos = pd.DataFrame(count_pos.groupby(level=0).mean(), columns = ['AvCopy#_pos_pres'])
  num_hits_pos = pd.DataFrame(count_pos.groupby(level=0).sum(), columns = ['Hits_pos_pres'])

  count_neg = subset_neg.groupby(['kmer','gene']).size()
  grouped_neg = pd.DataFrame(count_neg.groupby(level=0).mean(), columns = ['AvCopy#_neg_pres'])
  num_hits_neg = pd.DataFrame(count_neg.groupby(level=0).sum(), columns = ['Hits_neg_pres'])

  df = pd.merge(df, grouped_pos, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, grouped_neg, how = 'left', right_index=True, left_index=True)

  # Total number of kmer hits in promoters of pos and neg genes - usefor for fisher's enrichement tests.
  hits = pd.merge(num_hits_pos, num_hits_neg, how = 'left', right_index=True, left_index=True)

  # Perform Welch's T-test to look for statistical difference in pCRE copy number between pos and neg genes (non-normal dist)
  df_temp = pd.DataFrame(columns = ['Copy#_Tstat', 'Copy#_pval'], index = kmers)
  #print(count_pos)
  for k in kmers:
    val1 = count_pos.loc[k]
    val2 = count_neg.loc[k]
    tstat, pval = stats.ttest_ind(val1, val2, nan_policy = 'omit', equal_var = False)
    df_temp.loc[k] = [tstat, pval]

  df = pd.merge(df, df_temp, how = 'left', right_index=True, left_index=True )

if DAP:
  feat = pd.read_csv(DAP, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  feat_hits = ovrp_enriched(feat, hits)
  feat_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', "DAP_odds", "DAP_pval"]
  df = pd.merge(df, feat_hits.loc[:, ['DAP_odds', "DAP_pval"]], how = 'left', right_index=True, left_index=True)

if DHS:
  feat = pd.read_csv(DHS, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  feat_hits = ovrp_enriched(feat, hits)
  feat_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', "DHS_odds", "DHS_pval"]
  df = pd.merge(df, feat_hits.loc[:, ['DHS_odds', "DHS_pval"]], how = 'left', right_index=True, left_index=True)

if CNS:
  feat = pd.read_csv(CNS, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
  feat_hits = ovrp_enriched(feat, hits)
  feat_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', "CNS_odds", "CNS_pval"]
  df = pd.merge(df, feat_hits.loc[:, ['CNS_odds', "CNS_pval"]], how = 'left', right_index=True, left_index=True)

if WSC:
  wsc = pd.read_csv(WSC, sep='\t', header = 0)
  wsc['gene_and_region'], wsc['region'], wsc['kmer'] = wsc['Region'].str.split('|').str
  wsc['chr'], wsc['start'], wsc['stop'], wsc['gene'] = wsc['gene_and_region'].str.split('_').str
  subset = wsc[wsc['kmer'].isin(kmers)]
  subset['NtDiversity'] = pd.to_numeric(subset['NtDiversity'], errors = 'coerce') 
  #print(subset['NtDiversity'].unique())
  
  # Get just pos/neg example genes and convert NtDiv to numeric
  subset_pos = subset[subset['gene'].isin(pos_genes)]
  subset_neg = subset[subset['gene'].isin(neg_genes)]

  # Calculate means NtDiv for each kmer for genes in pos and neg datasets
  av_pos_mean = pd.DataFrame(subset_pos.groupby(['kmer'])['NtDiversity'].mean())
  av_pos_mean.columns = ['Mean_WSC_pos']
  av_neg_mean = pd.DataFrame(subset_neg.groupby(['kmer'])['NtDiversity'].mean())
  av_neg_mean.columns = ['Mean_WSC_neg']

  df = pd.merge(df, av_pos_mean, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, av_neg_mean, how = 'left', right_index=True, left_index=True)

  # Run Welch's T-test to determine if NtDiv of a kmer is statistically different between pos and neg genes.
  df_temp = pd.DataFrame(columns = ['WSC_tstat', 'WSC_pval'], index = kmers)
  for k in kmers:
    val1 = subset_pos[subset_pos['kmer'] == k]
    val2 = subset_neg[subset_neg['kmer'] == k]
    tstat, pval = stats.ttest_ind(val1['NtDiversity'], val2['NtDiversity'], nan_policy = 'omit', equal_var = False)
    df_temp.loc[k] = [tstat, pval]

  df = pd.merge(df, df_temp, how = 'left', right_index=True, left_index=True )


if HIST:
  with open(HIST, 'r') as hist_files:
    for l in hist_files:
      hist_name, act_repress, hist_file = l.strip().split('\t')
      his = pd.read_csv(hist_file, sep='\t', header = None, index_col = None, names = ['chr', 'gene', 'kmer', 'start','end','x','dir','y', 'details'], usecols = ['gene', 'kmer', 'details'])
      his_hits = ovrp_enriched(his, hits)
      odds = hist_name + "_" + act_repress[0] + "_Odds"
      pval = hist_name + "_" + act_repress[0] + "_pval"
      his_hits.columns = ["Hits_pos_pres", "Hits_neg_pres", "t1", "t3", "t2", "t4", 'temp', odds, pval]
      df = pd.merge(df, his_hits.loc[:, [odds, pval]], how = 'left', right_index=True, left_index=True)
      odds_direction = 'pos_enriched'

if FS:
  def overlap_function(x):
    """ Checks to see if there is overlap between a 6mer and the kmer identified by feature selection """
    for six_mer in fs_kmers:
      if six_mer in x:
        return 1
    return 0
  fs = pd.read_csv(FS, sep = "\t", header = 0, index_col=0)
  fs_kmers = np.delete(fs.columns.values, 0)
  df['FS_ovlp'] = df.index.values
  df['FS_ovlp'] = df['FS_ovlp'].apply(overlap_function)

if ARULES:
  ar_df = pd.read_csv(ARULES, sep="\t", header = 0, index_col=0)
  features = np.delete(ar_df.columns.values, 0)
  rule_kmers = []
  for f in features:
    for kmer in f.strip().split(":"):
      rule_kmers.append(kmer)
  rule_kmers = list(set(rule_kmers))
  df['arule'] = 0
  for k in rule_kmers:
    df.set_value(k, 'arule', 1)
  print(rule_kmers)

if GO:
  # Note that because the pCRE finding pipeline only looks at genes with promoters that don't 
  # overlap neighboring genes, some GO term genes might not be in the dataframe...
  
  go = pd.read_csv(GO, sep='\t', header=0, index_col=1, skiprows = 9)
  sig_go = go[go['CORRECTED_PVALUE']<=0.05]
  subset = sig_go['ANNOTATED_GENES']
  pos_hits_df = df_used[df_used.Class==1]

  terms = subset.to_dict()

  for t in terms:
    tx = t.replace(' ', '_')
    name1 = "GO_" + tx
    count_array = []
    for k in kmers:
      count = 0
      for term_genes in terms[t].split(", "):
        try:
          if pos_hits_df[k].loc[term_genes] == 1:
            count += 1
        except: 
          pass #print(t,k)
      count_array.append(count)
    df[name1] = count_array

  for t in terms:
    tx = t.replace(' ', '_')
    name = "GO_" + tx + "_%"
    name1 = "GO_" + tx
    df[name] = df[name1]/len(terms[t].split(", "))

if DIST:
  dist = pd.read_csv(DIST, header=None, sep='\t', index_col=0, names = ['1kb-900','900-800','800-700','700-600','600-500','500-400','400-300','300-200','200-100','100-TSS','drop'])
  dist = dist.drop('drop',1)
  dist['dist_max'] = dist.max(numeric_only=True, axis=1)
  dist['dist_ave'] = dist.mean(numeric_only=True, axis=1)
  df = pd.merge(df, dist, how = 'left', right_index=True, left_index=True)

if TFBM_CISBP:
  cisbp = pd.read_csv(TFBM_CISBP, sep='\t', header=0, index_col=0)
  min_val =  pd.DataFrame(cisbp.min(axis=1), columns = ['CISBP_TFBM_PCCD'])
  min_id = pd.DataFrame(cisbp.idxmin(axis=1), columns = ['CISBP_TFBM'])
  df = pd.merge(df, min_val, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, min_id, how = 'left', right_index=True, left_index=True)
  #df = df.drop_duplicates()

if TFBM_DAPSeq:
  dapseq = pd.read_csv(TFBM_DAPSeq, sep='\t', header=0, index_col=0)
  min_val = pd.DataFrame(dapseq.min(axis=1), columns = ['DAPSeq_TFBM_PCCD'])
  min_id = pd.DataFrame(dapseq.idxmin(axis=1), columns = ['DAPSeq_TFBM'])
  df = pd.merge(df, min_val, how = 'left', right_index=True, left_index=True)
  df = pd.merge(df, min_id, how = 'left', right_index=True, left_index=True)
  #df = df.drop_duplicates()

# For some reason the TFBM DISBP and DAPSeq sections duplicate the rows a bunch - this gets rid of the duplicates....
df = df.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')
print(df.head(10))
name = DF + '_pCRE_Criteria'
df.to_csv(name, sep = '\t')
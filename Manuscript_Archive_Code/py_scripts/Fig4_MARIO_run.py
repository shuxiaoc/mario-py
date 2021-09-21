#### Bokai Zhu
#### Script used to perform MARIO with the data presented in Figure 2
#### bone marrow data from citeseq and cytof


###############
### NOTICE ####
###############

### the script is run with devo version of MARIO, and is submitted as a record
### current version of MARIO should not run this script
### MARIO devo used here has same core-algorithm, though the grammar and format/structure of the code is legacy type

##### import necessary packages
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.sparse.csgraph import maximum_bipartite_matching
import pandas as pd
import timeit
import sys
sys.path.append('../MARIO_devo_version/code_0316/') # this is devo version of MARIO
from utils import *
from clustering import *
from matching import *
import matplotlib
import matplotlib.pyplot as plt

# read in protein expression matrices, with the annotations etc

df1_full = pd.read_csv("../data/murine_spleen/clean_cells.csv") # codex matrix
df2_full = pd.read_csv("../data/murine_spleen/citeMurine-206ForMatch-0423_whole2x.csv") # citeseq matrix duplicated version
# some protein are same but with different names, change them
df2_full.rename(columns={'Ly.6G':'Ly6G'}, inplace=True)
df2_full.rename(columns={'Ly.6C':'Ly6C'}, inplace=True)
df2_full.rename(columns={'CD90.1.Thy.1.1.':'CD90'}, inplace=True)
df2_full.rename(columns={'F4.80':'F480'}, inplace=True)
df2_full.rename(columns={'CD79b.Igb.':'CD79b'}, inplace=True)
df2_full.rename(columns={'CD45R.B220':'B220'}, inplace=True)
df2_full.rename(columns={'CD21.CD35.CR2.CR1.':'CD2135'}, inplace=True)
df2_full.rename(columns={'CD335.NKp46.':'NKp46'}, inplace=True)
df2_full.rename(columns={'CD169.Siglec.1.':'CD169'}, inplace=True)
df2_full.rename(columns={'CD16.32':'CD1632'}, inplace=True)
df2_full.rename(columns={'I.A.I.E':'MHCII'}, inplace=True)
df2_full.rename(columns={'TER.119.ErythroidCells':'Ter119'}, inplace=True)
df2_full.rename(columns={'TCRbchain':'TCR'}, inplace=True)

# for codex dataset we do not need to use all the channels
df1columnsUse=['Unnamed: 0', 'Ly6C', 'TCR', 'Ly6G', 'CD19',
        'CD169','CD106', 'CD3', 'CD1632', 
       'CD8a', 'CD90',  'F480', 'CD11c',
       'CD11b', 'IgD','CD27', 'CD5', 'CD79b',
       'CD4', 'IgM',
       'B220', 'MHCII', 'CD35', 'CD2135',
       'nucl', 'NKp46', 'cluster.sr',
       'cluster.term', 'cellLabelInImage', 'PointNum']
df1_full=df1_full[df1columnsUse]

### randomnize data, because orignal data is spatially slotted
### this step is automatic in current MARIO
n, _ = df1_full.shape
df1_full_rand=df1_full.iloc[np.random.choice(n, n, replace=False),:]

### here we MARIO match with batches, the current MARIO has integrated in the function already
### only legacy version use manual
loop=31 # 32 batches
idx=0
size = round(n1/loop) # in each batch how many cells in dataX used for matching
for i in range(loop):
    idx=idx+1
    print(idx)
    if idx<loop:
        df1=df1_full_rand.iloc[size*(idx-1):size*idx,:]
    else:
        df1=df1_full_rand.iloc[size*(idx-1):,:]
    df2=df2_full # match to the full citeseq dataset
    df1_save=df1
    df2_save=df2
    X_celcelID=df1['Unnamed: 0'] # not actually used
    Y_celcelID=df2['cellID'] # not actually used
    X_labels = df1['cluster.term'].to_numpy()
    Y_labels = df2['cluster.info'].to_numpy()
    
    n1, p1 = df1.shape
    n2, p2 = df2.shape
    # remove un used columns for mario
    df1 = df1.drop(['Unnamed: 0','cluster.sr','cluster.term','cellLabelInImage','PointNum','nucl'], 1)
    df2 = df2.drop(['Unnamed: 0','cellID','cluster.orig','cluster.info'], 1)
   
    # contruct object
    match = CellMatching(df1, df2, normalization=True)
    del df1, df2
    # ovlp
    dist_ovlp, s = match.compute_dist_ovlp(n_components=20)
    match.specify_matching_params(m_min=1, m_max=1, num_cells_to_use=match.n1)
    match_ovlp = match.match_cells(dist_ovlp, sparsity=1000, mode='auto')
    print("ovlp finished")
    # non_ovlp
    dist_non_ovlp = match.compute_dist_all(match_ovlp, n_components=15)
    match_non_ovlp = match.match_cells(dist_non_ovlp, mode='dense')
    print("nonovlp finished")
    # both
    wt = 0.9 # this value is preslected based on hyperparamter screening
    match_all = match.match_cells(wt*dist_ovlp+(1-wt)*dist_non_ovlp, sparsity=1000, mode='auto')
    eval_matching_accuracy(X_labels, Y_labels, match_all, 'maj')
    # filter
    match_final = match.filter_bad_matches(match_all, n_clusters=10, n_components=15, bad_prop=0.05,
                         max_iter=10, tol=1e-4, verbose=True)
    print('Filtering finished!')
    print('Matching accuracy is {}'.format(eval_matching_accuracy(X_labels, Y_labels, match_final, 'maj')))
    print("{}% cells remain.".format(np.sum([1 for ii in range(match.n1) if len(match_final[ii]) != 0])/match.n1*100))
    # get iloc for df1 and df2 after matching and filtering
    match_final_df1=list(range(len(match_final)))
    filtered_out=[i for i,x in enumerate(match_final) if not x]
    match_final_df1_filterd=[e for e in match_final_df1 if not e in filtered_out]
    match_final_df2 = [item for sublist in match_final for item in sublist]
    # actual data
    if idx<loop:
        df1_1_orig=df1_full_rand.iloc[size*(idx-1):size*idx,:]
    else:
        df1_1_orig=df1_full_rand.iloc[size*(idx-1):,:]
        
    df1_m=df1_1_orig.iloc[match_final_df1_filterd,]
    df2_m=df2_full.iloc[match_final_df2,]

    
    path_core="../temp/location/"
    print("saving files")
    path1="codex_matched"+str(idx)+".csv"
    path2="murineCite_matched"+str(idx)+".csv"
    df1_m.to_csv(path_core+path1)
    df2_m.to_csv(path_core+path2)
    
##### manual batch loop finished 
## gather maual batched, all these are automatic in current MARIO
df1 = pd.DataFrame()
# read and gather all the codex matched cells
for i in range(31):
    f="../temp/location/codex_matched"
    f_full=f+str(i+1)+".csv"
    csv = pd.read_csv(f_full)
    df1 = df1.append(csv)
    
# read and gather all the citeseq matched cells
df2 = pd.DataFrame()
for i in range(31):
    f="../temp/location/murineCite_matched"
    f_full=f+str(i+1)+".csv"
    csv = pd.read_csv(f_full)
    df2 = df2.append(csv)
# save out, these are the matched files used for downstream analysis
df1.to_csv("../data/murine_spleen/MurineCodex_mario_matched.csv")
df2.to_csv("../data/murine_spleen/MurineCiteseq_mario_matched.csv")
# produce cca scores
df1_value=normalize(df1.iloc[:,1:26])
df2_value=normalize(df2.iloc[:,2:213])
cancor, cca = get_cancor(df1_value, df2_value, 20, max_iter=1000)
df1_cca, df2_cca = cca.transform(df1_value, df2_value)
np.savetxt("../data/murine_spleen/codex_cca.csv", df1_cca, delimiter=",")
np.savetxt("../data/murine_spleen/cite_cca.csv", df2_cca, delimiter=",")
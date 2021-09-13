#### Bokai Zhu
#### Script used to perform MARIO with the data presented in Figure 5
#### macrophage cells from CODEX and CITE-seq covid-19 patients


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

df1 = pd.read_csv("../data/COVID-19/CODEX_covid19_macrophage.csv")
df2 = pd.read_csv("../data/COVID-19/citeseq_covid19_macrophage.csv")
df1_full=df1
df2_full=df2
### randomnize data
n, _ = df1_full.shape
df1_full_rand=df1_full.iloc[np.random.choice(n, n, replace=False),:] # in-built in current MARIO
###
loop=20 # currently called batch
idx=0
size = round(n1/loop) # in each batch how many cells in dataX used for matching

for i in range(loop):
    idx=idx+1
    print(idx)
    df1=df1_full_rand.iloc[size*(idx-1):size*idx,:]
    df2=df2_full # match to the full citeseq dataset
    df1_save=df1
    df2_save=df2
    X_celcelID=df1['Unnamed: 0']
    Y_celcelID=df2['X']
    X_labels = df1['cluster.term'].to_numpy()
    Y_labels = df2['cluster.sr'].to_numpy()
    n1, p1 = df1.shape
    n2, p2 = df2.shape
    df1 = df1.drop(['Unnamed: 0','ClusterID',"EventID","File.Name", # unused columns
                    "Index.in.File","Event.Name","Comment","annotation",
                    'cluster.term','tile_num','x','y','z','x_tile','y_tile','size',
                    'cell_id','tma'], 1)
    df2 = df2.drop(['Unnamed: 0','X','X.1', 'celltype','cluster.term','cluster.sr'], 1)
    # remove non-variable conlumns
    df1 = df1.loc[:, df1.std() > 0]
    df2 = df2.loc[:, df2.std() > 0]
    # contruct object
    match = CellMatching(df1, df2, normalization=True)
    del df1, df2
    # ovlp
    dist_ovlp, s = match.compute_dist_ovlp(n_components=25)
    match.specify_matching_params(m_min=1, m_max=1, num_cells_to_use=match.n1)
    match_ovlp = match.match_cells(dist_ovlp, sparsity=1000, mode='auto')
    eval_matching_accuracy(X_labels, Y_labels, match_ovlp, 'maj')
    print("ovlp finished")
    # non_ovlp
    dist_non_ovlp = match.compute_dist_all(match_ovlp, n_components=25)
    match_non_ovlp = match.match_cells(dist_non_ovlp, sparsity=1000, mode='auto')
    print("nonovlp finished")
    # both
    wt = 1 # value preselected from pre-screenining
    match_all = match.match_cells(wt*dist_ovlp+(1-wt)*dist_non_ovlp, sparsity=1000, mode='auto')
    eval_matching_accuracy(X_labels, Y_labels, match_all, 'maj')
    # filter
    match_final = match.filter_bad_matches(match_all, n_clusters=3, n_components=5, bad_prop=0.1,
                         max_iter=15, tol=1e-4, verbose=True)
    print('Filtering finished!')
    print('Matching accuracy is {}'.format(eval_matching_accuracy(X_labels, Y_labels, match_final, 'maj')))
    print("{}% cells remain.".format(np.sum([1 for ii in range(match.n1) if len(match_final[ii]) != 0])/match.n1*100))
    # get iloc for df1 and df2 after matching and filtering
    match_final_df1=list(range(len(match_final)))
    filtered_out=[i for i,x in enumerate(match_final) if not x]
    match_final_df1_filterd=[e for e in match_final_df1 if not e in filtered_out]
    match_final_df2 = [item for sublist in match_final for item in sublist]
    # actual data
    df1_1_orig=df1_full_rand.iloc[size*(idx-1):size*idx,:]
    df1_m=df1_1_orig.iloc[match_final_df1_filterd,]
    df2_m=df2_full.iloc[match_final_df2,]    
    path_core="../temp/location/"
    print("saving files")
    path1="codex_matched"+str(idx)+".csv"
    path2="citebleg"+str(idx)+".csv"
    df1_m.to_csv(path_core+path1)
    df2_m.to_csv(path_core+path2)
    
## now gather the matched cells and calculate the cca scores
# 20 batches remerge, this is avoided in current MARIO (automatic)
df1 = pd.DataFrame()
for i in range(20):
    f="../temp/location/codex_matched"
    f_full=f+str(i+1)+".csv"
    csv = pd.read_csv(f_full)
    df1 = df1.append(csv)
    
df2 = pd.DataFrame()
for i in range(20):
    f="../temp/location/citebleg"
    f_full=f+str(i+1)+".csv"
    csv = pd.read_csv(f_full)
    df2 = df2.append(csv)
    
df1.to_csv("../temp/location/lung_mph_mario_matched.csv")
df2.to_csv("../temp/location/balf_mph_mario_matched.csv")

### no cca is produced since here we only using matching information

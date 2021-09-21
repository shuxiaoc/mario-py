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
df1_full = pd.read_csv("../data/bonemarrow/bmcite_29000.csv")
df2_full = pd.read_csv("../data/bonemarrow/levine32_102977.csv")
df1=df1_full
df2=df2_full
df1.rename(columns={'HLA.DR':'HLA-DR'}, inplace=True)
#df1.rename(columns={'CD8a':'CD8'}, inplace=True) # changed already
X_celcelID=df1['Unnamed: 0']
Y_celcelID=df2['Unnamed: 0']
X_labels = df1['cluster.term'].to_numpy()
Y_labels = df2['cluster.term'].to_numpy()
X_labels_orig = df1['cluster.orig'].to_numpy()
Y_labels_orig = df2['cluster.orig'].to_numpy()
n1, p1 = df1.shape
n2, p2 = df2.shape
# dont need several columns for matching
df1 = df1.drop(['Unnamed: 0', 'cluster.orig','cluster.term'], 1)
df2 = df2.drop(['Unnamed: 0', 'cluster.orig','cluster.term'], 1)
shared_features = [x for x in df1.columns if x in df2.columns and x != 'cluster.orig'] # not used actually


### here we MARIO match with batches, the current MARIO has integrated in the function already
### only legacy version
loop=4
idx=0
size = round(n1/loop) # in each batch how many cells in dataX used for matching
for i in range(loop):
    idx=idx+1
    df1_1=df1.iloc[size*(idx-1):size*idx,:]
    df2=df2_full.drop(['Unnamed: 0', 'cluster.orig','cluster.term'], 1)
    # contruct object
    match = CellMatching(df1_1, df2, normalization=True) # this is legacy MARIO
    del df1_1, df2
    # ovlp
    dist_ovlp, s = match.compute_dist_ovlp(n_components=10)
    match.specify_matching_params(m_min=1, m_max=1, num_cells_to_use=match.n1)
    match_ovlp = match.match_cells(dist_ovlp, sparsity=1000, mode='auto')
    eval_matching_accuracy(X_labels[size*(idx-1):size*idx], Y_labels, match_ovlp, 'maj')
    print("ovlp finished")
    # non_ovlp
    dist_non_ovlp = match.compute_dist_all(match_ovlp, n_components=20)
    match_non_ovlp = match.match_cells(dist_non_ovlp, sparsity=1000, mode='auto')
    # both
    wt = 0.5 # this value is preslected based on hyperparamter screening
    match_all = match.match_cells(wt*dist_ovlp+(1-wt)*dist_non_ovlp, sparsity=None, mode='auto')
    eval_matching_accuracy(X_labels[size*(idx-1):size*idx], Y_labels, match_all, 'maj')
    # filter
    match_final = match.filter_bad_matches(match_all, n_clusters=15, n_components=20, bad_prop=0.2,
                         max_iter=10, tol=1e-4, verbose=True)
    print('Filtering finished!')
    print('Matching accuracy is {}'.format(eval_matching_accuracy(X_labels[size*(idx-1):size*idx],
                                                                  Y_labels, match_final, 'maj')))
    print("{}% cells remain.".format(np.sum([1 for ii in range(match.n1) if len(match_final[ii]) != 0])/match.n1*100))
    # matching finished
    # get iloc for df1 and df2 after matching and filtering
    match_final_df1=list(range(len(match_final)))
    filtered_out=[i for i,x in enumerate(match_final) if not x]
    match_final_df1_filterd=[e for e in match_final_df1 if not e in filtered_out]
    match_final_df2 = [item for sublist in match_final for item in sublist]
    # actual data
    df1_1_orig=df1_full.iloc[size*(idx-1):size*idx,:]
    df1_m=df1_1_orig.iloc[match_final_df1_filterd,]
    df2_m=df2_full.iloc[match_final_df2,]
    
    path_core="../temp/location/"
    print("saving files")
    path1= path_core + "cite"+str(idx)+".csv"
    path2=path_core + "levine"+str(idx)+".csv"
    df1_m2.to_csv(path1)
    df2_m2.to_csv(path2)
    
    
## now gather the matched cells and calculate the cca scores
# four batches remerge, this is avoided in current MARIO
bmc1 = pd.read_csv("../temp/location/cite1.csv")
bmc2 = pd.read_csv("../temp/location/cite2.csv")
bmc3 = pd.read_csv("../temp/location/cite3.csv")
bmc4 = pd.read_csv("../temp/location/cite4.csv")
bmcfull=bmc1.copy()
bmcfull = bmcfull.append([bmc2, bmc3 ,bmc4])
lvc1 = pd.read_csv("../temp/location/levine1.csv")
lvc2 = pd.read_csv("../temp/location/levine2.csv")
lvc3 = pd.read_csv("../temp/location/levine3.csv")
lvc4 = pd.read_csv("../temp/location/levine4.csv")
lvcfull=lvc1.copy()
lvcfull = lvcfull.append([lvc2, lvc3 ,lvc4])
# scale the values before calculation, automation in current MARIO
lvcfull_value=normalize(lvcfull.iloc[:,2:34])
bmcfull_value=normalize(bmcfull.iloc[:,2:27])
# calculate CCA scores
cancor, cca = get_cancor(bmcfull_value, lvcfull_value, 20, max_iter=1000)
df1_cca, df2_cca = cca.transform(bmcfull_value, lvcfull_value)
# save the CCA scores
np.savetxt("../temp/location/bm_cca.csv", df1_cca, delimiter=",")
np.savetxt("../temp/location/lv_cca.csv", df2_cca, delimiter=",")

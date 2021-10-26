#### Bokai Zhu
#### Script used to perform MARIO with the data presented in Figure S1.2
#### pbmc data from citeseq and cytof

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
df1_full = pd.read_csv("../data/pbmc/10xciteseq_5107_clean.csv")
df2_full = pd.read_csv("../data/pbmc/felixcytof_38866.csv")
# change one feature name the same
df2_full.rename(columns={'CCR7':'CD197'}, inplace=True)
df1=df1_full
df2=df2_full
X_celcelID=df1['X'] # id not used 
Y_celcelID=df2['X']
X_labels = df1['cluster.info'].to_numpy() # for testing of performance
Y_labels = df2['cluster.info'].to_numpy()
n1, p1 = df1.shape
n2, p2 = df2.shape
df1 = df1.drop(['Unnamed: 0', 'cluster.info','X'], 1) # drop unused columns during MARIO
df2 = df2.drop(['Unnamed: 0', 'cluster.info','X'], 1)
shared_features = [x for x in df1.columns if x in df2.columns and x != 'cluster.info'] # not used actually
match = CellMatching(df1, df2)
del df1, df2
dist_ovlp, s = match.compute_dist_ovlp(n_components=15)
# quick check whtat n_components to use ovlp matching
plt.scatter([i for i in range(len(s))], np.sqrt(s))
plt.title('eigenvalues')
# ovlp mathcing
dist_ovlp, s = match.compute_dist_ovlp(n_components=12)
match.specify_matching_params(m_min=1, m_max=1, num_cells_to_use=match.n1) # legacy version of MARIO
match_ovlp = match.match_cells(dist_ovlp, sparsity=1000, mode='dense')
# quick check for n_components for all matching
cancor, cca = match.fit_cca(match_ovlp, 28)
plt.scatter(np.arange(28), cancor)
plt.title("Canonical correlations")
# all matching 
dist_all = match.compute_dist_all(match_ovlp, n_components=17)
match_all = match.match_cells(dist_all, mode='auto')
# since this dataset is very small just interporlate a list of weights and find best one

## interprolation looping starts 
num_wts = 8
wt_vec = np.linspace(0, 1, num_wts)
num_cancors_to_look_at = 5
sparsity_level = 1000
dist_list = []
matching_list = []
cancor_list = []
acc_list = []
for ii in range(num_wts):
    print(ii)
    if ii == 0:
        dist_list.append(dist_ovlp)
        matching_list.append(match_ovlp)
    elif ii == num_wts-1:
        dist_list.append(dist_all)
        matching_list.append(match_all)
    else:
        # ii small --> more close to dist_ovlp
        dist_list.append((1 - wt_vec[ii])*dist_ovlp + wt_vec[ii]*dist_all)
        matching_list.append(match.match_cells(dist_list[ii], sparsity=sparsity_level, mode='auto'))
    # compute median of top  cancors
    cancor_list.append(np.median(match.fit_cca(matching_list[ii], n_components=num_cancors_to_look_at)[0]))
    acc_list.append(eval_matching_accuracy(X_labels, Y_labels, matching_list[ii], 'maj'))
best_wt_idx = np.argmax(cancor_list)
# plot out the median of cancor scores
plt.scatter(np.arange(8), cancor_list)
plt.title("Canonical correlations")
# best interpolation based on plotting
match_wted = matching_list[3]
dist_wted = dist_list[3]
match_all = match.match_cells(dist_wted, sparsity=None, mode='auto')
# then regularized filtering
match_final = match.filter_bad_matches(match_all, n_clusters=15,
                                       n_components=20, bad_prop=0.2,
                                       max_iter=10, tol=1e-4, verbose=True)
# quick check of matching accuracy
print("{} many cells remain.".format(np.sum([1 for ii in range(match.n1) if len(match_final[ii]) != 0])/match.n1))
# 94%

## start prepping the final files for downstream analysis
match_final_df1=list(range(len(match_final)))
filtered_out=[i for i,x in enumerate(match_final) if not x]
match_final_df1_filterd=[e for e in match_final_df1 if not e in filtered_out]
match_final_df2 = [item for sublist in match_final for item in sublist]
# actual data
df1_1_orig=df1_full
df1_m=df1_1_orig.iloc[match_final_df1_filterd,]
df2_m=df2_full.iloc[match_final_df2,]
# also directly save the cca scores in the same file
concor,cca=match.fit_cca(match_final, n_components=20, max_iter=1000)
df1_aligned, df2_aligned = match._align_modalities(match_final)
df1_cca, df2_cca = cca.transform(df1_aligned, df2_aligned)
n_components=10
df1CCA_frame = pd.DataFrame(data=df1_cca, index=np.arange(df1_cca.shape[0]), 
             columns=['CCA'+str(i) for i in range(df1_cca.shape[1])])
df2CCA_frame = pd.DataFrame(data=df2_cca, index=np.arange(df2_cca.shape[0]), 
             columns=['CCA'+str(i) for i in range(df2_cca.shape[1])])
df1CCA_frame.index=df1_m.index
df2CCA_frame.index=df2_m.index
df1_m2 = pd.concat([df1_m, df1CCA_frame], axis=1)
df2_m2 = pd.concat([df2_m, df2CCA_frame], axis=1)
path_core="../temp/location/"
print("saving files")
path1=path_core + "tenxGen"+str(1)+".csv"
path2=path_core + "Felix"+str(1)+".csv"
df1_m2.to_csv(path_core + path1) # this file contains the origianl protein features rows matched, with CCA scores 
df2_m2.to_csv(path_core + path2)


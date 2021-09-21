#### Bokai Zhu
#### Script used to perform MARIO with the data presented in Figure S3.2
#### human whole blood sample challenged with H1N1
#### human and non-human primate whole blood sample challenged with IL4


###############
### NOTICE ####
###############

### the script is run with devo version of MARIO, and is submitted as a record
### current version of MARIO should not run this script
### MARIO devo used here has same core-algorithm, though the grammar and format/structure of the code is legacy type

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

##### this script includes mathcing of more then two datasets
##### current version of MARIO can auto detect multiple datasets
##### legacy MARIO need to do this manually:

###################### first match wcct human -> xspecies human IL4

# read in data
df1_full = pd.read_csv("../data/Cytof-Xspecies-IFNG/wcct4d_term_srcluster120k.csv") # wcct human is still the same input
df2_full = pd.read_csv("../data/Cytof-Xspecies-IL4/zac_human_il4_120k.csv") # il4 stimulated human whole blood cell
n1, p1 = df1_full.shape

### here we MARIO match with batches, the current MARIO has integrated in the function already
### only legacy version
loop=4 # run with four batches
idx=0
size = round(n1/loop)
for i in range(loop):
    idx=idx+1
    df1=df1_full
    df2=df2_full # match to the full human ifng dataset
    df1 = df1.iloc[size*(idx-1):size*idx,:]
    X_celcelID=df1['index']
    Y_celcelID=df2['index']
    X_labels = df1['cluster.sr'].to_numpy()
    Y_labels = df2['cluster.sr'].to_numpy()
    n1, p1 = df1.shape
    n2, p2 = df2.shape
    df1 = df1.drop(['Unnamed: 0','index', 'cluster.orig','cluster.sr'], 1) # remove non-value columns
    df2 = df2.drop(['Unnamed: 0','index','cluster.orig','cluster.sr'], 1)
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
    wt = 0.8 # this value is preslected based on hyperparamter screening
    match_all = match.match_cells(wt*dist_ovlp+(1-wt)*dist_non_ovlp, sparsity=1000, mode='auto')
    # final filter
    match_final = match.filter_bad_matches(match_all, n_clusters=10, n_components=10, bad_prop=0.1,
                         max_iter=8, tol=1e-4, verbose=True)
    print('Filtering finished!')
    print('Matching accuracy is {}'.format(eval_matching_accuracy(X_labels, Y_labels, match_final, 'maj')))
    print("{}% cells remain.".format(np.sum([1 for ii in range(match.n1) if len(match_final[ii]) != 0])/match.n1*100))
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
    path1="wcct_sr_matched"+str(idx)+".csv"
    path2="Zac_human_il4_sr_matched"+str(idx)+".csv"
    df1_m.to_csv(path_core+path1)
    df2_m.to_csv(path_core+path2) # save out the row matched files

    
###################### second match wcct human -> xspecies rhesus interferon gamma    

df1_full = pd.read_csv("../data/Cytof-Xspecies-IFNG/wcct4d_term_srcluster120k.csv") # the same wcct human 
df2_full = pd.read_csv("../data/Cytof-Xspecies-IL4/zac_rhesus_il4_120k.csv") # rhesus data

######## legacy looping
loop=4 # run with four batches
idx=0
size = round(n1/loop)
for i in range(loop):
    df1=df1_full
    df2=df2_full # match to the full human ifng dataset
    df1 = df1.iloc[size*(idx-1):size*idx,:]
    X_celcelID=df1['index']
    Y_celcelID=df2['index']
    X_labels = df1['cluster.sr'].to_numpy()
    Y_labels = df2['cluster.sr'].to_numpy()
    n1, p1 = df1.shape
    n2, p2 = df2.shape
    df1 = df1.drop(['Unnamed: 0','index', 'cluster.orig','cluster.sr'], 1)
    df2 = df2.drop(['Unnamed: 0','index','cluster.orig','cluster.sr'], 1)
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
    wt = 0.8 # this value is preslected based on hyperparamter screening
    match_all = match.match_cells(wt*dist_ovlp+(1-wt)*dist_non_ovlp, sparsity=2000, mode='auto')
    eval_matching_accuracy(X_labels, Y_labels, match_all, 'maj')
    # filter
    match_final = match.filter_bad_matches(match_all, n_clusters=10, n_components=10, bad_prop=0.1,
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
    df1_1_orig=df1_full.iloc[size*(idx-1):size*idx,:]
    df1_m=df1_1_orig.iloc[match_final_df1_filterd,]
    df2_m=df2_full.iloc[match_final_df2,]
    path_core="../temp/location/"
    print("saving files")
    path1="wcct_sr_matched"+str(idx)+".csv"
    path2="Zac_rheus_il4_sr_matched"+str(idx)+".csv"
    df1_m.to_csv(path_core+path1)
    df2_m.to_csv(path_core+path2)

    
    
###################### thrid match wcct human -> xspecies cyno-crab eating interferon gamma    

df1_full = pd.read_csv("../data/Cytof-Xspecies-IFNG/wcct4d_term_srcluster120k.csv") # the same wcct human 
df2_full = pd.read_csv("../data/Cytof-Xspecies-IL4/zac_cyno_il4_120k.csv") # cynomolgus data

######## legacy looping
loop=4 # run with four batches
idx=0
size = round(n1/loop)
for i in range(loop):
    df1=df1_full
    df2=df2_full # match to the full human ifng dataset
    df1 = df1.iloc[size*(idx-1):size*idx,:]
    X_celcelID=df1['index']
    Y_celcelID=df2['index']
    X_labels = df1['cluster.sr'].to_numpy()
    Y_labels = df2['cluster.sr'].to_numpy()
    n1, p1 = df1.shape
    n2, p2 = df2.shape
    df1 = df1.drop(['Unnamed: 0','index', 'cluster.orig','cluster.sr'], 1)
    df2 = df2.drop(['Unnamed: 0','index','cluster.orig','cluster.sr'], 1)
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
    wt = 0.8 # this value is preslected based on hyperparamter screening
    match_all = match.match_cells(wt*dist_ovlp+(1-wt)*dist_non_ovlp, sparsity=2000, mode='auto')
    eval_matching_accuracy(X_labels, Y_labels, match_all, 'maj')
    # filter
    match_final = match.filter_bad_matches(match_all, n_clusters=10, n_components=10, bad_prop=0.1,
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
    df1_1_orig=df1_full.iloc[size*(idx-1):size*idx,:]
    df1_m=df1_1_orig.iloc[match_final_df1_filterd,]
    df2_m=df2_full.iloc[match_final_df2,]
    path_core="../temp/location/"
    print("saving files")
    path1="wcct_sr_matched"+str(idx)+".csv"
    path2="Zac_cyno_il4_sr_matched"+str(idx)+".csv"
    df1_m.to_csv(path_core+path1)
    df2_m.to_csv(path_core+path2)
    
    
####################   now start the generalize cca production


### functions for legacy only, current MARIO integrated these functions
from sklearn.cross_decomposition import CCA
def gcca_init(data_list, n_components=10, normalization=True):
    """
    Parameters
    ----------
    data_list: a list of data matrices, e.g., [X1, X2, X3, X4]
    Returns
    ----------
    score_list: a list of projected data
    """
    # normalize
    for data in data_list:
        data = normalize(data)
    # run CCA on the first two data matrices
    cca = CCA(n_components=n_components, max_iter=1000)
    score_list = list(cca.fit_transform(data_list[0], data_list[1]))
    avg_score = (score_list[0] + score_list[1]) / 2
    # run OLS regressions for the rest of data matrices
    for data in data_list[2:]:
        B = np.linalg.lstsq(data, avg_score, rcond=None)[0]
        score_list.append(data @ B)
    
    return score_list

def gcca_refine(data_list, init_scores, max_iter=500, tol=1e-3, verbose=True):
    """
    Parameters
    ---------
    data_list: a list of data matrices
    init_scores: a list of initial scores
    """
     
    m = len(init_scores)
    if len(init_scores[0]) == 0:
        n_components = 1
        n = len(init_scores[0])
    else: 
        n, n_components = init_scores[0].shape
        
    all_scores = np.array([np.empty_like(init_scores[0])]*m)
    for ii in range(n_components):
        if verbose:
            print("Computing the {}-th canonical score...".format(ii))
        # regress out the influence of all_score[:, :, :ii]
        if ii > 0:
            curr_data_list = []
            for idx, data in enumerate(data_list):
                coef = np.linalg.lstsq(all_scores[idx, :, :ii], data, rcond=None)[0]
                curr_data_list.append(data - all_scores[idx, :, :ii] @ coef)
        else:
            curr_data_list = data_list
        
        
        # preparations for entering main iteration 
        if n_components == 1:
            curr_scores = [score for score in init_scores]
        else:
            curr_scores = [score[:, ii] for score in init_scores]
        score_sum = sum(score for score in curr_scores)
        prev_obj = np.inf
        
        # alternating minimization
        for iter_idx in range(max_iter):
            for jj in range(m):
                y = (score_sum - curr_scores[jj]) / (m-1)
                X = curr_data_list[jj]
                coef = np.linalg.lstsq(X, y, rcond=None)[0]
                coef = coef / np.sqrt(np.sum(coef**2))
                score = X @ coef
                score_sum = score_sum - curr_scores[jj] + score
                curr_scores[jj] = score
            # check convergence
            obj = 0
            for jj in range(m-1):
                obj += np.sqrt(np.sum((curr_scores[jj] - curr_scores[jj+1])**2))
            
            if verbose and iter_idx % 50 == 0:
                print("At iteration {}, the objective value is {}.".format(iter_idx, obj))
            
            if abs(obj - prev_obj) < tol:
                break
            else:
                prev_obj = obj
        
        if verbose:
            print("Finished computing the {}-th canonical score, the objective is {}.".format(ii, obj))
            
        for jj in range(m):
            if n_components == 1:
                all_scores[jj, :] = curr_scores[jj]
            else:
                all_scores[jj, :, ii] = curr_scores[jj]
                
        # curr_scores are linear combinations of the residual matrices
        # convert them into linear combinations of the original data matrices
#         for jj in range(m):
#             coef = np.linalg.lstsq(data_list[jj], curr_scores[jj], rcond=None)[0]
#             if n_components == 1:
#                 all_scores[jj, :] = data_list[jj] @ coef
#             else:
#                 all_scores[jj, :, ii] = data_list[jj] @ coef
        
    return all_scores

def gcca(data_list, n_components=10, normalization=True, max_iter=500, tol=1e-3, verbose=True):
    """
    Returns
    --------
    all_scores: a list of projected scores
    """
    init_scores = gcca_init(data_list, n_components, normalization)
    return gcca_refine(data_list, init_scores, max_iter, tol, verbose)

########## custome function finish defining

## legacy reading result, now automated in current MARIO
# read in match 1 wcct to human ifng
wcct1 = pd.read_csv("../temp/location/wcct_sr_matched1.csv")
wcct2 = pd.read_csv("../temp/location/wcct_sr_matched2.csv")
wcct3 = pd.read_csv("../temp/location/wcct_sr_matched3.csv")
wcct4 = pd.read_csv("../temp/location/wcct_sr_matched4.csv")
wcctfull_1=wcct1.copy()
wcctfull_1 = wcctfull_1.append([wcct2, wcct3 ,wcct4])
zh1 = pd.read_csv("../temp/location/Zac_human_il4_sr_matched1.csv")
zh2 = pd.read_csv("../temp/location/Zac_human_il4_sr_matched2.csv")
zh3 = pd.read_csv("../temp/location/Zac_human_il4_sr_matched3.csv")
zh4 = pd.read_csv("../temp/location/Zac_human_il4_sr_matched4.csv")
zhfull=zh1.copy()
zhfull = zhfull.append([zh2, zh3 ,zh4])

# read in match 2 wcct to rhesus ifng
wcct1 = pd.read_csv("../temp/location/wcct_sr_matched1.csv")
wcct2 = pd.read_csv("../temp/location/wcct_sr_matched2.csv")
wcct3 = pd.read_csv("../temp/location/wcct_sr_matched3.csv")
wcct4 = pd.read_csv("../temp/location/wcct_sr_matched4.csv")
wcctfull_2=wcct1.copy()
wcctfull_2 = wcctfull_2.append([wcct2, wcct3 ,wcct4])

zh1 = pd.read_csv("../temp/location/Zac_rheus_sr_matched1.csv")
zh2 = pd.read_csv("../temp/location/Zac_rheus_sr_matched2.csv")
zh3 = pd.read_csv("../temp/location/Zac_rheus_sr_matched3.csv")
zh4 = pd.read_csv("../temp/location/Zac_rheus_sr_matched4.csv")
zrfull=zh1.copy()
zrfull = zrfull.append([zh2, zh3 ,zh4])

# read in match 3 wcct to cyno ifng
wcct1 = pd.read_csv("../temp/location/wcct_sr_matched1.csv")
wcct2 = pd.read_csv("../temp/location/wcct_sr_matched2.csv")
wcct3 = pd.read_csv("../temp/location/wcct_sr_matched3.csv")
wcct4 = pd.read_csv("../temp/location/wcct_sr_matched4.csv")
wcctfull_3=wcct1.copy()
wcctfull_3 = wcctfull_3.append([wcct2, wcct3 ,wcct4])

zh1 = pd.read_csv("../temp/location/Zac_cyno_sr_matched1.csv")
zh2 = pd.read_csv("../temp/location/Zac_cyno_sr_matched2.csv")
zh3 = pd.read_csv("../temp/location/Zac_cyno_sr_matched3.csv")
zh4 = pd.read_csv("../temp/location/Zac_cyno_sr_matched4.csv")
zcfull=zh1.copy()
zcfull = zcfull.append([zh2, zh3 ,zh4])
# need cells (wcct) been matched to all three datasets
wc1_2=np.intersect1d(wcctfull_1["Unnamed: 0.1"], wcctfull_2["Unnamed: 0.1"],return_indices=True)
wc1_3=np.intersect1d(wcctfull_1["Unnamed: 0.1"], wcctfull_3["Unnamed: 0.1"],return_indices=True)
shared = np.intersect1d(wc1_2[0],wc1_3[0]) # wcct cells matched to all other three datasets
# get wcct sharerd values
wcctfull_1_shared=wcctfull_1.loc[wcctfull_1["Unnamed: 0.1"].isin(shared)]
wcctfull_2_shared=wcctfull_2.loc[wcctfull_2["Unnamed: 0.1"].isin(shared)]
wcctfull_3_shared=wcctfull_3.loc[wcctfull_3["Unnamed: 0.1"].isin(shared)]
# stupid: wcctfull1shared are all the same lololol
# get zh zr zc shared value
zh_shared=zhfull.loc[wcctfull_1["Unnamed: 0.1"].isin(shared)]
zr_shared=zrfull.loc[wcctfull_2["Unnamed: 0.1"].isin(shared)]
zc_shared=zcfull.loc[wcctfull_3["Unnamed: 0.1"].isin(shared)]
# need to get value and scale
wcctfull_1_shared_scale=normalize(wcctfull_1_shared.iloc[:,3:44]) # only get the values
wcctfull_2_shared_scale=normalize(wcctfull_2_shared.iloc[:,3:44])
wcctfull_3_shared_scale=normalize(wcctfull_3_shared.iloc[:,3:44])
# not necessary only need first line
zh_shared_scale=normalize(zh_shared.iloc[:,3:42]) # only get the values
zr_shared_scale=normalize(zr_shared.iloc[:,3:42])
zc_shared_scale=normalize(zc_shared.iloc[:,3:42])
# get the annotation for downstream analysis
all_labels = np.concatenate([wcctfull_1_shared['cluster.sr'], zh_shared['cluster.sr'], zr_shared['cluster.sr'],zc_shared['cluster.sr']])
# dont need that many cell other wise will kill the server when making tsne plots
idx=np.random.choice(65974, 20000, replace=False) # downsample to 20k cells each dataset
wcctfull_1_shared_scale_small=wcctfull_1_shared_scale.iloc[idx,:]
zh_shared_scale_small=zh_shared_scale.iloc[idx,:]
zr_shared_scale_small=zr_shared_scale.iloc[idx,:]
zc_shared_scale_small=zc_shared_scale.iloc[idx,:]
# make sure columns are not constant
wcctfull_1_shared_scale2 = wcctfull_1_shared_scale.loc[:, wcctfull_1_shared_scale.std() > 0]
zh_shared_scale2 = zh_shared_scale.loc[:, zh_shared_scale.std() > 0]
zr_shared_scale2 = zr_shared_scale.loc[:, zr_shared_scale.std() > 0]
zc_shared_scale2 = zc_shared_scale.loc[:, zc_shared_scale.std() > 0]
# legacy gcca steps, automatic in current MARIO
# gcca initiation
score_list = gcca_init([wcctfull_1_shared_scale_small,
                        zh_shared_scale_small, zr_shared_scale_small,zc_shared_scale_small],
                       normalization=False,n_components=10) # only 10 components
Xc, Yc, Zc, Uc = score_list
# gcca refine
all_scores = gcca([wcctfull_1_shared_scale_small.to_numpy(), zh_shared_scale_small.to_numpy(), zr_shared_scale_small.to_numpy(),zc_shared_scale_small.to_numpy()], n_components=10)
# save out the gcca scores
x,y,z,u=all_scores
all_cca=np.vstack([x,y,z,u])
np.savetxt("../temp/location/gccaREF_10.csv", all_cca, delimiter=",")
# also save out the matched files
wcctfull_1_shared_sub=wcctfull_1_shared.iloc[idx,:] # cells matched in all datasets
zh_shared_sub=zh_shared.iloc[idx,:]
zr_shared_sub=zr_shared.iloc[idx,:]
zc_shared_sub=zc_shared.iloc[idx,:]
wcctfull_1_shared_sub.to_csv("../temp/location/wcctfull_shared_sub.csv")
zh_shared_sub.to_csv("../temp/location/zh_shared_sub.csv")
zr_shared_sub.to_csv("../temp/location/zr_shared_sub.csv")
zc_shared_sub.to_csv("../temp/location/zc_shared_sub.csv")
## these files are later concatenate into the matched_files in the data folder
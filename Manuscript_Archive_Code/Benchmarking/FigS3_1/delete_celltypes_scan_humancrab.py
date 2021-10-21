import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanorama
# some parameters
path0 = "../Benchmarking/data/delete_cell_types/human_crab/sxc/"
path = "../Benchmarking/data/delete_cell_types/human_crab/"
x = pd.read_csv(path0+"orig_x.csv")
y = pd.read_csv(path0+"orig_y.csv")
## some parameters
celltype_2_delete = ['B','NK','Neutrophil','CD8 T','CD4 T','cMC']
#feature_2_delete = []
path_out = path + "scan/"
celltype_2_delete.insert(0, 'Full')
## record the labels
x_labels = x['label']
y_labels = y['label']
x = x.drop('label',1)
y = y.drop('label',1)
res = {'1v1_acc': [],'prop_remain': [], 'celltype_error' : []}
metrics_fname = path_out + "metrics"
embeddings_fname = path_out + "embedding"
for i, feature in enumerate(celltype_2_delete):
    x_use = x
    y_use = y
    x_label_use = x_labels
    y_label_use = y_labels
    if feature != 'Full':
        y_use = y[y_labels != feature]
        y_label_use = y_labels[y_labels!=feature]
    #break
    shared_features = [x for x in x.columns if x in y.columns]
    # need some scaling
    x_sclae = x[shared_features] - np.mean(x[shared_features], axis=0)
    x_sclae = x_sclae / np.std(x_sclae, axis=0)
    x_sclae = x_sclae.to_numpy()
    y_sclae = y_use[shared_features] - np.mean(y_use[shared_features], axis=0)
    y_sclae = y_sclae / np.std(y_sclae, axis=0)
    y_sclae = y_sclae.to_numpy()
    # scanorama matching process
    # reduction values
    integrated, genes=scanorama.integrate([x_sclae,y_sclae],[shared_features,shared_features])
    # matching information
    alignments, matches = scanorama.find_alignments([x_sclae,y_sclae])
    # match results
    df = pd.DataFrame(np.array(list(matches[(0, 1)])), columns = ['x','y'])
    acc = np.sum(x_labels[df['x']].values == y_label_use.iloc[df['y'],].values) / df.shape[0]
    res['1v1_acc'].append(acc)
    remain_prop = len(np.unique(df['x'])) / x.shape[0]
    res['prop_remain'].append(remain_prop)
    dropcell = len(x_label_use[x_label_use == feature])
    dropcell_match = sum(x_labels[df['x']].values == feature)
    error_match = dropcell_match/dropcell
    res['celltype_error'].append(error_match)
pd.DataFrame.from_dict(res).to_csv("{}.csv".format(metrics_fname), index_label='index')
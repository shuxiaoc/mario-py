import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanorama

# some parameters
path0 = "../Benchmarking/data/feature_deletion/felix/sxc/"
path = "../Benchmarking/data/feature_deletion/felix/"
x = pd.read_csv(path0+"orig_x.csv")
y = pd.read_csv(path0+"orig_y.csv")
## some parameters
feature_2_delete = ["CD11b",
        "CD127",
        "CD16",
        "CD19",
        "CD25",
        "CD27",
        "CD3",
        "CD4",
        "CD45RA"]
#feature_2_delete = []
path_out = path + "scan/"
feature_2_delete.insert(0, 'Full')
## record the labels
x_labels = x['label']
y_labels = y['label']
x = x.drop('label',1)
y = y.drop('label',1)

res = {'1v1_acc': [],'prop_remain': []}
metrics_fname = path_out + "metrics"
embeddings_fname = path_out + "embedding"
for i, feature in enumerate(feature_2_delete):
    if feature != 'Full':
            x = x.drop(feature, 1)
            y = y.drop(feature, 1)
    # get shared feauters in current loop
    shared_features = [x for x in x.columns if x in y.columns]
    # need some scaling
    x_sclae = x[shared_features] - np.mean(x[shared_features], axis=0)
    x_sclae = x_sclae / np.std(x_sclae, axis=0)
    x_sclae = x_sclae.to_numpy()
    y_sclae = y[shared_features] - np.mean(y[shared_features], axis=0)
    y_sclae = y_sclae / np.std(y_sclae, axis=0)
    y_sclae = y_sclae.to_numpy()
    # scanorama matching process
    # reduction values
    integrated, genes=scanorama.integrate([x_sclae,y_sclae],[shared_features,shared_features])
    # matching information
    alignments, matches = scanorama.find_alignments([x_sclae,y_sclae])
    # match results
    df = pd.DataFrame(np.array(list(matches[(0, 1)])), columns = ['x','y'])
    acc = np.sum(x_labels[df['x']].values == y_labels[df['y']].values) / df.shape[0]
    res['1v1_acc'].append(acc)
    remain_prop = len(np.unique(df['x'])) / x.shape[0]
    res['prop_remain'].append(remain_prop)
    # save the result csv
    #metrics_fname = path_out + "metrics"
    #pd.DataFrame.from_dict(res).to_csv("{}.csv".format(metrics_fname), index_label='index')
    # embedding 
    embedding1 = integrated[0]
    embedding2 = integrated[1]
    x_embed = pd.DataFrame(embedding1)
    y_embed = pd.DataFrame(embedding2)
    
    x_embed.to_csv("{}_x{}.csv".format(embeddings_fname, i), index=False)
    y_embed.to_csv("{}_y{}.csv".format(embeddings_fname, i), index=False)
    #print(acc)
    #print(remain_prop)
pd.DataFrame.from_dict(res).to_csv("{}.csv".format(metrics_fname), index_label='index')

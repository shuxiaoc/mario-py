import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanorama
# some parameters
path0 = "../Benchmarking/data/noise_addition/felix/sxc/"
path = "../Benchmarking/data/noise_addition/felix/"
x = pd.read_csv(path0+"orig_x.csv")
y = pd.read_csv(path0+"orig_y.csv")
## some parameters
noise_2_add = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]
#feature_2_delete = []
path_out = path + "scan/"
noise_2_add.insert(0, 0)
## record the labels
x_labels = x['label']
y_labels = y['label']
x = x.drop('label',1)
y = y.drop('label',1)
res = {'1v1_acc': [],'prop_remain': []}
metrics_fname = path_out + "metrics"
embeddings_fname = path_out + "embedding"
for i, feature in enumerate(noise_2_add):
    shared_features = [x for x in x.columns if x in y.columns]
    # need some scaling
    x_sclae = x[shared_features] - np.mean(x[shared_features], axis=0)
    x_sclae = x_sclae / np.std(x_sclae, axis=0)
    x_sclae = x_sclae.to_numpy()
    y_sclae = y[shared_features] - np.mean(y[shared_features], axis=0)
    y_sclae = y_sclae / np.std(y_sclae, axis=0)
    y_sclae = y_sclae.to_numpy()
    ## add noise
    x_sclae = x_sclae + np.random.normal(scale=feature, size=x_sclae.shape)
    y_sclae = y_sclae + np.random.normal(scale=feature, size=y_sclae.shape)
    # scanorama matching process
    # reduction values
    integrated, genes=scanorama.integrate([x_sclae,y_sclae],[shared_features,shared_features])
    # matching information
    alignments, matches = scanorama.find_alignments([x_sclae,y_sclae])
    # match results
    df = pd.DataFrame(np.array(list(matches[(0, 1)])), columns = ['x','y'])
    acc = np.sum(x_labels[df['x']].values == y_labels.iloc[df['y'],].values) / df.shape[0]
    res['1v1_acc'].append(acc)
    remain_prop = len(np.unique(df['x'])) / x.shape[0]
    res['prop_remain'].append(remain_prop)
    print(acc)
    print(remain_prop)

pd.DataFrame.from_dict(res).to_csv("{}.csv".format(metrics_fname), index_label='index')
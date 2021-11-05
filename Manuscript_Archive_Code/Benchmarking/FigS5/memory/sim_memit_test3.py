import numpy as np
import pandas as pd
import logging
import os
from memory_profiler import profile

from matching import CellMatching

@profile
def run_menit(df1, df2, subsample_prop=0.2, sparsity=1000, mode='sparse'):
    n1 = int(df1.shape[0] * subsample_prop)
    n2 = int(df2.shape[0] * subsample_prop)
    df1_subsampled = df1.iloc[np.random.choice(df1.shape[0], n1, replace=False), :]
    df2_subsampled = df2.iloc[np.random.choice(df2.shape[0], n2, replace=False), :]

    # construct the object
    match = CellMatching(df1_subsampled, df2_subsampled, normalization=True)
    match.specify_matching_params(1)

    # ovlp mathing
    _ = match.compute_dist_ovlp(20)
    _, min_sparsity = match.search_minimum_sparsity(
        match.dist['ovlp'], slackness=100
    )
    logging.info("min_sparsity={}".format(min_sparsity))
    _ = match.compute_dist_ovlp(20)
    match.match_cells("ovlp", sparsity=sparsity, mode=mode)

    # all-feature matching
    _ = match.compute_dist_all('ovlp', 20)
    match.match_cells("all", sparsity=sparsity, mode=mode)


if __name__ == "__main__":
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    np.random.seed(2667)
    df1 = pd.read_csv('data/murineCodex_forSIM.csv')
    df2 = pd.read_csv('data/murineCite_forSIM_80k.csv')
    n1, _ = df1.shape
    n2, _ = df2.shape
    df1 = df1.iloc[np.random.choice(n1, 20000, replace=False), :]
    df2 = df2.iloc[np.random.choice(n2, 80000, replace=False), :]

    # preprocessing to make two dfs in standard format: they contain cells, features and their labels
    df1 = df1.drop(['Unnamed: 0', 'cluster.info'], 1)
    df2 = df2.drop(['Unnamed: 0', 'cluster.info'], 1)

    run_menit(df1, df2, subsample_prop=1, sparsity=200, mode='sparse')
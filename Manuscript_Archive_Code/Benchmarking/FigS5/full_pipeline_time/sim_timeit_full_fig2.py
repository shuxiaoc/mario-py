import numpy as np
import pandas as pd
import logging
import os
import timeit

from matching import CellMatching

def run_timeit_full(df1, df2, fname, n1_vec, n2_vec, nsim=20,
                    n_components_ovlp=20, n_components_all=20,
                    top_k=10, n_wts=4, filter_iter = 20,
                    n_clusters=10, n_components_filter=10, bad_prop=0.1):
    """
    Construct a pandas dataframe with the following columns:
        'n1': sample size of df1
        'n2': sample size of df2
        'time': average time (over nsim runs) taken to the whole pipeline

    Save this dataframe to fname/timeit_full.csv

    Parameters
    ----------
    df1: a pandas dataframe with cells and their features
    df2: another pandas dataframe with cells and their features
    fname: results are saved at fname/timeit_full.csv
    n1_vec: each element is a specific number of cells in df1
    n2_vec: each element is a specific number of cells in df2
    nsim: number of repetitive runs
    n_components_ovlp: number of SVD components when doing initial matching
    n_components_all: number of CCA components when doing refined matching
    top_k: the number of top CCA scores, whose avg will be taken as a proxy for matching quality
    n_wts: number of different weights to try when doing interpolation
    filter_iter: number of iterations doing filtering
    n_clusters: number of clusters when doing filtering
    n_components_filter: number of SVD components when doing filtering
    bad_prop: approximate bad proportion when doing filtering
    """
    res = {'n1': [], 'n2': [], 'time': []}

    def record_time(n1, n2, t):
        res['n1'].append(n1)
        res['n2'].append(n2)
        res['time'].append(t)

    for n1, n2 in zip(n1_vec, n2_vec):
        logging.info('n1 = {}, n2 = {}...'.format(n1, n2))
        df1_subsampled = df1.iloc[np.random.choice(df1.shape[0], n1, replace=False), :]
        df2_subsampled = df2.iloc[np.random.choice(df2.shape[0], n2, replace=False), :]

        ##### start time  #####
        start = timeit.default_timer()
        for _ in range(nsim):
            # construct the object
            match = CellMatching(df1_subsampled, df2_subsampled, normalization=True)
            match.specify_matching_params(1)
            # match using overlapping features
            logging.info('Matching using overlapping features...')
            _ = match.compute_dist_ovlp(n_components_ovlp)
            # search for minimum sparsity
            _, min_sparsity = match.search_minimum_sparsity(
                match.dist['ovlp'], slackness=100
            )
            # match using overlapping features
            logging.info('Matching using overlapping features...')
            # use min sparsity + 100 for this simulation
            _ = match.match_cells('ovlp', sparsity=min_sparsity+100)
            # match using all the features
            logging.info('Matching using all features...')
            _ = match.compute_dist_all('ovlp', n_components_all)
            _ = match.match_cells('all', sparsity=min_sparsity+100)
            # no need to run matchability test
            # find the best interpolated matching
            logging.info('Finding the best interpolated matching...')
            wt, _ = match.interpolate(n_wts, top_k, verbose=False)
            logging.info('The best weight is {}.'.format(wt))
            # filter bad matches
            logging.info('Filtering bad matched pairs...')
            _ = match.filter_bad_matches('wted', n_clusters, n_components_filter, bad_prop,
                                     max_iter=filter_iter, tol=1e-4, verbose=False)

            # calculate embeddings and save them
            _, cca = match.fit_cca(match.matching['final'], n_components_all)
            _ = cca.transform(match.df1, match.df2)

        ######## time end ########
        full_time = (timeit.default_timer() - start) / nsim
        
        # save the result
        record_time(n1, n2, full_time)
    # save the metrics
    pd.DataFrame.from_dict(res).to_csv("{}/timeit_full.csv".format(fname), index_label='index')
    logging.info('Done!')


if __name__ == "__main__":
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    # path to save all simulation results
    # path = '.'
    # os.makedirs(path, exist_ok=True)

    np.random.seed(2667)
    df1 = pd.read_csv('data/bmcite_forSIM.csv')
    df2 = pd.read_csv('data/levine32_forSIM.csv')

    # preprocessing to make two dfs in standard format: they contain cells, features and their labels
    df1 = df1.drop(['Unnamed: 0', 'cluster.info'], 1)
    df2 = df2.drop(['Unnamed: 0', 'cluster.info'], 1)

    logging.info('Start simulation')

    n1_vec = [2000, 4000, 6000, 8000, 10000]
    n2_vec = [8000, 16000, 24000, 32000, 40000]
    run_timeit_full(df1, df2, fname='.', n1_vec=n1_vec, n2_vec=n2_vec, nsim=2,
                    n_components_ovlp=20, n_components_all=20,
                    top_k=10, n_wts=4, filter_iter=10,
                    n_clusters=10, n_components_filter=10, bad_prop=0.1)
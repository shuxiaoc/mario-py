import numpy as np
import pandas as pd
import logging
import os
import timeit

from matching import CellMatching

def run_timeit(df1, df2, fname, subsample_prop=[0.1, 0.5, 1], sparsity_trials=10, nsim=50):
    """
    Construct a pandas dataframe with the following columns:
        'n1': sample size of df1
        'n2': sample size of df2
        'sparsity': sparsity level (smaller, more sparse)
        'mode': sparse or dense matching used
        'time_ovlp': average time (over nsim runs) taken to do overlapping matching
        'time_all': average time (over nsim runs) taken to do all-feature matching

    Save this dataframe to fname/timeit.csv

    Parameters
    ----------
    df1: a pandas dataframe with cells and their features
    df2: another pandas dataframe with cells and their features
    fname: results are saved at fname/timeit.csv
    subsample_prop: sub-sample proportion for df1 and df2
    sparsity_trials: evenly spaced integers from min_sparsity to n2 will be tried
    nsim: number of repetitive runs
    """
    res = {'n1': [], 'n2': [], 'sparsity': [], 'mode': [], 'time_ovlp': [], 'time_all': []}

    def record_time(n1, n2, sparsity, mode, time_ovlp, time_all):
        res['n1'].append(n1)
        res['n2'].append(n2)
        res['sparsity'].append(sparsity)
        res['mode'].append(mode)
        res['time_ovlp'].append(time_ovlp)
        res['time_all'].append(time_all)

    for prop in subsample_prop:
        logging.info('Subsample proportion = {}...'.format(prop))
        # drop a feature
        n1 = int(df1.shape[0] * prop)
        n2 = int(df2.shape[0] * prop)
        df1_subsampled = df1.iloc[np.random.choice(df1.shape[0], n1, replace=False), :]
        df2_subsampled = df2.iloc[np.random.choice(df2.shape[0], n2, replace=False), :]

        # construct the object
        match = CellMatching(df1_subsampled, df2_subsampled, normalization=True)
        match.specify_matching_params(1)

        # match using overlapping features
        logging.info('Matching using overlapping features...')
        _ = match.compute_dist_ovlp(20)
        # search for minimum sparsity
        _, min_sparsity = match.search_minimum_sparsity(
            match.dist['ovlp'], slackness=100
        )

        sparsity_lst = np.linspace(min_sparsity, n1, num=sparsity_trials, dtype=int) # sparsity larger then n1 is just dense

        for sparsity in sparsity_lst:
            logging.info('Sparsity = {}...'.format(sparsity))
            # log the time for ovlp matching
            _ = match.compute_dist_ovlp(20)
            # sparse matching
            start = timeit.default_timer()
            for _ in range(nsim):
                match.match_cells("ovlp", sparsity=sparsity, mode="sparse")
            time_ovlp_sparse = (timeit.default_timer() - start) / nsim
            # dense matching
            start = timeit.default_timer()
            for _ in range(nsim):
                match.match_cells("ovlp", sparsity=sparsity, mode="dense")
            time_ovlp_dense = (timeit.default_timer() - start) / nsim

            # log the time for all-feature matching
            _ = match.compute_dist_all('ovlp', 20)
            # note that min_sparsity for ovlp and all matching may be different
            # if current sparsity does not work, we find the minimum valid new sparsity
            try:
                match.match_cells("all", sparsity=sparsity, mode="sparse")
            except ValueError:
                logging.info("current sparsity is too small for all-feature matching,"
                             "do matching with a minimum valid new sparsity that is >= the old sparsity")
                _, sparsity = match.search_minimum_sparsity(
                    match.dist['all'], slackness=0, init_sparsity=sparsity
                )

            # sparse matching
            start = timeit.default_timer()
            for _ in range(nsim):
                match.match_cells("all", sparsity=sparsity, mode="dense") # no need sparse for all
            time_all_sparse = (timeit.default_timer() - start) / nsim
            # dense matching
            start = timeit.default_timer()
            for _ in range(nsim):
                match.match_cells("all", sparsity=sparsity, mode="dense") # no need sparse for all matching
            time_all_dense = (timeit.default_timer() - start) / nsim

            # save the result
            record_time(n1, n2, sparsity, "sparse", time_ovlp_sparse, time_all_sparse)
            record_time(n1, n2, sparsity, "dense", time_ovlp_dense, time_all_dense)

    # save the metrics
    pd.DataFrame.from_dict(res).to_csv("{}/timeit.csv".format(fname), index_label='index')
    logging.info('Done!')


if __name__ == "__main__":
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    # path to save all simulation results
    path = 'data/timeit_fig4'
    os.makedirs(path, exist_ok=True)

    np.random.seed(2667)
    df1 = pd.read_csv('data/murineCodex_forSIM.csv')
    df2 = pd.read_csv('data/murineCite_forSIM_40k.csv')
    n1, _ = df1.shape
    n2, _ = df2.shape
    df1 = df1.iloc[np.random.choice(n1, 10000, replace=False), :]
    df2 = df2.iloc[np.random.choice(n2, 40000, replace=False), :]

    # preprocessing to make two dfs in standard format: they contain cells, features and their labels
    df1 = df1.drop('Unnamed: 0', 1)
    df2 = df2.drop('Unnamed: 0', 1)
    df1 = df1.rename(columns={'cluster.info': 'label'})
    df2 = df2.rename(columns={'cluster.info': 'label'})

    # df1 and df2 should only contain cells and features
    df1 = df1.drop('label', 1)
    df2 = df2.drop('label', 1)

    logging.info('Start simulation')

    run_timeit(df1, df2, path, subsample_prop=[0.2,0.4,0.6,0.8,1], sparsity_trials=5, nsim=3)
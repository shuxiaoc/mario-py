import numpy as np
import pandas as pd
import logging
import os

from matching import CellMatching, eval_matching_accuracy
from utils import normalize


def run_add_noise(df1, df2, x_labels, y_labels, noise_levels,
                        metrics_fname, embeddings_fname,
                        n_components_ovlp=20, n_components_all=20,
                        top_k=10, n_test=20, n_wts=10,
                        n_clusters=10, n_components_filter=10, bad_prop=0.1, knn=5):
    """
    Construct a pandas dataframe with rows being noise_levels and with the following features:
        'pval_ovlp': p-value from the matchability test for overlapping features
        'pval_all': p-value from the matchability test for all the features
        'ovlp_acc': matching accuracy using only overlapping features
        'all_acc': matching accuracy using all the features
        'wted_acc': matching accuracy using interpolation
        'final_acc': matching accuracy after filtering
        'prop_remain': proportion of number of remaining cells
        'transfer_acc': k-NN label transfer accuracy
    Save this dataframe to metrics_faname.csv
    Then obtain embeddings of df1 and df2 and save them to embeddings_fname_xi.csv and embeddings_fname_yi.csv,
    where i ranges from 0, 1, ..., len(noise_levels). Here i=0 means no noises have been added.

    Parameters
    ----------
    df1: a pandas dataframe with cells and their features
    df2: another pandas dataframe with cells and their features
    x_labels: a list of cluster labels for df1
    y_labels: a list of cluster labels for df2
    noise_levels: a list of noises, will add i.i.d. N(0, noise) to both df1 and df2
    metrics_fname: file name for saving the metrics
    embeddings_fname: file name for saving the embedddings
    n_components_ovlp: number of components to keep when doing stacked SVD
    n_components_all: number of CCA components to keep when calculating distance using all the features,
                      this is also the n_component when calculating the embeddings
    top_k: whenever we evaluate the quality of a matching (matchability test, weight selection),
           we take the median of the top_k canonical correlations as a proxy
    n_test: number of iterations in machability test
    n_wts: number of weights to try when computing the interpolated matching
    n_clusters, n_components_filter, bad_prop: parameters for filtering
    knn: number of nearest neighbors to count when doing label transfer
    """
    noise_levels.insert(0, 0)
    res = {'pval_ovlp': [], 'pval_all': [], 'ovlp_acc': [], 'all_acc': [], 'wted_acc': [], 'final_acc': [],
           'prop_remain': [], 'transfer_acc': []}
    for i, noise in enumerate(noise_levels):
        logging.info('Noise level = {}...'.format(noise))

        # normalize then add noise
        df1 = normalize(df1)
        df2 = normalize(df2)
        df1 = df1 + np.random.normal(scale=noise, size=df1.shape)
        df2 = df2 + np.random.normal(scale=noise, size=df2.shape)

        # construct the object
        match = CellMatching(df1, df2, normalization=False)
        match.specify_matching_params(1)

        # match using overlapping features
        logging.info('Matching using overlapping features...')
        _ = match.compute_dist_ovlp(n_components_ovlp)
        _ = match.match_cells('ovlp', sparsity=None)
        res['ovlp_acc'].append(eval_matching_accuracy(x_labels, y_labels, match.matching['ovlp'], 'maj'))

        # match using all the features
        logging.info('Matching using all features...')
        _ = match.compute_dist_all('ovlp', n_components_all)
        _ = match.match_cells('all', sparsity=None)
        res['all_acc'].append(eval_matching_accuracy(x_labels, y_labels, match.matching['all'], 'maj'))

        # run matchability test
        logging.info('Running matchability test...')
        # ignore old version matchable test
        #pval_ovlp, pval_all = match.matchable(n_test, top_k, verbose=False)
        pval_ovlp, pval_all = match.matchable(n_test, top_k, flip_prob = 0.2,
                                              subsample_prop =1,
                                              subsample_rounds =1, verbose=False)
        logging.info('p-value for matching using overlapping features is {}.'.format(pval_ovlp))
        logging.info('p-value for matching using all features is {}.'.format(pval_all))
        res['pval_ovlp'].append(pval_ovlp)
        res['pval_all'].append(pval_all)

        # fin the best interpolated matching
        logging.info('Finding the best interpolated matching...')
        wt, _ = match.interpolate(n_wts, top_k, verbose=False)
        logging.info('The best weight is {}.'.format(wt))
        res['wted_acc'].append(eval_matching_accuracy(x_labels, y_labels, match.matching['wted'], 'maj'))

        # filter bad matches
        logging.info('Filtering bad matched pairs...')
        _ = match.filter_bad_matches('wted', n_clusters, n_components_filter, bad_prop,
                                     max_iter=20, tol=1e-4, verbose=False)

        # calculate metrics
        logging.info('Calculating metrics...')
        res['final_acc'].append(eval_matching_accuracy(x_labels, y_labels, match.matching['final'], 'maj'))
        prop_remain = np.sum([1 for ii in range(match.n1) if len(match.matching['final'][ii]) != 0])
        prop_remain /= match.n1
        res['prop_remain'].append(prop_remain)

        # do label transfer
        wt = match.best_wt
        dist_argsort = np.argsort((1 - wt) * match.dist['ovlp'] + wt * match.dist['all'], axis=1)
        match_maj_vote = dist_argsort[:, :knn]
        match_maj_vote = [match_maj_vote[ii, :].tolist() for ii in range(match.n1)]
        res['transfer_acc'].append(eval_matching_accuracy(x_labels, y_labels, match_maj_vote, 'maj'))

        logging.info('ovlp_acc={}, all_acc={}, wted_acc={}, final_acc={}, prop_remain={}, transfer_acc={}'.format(
            res['ovlp_acc'][-1], res['all_acc'][-1], res['wted_acc'][-1], res['final_acc'][-1],
            res['prop_remain'][-1], res['transfer_acc'][-1]
        ))

        # calculate embeddings and save them
        _, cca = match.fit_cca(match.matching['final'], n_components_all)
        x_embed, y_embed = cca.transform(match.df1, match.df2)
        x_embed = pd.DataFrame(x_embed)
        # x_embed['label'] = x_labels
        y_embed = pd.DataFrame(y_embed)
        # y_embed['label'] = y_labels
        x_embed.to_csv("{}_x{}.csv".format(embeddings_fname, i), index=False)
        y_embed.to_csv("{}_y{}.csv".format(embeddings_fname, i), index=False)

    # save the metrics
    pd.DataFrame.from_dict(res).to_csv("{}.csv".format(metrics_fname), index_label='index')
    logging.info('Done!')


if __name__ == "__main__":
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    # path to save all simulation results
    path = 'data/noise_addition/sxc_bm'
    os.makedirs(path, exist_ok=True)

    np.random.seed(2667)
    df1 = pd.read_csv('data/bmcite_forSIM.csv')
    df2 = pd.read_csv('data/levine32_forSIM.csv')
    n1, _ = df1.shape
    n2, _ = df2.shape
    df1 = df1.iloc[np.random.choice(n1, 5000, replace=False), :]
    df2 = df2.iloc[np.random.choice(n2, 20000, replace=False), :]

    # preprocessing to make two dfs in standard format: they contain cells, features and their labels
    df1 = df1.drop('Unnamed: 0', 1)
    df2 = df2.drop('Unnamed: 0', 1)
    df1 = df1.rename(columns={'cluster.info': 'label'})
    df2 = df2.rename(columns={'cluster.info': 'label'})

    # save preprocessed data for later usage
    df1.to_csv('{}/orig_x.csv'.format(path), index=False)
    df2.to_csv('{}/orig_y.csv'.format(path), index=False)

    normalize(df1).to_csv('{}/orig_normalized_x.csv'.format(path), index=False)
    normalize(df2).to_csv('{}/orig_normalized_y.csv'.format(path), index=False)

    x_labels = df1['label'].to_list()
    y_labels = df2['label'].to_list()

    # df1 and df2 should only contain cells and features
    df1 = df1.drop('label', 1)
    df2 = df2.drop('label', 1)

    noise_levels = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

    logging.info('Start simulation, will add noises with standard deviation {} sequentially.'.format(noise_levels))

    metrics_fname = '{}/metrics'.format(path)
    embeddings_fname = '{}/embedding'.format(path)

    run_add_noise(df1, df2, x_labels, y_labels, noise_levels,
                        metrics_fname, embeddings_fname,
                        n_components_ovlp=10, n_components_all=20,
                        top_k=3, n_test=50, n_wts=10,
                        n_clusters=10, n_components_filter=10, bad_prop=0.2, knn=5)
import numpy as np
import pandas as pd
import logging
import os

from matching import CellMatching, eval_matching_accuracy


def run_sparsify(df1, df2, x_labels, y_labels,
                        metrics_fname, embeddings_fname,
                        sparsity_trials=10,
                        n_components_ovlp=20, n_components_all=20,
                        top_k=10, n_test=20, flip_prob=0.4, subsample_prop=1, subsample_rounds=1, n_wts=10,
                        n_clusters=10, n_components_filter=10, bad_prop=0.1, knn=5):
    """
    Construct a pandas dataframe with the following features:
        'sparsity': sparsity level
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
    where i ranges from 0, 1, ..., len(features_to_delete). Here i=0 means no features have been deleted.

    Parameters
    ----------
    df1: a pandas dataframe with cells and their features
    df2: another pandas dataframe with cells and their features
    x_labels: a list of cluster labels for df1
    y_labels: a list of cluster labels for df2
    metrics_fname: file name for saving the metrics
    embeddings_fname: file name for saving the embedddings
    sparsity_trials: try evenly spaced integers from minimum sparsity to full sparsity (i.e., dense)
    n_components_ovlp: number of components to keep when doing stacked SVD
    n_components_all: number of CCA components to keep when calculating distance using all the features,
                      this is also the n_component when calculating the embeddings
    top_k: whenever we evaluate the quality of a matching (matchability test, weight selection),
           we take the median of the top_k canonical correlations as a proxy
    n_test, flip_prob, subsample_prop, subsample_rounds: hyperparameters for the matchability test
    n_wts: number of weights to try when computing the interpolated matching
    n_clusters, n_components_filter, bad_prop: parameters for filtering
    knn: number of nearest neighbors to count when doing label transfer
    """
    res = {'sparsity': [], 'pval_ovlp': [], 'pval_all': [],
           'ovlp_acc': [], 'all_acc': [], 'wted_acc': [], 'final_acc': [],
           'prop_remain': [], 'transfer_acc': []}

    match = CellMatching(df1, df2, normalization=True)
    match.specify_matching_params(1)
    _ = match.compute_dist_ovlp(n_components_ovlp)
    _, min_sparsity = match.search_minimum_sparsity(
        match.dist['ovlp'], slackness=100
    )
    sparsity_lst = np.linspace(min_sparsity, match.n1 // 2, num=sparsity_trials, dtype=int)
    sparsity_lst = np.append(sparsity_lst, match.n1 // 2 + 1000) # add dense??
    
    for i, sparsity in enumerate(sparsity_lst):
        logging.info('sparsity={}'.format(sparsity))
        res['sparsity'].append(sparsity)
        # match using overlapping features
        logging.info('Matching using overlapping features...')
        _ = match.match_cells('ovlp', sparsity=sparsity)
        res['ovlp_acc'].append(eval_matching_accuracy(x_labels, y_labels, match.matching['ovlp'], 'maj'))

        # match using all the features
        logging.info('Matching using all features...')
        _ = match.compute_dist_all('ovlp', n_components_all)
        # note that min_sparsity for ovlp and all matching may be different
        # if current sparsity does not work, we find the minimum valid new sparsity
        try:
            match.match_cells("all", sparsity=sparsity)
        except ValueError:
            logging.info("current sparsity is too small for all-feature matching,"
                         "do matching with a minimum valid new sparsity that is >= the old sparsity")
            _, sparsity = match.search_minimum_sparsity(
                match.dist['all'], slackness=0, init_sparsity=sparsity
            )
            match.match_cells("all", sparsity=sparsity)
        res['all_acc'].append(eval_matching_accuracy(x_labels, y_labels, match.matching['all'], 'maj'))

        # run matchability test
        logging.info('Running matchability test...')
        #pval_ovlp, pval_all = match.matchable(n_test, top_k, flip_prob, subsample_prop, subsample_rounds, verbose=False)
        pval_ovlp = 0
        pval_all = 0
        
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
                                     max_iter=10, tol=1e-4, verbose=False)

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
    path = 'data/sparsify_fig2'
    os.makedirs(path, exist_ok=True)

    np.random.seed(2667)
    df1 = pd.read_csv('data/bmcite_forSIM.csv')
    df2 = pd.read_csv('data/levine32_forSIM.csv')
    n1, _ = df1.shape
    n2, _ = df2.shape
    df1 = df1.iloc[np.random.choice(n1, 10000, replace=False), :]
    df2 = df2.iloc[np.random.choice(n2, 40000, replace=False), :]

    # preprocessing to make two dfs in standard format: they contain cells, features and their labels
    df1 = df1.drop('Unnamed: 0', 1)
    df2 = df2.drop('Unnamed: 0', 1)
    df1 = df1.rename(columns={'cluster.info': 'label'})
    df2 = df2.rename(columns={'cluster.info': 'label'})

    # save preprocessed data for later usage
    df1.to_csv('{}/orig_x.csv'.format(path), index=False)
    df2.to_csv('{}/orig_y.csv'.format(path), index=False)

    x_labels = df1['label'].to_list()
    y_labels = df2['label'].to_list()

    # df1 and df2 should only contain cells and features
    df1 = df1.drop('label', 1)
    df2 = df2.drop('label', 1)

    logging.info('Start simulation.')

    metrics_fname = '{}/metrics'.format(path)
    embeddings_fname = '{}/embedding'.format(path)

    run_sparsify(df1, df2, x_labels, y_labels,
                 metrics_fname, embeddings_fname,
                 sparsity_trials=5,
                 n_components_ovlp=20, n_components_all=20,
                 top_k=10, n_test=20, flip_prob=0.4, subsample_prop=1, subsample_rounds=1, n_wts=10,
                 n_clusters=10, n_components_filter=10, bad_prop=0.1, knn=5)
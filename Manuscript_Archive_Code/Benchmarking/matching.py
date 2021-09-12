import numpy as np
import logging
from collections import Counter

from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.sparse.csgraph import maximum_bipartite_matching
from sklearn.utils.extmath import randomized_svd

from utils import cdist_correlation, normalize
from cca import get_cancor
from clustering import spectral_clustering, jr_kmeans


def _check_min_dist(dist_mat, min_dist):
    """
    Make sure min(dist_mat) >= min_dist

    Parameters
    ----------
    dist_mat: distance matrix
    min_dist: required minimum distance

    Returns
    -------
    dist_mat: a distance matrix with min(dist_mat) >= min_dist

    """
    current_min_dist = np.min(dist_mat)
    if current_min_dist < min_dist:
        dist_mat = dist_mat - current_min_dist + min_dist

    return dist_mat


def _search_minimum_sparsity(dist_mat, slackness, init_sparsity,
                             m_min, m_max, num_cells_to_use,
                             min_dist, verbose=True):
    """
    Use binary search to search for the minimum sparsity level k such that
        if dist_mat is trimmed to be a k-NN graph,
        a valid matching still exists.
    The search starts with k_left=1 and k_right=dist_mat.shape[0],
        and it is always true that:
        1) any k < k_left doesn't give a valid matching
        2) any k >= k_right gives a valid matching

    Parameters
    ----------
    dist_mat: a numpy array of distance matrix
    slackness: binary search terminates when k_right - k_left <= slackness,
               an exact binary search corresponds to slackness = 0
    init_sparsity: binary search starts from k=init_sparsity
    m_min: each cell in df1 is matched to at least m_min many cells in df2
    m_max: each cell in df1 is matched to at most m_max many cells in df2
    num_cells_to_use: total number of cells to use in df2
    min_dist: it must be true that minimum entry in dist_mat is >= min_dist

    Returns
    -------
    k_left: if sparsity<k_left, then there is no valid matching
    k_right: if sparsity>=k_right, then there is a valid matching
    """

    assert np.min(dist_mat) >= min_dist
    n1, n2 = dist_mat.shape

    num_sinks = max(n1 * m_max - num_cells_to_use, 0)

    # construct distance matrix that's ready for matching
    if m_max > 1:
        dist_mat = np.tile(dist_mat, (m_max, 1))

    argsort_res = np.argsort(dist_mat, axis=1)

    k_left = 1
    k_right = n2
    cnt = 0

    # start binary search
    while k_left < k_right - slackness:

        if verbose:
            logging.info('If sparsity>={}, then there is a valid matching; '
                         'if sparsity<{}, then there is no valid matching.'.format(k_right, k_left))
        if cnt == 0 and init_sparsity is not None:
            k = init_sparsity
        else:
            k = (k_left + k_right) // 2

        # construct k-NN graph from dist_mat
        # indices for the largest k entries
        largest_k = argsort_res[:, -(n2 - k):]
        # one means there is an edge and zero means there's no edge
        dist_bin = np.ones_like(dist_mat)
        dist_bin[np.arange(dist_mat.shape[0])[:, None], largest_k] = 0

        # add sinks if necessary
        if num_sinks > 0:
            # again, one means there is an edge,
            # and zero means there is no edge (i.e., the distance is np.inf)
            dist_sinks_bin = np.ones((dist_mat.shape[0], num_sinks))
            dist_sinks_bin[:m_min * n1, :] = 0
            dist_bin = np.concatenate((dist_bin, dist_sinks_bin), axis=1)

        dist_bin = csr_matrix(dist_bin)

        col_idx = maximum_bipartite_matching(dist_bin, perm_type='column')
        n_matched = np.sum(col_idx != -1)

        if n_matched == dist_bin.shape[0]:
            # the current k gives a valid matching
            # can trim more aggresively --> try a smaller k
            k_right = k
        else:
            # the current k doesn't give a valid matching
            # try a larger k
            k_left = k + 1

        cnt = cnt + 1

    # we know that k_right must give a valid matching
    if verbose:
        logging.info('If sparsity>={}, then there is a valid matching; '
                     'if sparsity<{}, then there is no valid matching.'.format(k_right, k_left))
    return k_left, k_right


def _get_matching_from_indices(col_idx, n1, n2, m_max):
    """
    Assume col_idx is obtained from min_weight_full_bipartite_matching
    with some dist_mat as input.
    And this dist_mat is organized such that the sinks are in dist_mat[:, :n2].
    This function calculates the matching from col_idx.

    Args:
        col_idx: output from min_weight_full_bipartite_matching
        n1: number of cells in df1
        n2: number of cells in df2
        m_max: each cell in df1 is matched to at most m_max cells in df2

    Returns:
       res: a list of (potentially variable length) lists;
            it holds that cell ii in df1 is matched to res[ii] in df2
    """
    if m_max == 1:
        # make sure the result is a list of (length-one) lists
        return [[col_idx[ii]] for ii in range(n1)]

    res = {ii: [] for ii in range(n1)}

    for kk in range(m_max):
        for ii in range(n1):
            candidate = col_idx[ii + (kk - 1) * n1]
            if candidate < n2:
                res[ii].append(candidate)

    return [res[ii] for ii in range(n1)]


def _match_cells(dist_mat, sparsity, m_min, m_max, num_cells_to_use,
                 min_dist, mode='auto'):
    """
    Given dist_mat, first trim it to a k-NN graph according to
    sparsity leve, then do a matching according to the specified parameters.

    Args:
        dist_mat: a numpy array of distance matrix
        sparsity: an integer k such that dist_mat will be trimmed into a k-NN graph
        m_min: each cell in df1 is matched to at least m_min many cells in df2
        m_max: each cell in df1 is matched to at most m_max many cells in df2
        num_cells_to_use: total number of cells to use in df2
        min_dist: it must be true that minimum entry in dist_mat is >= min_dist
        mode: 'sparse' means using min_weight_full_bipartite_matching;
              'dense' means using linear_sum_assignment;
              'auto': if sparsity<=n//2, use 'sparse', else use 'dense'

    Returns:
       res: a list of (potentially variable length) lists;
            it holds that cell ii in df1 is matched to res[ii] in df2
    """
    assert np.min(dist_mat) >= min_dist
    n1, n2 = dist_mat.shape

    if m_max > 1:
        dist_mat = np.tile(dist_mat, (m_max, 1))

    num_sinks = max(n1 * m_max - num_cells_to_use, 0)

    if sparsity is None:
        mode = 'dense'
    elif mode == 'auto':
        mode = 'sparse' if sparsity <= n2 // 2 else 'dense'

    infinity_placeholder = 0 if mode == 'sparse' else np.Inf

    # trim nodes if necessary
    if sparsity is not None and sparsity != n2:
        argsort_res = np.argsort(dist_mat, axis=1)
        largest_k = argsort_res[:, -(n2 - sparsity):]
        # make a copy because some operations are in-place
        # and may modify the original dist_mat
        dist_mat_cp = np.copy(dist_mat)
        dist_mat_cp[np.arange(dist_mat.shape[0])[:, None], largest_k] = infinity_placeholder
    else:
        dist_mat_cp = dist_mat

    # add sinks if necessary
    if num_sinks > 0:
        # we need make sure that
        # those sinks are favored compared to all other nodes
        dist_sinks = np.zeros((dist_mat.shape[0], num_sinks)) + min_dist / 100
        dist_sinks[:m_min * n1, :] = infinity_placeholder
        dist_mat_cp = np.concatenate((dist_mat_cp, dist_sinks), axis=1)

    if mode == 'sparse':
        dist_mat_cp = csr_matrix(dist_mat_cp)
        _, col_idx = min_weight_full_bipartite_matching(dist_mat_cp)
    elif mode == 'dense':
        _, col_idx = linear_sum_assignment(dist_mat_cp)
    else:
        raise NotImplementedError

    return _get_matching_from_indices(col_idx, n1, n2, m_max)


def _is_perfect_matching(matching):
    """
    Check if matching is perfect (i.e., all nodes are matched) or not
    Parameters
    ----------
    matching: if matching[ii] is an empty list, then ii is not matched to anyone

    Returns
    -------
    True if matching is perfect

    """
    matching_len = np.array([len(matching[ii]) for ii in range(len(matching))])
    matched_flag = (matching_len > 0)
    return all(matched_flag)


def _calculate_duplications(matching):
    """
    Calculate the number of i such that matching[i] appears more than one times.
    Then among those i's, calculate the average of times that they appear.
    Return those two numbers.
    """
    cells = []
    for matched in matching:
        cells.extend(matched)
    counter = Counter(cells)
    duplicated_cells = [cell for cell, freq in counter.items() if freq > 1]
    dup_prop = len(duplicated_cells) / len(cells)
    dup_freq = np.mean([counter[cell] for cell in duplicated_cells]) if len(duplicated_cells) > 0 else 0
    return dup_prop, dup_freq


def _batching(df, num_batches):
    """
    Given a pandas dataframe, randomly cut it into approximately equal-sized batches and return them as a list.
    Also returns the random permutations as a list.
    """
    n = df.shape[0]
    perm = np.random.permutation(n)
    size_each, remainder = divmod(n, num_batches)
    size_last_batch = size_each + remainder
    start_idx = 0
    df_lst = []
    perm_lst = []
    for i in range(num_batches):
        end_idx = start_idx + size_last_batch if i == num_batches-1 else start_idx + size_each
        indices = perm[start_idx:end_idx]
        perm_lst.append(indices)
        df_lst.append(df.iloc[indices])
        start_idx = end_idx
    return df_lst, perm_lst


def _stitch(perm_lst, matching_lst, n):
    """
    Assume perm_lst[i][j] should be matched to matching_lst[i][j].
    Return stitched_matching so that i should be matched to stitched_matching[i].
    """
    stitched_dict = {} # store perm_lst[i][j]: matching_lst[i][j]
    for i in range(len(perm_lst)):
        for j in range(len(perm_lst[i])):
            stitched_dict[perm_lst[i][j]] = matching_lst[i][j]
    return [stitched_dict[i] for i in range(n)]


class CellMatching(object):
    def __init__(self, df1, df2, normalization=True):
        self.min_dist = 1e-5

        # parameters related to datasets
        if normalization:
            self.df1 = normalize(df1)
            self.df2 = normalize(df2)
        else:
            self.df1 = df1
            self.df2 = df2
        self.n1, self.p1 = df1.shape
        self.n2, self.p2 = df2.shape
        assert self.n1 <= self.n2
        self.ovlp_features = [x for x in df1.columns if x in df2.columns]

        # hyper-parameters
        self.n_matched_per_cell = None
        self.m_min = None
        self.m_max = None
        self.num_cells_to_use = None
        self.num_sinks = None
        self.sparsity = {'ovlp': None, 'all': None}
        self.n_components = {'ovlp': None, 'all': None}

        # cache some results
        self.dist = {'ovlp': None, 'all': None}
        self.matching = {'ovlp': None, 'all': None, 'wted': None, 'final': None}
        self.best_wt = None
        self.stacked_svd = {'U': None, 's': None, 'Vh': None}
        self.ovlp_cancor = None
        self.ovlp_scores = {'x': None, 'y': None}

    def compute_dist_ovlp(self, n_components=10):
        """
        Compute distance matrix based on overlapping features.

        Parameters
        ----------
        n_components: number of SVD components to keep

        Returns
        -------
        dist_ovlp: distance matrix
        s: vector of singular values

        """
        if n_components > len(self.ovlp_features):
            logging.warning("n_components exceed the number of overlapping features,"
                            " set it to be the number of overlapping features.")
            n_components = len(self.ovlp_features)

        self.n_components['ovlp'] = n_components

        if not (
                self.stacked_svd['U'] is not None and self.stacked_svd['s'] is not None
                and self.stacked_svd['Vh'] is not None and len(self.stacked_svd['s']) >= n_components
        ):
            # Cached results are not valid, do SVD
            arr1 = normalize(self.df1[self.ovlp_features]).to_numpy()
            arr2 = normalize(self.df2[self.ovlp_features]).to_numpy()

            self.stacked_svd['U'], self.stacked_svd['s'], self.stacked_svd['Vh'] = \
                randomized_svd(np.concatenate((arr1, arr2), axis=0), n_components=n_components)
            if n_components == len(self.ovlp_features):
                dist_mat = cdist_correlation(arr1, arr2)
            else:
                svd1 = self.stacked_svd['U'][:self.n1, :] @ np.diag(self.stacked_svd['s']) @ self.stacked_svd['Vh']
                svd2 = self.stacked_svd['U'][self.n1:, :] @ np.diag(self.stacked_svd['s']) @ self.stacked_svd['Vh']
                dist_mat = cdist_correlation(svd1, svd2)
        else:
            # use cached results
            svd1 = self.stacked_svd['U'][:self.n1, :n_components] @ np.diag(self.stacked_svd['s'][:n_components]) \
                   @ self.stacked_svd['Vh'][:n_components, :]
            svd2 = self.stacked_svd['U'][self.n1:, :n_components] @ np.diag(self.stacked_svd['s'][:n_components]) \
                   @ self.stacked_svd['Vh'][:n_components, :]
            dist_mat = cdist_correlation(svd1, svd2)

        # make sure min_dist is at least self.min_dist
        dist_mat = _check_min_dist(dist_mat, self.min_dist)
        self.dist['ovlp'] = dist_mat
        return dist_mat, self.stacked_svd['s'][:n_components]

    def search_minimum_sparsity(self, dist_mat, slackness=200, init_sparsity=None, verbose=True):
        """
        Search for the minimum sparsity level so that we still get a valid matching.
        """
        return _search_minimum_sparsity(dist_mat, slackness, init_sparsity, self.m_min,
                                        self.m_max, self.num_cells_to_use, self.min_dist, verbose)

    def specify_matching_params(self, n_matched_per_cell):
        """
        Specify how many cells in Y data are to be matched with one cell in X data.

        Parameters
        ----------
        n_matched_per_cell: how many cells in Y data are to be matched with one cell in X data.
        """
        self.n_matched_per_cell = n_matched_per_cell
        if self.n1 * n_matched_per_cell > self.n2:
            raise ValueError("Not enough cells in Y data!")
        self._specify_matching_params(1, n_matched_per_cell, self.n1 * n_matched_per_cell)

    def _specify_matching_params(self, m_min, m_max, num_cells_to_use):
        assert m_min >= 1
        assert self.n1 <= num_cells_to_use
        assert num_cells_to_use <= self.n1 * m_max
        assert num_cells_to_use <= self.n2
        # m_max cannot be too large, otherwise a matching does not exist
        assert m_max <= num_cells_to_use - (self.n1 - 1)

        self.m_min = m_min
        self.m_max = m_max
        self.num_cells_to_use = num_cells_to_use
        self.num_sinks = max(self.n1 * m_max - num_cells_to_use, 0)

    def match_cells(self, dist_mat='ovlp', sparsity=None, mode='auto'):
        """
        Do cell matching.

        Parameters
        ----------
        dist_mat: either 'ovlp' or 'all', or a user-specified distance matrix.
        sparsity: number of nearest neighbors to keep in the distance matrix
        mode: 'sparse' means using min_weight_full_bipartite_matching;
              'dense' means using linear_sum_assignment;
              'auto': if sparsity<=n//2, use 'sparse', else use 'dense'

        Returns
        -------
        a list of (potentially variable length) lists; it holds that cell ii in df1 is matched to res[ii] in df2
        """

        if isinstance(dist_mat, str):
            if dist_mat not in self.dist or self.dist[dist_mat] is None:
                raise ValueError("Distance not found!")
            self.sparsity[dist_mat] = sparsity
            self.matching[dist_mat] = _match_cells(
                self.dist[dist_mat], sparsity, self.m_min, self.m_max, self.num_cells_to_use, self.min_dist, mode
            )
            return self.matching[dist_mat]
        else:
            return _match_cells(
                dist_mat, sparsity, self.m_min, self.m_max, self.num_cells_to_use, self.min_dist, mode
            )

    def _align_modalities(self, matching):
        """
        Args:
            matching: a list of (potentially variable length) lists;
        cell ii in df1 is matched to cell matching[ii] in df2

        Returns:
            X, Y: numpy arrays so that X[ii, :] is aligned with Y[ii, :]
        """
        assert len(matching) == self.n1

        # if cell ii in df1 is filtered out, then matching[ii] is an empty list
        X = []
        Y = []
        for ii in range(self.n1):
            if len(matching[ii]) == 0:
                continue

            X.append(self.df1.iloc[ii, :])
            if len(matching[ii]) == 1:
                Y.append(self.df2.iloc[matching[ii][0], :])
            else:
                Y.append(self.df2.iloc[matching[ii], :].mean(axis=0))

        X = np.array(X)
        Y = np.array(Y)

        return X, Y

    def fit_cca(self, matching, n_components=20, max_iter=3000):
        """
        Align df1 and df2 using matching, then fit a CCA.

        Parameters
        ----------
        matching: a user-specified list of (potentially variable length) lists,
                  cell ii in df1 is matched to cell matching[ii] in df2.
        n_components: number of components for CCA.
        max_iter: maximum iteration for CCA.

        Returns
        -------
        cancor: canonical correlations
        cca: CCA object
        """
        if n_components > min(self.p1, self.p2):
            logging.warning('n_components must be <= the dimensions of the two datasets, '
                            'set it to be equal to the minimum of the dimensions of the two datasets')
            n_components = min(self.p1, self.p2)

        X, Y = self._align_modalities(matching)
        cancor, cca = get_cancor(X, Y, n_components, max_iter)

        return cancor, cca

    def compute_dist_all(self, matching='ovlp', n_components=20, max_iter=5000):
        """
        Given matching, align df1 and df2, fit a CCA,
        then use CCA scores to get the distance matrix.

        Parameters
        ----------
        matching: either 'ovlp' or a list of (potentially variable length) lists;
                  cell ii in df1 is matched to cell matching[ii] in df2
        n_components: number of CCA components
        max_iter: max number of iterations when fitting CCA

        Returns
        -------
        dist_mat: distance matrix
        cancor: cononical correlations
        """
        if n_components <= 1:
            n_components = 2
            logging.warning('n_components must be at least 2, '
                            'set it to 2')

        if n_components > min(self.p1, self.p2):
            logging.warning('n_components must be <= the dimensions of the two datasets, '
                            'set it to be equal to the minimum of the dimensions of the two datasets')
            n_components = min(self.p1, self.p2)

        self.n_components['all'] = n_components

        if isinstance(matching, str) and matching == 'ovlp':
            if self.matching['ovlp'] is None:
                raise ValueError('Initial matching not found!')
            if not (
                    self.ovlp_scores['x'] is not None and self.ovlp_scores['y'] is not None
                    and self.ovlp_scores['x'].shape[1] >= n_components
            ):
                # cached results are not valid, do CCA
                cancor, cca = self.fit_cca(self.matching['ovlp'], n_components, max_iter)
                self.ovlp_cancor = cancor
                self.ovlp_scores['x'], self.ovlp_scores['y'] = cca.transform(self.df1, self.df2)
                dist_mat = cdist_correlation(self.ovlp_scores['x'], self.ovlp_scores['y'])
            else:
                # use cached results
                dist_mat = cdist_correlation(self.ovlp_scores['x'][:, :n_components],
                                             self.ovlp_scores['y'][:, :n_components])
                cancor = self.ovlp_cancor[:n_components]
        else:
            # use user-specified matching
            cancor, cca = self.fit_cca(matching, n_components, max_iter)
            df1_cca, df2_cca = cca.transform(self.df1, self.df2)
            dist_mat = cdist_correlation(df1_cca, df2_cca)

        self.dist['all'] = _check_min_dist(dist_mat, self.min_dist)
        return self.dist['all'], cancor

    def _matchable(self, n_sim=20, top_k=10, flip_prob=0.4, subsample_prop=0.3, verbose=True):
        """
        Calculate the p-value for the hypothesis that two datasets are matchable
        Parameters
        ----------
        n_sim: number of simulation rounds
        top_k: the median/mean of the top_k canonical correlations will be taken as the test statistic
        flip_prob: probability of sign flip
        subsample_prop: random subsample df1 and df2 to speed up the computation.

        Returns
        -------
        (pval_ovlp, pval_all): the p value for match_ovlp and match_all
        """
        cancor_ovlp_list = []
        cancor_all_list = []
        if self.n_components['ovlp'] is None or self.n_components['all'] is None:
            raise ValueError("Please do both initial and refined matching first before calling this function.")

        df1_subsampled = self.df1.iloc[np.random.choice(self.n1, int(self.n1 * subsample_prop), replace=False)]
        df2_subsampled = self.df2.iloc[np.random.choice(self.n2, int(self.n2 * subsample_prop), replace=False)]
        match_subsampled = CellMatching(df1_subsampled, df2_subsampled, normalization=True)
        match_subsampled.specify_matching_params(self.n_matched_per_cell)
        # do matching
        _ = match_subsampled.compute_dist_ovlp(self.n_components['ovlp'])
        try:
            _ = match_subsampled.match_cells('ovlp', self.sparsity['ovlp'], mode='auto')
        except ValueError:
            # too sparse, find the suitable sparsity level
            if verbose:
                logging.info(
                    'When using ovlp features, current sparsity config '
                    'is too sparse, finding a suitable sparsity level...'
                )
            _, new_sparsity = match_subsampled.search_minimum_sparsity(
                match_subsampled.dist['ovlp'], slackness=200, init_sparsity=self.sparsity['ovlp'] + 1, verbose=verbose
            )
            _ = match_subsampled.match_cells('ovlp', new_sparsity, mode='auto')

        # use all features
        _ = match_subsampled.compute_dist_all('ovlp', n_components=self.n_components['all'])
        try:
            _ = match_subsampled.match_cells('all', self.sparsity['all'], mode='auto')
        except ValueError:
            if verbose:
                logging.info(
                    'When using all features, current sparsity config '
                    'is too sparse, finding a suitable sparsity level...'
                )
            # too sparse, find the suitable sparsity level
            _, new_sparsity = match_subsampled.search_minimum_sparsity(
                match_subsampled.dist['all'], slackness=200, init_sparsity=self.sparsity['all'] + 1, verbose=verbose
            )
            _ = match_subsampled.match_cells('all', new_sparsity, mode='auto')

        # calculate test statistics
        cancor_ovlp_obs = np.mean(match_subsampled.fit_cca(match_subsampled.matching['ovlp'], n_components=top_k)[0])
        cancor_all_obs = np.mean(match_subsampled.fit_cca(match_subsampled.matching['all'], n_components=top_k)[0])

        for ii in range(n_sim):
            max_trial = 100
            # looks like we may get SVD not converge error
            # and it is due to linear algebra libraries in the linux server
            # so try max_trial times until success
            trial_idx = 0
            while trial_idx < max_trial:
                trial_idx += 1
                logging.warning("Now at the {}-th trial".format(trial_idx))
                try:
                    rand_signs_1 = 2 * np.random.binomial(1, 1-flip_prob, match_subsampled.n1) - 1
                    rand_signs_2 = 2 * np.random.binomial(1, 1-flip_prob, match_subsampled.n2) - 1
                    df1_flipped = (match_subsampled.df1.T * rand_signs_1).T
                    df2_flipped = (match_subsampled.df2.T * rand_signs_2).T

                    match_flipped = CellMatching(df1_flipped, df2_flipped, normalization=False)
                    match_flipped.specify_matching_params(self.n_matched_per_cell)

                    # use ovlp features
                    _ = match_flipped.compute_dist_ovlp(self.n_components['ovlp'])
                    try:
                        _ = match_flipped.match_cells('ovlp', self.sparsity['ovlp'], mode='auto')
                    except ValueError:
                        # too sparse, find the suitable sparsity level
                        if verbose:
                            logging.info(
                                'When using ovlp features, current sparsity config '
                                'is too sparse, finding a suitable sparsity level...'
                            )
                        _, new_sparsity = match_flipped.search_minimum_sparsity(
                            match_flipped.dist['ovlp'], slackness=200,
                            init_sparsity=self.sparsity['ovlp'] + 1, verbose=verbose
                        )
                        _ = match_flipped.match_cells('ovlp', new_sparsity, mode='auto')

                    # calculate the median/mean of top_k cancors
                    cancor_ovlp_list.append(np.mean(
                        match_flipped.fit_cca(match_flipped.matching['ovlp'], n_components=top_k)[0]
                    ))

                    # use all features
                    _ = match_flipped.compute_dist_all('ovlp', n_components=self.n_components['all'])
                    try:
                        _ = match_flipped.match_cells('all', self.sparsity['all'], mode='auto')
                    except ValueError:
                        if verbose:
                            logging.info(
                                'When using all features, current sparsity config '
                                'is too sparse, finding a suitable sparsity level...'
                            )
                        # too sparse, find the suitable sparsity level
                        _, new_sparsity = match_flipped.search_minimum_sparsity(
                            match_flipped.dist['all'], slackness=200,
                            init_sparsity=self.sparsity['all'] + 1, verbose=verbose
                        )
                        _ = match_flipped.match_cells('all', new_sparsity, mode='auto')

                    # calculate the median/mean of top_k cancors
                    cancor_all_list.append(np.mean(
                        match_flipped.fit_cca(match_flipped.matching['all'], n_components=top_k)[0]
                    ))
                except Exception as exception:
                    logging.warning("{} occurs at the {}/{}-th trial.".format(
                        type(exception).__name__, trial_idx, max_trial)
                    )
                    continue
                else:
                    # the code will arrive here only if try has succeeded
                    break
            if trial_idx == max_trial:
                raise ValueError("{} consecutive trials have all failed.".format(max_trial))
            else:
                logging.warning("Success at {}/{}-th trial".format(trial_idx, max_trial))

        # calculate p-values
        pval_ovlp = np.mean([cancor_ovlp_list[ii] >= cancor_ovlp_obs for ii in range(n_sim)])
        pval_all = np.mean([cancor_all_list[ii] >= cancor_all_obs for ii in range(n_sim)])
        return pval_ovlp, pval_all

    def matchable(self, n_sim=20, top_k=10, flip_prob=0.3, subsample_prop=1, subsample_rounds=1, verbose=True):
        """
        Calculate the p-value for the hypothesis that two datasets are matchable

        Parameters
        ----------
        n_sim: number of simulation rounds
        top_k: the median/mean of the top_k canonical correlations will be taken as the test statistic
        flip_prob: probability of sign flip
        subsample_prop: subsample to speed up computation
        subsample_rounds: the subsampling operation is done this many rounds
        verbose: print details if True

        Returns
        -------
        (pval_ovlp, pval_all): the p value for match_ovlp and match_all
        """
        pval_ovlp = 0
        pval_all = 0
        for i in range(subsample_rounds):
            if verbose:
                logging.info(f"Now at subsample round {i}...")
            curr_pval_ovlp, curr_pval_all = self._matchable(n_sim, top_k, flip_prob, subsample_prop, verbose)
            pval_ovlp += curr_pval_ovlp
            pval_all += curr_pval_all

        return pval_ovlp/subsample_rounds, pval_all/subsample_rounds

    def interpolate(self, n_wts=10, top_k=10, verbose=True):
        """
        Let wt_vec be an evenly spaced list from 0 to 1 with length n_wts.
        For each wt in wt_vec, do matching on (1-wt)*dist_ovlp + wt*dist_all,
        and select the best wt according to the median/mean of top_k canonical correlations.

        Parameters
        ----------
        n_wts: wt_vec is a evenly spaced list from 0 to 1 with length n_wts
        top_k: the median/mean of top_k canonical correlations is taken as the quality measure

        Returns
        -------
        best_wt: the best wt in wt_vec
        best_matching: the matching corresponds to best_wt
        """
        wt_vec = np.linspace(0, 1, n_wts)
        max_cancor = float('-inf')
        for ii in range(n_wts):
            if verbose:
                logging.info('Now at iteration {}, wt={}'.format(ii, wt_vec[ii]))
            if ii == 0:
                # curr_dist = self.dist['ovlp']
                curr_matching = self.matching['ovlp']
            elif ii == n_wts - 1:
                # curr_dist = self.dist['all']
                curr_matching = self.matching['all']
            else:
                # ii small --> more close to dist_ovlp
                curr_dist = (1 - wt_vec[ii]) * self.dist['ovlp'] + wt_vec[ii] * self.dist['all']
                if self.sparsity['ovlp'] is None or self.sparsity['all'] is None:
                    curr_sparsity = self.sparsity['ovlp'] or self.sparsity['all']
                else:
                    curr_sparsity = max(self.sparsity['ovlp'], self.sparsity['all'])
                try:
                    curr_matching = self.match_cells(curr_dist, curr_sparsity, 'auto')
                except ValueError:
                    # too sparse, find the suitable sparsity level
                    if verbose:
                        logging.info(
                            'Current sparsity config '
                            'is too sparse, finding a suitable sparsity level...'
                        )
                    _, curr_sparsity = self.search_minimum_sparsity(
                        curr_dist, slackness=200, init_sparsity=curr_sparsity + 1, verbose=verbose
                    )
                    curr_matching = self.match_cells(curr_dist, curr_sparsity, mode='auto')

            # compute median/mean of top_k cancors
            curr_cancor = np.mean(self.fit_cca(curr_matching, n_components=top_k)[0])
            if curr_cancor > max_cancor:
                max_cancor = curr_cancor
                self.matching['wted'] = curr_matching
                self.best_wt = wt_vec[ii]

        return self.best_wt, self.matching['wted']

    def _jr_clustering(self, matching, n_clusters, n_components=20, mismatch_prob=0.1, max_iter=50, tol=1e-5,
                       verbose=True):
        """
        Given matching:
        1) align df1 and df2, get df1_aligned and df2_aligned
        2) do CCA on df1_aligned and df2_aligned to get CCA scores, call them df1_cca and df2_cca
        3) do spectral clustering on the averaged df, i.e., (df1_cca+df2_cca)/2 to get init_labels
        4) given init_labels, run jr_kmeans on df1_aligned and df2_aligned

        Parameters
        ----------
        matching: matching used to align df1 and df2, must be a perfect matching
        n_clusters: number of clusters desired
        n_components: number of CCA components for initial clustering
        mismatch_prob: approximately mismatch_prob proportion of cells have different cluster labels
        max_iter: max iteration desired
        tol: if the objective function changes less than tol, terminates jr_kmeans
        verbose: if True, print details when running jr_kmeans

        Returns
        -------
        cluster_labels: cluster_labels[0] (resp. [1]) is the cluster labels for df1_aligned (resp. df2_aligned)
        centroids: centroids[0] (resp. [1]) is a (n_clusters, df1.shape[1])
                   (resp. df2.shape[1]) array representing the centroids for df1_aligned (resp. df2_aligned)
        """

        assert _is_perfect_matching(matching)

        df1_aligned, df2_aligned = self._align_modalities(matching)
        _, cca = get_cancor(df1_aligned, df2_aligned, n_components)
        df1_cca, df2_cca = cca.transform(df1_aligned, df2_aligned)
        # df1_cca = df1_cca / np.std(df1_cca, axis=0)
        # df2_cca = df2_cca / np.std(df2_cca, axis=0)
        df1_cca = normalize(df1_cca)
        df2_cca = normalize(df2_cca)
        init_cluster_labels = spectral_clustering(
            (df1_cca + df2_cca) / 2,
            n_clusters=n_clusters,
            need_svd=True
        )

        cluster_labels, _, centroids = jr_kmeans(
            [df1_aligned, df2_aligned],
            [init_cluster_labels, init_cluster_labels],
            n_clusters, mismatch_prob, max_iter, tol, verbose
        )

        return cluster_labels, centroids

    def filter_bad_matches(self, matching='wted', n_clusters=15,
                           n_components=20, bad_prop=0.1, max_iter=50,
                           tol=1e-5, verbose=True):
        """
        Given matching:
        1) run self.jr_clustering to get cluster_labels,
           note that cluster_labels[0] is the cluster labels for df1_aligned,
           and cluster_label[1] is the cluster labels for df2_aligned
        2) matching[ii] is bad if cluster_labels[0][ii] != cluster_labels[1][ii],
           then set matching[ii] to be an empty list

        Parameters
        ----------
        matching: matching to be filtered, must be a perfect matching
        n_clusters: how many clusters desired in self.jr_clustering
        n_components: number of CCA components in the initialization step in self.jr_clustering
        bad_prop: approximate proportion of bad matches to be filtered out
        max_iter: run self.jr_clustering for how many iterations
        tol: self.jr_clustering terminates when the change of objective function is <= tol
        verbose: print all the details when True

        Returns
        -------
        good_matching: for some ii, matching[ii] is an empty list, meaning that ii is bad and
                       has been filtered out
        """

        # proportion of good cells is approximately (1-mismatch_prob)^2 + mismatch_prob^2,
        # this is because the model says that a cluster label flips w.p. mismatch_prob
        # and we declare a pair to be good if either (1) both not flipped, or (2) both flipped
        # so given bad_prop, can solve for the corresponding mismatch_prob

        assert 1e-10 < bad_prop < 0.5 - 1e-10

        if n_components > min(self.p1, self.p2):
            logging.warning('n_components must be <= the dimensions of the two datasets, '
                            'set it to be equal to the minimum of the dimensions of the two datasets')
            n_components = min(self.p1, self.p2)

        mismatch_prob = 0.5 - np.sqrt(4 - 8 * bad_prop) / 4

        if isinstance(matching, str):
            if matching not in self.matching or self.matching[matching] is None:
                raise ValueError("Matching not found.")
            matching = self.matching[matching]

        cluster_labels, _ = self._jr_clustering(
            matching, n_clusters, n_components, mismatch_prob,
            max_iter, tol, verbose
        )

        good_matching = [x for x in matching]

        for ii in range(self.n1):
            if cluster_labels[0][ii] != cluster_labels[1][ii]:
                good_matching[ii] = []

        self.matching['final'] = good_matching

        return self.matching['final']


def eval_matching_accuracy(X_labels, Y_labels, matching, criterion='maj'):
    """
    'any': correct if X_labels[i] == any of Y_labels[i]
    'all': correct if X_labels[i] == all of Y_labels[i]
    'maj': correct if X_labels[i] == the mode in Y_labels[i], with
           ties broken in an arbitrary fashion
    """
    X_labels = np.array(X_labels)
    Y_labels = np.array(Y_labels)
    n = len(X_labels)
    assert n == len(matching)

    # if np.sum([len(perm[ii]) for ii in range(n)]) == n:
    #     # we know it's a 1-1 matching
    #     flags = [X_labels[ii] == Y_labels[perm][ii] for ii in range(n)]
    #     return np.mean(flags)

    # the matching is of variable length
    acc = 0
    cnt = 0
    for ii in range(n):
        if len(matching[ii]) == 0:
            continue
        elif len(matching[ii]) == 1:
            cnt = cnt + 1
            acc = acc + int(X_labels[ii] == Y_labels[matching[ii][0]])
        else:
            cnt = cnt + 1
            flag = [X_labels[ii] == y for y in Y_labels[matching[ii]]]
            if criterion == 'any':
                acc = acc + any(flag)
            elif criterion == 'all':
                acc = acc + all(flag)
            elif criterion == 'maj':
                acc = acc + (np.sum(flag) > (len(flag) // 2))
            else:
                raise NotImplementedError

    return acc / cnt
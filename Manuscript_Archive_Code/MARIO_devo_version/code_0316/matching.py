import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.sparse.csgraph import maximum_bipartite_matching
from sklearn.utils.extmath import randomized_svd
from utils import cdist_correlation, normalize, get_cancor
from clustering import spectral_clustering, jr_kmeans


def _compute_dist_ovlp(arr1, arr2, n_components, min_dist):
    """
    Do stacked SVD for arr1 and arr2,
    then compute their cross distance

    Parameters
    ----------
    arr1: a numpy array
    arr2: another numpy array
    min_dist: in the final distance matrix, the minimum value must be >= min_dist

    Returns
    -------
    dist: distance matrix

    """
    n1, p1 = arr1.shape
    n2, p2 = arr2.shape
    assert p1 == p2
    assert n_components <= p1

    U, s, Vh = randomized_svd(np.concatenate((arr1, arr2), axis=0),
                              n_components=n_components)

    if n_components == p1:
        dist_mat = cdist_correlation(arr1, arr2)
    else:
        arr1_svd = U[:n1, :] @ np.diag(s) @ Vh
        arr2_svd = U[n1:, :] @ np.diag(s) @ Vh
        dist_mat = cdist_correlation(arr1_svd, arr2_svd)

    # make sure min_dist is at least self.min_dist
    dist_mat = _check_min_dist(dist_mat, min_dist)

    return dist_mat, s


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
                             min_dist):
    """
    Use binary search to search for the minimum sparsity level k such that
        if dist_mat is trimmed to be a k-NN graph,
        a valid matching still exists.
    The search starts with k_left=1 and k_right=dist_mat.shape[0],
        and it is always true that:
        1) any k < k_left doesn't give a valid matching
        2) any k >= k_right gives a valid matching

    Args:
        dist_mat: a numpy array of distance matrix
        slackness: binary search terminates when k_right - k_left <= slackness,
                   an exact binary search corresponds to slackness = 0
        init_sparsity: binary search starts from k=init_sparsity
        m_min: each cell in df1 is matched to at least m_min many cells in df2
        m_max: each cell in df1 is matched to at most m_max many cells in df2
        num_cells_to_use: total number of cells to use in df2
        min_dist: it must be true that minimum entry in dist_mat is >= min_dist

    TODO:
        Currently I'm using print(), use fancier logging tools
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

        print('If sparsity>={}, then there is a valid matching; '
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
    print('If sparsity>={}, then there is a valid matching; '
          'if sparsity<{}, then there is no valid matching.'.format(k_right, k_left))

    # return argsort_res


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
    if sparsity is not None:
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


class CellMatching(object):
    def __init__(self, df1, df2, normalization=True):
        self.min_dist = 1e-5

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

        # parameters related to matching
        self.m_min = None
        self.m_max = None
        self.num_cells_to_use = None
        self.num_sinks = None

        # self.df1_cca = None
        # self.df2_cca = None

    def specify_matching_params(self, m_min, m_max, num_cells_to_use):
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

    def compute_dist_ovlp(self, n_components=10):
        """
        Compute distance matrix based on overlapping features.
        Args:
            n_components: number of SVD components to keep
        Returns:
            dist_ovlp: distance matrix
            s: vector of singular values
        """
        assert n_components <= len(self.ovlp_features)
        df1_ovlp = self.df1[self.ovlp_features]
        df2_ovlp = self.df2[self.ovlp_features]
        df1_ovlp = normalize(df1_ovlp).to_numpy()
        df2_ovlp = normalize(df2_ovlp).to_numpy()

        return _compute_dist_ovlp(df1_ovlp, df2_ovlp, n_components, self.min_dist)

    def search_minimum_sparsity(self, dist_mat, slackness=200, init_sparsity=None):
        """
        Wrapper for _search_minum_sparsity()
        """
        _search_minimum_sparsity(dist_mat, slackness, init_sparsity, self.m_min,
                                 self.m_max, self.num_cells_to_use, self.min_dist)

    def match_cells(self, dist_mat, sparsity=None, mode='auto'):
        """
        Wrapper for _match_cells()
        """
        return _match_cells(dist_mat, sparsity, self.m_min, self.m_max, self.num_cells_to_use,
                            self.min_dist, mode)

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

    def fit_cca(self, matching, n_components=20, max_iter=1000):
        """
        Align df1 and df2 using matching, then fit a CCA.
        Args:
            matching: a list of (potentially variable length) lists;
                      cell ii in df1 is matched to cell matching[ii] in df2
            n_components: number of components for CCA
            max_iter: maximum iteration for CCA
        Returns:
            cancor: canonical correlations
            cca: CCA object
        """

        X, Y = self._align_modalities(matching)
        cancor, cca = get_cancor(X, Y, n_components, max_iter)

        return cancor, cca

    # def align_features(self, matching, n_components=20, max_iter=1000):
    #     """
    #     Use matching to fit a CCA, then use CCA to embed df1 and df2 to a common subspace,
    #
    #     Parameters
    #     ----------
    #     matching: a list of (potentially variable length) lists;
    #               cell ii in df1 is matched to cell matching[ii] in df2
    #     n_components: number of CCA components
    #     max_iter: number of max CCA iterations
    #
    #     Returns
    #     -------
    #     df1_cca, df2_cca: arrays of CCA scores
    #     """
    #     _, cca = self.fit_cca(matching, n_components, max_iter)
    #     self.df1_cca, self.df2_cca = cca.transform(self.df1, self.df2)
    #     self.df1_cca = pd.DataFrame(data=self.df1_cca, index=np.arange(self.n1),
    #                                 columns=['CCA' + str(ii) for ii in range(n_components)])
    #     self.df2_cca = pd.DataFrame(data=self.df2_cca, index=np.arange(self.n2),
    #                                 columns=['CCA' + str(ii) for ii in range(n_components)])

    def compute_dist_all(self, matching, n_components=20, max_iter=1000):
        """
        Given matching, align df1 and df2, fit a CCA,
        then use CCA scores to get the distance matrix.

        Parameters
        ----------
        matching: a list of (potentially variable length) lists;
                  cell ii in df1 is matched to cell matching[ii] in df2
        n_components: number of CCA components
        max_iter: max number of iterations when fitting CCA

        Returns
        -------
        dist_mat: distance matrix
        """
        _, cca = self.fit_cca(matching, n_components, max_iter)
        df1_cca, df2_cca = cca.transform(self.df1, self.df2)
        dist_mat = cdist_correlation(df1_cca, df2_cca)
        dist_mat = _check_min_dist(dist_mat, self.min_dist)
        return dist_mat

        # assert self.df1_cca is not None
        # assert self.df2_cca is not None
        # return _compute_dist_ovlp(self.df1_cca.to_numpy(), self.df2_cca.to_numpy(), n_components, self.min_dist)

    def _jr_clustering(self, matching, n_clusters, n_components=20, mismatch_prob=0.1, max_iter=50, tol=1e-5,
                       verbose=True):
        """
        Given matching:
        1) align df1 and df2, get df1_aligned and df2_aligned
        2) do CCA on df1_aligned and df2_aligned to get CCA scores, call them df1_cca and df2_cca
        3) do spectral clustering on the averaged df, i.e., (df1_cca+df2_cca) to get init_labels
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

    def filter_bad_matches(self, matching, n_clusters=15,
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

        mismatch_prob = 0.5 - np.sqrt(4 - 8 * bad_prop) / 4
        cluster_labels, _ = self._jr_clustering(
            matching, n_clusters, n_components, mismatch_prob,
            max_iter, tol, verbose
        )

        good_matching = [x for x in matching]

        for ii in range(self.n1):
            if cluster_labels[0][ii] != cluster_labels[1][ii]:
                good_matching[ii] = []

        return good_matching

    # def _filter_bad_matches_debug(self, matching, n_clusters=15,
    #                        n_components=20, bad_prop=0.1, max_iter=50,
    #                        tol=1e-5, verbose=True):
    #     mismatch_prob = 0.5 - np.sqrt(4-8*bad_prop)/4
    #
    #     df1_aligned, df2_aligned = self._align_modalities(matching)
    #     _, cca = get_cancor(df1_aligned, df2_aligned, n_components)
    #     df1_cca, df2_cca = cca.transform(df1_aligned, df2_aligned)
    #
    #     # setting df1_cca and df2_cca will let the data be "too uniform",
    #     # this can hurt the clustering performance
    #     init_cluster_labels = spectral_clustering(
    #         (df1_cca + df2_cca) / 2,
    #         n_clusters=n_clusters,
    #         need_svd=True
    #     )
    # cluster_labels, _, centroids = jr_kmeans(
    #     [df1_cca, df2_cca],
    #     [init_cluster_labels, init_cluster_labels],
    #     n_clusters, mismatch_prob, max_iter, tol, verbose
    # )
    #
    # good_matching = [x for x in matching]
    #
    # for ii in range(self.n1):
    #     if cluster_labels[0][ii] != cluster_labels[1][ii]:
    #         good_matching[ii] = []
    #
    # return good_matching


def eval_matching_accuracy(X_labels, Y_labels, matching, criterion='any'):
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


def some_function():
    # TODO: cut df1 in multiple pieces then to matching separately
    # this will be more accurate in terms of cluster-level accuracy
    pass

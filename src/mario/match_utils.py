import numpy as np
from collections import Counter
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import min_weight_full_bipartite_matching


def search_minimum_sparsity(dist_mat, slackness, init_sparsity,
                            m_min, m_max, num_cells_to_use,
                            min_dist, verbose=True):
    """
    Use binary search to search for the minimum sparsity level k such that
    if dist_mat is trimmed to be a k-NN graph, a valid matching still exists.
    The search starts with k_left=1 and k_right=dist_mat.shape[0], and it is always true that:
    1) any k < k_left doesn't give a valid matching;
    2) any k >= k_right gives a valid matching.

    Parameters
    ----------
    dist_mat : array-like of shape (n_samples_1, n_samples_2)
        A two-dimensional distance matrix.
    slackness : int
        Binary search terminates when k_right - k_left <= slackness;
        an exact binary search corresponds to slackness = 0
    init_sparsity : int
        Binary search starts from k=init_sparsity.
    m_min : int
        Each row in the first dataset is matched to at least m_min many rows in the second dataset.
    m_max : int
        Each row in the first dataset is matched to at most m_max many rows in the second dataset.
    num_cells_to_use : int
        Total number of rows to use in the second dataset.
    min_dist : float
        It must be true that minimum entry in dist_mat is >= min_dist.

    Returns
    -------
    k_left : int
        If sparsity<k_left, then there is no valid matching.
    k_right : int
        If sparsity>=k_right, then there is a valid matching.
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
            print(
                'If sparsity>={}, then there is a valid matching; '
                'if sparsity<{}, then there is no valid matching.'.format(k_right, k_left), flush=True
            )
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
            # can trim more aggressively --> try a smaller k
            k_right = k
        else:
            # the current k doesn't give a valid matching
            # try a larger k
            k_left = k + 1

        cnt = cnt + 1

    # we know that k_right must give a valid matching
    if verbose:
        print(
            'If sparsity>={}, then there is a valid matching; '
            'if sparsity<{}, then there is no valid matching.'.format(k_right, k_left), flush=True
        )
    return k_left, k_right


def get_matching_from_indices(col_idx, n1, n2, m_max):
    """
    Assume col_idx is obtained from min_weight_full_bipartite_matching or linear_sum_assignment
    with some dist_mat as input.
    And this dist_mat is organized such that the sinks are in dist_mat[:, :n2].
    This function calculates the matching from col_idx.

    Parameters
    ----------
    col_idx : list
        output from min_weight_full_bipartite_matching or linear_sum_assignment.
    n1 : int
        Sample size of the first dataset.
    n2 : int
        Sample size of the second dataset.
    m_max : int
        Each row in the first dataset is matched to at most m_max many rows in the second dataset.

    Returns
    -------
    list
        A list of (potentially variable length) lists;
        it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.
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


def match_cells(dist_mat, sparsity, m_min, m_max, num_cells_to_use,
                min_dist, mode='auto'):
    """
    Given dist_mat, first trim it to a k-NN graph according to the desired sparsity
    sparsity level, then do a matching according to the specified parameters.

    Parameters
    ----------
    dist_mat : array-like of shape (n_samples_1, n_samples_2)
        A two-dimensional distance matrix.
    sparsity : int
        An integer k such that dist_mat will be trimmed into a k-NN graph.
    m_min : int
        Each row in the first dataset is matched to at least m_min many rows in the second dataset.
    m_max : int
        Each row in the first dataset is matched to at most m_max many rows in the second dataset.
    num_cells_to_use : int
        Total number of samples to use in the second dataset.
    min_dist : float
        It must be true that minimum entry in dist_mat is >= min_dist.
    mode : str, default='auto'
        If 'sparse', use min_weight_full_bipartite_matching;
        if 'dense', use linear_sum_assignment;
        if 'auto': when sparsity<=n//2, use 'sparse', else use 'dense'.

    Returns
    -------
    list
        A list of (potentially variable length) lists;
        it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.
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

    return get_matching_from_indices(col_idx, n1, n2, m_max)


def is_perfect_matching(matching):
    """Check if matching is perfect (i.e., all nodes are matched) or not.
    Parameters
    ----------
    matching : list
        If matching[ii] is an empty list, then ii is not matched to anyone.

    Returns
    -------
    bool
        True if matching is perfect.
    """
    matching_len = np.array([len(matching[ii]) for ii in range(len(matching))])
    matched_flag = (matching_len > 0)
    return all(matched_flag)


def calculate_duplications(matching):
    """
    Calculate the number of i such that matching[i] appears more than one times.
    Also calculate the proportion of duplicated i's among all matched samples.

    Parameters
    ----------
    matching : list
        A list of (potentially variable length) lists;
        it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.

    Returns
    -------
    dup_prop : float
         Proportion of i (among all matched samples) such that matching[i] appears more than one times.
    dup_freq : int
        Number of i such that matching[i] appears more than one times.
    """
    cells = []
    for matched in matching:
        cells.extend(matched)
    counter = Counter(cells)
    duplicated_cells = [cell for cell, freq in counter.items() if freq > 1]
    dup_prop = len(duplicated_cells) / len(cells)
    dup_freq = np.mean([counter[cell] for cell in duplicated_cells]) if len(duplicated_cells) > 0 else 0
    return dup_prop, dup_freq


def batching(df, num_batches):
    """
    Given a pandas dataframe, shuffle it and then cut it into approximately equal-sized batches.

    Parameters
    ----------
    df : pandas.DataFrame
        A pandas dataframe.
    num_batches : int
        Number of batches.

    Returns
    -------
    df_lst : list
        A list of pandas dataframe; each entry corresponds to one batch.
    perm_lst : list
        A list of nested lists of indices; perm_lst[i] corresponds to the original row indices of df_lst[i].
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


def stitch(perm_lst, matching_lst, n):
    """Stitch multiple matching into one.

    Parameters
    ----------
    perm_lst : list
        Output from batching.
    matching_lst : list
        perm_lst[i][j] should be matched to matching_lst[i][j].
    n : int
        Sample size.

    Returns
    -------
    list
        A list of (potentially variable length) lists;
        it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.
    """
    stitched_dict = {} # store perm_lst[i][j]: matching_lst[i][j]
    for i in range(len(perm_lst)):
        for j in range(len(perm_lst[i])):
            stitched_dict[perm_lst[i][j]] = matching_lst[i][j]
    return [stitched_dict[i] for i in range(n)]


def eval_matching_accuracy(X_labels, Y_labels, matching, criterion='maj'):
    """Evaluate matching accuracy.

    Parameters
    ----------
    X_labels : array-like of shape (n_samples_of_X)
        The cluster labels of the first dataset.
    Y_labels : array-like of shape (n_samples_of_Y)
        The cluster labels of the second dataset.
    matching : list
        A list of (potentially variable length) lists;
        it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.
    criterion : str, default='maj'
        If 'any': correct if X_labels[i] == any of Y_labels[i];
        if 'all': correct if X_labels[i] == all of Y_labels[i];
        if 'maj': correct if X_labels[i] == the mode in Y_labels[i], with ties broken in an arbitrary fashion.

    Returns
    -------
    float
        The proportion of correctly matched cells in X.
    """
    X_labels = np.array(X_labels)
    Y_labels = np.array(Y_labels)
    n = len(X_labels)
    assert n == len(matching)

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
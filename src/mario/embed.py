import numpy as np
from sklearn.cross_decomposition import CCA
from . import utils


def get_cancor(X, Y, n_components=10, max_iter=2000):
    """Fit CCA and calculate the canonical correlations.

    Parameters
    ----------
    X: array-like of shape (n_samples, n_features_of_X)
        First dataset.
    Y: array-like of shape (n_samples, n_features_of_Y)
        Second dataset.
    n_components: int, default=2
        Number of CCA components to calculate; must be <= min(n_features_of_X, n_features_of_Y)
    max_iter: int, default=1000
        Maximum number of iterations.

    Returns
    -------
    cancor: array-like of shape (n_components, )
        Vector of canonical components.
    cca: CCA
        CCA object.
    """
    assert X.shape[1] >= n_components
    assert Y.shape[1] >= n_components

    cca = CCA(n_components=n_components, max_iter=max_iter)
    cca.fit(X, Y)
    X_c, Y_c = cca.transform(X, Y)
    cancor = np.corrcoef(
        X_c, Y_c, rowvar=False).diagonal(
        offset=cca.n_components)
    return cancor, cca


def gcca_init(data_list, n_components=10, normalization=True):
    """Initialization of GCCA.
    Parameters
    ----------
    data_list: list
        A list of data matrices, e.g., [X1, X2, X3, X4]
    n_components : int, default=10
        Number of GCCA components.
    normalization : bool, default=True
        If true, center each column and scale each column to have unit standard deviation.

    Returns
    ----------
    score_list: list
        A list of projected data.
    """

    # normalize
    if normalization:
        for data in data_list:
            data = utils.normalize(data)
    # run CCA on the first two data matrices
    _, cca = get_cancor(data_list[0], data_list[1], n_components=n_components)
    score_list = list(cca.fit_transform(data_list[0], data_list[1]))
    avg_score = (score_list[0] + score_list[1]) / 2
    # run OLS regressions for the rest of data matrices
    for data in data_list[2:]:
        B = np.linalg.lstsq(data, avg_score, rcond=None)[0]
        score_list.append(data @ B)

    return score_list


def gcca_refine(data_list, init_scores, max_iter=500, tol=1e-3, verbose=True):
    """Refinement of GCCA.

    Parameters
    ---------
    data_list : list
        A list of data matrices.
    init_scores : list
        A list of initial scores.
    max_iter : int, default=500
        Number of maximum iterations.
    tol : float, default=1e-3
        Algorithm terminates when the change of objective function value is <=tol.
    verbose : bool, default=True
        Print details if True.

    Returns
    ----------
    all_scores : list
        A list of projected data.
    """

    m = len(init_scores)
    if len(init_scores[0]) == 0:
        n_components = 1
        n = len(init_scores[0])
    else:
        n, n_components = init_scores[0].shape

    all_scores = np.array([np.empty_like(init_scores[0])] * m)
    for ii in range(n_components):
        if verbose:
            print("Computing the {}-th canonical score...".format(ii), flush=True)
        # regress out the influence of all_score[:, :, :ii]
        if ii > 0:
            curr_data_list = []
            for idx, data in enumerate(data_list):
                coef = np.linalg.lstsq(all_scores[idx, :, :ii], data, rcond=None)[0]
                curr_data_list.append(data - all_scores[idx, :, :ii] @ coef)
        else:
            curr_data_list = data_list

        # preparations for entering main iteration
        if n_components == 1:
            curr_scores = [score for score in init_scores]
        else:
            curr_scores = [score[:, ii] for score in init_scores]
        score_sum = sum(score for score in curr_scores)
        prev_obj = np.inf

        # alternating minimization
        for iter_idx in range(max_iter):
            for jj in range(m):
                y = (score_sum - curr_scores[jj]) / (m - 1)
                X = curr_data_list[jj]
                coef = np.linalg.lstsq(X, y, rcond=None)[0]
                coef = coef / np.sqrt(np.sum(coef ** 2))
                score = X @ coef
                score_sum = score_sum - curr_scores[jj] + score
                curr_scores[jj] = score
            # check convergence
            obj = 0
            for jj in range(m - 1):
                obj += np.sqrt(np.sum((curr_scores[jj] - curr_scores[jj + 1]) ** 2))

            if verbose and iter_idx % 50 == 0:
                print("At iteration {}, the objective value is {}.".format(iter_idx, obj), flush=True)

            if abs(obj - prev_obj) < tol:
                break
            else:
                prev_obj = obj

        if verbose:
            print("Finished computing the {}-th canonical score, the objective is {}.".format(ii, obj), flush=True)

        for jj in range(m):
            if n_components == 1:
                all_scores[jj, :] = curr_scores[jj]
            else:
                all_scores[jj, :, ii] = curr_scores[jj]

        # curr_scores are linear combinations of the residual matrices
        # convert them into linear combinations of the original data matrices
    #         for jj in range(m):
    #             coef = np.linalg.lstsq(data_list[jj], curr_scores[jj], rcond=None)[0]
    #             if n_components == 1:
    #                 all_scores[jj, :] = data_list[jj] @ coef
    #             else:
    #                 all_scores[jj, :, ii] = data_list[jj] @ coef

    return [all_scores[jj, :, :] for jj in range(m)]


def gcca(data_list, n_components=10, normalization=True, max_iter=500, tol=1e-3, verbose=True):
    """Embedding multiple datasets via generalized canonical correlation analysis.
    Parameters
    ---------
    data_list : list
        A list of data matrices.
    n_components : int, default=10
        Number of GCCA components.
    normalization : bool, default=True
        If true, center each column and scale each column to have unit standard deviation.
    max_iter : int, default=500
        Number of maximum iterations.
    tol : float, default=1e-3
        Terminate the algorithm when the change in objective value is <=tol.
    verbose : bool, default=True
        Print details if True.

    Returns
    ----------
    list
        A list of projected data.
    """
    init_scores = gcca_init(data_list, n_components, normalization)
    return gcca_refine(data_list, init_scores, max_iter, tol, verbose)
import numpy as np
import itertools
from sklearn.utils.extmath import randomized_svd
from sklearn.cluster import KMeans


def align_clusters(z1, z2, K, comb_search=False):
    """Align two cluster labels so that the Hamming distance is minimized.

    Parameters
    ----------
    z1 : array-like of shape (n_samples, )
        First cluster label vector.
    z2 : array-like of shape (n_samples, )
        Second cluster label vector.
    K : int
        Number of clusters, cluster labels are encodes from 0, 1, ..., K-1
    comb_search : bool, default=False
        If False, the function is implemented according
        to Lemma 19 of Gao, C., Ma, Z., Zhang, A. Y., & Zhou, H. H. (2017).
        "Achieving optimal misclassification proportion in
        stochastic block models." JMLR.
        In this case, the function will succeed only if
        there exists a constant C>=1 such that
        1. the minimum cluster size of z1 and z2 are both >= n/(C*K);
        2. after alignment, the Hamming distance between z1 and z2
            is < 1/(C*K).
        Otherwise there is no guarantee that such an alignment can be done
        If True, then this function tries out
        all possible permutations of size K and pick the best one.

    Returns
    -------
    z1 : array-like of shape (n_samples, )
        First cluster label vector.
    z2_aligned : array-like of shape (n_samples, )
        Second cluster label vector after alignment.
    perm: list
        A permutation so that z2_aligned=[perm[k] for k in z2].
    """
    n = len(z1)
    assert n == len(z2)

    z1 = np.array(z1, dtype=int)
    z2 = np.array(z2, dtype=int)

    if comb_search:
        best_perm = None
        best_loss = 2
        for perm in itertools.permutations(np.arange(K)):
            z2_permuted = [perm[k] for k in z2]
            loss = np.mean(z1 != z2_permuted)
            if loss < best_loss:
                best_perm = perm
                best_loss = loss
        if best_perm is not None:
            z2_aligned = [best_perm[k] for k in z2]
            return z1, z2_aligned, best_perm
        else:
            raise(ValueError("Alignment failed!"))

    # don't do combinatorial search

    nbhd_z1 = [np.nonzero(z1 == k) for k in range(K)]
    nbhd_z2 = [np.nonzero(z2 == k) for k in range(K)]

    perm = []

    for i in range(K):
        intersec_sizes = [len(np.intersect1d(nbhd_z1[j], nbhd_z2[i]))
                          for j in range(K)]
        perm.append(np.argmax(intersec_sizes))

    z2_aligned = np.array([perm[k] for k in z2], dtype=int)

    return z1, z2_aligned, perm


def spectral_clustering(X, n_clusters=10, need_svd=True):
    """Do spectral clustering on X.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Dataset to do clustering.
    n_clusters: int, default=10
        Number of clusters.
    need_svd: bool, default=True
        If False, then directly do a k-means on rows of X.
        If True, then this function
        (1) computes the truncated SVD of X = U @ diag(s) @ Vh;
        (2) do a k-means clustering on rows of U @ diag(s).

    Returns
    -------
    cluster_label : array-like of shape (n_samples, )
        Cluster label vector.
    """
    if need_svd:
        U, s, _ = randomized_svd(X, n_components=n_clusters)
        X = U @ np.diag(s)

    kmeans = KMeans(n_clusters=n_clusters).fit(X)
    return kmeans.labels_


def compute_centroids(X, cluster_labels, n_clusters):
    """Given cluster labels, compute the centroids of each cluster.

    Parameters
    ----------
    X : array-like of shape(n_samples, n_features)
        Dataset for computation of centroids.
    cluster_labels : array-like of shape (n_samples, )
        Cluster label vector, each entry takes value in {0, ..., n_clusters-1}
    n_clusters : int
        Number of clusters.

    Returns
    -------
    array-like of shape (n_clusters, n_features)
        Each row is a centroid vector.
    """
    cluster_idx = [np.nonzero(cluster_labels == k) for k in range(n_clusters)]
    centroids = [np.mean(X[cluster_idx[k], :], axis=1).squeeze() for k in range(n_clusters)]
    return np.array(centroids)


def jr_kmeans_one_step(X, init_cluster_labels, n_clusters, mismatch_prob=0.1):
    """Do one step of jointly regularized k-means update.

    Given initial cluster labels, say z1, ..., zL, this function computes the initial centroids
    and solves for the next z1, ..., zL, zstar such that the objective function
    sum_{l} sum_{i} 0.5 * || X[l, i, :] - centroids[l, zl[i], :] ||^2
    + log[(1-mismatch_prob)/(mismatch_prob/(n_clusters-1))] * indicator(zl[i] != zstar[i])
    is minimized.

    Parameters
    ----------
    X : list
        A length L list of numpy arrays, where L is the number of modalities and the l-th element in this list
        is of shape (n_samples, n_features_of_modality_l).
    init_cluster_labels : list
        A length L list of integer arrays, encoding the preliminary cluster labels for each modality.
    n_clusters : int
        Number of clusters desired.
    mismatch_prob : float, default=0.1
        A number strictly in (0, 0.5) that determines the strength of regularization.

    Returns
    -------
    ind_cluster_labels : array-like of shape (L, n_samples)
        Each row is an integer vector encoding the clustering label for one modality.
    glob_cluster_labels : array-like of shape (n_samples, )
        An integer vector encoding the clustering label for the global structure.
    final_loss: float
        Value of the loss function.
    """
    if mismatch_prob < 1e-8:
        raise NotImplementedError("mismatch_prob=0, not implemented")

    L = len(X)
    n = X[0].shape[0]

    centroids = [compute_centroids(X[layer],
                                   init_cluster_labels[layer],
                                   n_clusters)
                 for layer in range(L)]

    ind_cluster_labels = np.zeros((L, n), dtype=np.int64)
    glob_cluster_labels = np.zeros(n, dtype=np.int64)

    reg_param = np.log((mismatch_prob/(n_clusters-1))/(1-mismatch_prob))

    final_loss = 0
    for i in range(n):
        best_loss = np.Inf
        best_Z = None
        best_zstar = None

        # suppose zstar_i = k
        for k in range(n_clusters):

            loss_k = np.zeros((n_clusters, L))
            # loss_k is organized as follows
            # -------------------------------
            # l         | 1  2  3           L
            # zl[i]=0   |
            # zl[i]=1   |
            # ...
            # zl[i]=    |

            for layer in range(L):
                loss_k[:, layer] = [0.5 * np.sum((X[layer][i, :]
                                                  - centroids[layer][k, :])**2)
                                    for k in range(n_clusters)]
                # less loss when zl[i] = k
                loss_k[k, layer] = loss_k[k, layer] + reg_param

            # Z is a length L list encoding
            # the estimator for zl[i] when zstar[i] = k
            Z = np.argmin(loss_k, axis=0)

            # calculate the loss w.r.t. Z
            loss_k = [loss_k[Z[layer], layer] for layer in range(L)]
            loss_k = np.sum(loss_k)

            if loss_k < best_loss:
                best_loss = loss_k
                best_Z = Z
                best_zstar = k

        assert best_Z is not None
        assert best_zstar is not None

        ind_cluster_labels[:, i] = np.array(best_Z)
        glob_cluster_labels[i] = best_zstar
        final_loss = final_loss + best_loss / n

    return ind_cluster_labels, glob_cluster_labels, final_loss


def jr_kmeans(X, init_cluster_labels, n_clusters,
              mismatch_prob=0.1, max_iter=20, tol=1e-5, verbose=False):
    """Training jointly regularized k-means clustering using alternating minimization.

    The objective function is
    sum_{l} sum_{i} 0.5* || X[l, i, :] - centroids[l, zl[i], :] ||^2
    + log[(1-mismatch_prob/(n_clusters-1))/mismatch_prob] * indicator(zl[i] != zstar[i]).

    Parameters
    ----------
    X : list
        A length L list of numpy arrays, where L is the number of modalities and the l-th element in this list
        is of shape (n_samples, n_features_of_modality_l).
    init_cluster_labels : list
        A length L list of integer arrays, encoding the preliminary cluster labels for each modality.
    n_clusters : int
        Number of clusters desired.
    mismatch_prob : float, default=0.1
        A number strictly in (0, 0.5) that determines the strength of regularization.
    max_iter : float, default=20
        Maximum number of iterations.
    tol : float, default=1e-5
        Convergence criterion.
    verbose : bool, default=True
        If True, print the progress of each iteration.

    Returns
    -------
    ind_cluster_labels : array-like of shape (L, n_samples)
        Each row is an integer vector encoding the clustering label for one modality.
    glob_cluster_labels : array-like of shape (n_samples, )
        An integer vector encoding the clustering label for the global structure.
    centroids : list
        A length L list, each is an array of shape (n_clusters, n_features_of_modality_l)
        representing the centroids for the l-th modality.
    """
    if mismatch_prob < 1e-8:
        raise NotImplementedError("mismatch_prob=0, not implemented")

    L = len(X)
    assert L == len(init_cluster_labels)
    n_vec = [X[layer].shape[0] for layer in range(L)]
    n = n_vec[0]
    assert np.sum(np.array([(n_each == n) for n_each in n_vec])) == len(n_vec)

    # prev_glob_cluster_labels = init_cluster_labels[0]
    prev_loss = np.Inf
    for ii in range(max_iter):
        ind_cluster_labels, glob_cluster_labels, loss = \
            jr_kmeans_one_step(X, init_cluster_labels,
                               n_clusters, mismatch_prob)

        if verbose:
            print("Now at iteration {}, "
                  "current loss is {}.".format(ii, loss), flush=True)

        if np.abs(prev_loss - loss) < tol:
            centroids = [compute_centroids(
                X[layer], init_cluster_labels[layer], n_clusters
            ) for layer in range(L)]

            return ind_cluster_labels, glob_cluster_labels, centroids

        prev_loss = loss
        init_cluster_labels = ind_cluster_labels

    centroids = [compute_centroids(
        X[layer], init_cluster_labels[layer], n_clusters
    ) for layer in range(L)]

    return ind_cluster_labels, glob_cluster_labels, centroids
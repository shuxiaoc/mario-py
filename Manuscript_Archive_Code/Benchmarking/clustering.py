import numpy as np
import logging
import itertools
from sklearn.utils.extmath import randomized_svd
from sklearn.cluster import KMeans
# from numba import jit


def align_clusters(z1, z2, K, comb_search=False):
    """
    Given two cluster label vectors of size n,
    each assumed to only take values in {0, 1, ..., K-1},
    this function compute a permutation of size K, called perm,
    such that z1 is aligned with z2' = [perm[k] for k in z2],
    in the sense that the Hamming distance between z1 and z2'
    is minimized over all possible permutations of size K.

    The default is comb_search=False, and is implemented according
    to Lemma 19 of Gao, C., Ma, Z., Zhang, A. Y., & Zhou, H. H. (2017).
    "Achieving optimal misclassification proportion in
    stochastic block models." JMLR.

    In this case, the function will succeed only if
    there exists a constant C>=1 such that
        1. the minimum cluster size of z1 and z2 are both >= n/(C*K);
        2. after alignment, the Hamming distance between z1 and z2
            is < 1/(C*K).
    Otherwise there is no guarantee that such an alignment can be done.

    When comb_search=True, then this function tries out
    all possible permutations of size K and pick the best one.

    Returns:
        z1: original z1
        z2: aligned version of z2
        perm: the permutation that aligns z1 and z2
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

    z2_aligned = [perm[k] for k in z2]

    return z1, z2_aligned, perm


def spectral_clustering(X, n_clusters=10, need_svd=True):
    """
    If need_svd=True, then this function
        (1) computes the truncated SVD of X = U @ diag(s) @ Vh;
        (2) do a k-means clustering on rows of U @ diag(s).
    This is the spectral clustering algorithm analyzed in
        Löffler, M., Zhang, A. Y., & Zhou, H. H. (2019).
        "Optimality of Spectral Clustering in the Gaussian Mixture Model".
        arXiv preprint arXiv:1911.00538.

    If need_svd=False, then directly do a k-means on rows of X

    Returns:
        cluster_label: an integer array of shape (n_clusters, )
    """

    if need_svd:
        U, s, _ = randomized_svd(X, n_components=n_clusters)
        X = U @ np.diag(s)

    kmeans = KMeans(n_clusters=n_clusters).fit(X)
    return kmeans.labels_


def compute_centroids(X, cluster_labels, n_clusters):
    """
    Compute the centroids given the cluster_labels.

    Inputs:
        X: (n, p) data matrix
        cluster_labels: integer array taking values in {0, ..., n_clusters-1}
            encoding cluster labels

    Returns:
        centroids: (n_clusters, p) array, each row is a centroid
    """
    cluster_idx = [np.nonzero(cluster_labels == k) for k in range(n_clusters)]
    centroids = [np.mean(X[cluster_idx[k], :], axis=1).squeeze() for k in range(n_clusters)]
    return np.array(centroids)


def jr_kmeans_one_step(X, init_cluster_labels, n_clusters, mismatch_prob=0.1):
    """
    Inputs:
        X: a length L list of numpy arrays, where L is
            the number of modalities and the l-th element in this list
            is of shape (n, pl).
        init_cluster_labels: a length L list of integer arrays, encoding
            the preliminary cluster labels for each modality.
        n_clusters: number of clusters required
        mismatch_prob: a number strictly in (0, 0.5) that determines
            the strength of regularization
            TO-DO: implement the mismatch_prob=0 case

    Given initial cluster labels, say z1, ..., zL, this function
    computes the initial centroids
    (a length L list, each is a n_clusters*pl array)
    and solves for the next z1, ..., zL, zstar such that the objective function
        sum_{l} sum_{i} 0.5* || X[l, i, :] - centroids[l, zl[i], :] ||^2
          + log[(1-mismatch_prob)/mismatch_prob] * indicator(zl[i] != zstar[i])
    is minimized.

    Returns:
        ind_cluster_labels: a lenght L list, each is an integer array
            encoding the clustering label for a modality
        glob_cluster_labels: a length n integer array encoding
            the clustering label for the global structure
        final_loss: value of the loss function
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

    reg_param = np.log(mismatch_prob/(1-mismatch_prob))

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
        # for layer in range(L):
        #     ind_cluster_labels[layer, i] = best_Z[layer]
        glob_cluster_labels[i] = best_zstar
        final_loss = final_loss + best_loss / n

    return ind_cluster_labels, glob_cluster_labels, final_loss


def jr_kmeans(X, init_cluster_labels, n_clusters,
              mismatch_prob=0.1, max_iter=20, tol=1e-5, verbose=False):
    """
    Use alternating minimization to minimize the jointly regularized
    kmeans objective:
        sum_{l} sum_{i} 0.5* || X[l, i, :] - centroids[l, zl[i], :] ||^2
          + log[(1-mismatch_prob)/mismatch_prob] * indicator(zl[i] != zstar[i])

    Inputs:
        X: a length L list of numpy arrays, where L is
            the number of modalities and the l-th element in this list
            is of shape (n, pl).
        init_cluster_labels: a length L list of integer arrays, encoding
            the preliminary cluster labels for each modality.
        n_clusters: number of clusters required
        mismatch_prob: a number strictly in (0, 0.5) that determines
            the strength of regularization
        maxiter: maximum number of iterations

    Returns:
        ind_cluster_labels: a length L list, each is an integer array
            encoding the clustering label for a modality
        glob_cluster_labels: a length n integer array encoding
            the clustering label for the global structure
        centroids: a length L list, each is a (n_clusters, pl) array
            representing the centroids for the l-th modality
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
    # if verbose:
        # print("Start training...")
    for ii in range(max_iter):
        ind_cluster_labels, glob_cluster_labels, loss = \
            jr_kmeans_one_step(X, init_cluster_labels,
                               n_clusters, mismatch_prob)

        if verbose:
            logging.info("Now at iteration {}, "
                  "current loss is {}.".format(ii, loss))
            # print("Current global cluster label is {}"
            #       .format(glob_cluster_labels))
            # print("Sampling a random modality...")
            # layer = np.random.randint(L)
            # print("Current cluster label for the {}-th modality is {}"
            #       .format(layer, ind_cluster_labels[layer]))

        if np.abs(prev_loss - loss) < tol:
            centroids = [compute_centroids(
                X[layer], init_cluster_labels[layer], n_clusters
            ) for layer in range(L)]

            return ind_cluster_labels, glob_cluster_labels, centroids

        # prev_glob_cluster_labels = glob_cluster_labels
        prev_loss = loss
        init_cluster_labels = ind_cluster_labels

    centroids = [compute_centroids(
        X[layer], init_cluster_labels[layer], n_clusters
    ) for layer in range(L)]

    return ind_cluster_labels, glob_cluster_labels, centroids



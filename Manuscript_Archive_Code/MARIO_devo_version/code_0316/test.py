import numpy as np
from numba import jit

@jit(nopython=True)
def compute_centroids(X, cluster_labels, n_clusters):
    centroids = np.zeros((n_clusters, X.shape[1]))
    for kk in range(n_clusters):
        # current_idx = np.nonzero(np.array(cluster_labels == kk))
        indices = np.nonzero(cluster_labels==kk)[0]
        print(indices)
        current_centroid = np.zeros(X.shape[1])
        for ii in indices:
            current_centroid = current_centroid + X[ii, :]/len(indices)
        centroids[kk, :] = current_centroid
    return centroids


if __name__ == "__main__":
    X = np.array([[1, 3, 2], [1, 2, 2], [1, 1, 1], [2, 2, 2]])
    cluster_labels = np.array([0, 0, 1, 1])
    n_clusters = 2
    compute_centroids(X, cluster_labels, n_clusters)

import numpy as np


def normalize(df):
    """Center each column and scale each column to have unit standard deviation.

    Parameters
    ----------
    df : array-like of shape (n_samples, n_features)
        Two dimensional array to be normalized.
    Returns
    -------
    df : array-like of shape (n_samples, n_features)
        Two dimensional array after normalization.
    """
    df = df - np.mean(df, axis=0)
    df = df / np.std(df, axis=0)
    return df


def log_normalize(X, scale_factor=1e6):
    """Log-normalize the data.

    Parameters
    ----------
    X: array-like of shape (n_samples, n_features)
        Two dimensional array to normalize.
    scale_factor: float, default=1e6
        Multiple all entries of X by scale_factor before taking log.

    Returns
    -------
    X: array-like of shape (n_samples, n_features)
        Two dimensional array after normalization.
    """
    row_cnts = np.sum(X, axis=1)
    X = ((X * scale_factor).T / row_cnts).T
    X = np.log(X + 1)
    return X


def cdist_correlation(X, Y):
    """Calculate pair-wise Pearson correlation between X and Y.

    Parameters
    ----------
    X: array-like of shape (n_samples_of_X, n_features)
        First dataset.
    Y: array-like of shape (n_samples_of_Y, n_features)
        Second dataset.

    Returns
    -------
    array-like of shape (n_samples_of_X, n_samples_of_Y)
        The (i, j)-th entry is the Pearson correlation between i-th row of X and j-th row of Y.
    """
    n, p = X.shape
    m, p2 = Y.shape
    assert p2 == p

    X = (X.T - np.mean(X, axis=1)).T
    Y = (Y.T - np.mean(Y, axis=1)).T

    X = (X.T / np.sqrt(np.sum(X ** 2, axis=1))).T
    Y = (Y.T / np.sqrt(np.sum(Y ** 2, axis=1))).T

    return 1 - X @ Y.T


def check_min_dist(dist_mat, min_dist):
    """Make sure min(dist_mat) >= min_dist.

    Parameters
    ----------
    dist_mat : array-like of shape (n_samples_1, n_samples_2)
        A two-dimensional distance matrix.
    min_dist : flaot
        Desired minimum distance.

    Returns
    -------
    dist_mat : array-like of shape (n_samples_1, n_samples_2)
        A two-dimensional distance matrix with min(dist_mat) >= min_dist.
    """
    current_min_dist = np.min(dist_mat)
    if current_min_dist < min_dist:
        dist_mat = dist_mat - current_min_dist + min_dist

    return dist_mat
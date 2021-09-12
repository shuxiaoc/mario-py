import numpy as np
from sklearn.cross_decomposition import CCA
from collections import Counter


def normalize(df):
    """
    Given a pandas dataframe df,
        1) center each column;
        2) scale each column to have unit standard deviation
    """
    df = df - np.mean(df, axis=0)
    df = df / np.std(df, axis=0)
    return df


def get_cancor(X, Y, n_components=10, max_iter=1000):
    assert X.shape[1] >= n_components
    assert Y.shape[1] >= n_components
    ###### updated on 0710 debug for lingalg
    #try:
    #    print("in try")
    cca = CCA(n_components=n_components, max_iter=max_iter)
    cca.fit(X, Y)
    X_c, Y_c = cca.transform(X, Y)
    #except:
    #    print("in except")
    #    pd.DataFrame(X).to_csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/debug/tempX.csv")
    #    pd.DataFrame(Y).to_csv("/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/debug/tempY.csv")
    
    ####### update end 
    cancor = np.corrcoef(
        X_c, Y_c, rowvar=False).diagonal(
        offset=cca.n_components)
    return cancor, cca


def log_normalize(X, scale_factor=1e6):
    n, p = X.shape
    row_cnts = np.sum(X, axis=1)
    X = ((X * scale_factor).T / row_cnts).T
    X = np.log(X + 1)
    return X


def cdist_correlation(X, Y):
    n, p = X.shape
    m, p2 = Y.shape
    assert p2 == p

    X = (X.T - np.mean(X, axis=1)).T
    Y = (Y.T - np.mean(Y, axis=1)).T

    X = (X.T / np.sqrt(np.sum(X ** 2, axis=1))).T
    Y = (Y.T / np.sqrt(np.sum(Y ** 2, axis=1))).T

    return (1 - X @ Y.T)



#
# def split_into_half(x):
#     assert len(x) > 1
#     x1 = np.random.choice(x, len(x) // 2, replace=False)
#     x2 = np.array([ele for ele in x if ele not in x1])
#     return x1, x2
#
#
# def convert_to_rank(D):
#     D = np.argsort(D, axis=1)
#     D = np.argsort(D, axis=1)
#     return D
#
#
#
# def align_columns(U1, U2):
#     n, p = U1.shape
#     assert n == U2.shape[0]
#     assert p == U2.shape[1]
#
#     U1 = np.sort(U1, axis=0)
#     U2 = np.sort(U2, axis=0)
#
#     sign_vec = []
#
#     for i in range(p):
#         err1 = np.sum((U1[:, i] - U2[:, i]) ** 2)
#         err2 = np.sum((U1[:, i] + U2[:, i]) ** 2)
#         if err1 <= err2:
#             sign_vec.append(1)
#         else:
#             sign_vec.append(-1)
#
#     sign_vec = np.array(sign_vec, dtype=int)
#     U2 = U2 * sign_vec
#     return U1, U2, sign_vec
#
#
# def iter_refine(A, B, perm_init, max_iter, method='assignemnt',
#                 tol=1e-3, verbose=True):
#     """
#     A, B must be symmetric matrices
#     attemps to solve the quadratic assignment problem
#     min <A, PI B PI.T> via a coordinate descent procedure:
#     PI_next = argmin <A, B PI_current.T>
#
#     note that effectively, the distance matrix is
#     dist(A, B@PI.T), i.e., dist(A, B[:, pi])
#
#     todo: output all estimated permutations and choose the best one
#     """
#
#     obj_prev = -np.Inf
#     for i in range(max_iter):
#         obj_val = np.sum(A * B[perm_init, :][:, perm_init])
#         if verbose:
#             print("At iteration {}, current objective value is {}.".format(
#                 i, obj_val
#             ))
#         if np.abs(obj_val - obj_prev) < tol:
#             return perm_init
#
#         obj_prev = obj_val
#         dist_mat = cdist_sqeuclidean(A, B[:, perm_init])
#         perm_init = match_cells(dist_mat, 'assignment')
#
#     return perm_init
#
#
# def cdist_sqeuclidean(X, Y):
#     # this is faster than scipy distance.cdist
#     # when the dimension is large
#     n, p = X.shape
#     m, p2 = Y.shape
#     assert p2 == p
#     X_row_norms = np.sum(X ** 2, axis=1)
#     Y_row_norms = np.sum(Y ** 2, axis=1)
#
#     # dist_mat = np.tile((X_row_norms[np.newaxis, :]).T, [1, m])
#     # dist_mat = dist_mat + np.tile(Y_row_norms, [n, 1])
#     # dist_mat = dist_mat - 2 * X @ Y.T
#
#     dist_mat = -2 * X @ Y.T
#     dist_mat = (dist_mat.T + X_row_norms).T
#     dist_mat = dist_mat + Y_row_norms
#
#     return dist_mat - np.min(dist_mat)
#
#
# def cdist_w1(X, Y, method='euclidean'):
#     n, p = X.shape
#     m, p2 = Y.shape
#     assert n == m
#     assert p == p2
#
#     X = np.sort(X, axis=1)
#     Y = np.sort(Y, axis=1)
#
#     if method == 'sqeuclidean':
#         return (cdist_sqeuclidean(X, Y))
#     elif method == 'correlation':
#         return (cdist_correlation(X, Y))
#     elif method == 'euclidean':
#         return (np.sqrt(cdist_sqeuclidean(X, Y)))
#     else:
#         raise NotImplementedError
#
#
#
# def align_eigenvecs(U):
#     """
#     The signs of each column of U is set to
#         +1 if originally positive entries are dominant
#         -1 if originally negative entries are dominant
#     """
#     pos_num = np.sum(U > 0, axis=0)
#     s = 2 * (pos_num > (U.shape[0] // 2)) - 1
#     return U * s
#
#
# def get_cell_scores(U, V):
#     """
#     Given the latent embeddings of two datasets,
#     give quality score to each row of U and V,
#     the higer the score, the more certein we are
#     that this pair should be matched together.
#     """
#
#     U = U / np.sqrt(np.sum(U ** 2, axis=0))
#     V = V / np.sqrt(np.sum(V ** 2, axis=0))
#     return (U @ V.T).diagonal()
#
#
# def invert_permutation(perm):
#     return np.argsort(perm)
#
#
# def select_top_pairs(dist_mat, perm_hat, topk=1000):
#     good_cell_idx = (
#                         np.argsort(dist_mat[:, perm_hat].diagonal())
#                     )[:topk]
#     return good_cell_idx[:topk]
#
# # if __name__ == '__main__':
# # X = np.random.normal(0, 1, [1000, 300])
# # Y = np.random.normal(0, 1, [1000, 300])
# # d1 = cdist_sqeuclidean_old(X, Y)
# # d2 = cdist_sqeuclidean(X, Y)
# # print(np.linalg.norm(d1-d2))

import numpy as np
import warnings
from sklearn.utils.extmath import randomized_svd
from . import match_utils, utils
from .cluster import spectral_clustering, jr_kmeans


class Mario(object):
    def __init__(self, df1, df2, normalization=True):
        """Initialize the Mario object.

        Parameters
        ----------
        df1 : array-like of shape (n_samples_1, n_features_1)
            The first dataset.
        df2 : array-like of shape (n_samples_2, n_features_2)
            The second dataset.
        normalization : bool, default=True
            If true, center each column and scale each column to have unit standard deviation.
        """
        self.min_dist = 1e-5

        # parameters related to datasets
        if normalization:
            self.df1 = utils.normalize(df1)
            self.df2 = utils.normalize(df2)
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
        """Compute distance matrix based on overlapping features.

        Parameters
        ----------
        n_components : int
            Number of SVD components to keep

        Returns
        -------
        dist_ovlp : array-like of shape (n1, n2)
            The distance matrix based on the overlapping features.
        s : array-like of shape (n_components, )
            Vector of singular values.
        """
        if n_components > len(self.ovlp_features):
            warnings.warn("n_components exceed the number of overlapping features,"
                            " set it to be the number of overlapping features.")
            n_components = len(self.ovlp_features)

        self.n_components['ovlp'] = n_components

        if not (
                self.stacked_svd['U'] is not None and self.stacked_svd['s'] is not None
                and self.stacked_svd['Vh'] is not None and len(self.stacked_svd['s']) >= n_components
        ):
            # Cached results are not valid, do SVD
            arr1 = utils.normalize(self.df1[self.ovlp_features]).to_numpy()
            arr2 = utils.normalize(self.df2[self.ovlp_features]).to_numpy()

            self.stacked_svd['U'], self.stacked_svd['s'], self.stacked_svd['Vh'] = \
                randomized_svd(np.concatenate((arr1, arr2), axis=0), n_components=n_components)
            if n_components == len(self.ovlp_features):
                dist_mat = utils.cdist_correlation(arr1, arr2)
            else:
                svd1 = self.stacked_svd['U'][:self.n1, :] @ np.diag(self.stacked_svd['s']) @ self.stacked_svd['Vh']
                svd2 = self.stacked_svd['U'][self.n1:, :] @ np.diag(self.stacked_svd['s']) @ self.stacked_svd['Vh']
                dist_mat = utils.cdist_correlation(svd1, svd2)
        else:
            # use cached results
            svd1 = self.stacked_svd['U'][:self.n1, :n_components] @ np.diag(self.stacked_svd['s'][:n_components]) \
                   @ self.stacked_svd['Vh'][:n_components, :]
            svd2 = self.stacked_svd['U'][self.n1:, :n_components] @ np.diag(self.stacked_svd['s'][:n_components]) \
                   @ self.stacked_svd['Vh'][:n_components, :]
            dist_mat = utils.cdist_correlation(svd1, svd2)

        # make sure min_dist is at least self.min_dist
        dist_mat = utils.check_min_dist(dist_mat, self.min_dist)
        self.dist['ovlp'] = dist_mat
        return dist_mat, self.stacked_svd['s'][:n_components]

    def search_minimum_sparsity(self, dist_mat, slackness=200, init_sparsity=None, verbose=True):
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
        slackness : int, default=200
            Binary search terminates when k_right - k_left <= slackness;
            an exact binary search corresponds to slackness = 0
        init_sparsity : int, default=None
            Binary search starts from k=init_sparsity. If None, start from the middle.
        verbose : bool, default=None
            If True, print the progress.

        Returns
        -------
        k_left : int
            If sparsity<k_left, then there is no valid matching.
        k_right : int
            If sparsity>=k_right, then there is a valid matching.
        """
        return match_utils.search_minimum_sparsity(
            dist_mat, slackness, init_sparsity, self.m_min,
            self.m_max, self.num_cells_to_use, self.min_dist, verbose
        )

    def specify_matching_params(self, n_matched_per_cell):
        """Specify how many cells in the second dataset are to be matched with one cell in the first dataset.

        Parameters
        ----------
        n_matched_per_cell : int
            How many cells in the second dataset are to be matched with one cell in the first dataset.
        """
        self.n_matched_per_cell = n_matched_per_cell
        if self.n1 * n_matched_per_cell > self.n2:
            raise ValueError("Not enough cells in Y data!")
        self._specify_matching_params(1, n_matched_per_cell, self.n1 * n_matched_per_cell)

    def _specify_matching_params(self, m_min, m_max, num_cells_to_use):
        """Specify the matching parameters.

        Parameters
        ----------
        m_min : int
            Each row in the first dataset is matched to at least m_min many rows in the second dataset.
        m_max : int
            Each row in the first dataset is matched to at most m_max many rows in the second dataset.
        num_cells_to_use : int
            Total number of rows to use in the second dataset.
        -------
        """
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
        """Do cell matching.

        Parameters
        ----------
        dist_mat : str or a user-specified distance matrix, default='ovlp'
            If 'ovlp', then match using the distance matrix computed from overlapping features;
            if 'all', then match using the distance matrix computed from all the features;
            if a user-specified array-like of shape (n1, n2), then match using this distance matrix.
        sparsity : int
            Number of nearest neighbors to keep in the distance matrix.
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
        if isinstance(dist_mat, str):
            if dist_mat not in self.dist or self.dist[dist_mat] is None:
                raise ValueError("Distance not found!")
            self.sparsity[dist_mat] = sparsity
            self.matching[dist_mat] = match_utils.match_cells(
                self.dist[dist_mat], sparsity, self.m_min, self.m_max, self.num_cells_to_use, self.min_dist, mode
            )
            return self.matching[dist_mat]
        else:
            return match_utils.match_cells(
                dist_mat, sparsity, self.m_min, self.m_max, self.num_cells_to_use, self.min_dist, mode
            )

    def _align_modalities(self, matching):
        """Align df1 so so that cell i in df1 is matched to the "averaged cell" in cells matching[i] in df2.

        Parameters
        ----------
        matching : list
            A list of (potentially variable length) lists;
            it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.

        Returns
        -------
        X : array-like of shape (n_samples, n_features_1)
            The first dataset.
        Y : array-like of shape (n_samples, n_features_2)
            The second dataset after alignment.
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
        """Align df1 and df2 using matching, then fit a CCA.

        Parameters
        ----------
        matching : list
            A list of (potentially variable length) lists;
            it holds that the i-th row in the first dataset is matched to the res[i]-th row in the second dataset.
        n_components : int
            Number of components for CCA.
        max_iter : int
            Maximum iteration for CCA.

        Returns
        -------
        cancor: array-like of shape (n_components, )
            Vector of canonical components.
        cca: CCA
            CCA object.
        """
        if n_components > min(self.p1, self.p2):
            warnings.warn('n_components must be <= the dimensions of the two datasets, '
                          'set it to be equal to the minimum of the dimensions of the two datasets')
            n_components = min(self.p1, self.p2)

        X, Y = self._align_modalities(matching)
        cancor, cca = utils.get_cancor(X, Y, n_components, max_iter)

        return cancor, cca

    def compute_dist_all(self, matching='ovlp', n_components=20, max_iter=5000):
        """Given matching, align df1 and df2, fit a CCA, then use CCA scores to get the distance matrix.

        Parameters
        ----------
        matching : str or list
            Either 'ovlp', meaning that we use the matching obtained from the overlapping features;
            or a list of (potentially variable length) lists, cell i in df1 is matched to cell matching[i] in df2.
        n_components : int
            Number of CCA components.
        max_iter : int
            Max number of iterations when fitting CCA.

        Returns
        -------
        dist_mat : array-like of shape (n1, n2)
            The distance matrix.
        cancor: array-like of shape (n_components, )
            Vector of canonical components.
        """
        if n_components <= 1:
            n_components = 2
            warnings.warn('n_components must be at least 2, '
                          'set it to 2')

        if n_components > min(self.p1, self.p2):
            warnings.warn('n_components must be <= the dimensions of the two datasets, '
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
                dist_mat = utils.cdist_correlation(self.ovlp_scores['x'], self.ovlp_scores['y'])
            else:
                # use cached results
                dist_mat = utils.cdist_correlation(self.ovlp_scores['x'][:, :n_components],
                                             self.ovlp_scores['y'][:, :n_components])
                cancor = self.ovlp_cancor[:n_components]
        else:
            # use user-specified matching
            cancor, cca = self.fit_cca(matching, n_components, max_iter)
            df1_cca, df2_cca = cca.transform(self.df1, self.df2)
            dist_mat = utils.cdist_correlation(df1_cca, df2_cca)

        self.dist['all'] = utils.check_min_dist(dist_mat, self.min_dist)
        return self.dist['all'], cancor

    def _matchable(self, n_sim=20, top_k=10, flip_prob=0.3, subsample_prop=1, verbose=True):
        """Calculate the p-value for the hypothesis that two datasets are matchable

        Parameters
        ----------
        n_sim : int, default=20
            Number of simulation rounds.
        top_k : int, default=10
            The mean of the top_k canonical correlations will be taken as the test statistic.
        flip_prob : float, default=0.3
            Probability of sign flip.
        subsample_prop : float, default=1
            Random subsample df1 and df2 to speed up the computation, 1 means using the full data.
        verbose : bool, default=True
            Print details if True.

        Returns
        -------
        pval_ovlp : float
            The p-value for the hypothesis that two datasets are matchable based on overlapping features.
        pval_all : float
            The p-value for the hypothesis that two datasets are matchable based on all the features.
        """
        cancor_ovlp_list = []
        cancor_all_list = []
        if self.n_components['ovlp'] is None or self.n_components['all'] is None:
            raise ValueError("Please do both initial and refined matching first before calling this function.")

        if subsample_prop < 1:
            df1_subsampled = self.df1.iloc[np.random.choice(self.n1, int(self.n1 * subsample_prop), replace=False)]
            df2_subsampled = self.df2.iloc[np.random.choice(self.n2, int(self.n2 * subsample_prop), replace=False)]
        else:
            df1_subsampled = self.df1
            df2_subsampled = self.df2

        match_subsampled = Mario(df1_subsampled, df2_subsampled, normalization=True)
        match_subsampled.specify_matching_params(self.n_matched_per_cell)
        # do matching
        _ = match_subsampled.compute_dist_ovlp(self.n_components['ovlp'])
        try:
            _ = match_subsampled.match_cells('ovlp', self.sparsity['ovlp'], mode='auto')
        except ValueError:
            # too sparse, find the suitable sparsity level
            warnings.warn(
                'When using ovlp features, current sparsity config '
                'is too sparse, finding a suitable sparsity level...'
            )
            _, new_sparsity = match_subsampled.search_minimum_sparsity(
                match_subsampled.dist['ovlp'], slackness=200, init_sparsity=self.sparsity['ovlp'] + 1, verbose=True
            )
            warnings.warn('The new sparsity level is {}'.format(new_sparsity))
            _ = match_subsampled.match_cells('ovlp', new_sparsity, mode='auto')

        # use all features
        _ = match_subsampled.compute_dist_all('ovlp', n_components=self.n_components['all'])
        try:
            _ = match_subsampled.match_cells('all', self.sparsity['all'], mode='auto')
        except ValueError:
            warnings.warn(
                'When using all features, current sparsity config '
                'is too sparse, finding a suitable sparsity level...'
            )
            # too sparse, find the suitable sparsity level
            _, new_sparsity = match_subsampled.search_minimum_sparsity(
                match_subsampled.dist['all'], slackness=200, init_sparsity=self.sparsity['all'] + 1, verbose=True
            )
            warnings.warn('The new sparsity level is {}'.format(new_sparsity))
            _ = match_subsampled.match_cells('all', new_sparsity, mode='auto')

        # calculate test statistics
        cancor_ovlp_obs = np.mean(match_subsampled.fit_cca(match_subsampled.matching['ovlp'], n_components=top_k)[0])
        cancor_all_obs = np.mean(match_subsampled.fit_cca(match_subsampled.matching['all'], n_components=top_k)[0])

        for ii in range(n_sim):
            if verbose:
                print("Random sign flip, round {}...".format(ii))
            max_trial = 100
            # looks like we may get SVD not converge error
            # and it is due to linear algebra libraries in the linux server
            # so try max_trial times until success
            trial_idx = 0
            while trial_idx < max_trial:
                trial_idx += 1
                try:
                    rand_signs_1 = 2 * np.random.binomial(1, 1-flip_prob, match_subsampled.n1) - 1
                    rand_signs_2 = 2 * np.random.binomial(1, 1-flip_prob, match_subsampled.n2) - 1
                    df1_flipped = (match_subsampled.df1.T * rand_signs_1).T
                    df2_flipped = (match_subsampled.df2.T * rand_signs_2).T

                    match_flipped = Mario(df1_flipped, df2_flipped, normalization=False)
                    match_flipped.specify_matching_params(self.n_matched_per_cell)

                    # use ovlp features
                    _ = match_flipped.compute_dist_ovlp(self.n_components['ovlp'])
                    try:
                        _ = match_flipped.match_cells('ovlp', self.sparsity['ovlp'], mode='auto')
                    except ValueError:
                        # too sparse, find the suitable sparsity level
                        warnings.warn(
                            'When using ovlp features, current sparsity config '
                            'is too sparse, finding a suitable sparsity level...'
                        )
                        _, new_sparsity = match_flipped.search_minimum_sparsity(
                            match_flipped.dist['ovlp'], slackness=200,
                            init_sparsity=self.sparsity['ovlp'] + 1, verbose=True
                        )
                        warnings.warn('The new sparsity level is {}'.format(new_sparsity))
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
                        warnings.warn(
                            'When using all features, current sparsity config '
                            'is too sparse, finding a suitable sparsity level...'
                        )
                        # too sparse, find the suitable sparsity level
                        _, new_sparsity = match_flipped.search_minimum_sparsity(
                            match_flipped.dist['all'], slackness=200,
                            init_sparsity=self.sparsity['all'] + 1, verbose=True
                        )
                        warnings.warn('The new sparsity level is {}'.format(new_sparsity))
                        _ = match_flipped.match_cells('all', new_sparsity, mode='auto')

                    # calculate the median/mean of top_k cancors
                    cancor_all_list.append(np.mean(
                        match_flipped.fit_cca(match_flipped.matching['all'], n_components=top_k)[0]
                    ))
                except Exception as exception:
                    Warning("{} occurs at the {}/{}-th trial, try flipping again...".format(
                        type(exception).__name__, trial_idx, max_trial)
                    )
                    continue
                else:
                    # the code will arrive here only if try has succeeded
                    break
            if trial_idx == max_trial:
                raise ValueError("{} consecutive trials have all failed.".format(max_trial))
            elif trial_idx > 1:
                warnings.warn("Success at {}/{}-th trial".format(trial_idx, max_trial))

        # calculate p-values
        pval_ovlp = np.mean([cancor_ovlp_list[ii] >= cancor_ovlp_obs for ii in range(n_sim)])
        pval_all = np.mean([cancor_all_list[ii] >= cancor_all_obs for ii in range(n_sim)])
        return pval_ovlp, pval_all

    def matchable(self, n_sim=20, top_k=10, flip_prob=0.3, subsample_prop=1, subsample_rounds=1, verbose=True):
        """Calculate the p-value for the hypothesis that two datasets are matchable.

        Parameters
        ----------
        n_sim : int, default=20
            Number of simulation rounds.
        top_k : int, default=10
            The mean of the top_k canonical correlations will be taken as the test statistic.
        flip_prob : float, default=0.3
            Probability of sign flip.
        subsample_prop : float, default=1
            Subsample to speed up computation, 1 means using the full data.
        subsample_rounds : int, default=1
            The subsampling operation is done this many rounds and the p-values are averaged.
        verbose : bool, default=True
            Print details if True.

        Returns
        -------
        pval_ovlp : float
            The p-value for the hypothesis that two datasets are matchable based on overlapping features.
        pval_all : float
            The p-value for the hypothesis that two datasets are matchable based on all the features.
        """
        pval_ovlp = 0
        pval_all = 0
        for i in range(subsample_rounds):
            if verbose and subsample_prop < 1:
                print(f"Now at subsample round {i}...")
            curr_pval_ovlp, curr_pval_all = self._matchable(n_sim, top_k, flip_prob, subsample_prop, verbose)
            pval_ovlp += curr_pval_ovlp
            pval_all += curr_pval_all

        return pval_ovlp/subsample_rounds, pval_all/subsample_rounds

    def interpolate(self, n_wts=10, top_k=10, verbose=True):
        """
        Let wt_vec be an evenly spaced list from 0 to 1 with length n_wts.
        For each wt in wt_vec, do matching on (1-wt)*dist_ovlp + wt*dist_all,
        and select the best wt according to the mean of top_k canonical correlations.

        Parameters
        ----------
        n_wts : int, default=10
            wt_vec is a evenly spaced list from 0 to 1 with length n_wts.
        top_k : int, default=10
            The mean of top_k canonical correlations is taken as the quality measure.
        verbose : bool, default=True
            Print details if True.

        Returns
        -------
        best_wt : float
            The best wt in wt_vec.
        best_matching : list
            The matching corresponds to best_wt.
        """
        wt_vec = np.linspace(0, 1, n_wts)
        max_cancor = float('-inf')
        for ii in range(n_wts):
            if verbose:
                print('Now at iteration {}, wt={}'.format(ii, wt_vec[ii]))
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
                        print(
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
        1) align df1 and df2, get df1_aligned and df2_aligned;
        2) do CCA on df1_aligned and df2_aligned to get CCA scores, call them df1_cca and df2_cca;
        3) do spectral clustering on the averaged df, i.e., (df1_cca+df2_cca)/2 to get init_labels;
        4) given init_labels, run jr_kmeans on df1_aligned and df2_aligned.

        Parameters
        ----------
        matching : list
            Matching used to align df1 and df2, must be a perfect matching
        n_clusters : int
            Number of clusters desired.
        n_components : int, default=20
            Number of CCA components for initial clustering.
        mismatch_prob : float, default=0.1
            Approximately mismatch_prob proportion of cells have different cluster labels.
        max_iter : int, default=50
            Max iteration desired.
        tol : float, default=1e-5.
            If the objective function changes less than tol, terminates jr_kmeans.
        verbose : bool, default=True
            If True, print details when running jr_kmeans.

        Returns
        -------
        cluster_labels : list
            cluster_labels[0] (resp. [1]) is the cluster labels for df1_aligned (resp. df2_aligned).
        centroids : list
            centroids[0] (resp. [1]) is an array of shape (n_clusters, df1.shape[1])
            (resp. of shape (n_clusters, df2.shape[1])) representing the centroids for df1_aligned (resp. df2_aligned).
        """

        assert match_utils.is_perfect_matching(matching)

        df1_aligned, df2_aligned = self._align_modalities(matching)
        _, cca = utils.get_cancor(df1_aligned, df2_aligned, n_components)
        df1_cca, df2_cca = cca.transform(df1_aligned, df2_aligned)
        # df1_cca = df1_cca / np.std(df1_cca, axis=0)
        # df2_cca = df2_cca / np.std(df2_cca, axis=0)
        df1_cca = utils.normalize(df1_cca)
        df2_cca = utils.normalize(df2_cca)
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
        matching : str or list, default='wted'
            A matching to be filtered, must be a perfect matching.
            If 'ovlp', 'all', or 'wted', then use the matching obtained from
            the overlapping features, all the features, and their interpolation, respectively.
            Otherwise, a list corresponds to user-specified matching.
        n_clusters : int, default=15,
            How many clusters desired in self.jr_clustering.
        n_components : int, default=20
            Number of CCA components in the initialization step in self.jr_clustering.
        bad_prop : float, default=0.1
            Approximate proportion of bad matches to be filtered out.
        max_iter : int, default=50
            Run self.jr_clustering for how many iterations.
        tol : int, default=1e-5
            Terminate self.jr_clustering when the change of objective function is <= tol.
        verbose : bool, default=True
            Print all the details when True.

        Returns
        -------
        good_matching : list
            If for some i, matching[i] is an empty list, it means that i is bad and has been filtered out.
        """

        # proportion of good cells is approximately (1-mismatch_prob)^2 + mismatch_prob^2,
        # this is because the model says that a cluster label flips w.p. mismatch_prob
        # and we declare a pair to be good if either (1) both not flipped, or (2) both flipped
        # so given bad_prop, can solve for the corresponding mismatch_prob

        assert 1e-10 < bad_prop < 0.5 - 1e-10

        if n_components > min(self.p1, self.p2):
            warnings.warn('n_components must be <= the dimensions of the two datasets, '
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
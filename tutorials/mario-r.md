Tutorial on Mario-R
================
someone
8/26/2021

**We can copy everything in mario-py tutorial here.**

We will use `reticulate` to call modules from `mario-py`.

``` r
library(reticulate)
library(tidyverse)
use_python("/Users/shuxiaochen/miniconda3/envs/cell/bin/python") # specify the python under which mario-py is installed
mario.match <- import("mario.match") # main mario-py module
Mario = mario.match$Mario # main object
pipelined_mario = mario.match$pipelined_mario # for running the overall pipeline
int = as.integer # need to manually convert numeric to integer when using reticulate
```

# Exploratory analysis with subsampled data

*INSERT MORE DETAILS.* We use subsampled data to select hyperparameters,
test for matchability, etc. Also, need to explain the meanings of each
argument

``` r
# import and subsample data
df1 <- read_csv("/Users/shuxiaochen/Dropbox/Research/ongoing/single-cell-integration/data-biology/drop_out_test/vaxaart-wb50k-cytof.csv") %>% select(-X1)
df2 <- read_csv("/Users/shuxiaochen/Dropbox/Research/ongoing/single-cell-integration/data-biology/drop_out_test/wcctg-wb50k-cytof.csv") %>% select(-X1)
set.seed(2667)
df1_sub <- df1 %>% slice_sample(prop = 0.01)
df2_sub <- df2 %>% slice_sample(prop = 0.02)

# extract cluster labels
df1_labels = df1 %>% select(cluster.info) %>% pull()
df2_labels = df2 %>% select(cluster.info) %>% pull()
df1_sub_labels = df1_sub %>% select(cluster.info) %>% pull()
df2_sub_labels = df2_sub %>% select(cluster.info) %>% pull()

# remove non-numerical columns
df1 <- df1 %>% select(-cluster.info)
df2 <- df2 %>% select(-cluster.info)
df1_sub <- df1_sub %>% select(-cluster.info)
df2_sub <- df2_sub %>% select(-cluster.info)
```

## Initial matching with overlapping features

``` r
mario = Mario(df1_sub, df2_sub, normalization = TRUE)
dist_ovlp = mario$compute_dist_ovlp(n_components = int(20))
```

``` r
mario$specify_matching_params(int(1))
```

``` r
mario$search_minimum_sparsity(mario$dist$ovlp, slackness=int(1), init_sparsity=NULL, verbose=TRUE)
```

    ## [[1]]
    ## [1] 13
    ## 
    ## [[2]]
    ## [1] 14

``` r
matching_ovlp = mario$match_cells('ovlp', sparsity=int(100), mode='auto')
```

## matching using all features

``` r
dist_all = mario$compute_dist_all(matching = 'ovlp', n_components = int(15))
```

``` r
matching_all = mario$match_cells(dist_mat = 'all', sparsity=NULL, mode='auto')
```

## Testing for matchability

``` r
mario$matchable(n_sim=int(5), top_k=int(5), flip_prob=0.3, subsample_prop=int(1), verbose=TRUE)
```

    ## [[1]]
    ## [1] 0
    ## 
    ## [[2]]
    ## [1] 0

## Finding the best interpolation

``` r
interpolate_res = mario$interpolate(n_wts=int(5), top_k=int(5), verbose=TRUE)
```

## Filtering low-quality matched pairs

``` r
matching_final  = mario$filter_bad_matches(
  matching='wted', n_clusters=int(10), n_components=int(10), bad_prop=0.2,
  max_iter=int(5), tol=1e-5, verbose=TRUE)
```

## k-NN matching

``` r
matching_knn = mario$knn_matching('wted', k=int(5))
```

## Joint embedding

``` r
cca_res = mario$fit_cca('final', n_components=as.integer(20), max_iter=as.integer(10000))
df1_cca = cca_res[[2]]$x_scores_
df2_cca = cca_res[[2]]$y_scores_
```

# Full pipeline

``` r
pipelined_mario = mario.match$pipelined_mario
pipelined_res = pipelined_mario(
  data_lst=list(df1, df2), normalization=TRUE, n_batches=int(4),
  n_matched_per_cell=int(3), sparsity_ovlp=NULL, sparsity_all=NULL,
  n_components_ovlp=int(20), n_components_all=int(20),
  n_cancor=int(5), n_wts=int(3),
  n_clusters=int(10), n_components_filter=int(10), bad_prop=0.2, max_iter_filter=int(5),
  knn=int(5), embed_dim=int(20), max_iter_embed=int(50), save_path='../code/data', verbose=FALSE
)
```

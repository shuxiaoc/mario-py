Tutorial on Mario-R
================
Bokai Zhu
09/20/2021

#MARIO matching and integration pipeline

This tutorial will guide you through a list of general steps to perform
matching and integration of single cell proteomic datasets. This is the
tutorial for using MARIO in R.

This script will include matching and integration of two or more then
two datasets. First we will start with analysis on bone marrow cells
cite-seq and cytof cells. The data is same as we presented in the
[paper](httpe%20whater). First we load up `mario-py` from python to R
with the :

Here we will be using the `reticulate` package to call the python
functions in R. First we need to tell R which python virtualenv to use:

``` r
library(reticulate)
library(tidyverse)
library(ggplot2) # for plotting
library(scales)

myenvs=conda_list() # get conda virtualenv list
envname=myenvs$name[12] # specify which virtualenv to use, should use the one for MARIO-py
use_condaenv(envname, required = TRUE)
mario.match <- import("mario.match") # import main mario-py module
Mario = mario.match$Mario # main object
pipelined_mario = mario.match$pipelined_mario # for running the overall pipeline
int = as.integer # need to manually convert numeric to integer when using reticulate
```

##1. Read in data

We read in the two datasets, both protein features, with row as cells
and column with features: overlapping or non-overlapping protein types
(remember to make sure overlapping proteins have the same name accross
datasets). Here the two datasets we used also contains a column of cell
IDs and a column of cell types. These are **optional** and will not
affect the **MARIO** mathcing and integration.

**Note**: *If there exists inner sturcture of rows and cell types in the
dataset (eg. an imaging data where row 1 - row 5000 are all B cells due
to same FOV of germinal center) the best practice is to randomnize the
rows before matching, as extremely unbalanced cell population when
running in certain parameters (high `batches`) can sometimes give
sub-optimal results.*

``` r
# import and subsample data
df1 <- read_csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/github_tutorial_relate/bmcite_forSIM.csv") %>% select(-X1)
df2 <- read_csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/github_tutorial_relate/levine32_forSIM.csv") %>% select(-X1)

head(df1)
```

    ## # A tibble: 6 x 26
    ##   CD11a CD11c  CD123 CD127.IL7Ra  CD14  CD16 CD161  CD19 CD197.CCR7  CD25  CD27
    ##   <dbl> <dbl>  <dbl>       <dbl> <dbl> <dbl> <dbl> <dbl>      <dbl> <dbl> <dbl>
    ## 1  1.44 0.912 0.343        0.460 0.363 0.900 0.363 0.256      0.423 0.210 1.13 
    ## 2  2.14 0.486 0.0607       1.31  0.223 0.446 0.223 0.172      0.223 0     3.25 
    ## 3  2.82 0.707 0.115        1.34  0.219 0.477 0.382 0.219      0.257 0.137 2.76 
    ## 4  3.82 2.83  0.501        0.281 2.19  0.397 0.240 0.281      0.225 0.103 0.871
    ## 5  1.84 1.19  0.536        0.336 0.357 0.570 0.224 2.37       0.611 0.200 0.673
    ## 6  2.66 0.682 0.184        1.74  0.265 0.699 0.495 0.156      0.453 0.239 3.00 
    ## # … with 15 more variables: CD278.ICOS <dbl>, CD28 <dbl>, CD3 <dbl>,
    ## #   CD34 <dbl>, CD38 <dbl>, CD4 <dbl>, CD45RA <dbl>, CD45RO <dbl>, CD56 <dbl>,
    ## #   CD57 <dbl>, CD69 <dbl>, CD79b <dbl>, CD8 <dbl>, HLA.DR <dbl>,
    ## #   cluster.info <chr>

## 2. Parameter Screening and step-by-step MARIO

**MARIO** involves using a series of hyper-parameters. You can use the
default numbers we suggested in the package for those hyperparameters,
although to achieve optimal performance, we suggest tailoring these
hyper-parameters for each dataset and do a step-by-step MARIO analysis.
Here we will screen these parameters in a subsampled version of our
dataset. You can also run this part with the full dataset but it could
take some time depending on the data size.

``` r
# subsample data for tutorial
set.seed(2667)
df1_sub <- df1 %>% slice_sample(prop = 0.1)
df2_sub <- df2 %>% slice_sample(prop = 0.1)
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

### Parameters for overlapping feature matching:

*Note: here reticulate will be calling the python functions in R, the
numeric values for hyperparameter need the * `int()`.

``` r
mario = Mario(df1_sub, df2_sub, normalization = TRUE)
dist_ovlp = mario$compute_dist_ovlp(n_components = int(12))

# check singular value by plotting
dftx = data.frame(idx = 1:length(dist_ovlp[[2]]), Singular_value =dist_ovlp[[2]])
ggplot(data=dftx, aes(x=idx, y=Singular_value)) +
  geom_line() + scale_x_continuous(breaks = pretty_breaks())
```

![](https://github.com/shuxiaoc/mario-py/blob/main/tutorials/plots/R_fig1.png)<!-- -->

We use this plot to decide how many components to use during mathcing
with overlapping features, similar to the concept used with Elbow plots
during PCA reduction. From this plot any value **above 6** looks good,
we choose 10 here.

``` r
dist_ovlp = mario$compute_dist_ovlp(n_components = int(10))
mario$specify_matching_params(int(1)) # specify that we want default 1v1 matching
```

**Optional**: We also provide the option to run the **MARIO** matching
with sparsity. Running with sparsity could reduce the run time if the
data size is large, and potentially be a denoising step. We provide a
function find the best sparsity number to use:

``` r
# [optional] check the minimum valid sparsity level
mario$search_minimum_sparsity(mario$dist$ovlp, slackness=int(1), init_sparsity=NULL, verbose=TRUE)
```

    ## [[1]]
    ## [1] 55
    ## 
    ## [[2]]
    ## [1] 56

Based on the result, any sparsity level \> 66 will give valid matching,
thus we will use 100 and perform matching with overlapping features:

*Note: the sparsity value is related (in ratio) to the amount of cells
in n1.*

``` r
matching_ovlp = mario$match_cells('ovlp', sparsity=int(100), mode='auto')
```

### Parameters for refined all feature matching:

Then we will start testing how many components to use in the refined
matching (which will use all features including non-overlapping
features):

``` r
# compute distance matrix using all the features
dist_all = mario$compute_dist_all(matching = 'ovlp', n_components = int(20))
# check singular value by plotting
dftx = data.frame(idx = 1:length(dist_all[[2]]), canonical_correlations = dist_all[[2]])
ggplot(data=dftx, aes(x=idx, y=canonical_correlations)) +
  geom_line() + scale_x_continuous(breaks = pretty_breaks())
```

![](mario-r-bk_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Based on the plot we can see a steep drop of canonical correlation at 9
components, therefore, we should choose **\> 9** components for refined
matching, here we choose 12:

``` r
dist_all = mario$compute_dist_all(matching = 'ovlp', n_components = int(12))
# perform the refined matching
matching_all = mario$match_cells(dist_mat = 'all', sparsity=NULL, mode='auto')
```

### \[Optional\] MARIO Matchability test

However, in most cases we don’t have an annotation file to quickly
validate if the matching of MARIO is solid. Moreover, many datasets
should not be forcefully match and integrated (due to underlying biology
or poor quality data). Therefore, we provide a rigorous statistical
testing to confirm the two dataset can be matched. Note this step could
be slow depending on the datasize and parameter selected:

``` r
mario$matchable(n_sim=int(20), top_k=int(5), flip_prob=0.2, subsample_prop=int(1), verbose=TRUE)
```

    ## [[1]]
    ## [1] 0
    ## 
    ## [[2]]
    ## [1] 0

The test will give you two p values (#1 for ovlp matching, #2 for all
matching). We suggest **not to perform** a matching/integration analysis
of any kind on datasets with an anverage *p* \> 0.05.

### Find the best interpolation of overlapping and all matching

Then we can find the best balance between the `ovlp` feature matching
and `all` feature mathcing by the interpolate function. The best weight
will be finded and the corresponding match will be stored inside the
object:

``` r
interpolate_res = mario$interpolate(n_wts=int(5), top_k=int(5), verbose=TRUE)
```

### Joint regularized filtering of sub-optimal matching

At the very last step, we will filter out the potential sub-optimal
matching pairs via joint regularized filtering. The filtering will be
formed on the matching based on weighted interpolation. The parameter
`n_clusters` should be approximately be the number of cell populations,
and `n_components` will be the number of components used during
filtering, suggest using number in between the components used during
`ovlp` and `all` matching. The parameter `bad_prop` controls the
filtering strength and will be more stringent if increased, we suggest
using 0.1 - 0.3. After filtering, matched cell pairs will be less rows
of original data `df1_sub`.

``` r
matching_final  = mario$filter_bad_matches(
  matching='wted', n_clusters=int(10), n_components=int(10), bad_prop=0.2,
  max_iter=int(5), tol=1e-5, verbose=TRUE)
```

### Matching result access

All the matching produced by MARIO are stored in object `mario$matching`
and different matching results can be accessed via keys eg. `ovlp`,
`all`,`wted`, `final`. The matching result contains a list of numbers
which corresponding of the matching rows (eg.
`mario$matching['final']$final[1]` is 10136, meaning 1th row of
`df1_sub` is matched to 10136th row of `df2_sub`. Empty means match pair
beeing removed after filtering.

``` r
match_final = mario$matching['final']$final # geting matching final list, first element is empty
# some reformating
df_idx = c(1:length(match_final))
filtered_out = c(which(lapply(match_final,length)==0))#
match_final_df1 = df_idx[-filtered_out]
match_final_df2 = unlist(match_final)
df = data.frame(idx_1 = match_final_df1, idx_2 = match_final_df2)
head(df)
```

    ##   idx_1 idx_2
    ## 1     1  7256
    ## 2     3  2501
    ## 3     4  8446
    ## 4     5  9827
    ## 5     6  1210
    ## 6     7  7511

**optional**: We also provides the option to do the standard k-NN
matching of the MARIO matched pairs:

``` r
matching_knn = mario$knn_matching('wted', k=int(5))
```

### Joint embedding

With the cells having high quality pairs, we can use CCA to produce the
sub-space scores, that can be used for downstream analysis eg.
visualization with t-sne/umap, joint clustering and more.

``` r
cca_res = mario$fit_cca('final', n_components=as.integer(20), max_iter=as.integer(10000))
df1_cca = cca_res[[2]]$x_scores_
df2_cca = cca_res[[2]]$y_scores_ # cca scores used to perform t-sne etc
```

## 3. Pipeline version MARIO with one line

Previously we have showcased a step-by-step style **MARIO** and
parameter screening. For easy impelmentation, we also provide a compact
function `pipelined_mario` to perform all the steps in **MARIO** in one
line. Please note the time and memory usage is directly linked to cell
numbers in datasets and should not run locally if matching a total of
more than 100k cells. In our tutorial case we are matching the full size
datasets and a total of **\~30k** cells matching towards **\~100k**
cells, thus the code chunk could take \~40mins and \~25GB of space
(peak). The time and memory usage can be reduced by increasing the
`n_batches` parameter (inputing smaller cell numbers each iteration and
process in stream-line).

Remember you can use the default parameters but **MARIO** will perform
better with customized ones. Please refer to our
[documentation](https://github.com/shuxiaoc/mario-py/blob/fef06cdfb74017456136741dfeaf61fe0b94b608/src/mario/match.py)
inside the package for detailed information.

*Note: although the minimal sparsity level we searched previously was
around 100, here the n1 of each batch is around 6k cells, so we roughly
made sparsity to 300-400.*

``` r
pipelined_mario = mario.match$pipelined_mario
pipelined_res = pipelined_mario(
  data_lst=list(df1, df2), normalization=TRUE, n_batches=int(5),
  n_matched_per_cell=int(1), sparsity_ovlp=int(400), sparsity_all=int(400),
  n_components_ovlp=int(10), n_components_all=int(12),
  n_cancor=int(5), n_wts=int(4),
  n_clusters=int(10), n_components_filter=int(10), bad_prop=0.2, max_iter_filter=int(20),
  knn=FALSE, embed_dim=int(10), max_iter_embed=int(500), save_path='../code/data', verbose=FALSE
)
```

## 4. Downstream Analysis

After the MARIO matching and integration pipeline, we can perform
various downstream analysis. Here we will just do some simple showcase
of t-sne visualization. These results are also saved out and are
compatible for furthur analysis with popular packages eg `seurat`

First we can make the t-sne plots:

``` r
library(Rtsne)
set.seed(42)
cca_values = rbind(pipelined_res[[2]][[1]],pipelined_res[[2]][[2]]) # the rbind cca values
all_tsne_10=Rtsne(cca_values, check_duplicates = FALSE, num_threads = 10) # run with cca
all_tsne_data_10 <- data.frame(x = all_tsne_10$Y[,1], y = all_tsne_10$Y[,2])
label=as.factor(c(rep("cite",dim(cca_values)[1]/2), rep("cytof",dim(cca_values)[1]/2))) # lazy input of modality information
# produce tsne colored by modality
p=ggplot(all_tsne_data_10)  + 
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("MARIO")
p
```

![](mario-r-bk_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
match_final = pipelined_res[[1]][[2]] # geting matching final list, first element is empty
# some reformating
df_idx = c(1:length(match_final))
filtered_out = c(which(lapply(match_final,length)==0))#
match_final_df1 = df_idx[-filtered_out]
match_final_df2 = unlist(match_final)

label = c(df1_labels[match_final_df1], df2_labels[match_final_df2])
p=ggplot(all_tsne_data_10)  + 
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("MARIO")
p
```

![](mario-r-bk_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

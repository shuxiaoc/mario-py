Inculdes all the dependencies to perform benchmarking.

Methods were implemented in their own language (Python: MARIO, Scanorama; R: Seurat, Fastmnn, Liger), and benchmarked in R (matching accuracy and multiple metrics etc).

General steps:

```
#### 1. Test feature dropping: sequentially drop shared protein features

step 1:
delete_feature1_[*dataset*].sh
# will run MARIO (python), Seurat, Fastmnn, Liger (R) with dropped features

step 2:
delete_feature_scan-[*dataset*].py
# will run Scanorama (python) with dropped features

step 3:
delete_feature3_[*dataset*].sh
# Parrallel running code for calculate_metrics.R, calculate related metrics from reduced scores (all methods)

#### 3. Poor quality data by adding random noise



```


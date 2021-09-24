Inculdes all the dependencies to perform benchmarking.

Methods were implemented in their own language (Python: MARIO, Scanorama; R: Seurat, Fastmnn, Liger), and benchmarked in R (matching accuracy and multiple metrics etc).

```data``` includes the files used as input for benchmarking.

General steps:

```
#### 1. Test feature dropping: sequentially drop shared protein features

step 1:
delete_feature1_[*dataset*].sh
# will run MARIO (python), Seurat, Fastmnn, Liger (R) with dropped features

step 2: (Run in scanorama virtualenv)
delete_feature_scan-[*dataset*].py
# will run Scanorama (python) with dropped features

step 3:
delete_feature3_[*dataset*].sh
# Parrallel running code for calculate_metrics.R, calculate related metrics from reduced scores (all methods)

#### 2. Change cell type composition by dropping single types in DataY

step 1:
delete_celltype1_[*dataset*].sh
# will run MARIO (python), Seurat, Fastmnn, Liger (R) with dropped cell types in dataY.

step 2:
delete_celltypes_scan-[*dataset*].py
# will run Scanorama (python) with dropped cell types.


#### 3. Poor quality data by adding random noise

step 1:
addNoise1_[*dataset*].sh
# will run MARIO (python), Seurat, Fastmnn, Liger (R) with increasing random noise on both datasets.

step 2:
add_noise_scan_[*dataset*].py
# will run Scanorama (python) with increasing random noise on both datasets.

### 4. figure production

plot_making_script_[*dataset*].R
# will make the figures as presented in the paper.

```


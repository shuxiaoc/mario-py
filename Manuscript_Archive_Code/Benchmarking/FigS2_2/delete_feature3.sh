/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/bm_levine/sxc/metrics.csv' 'data/feature_deletion/bm_levine/sxc/orig' 'data/feature_deletion/bm_levine/sxc/embedding' 9 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/bm_levine/seurat/metrics.csv' 'data/feature_deletion/bm_levine/sxc/orig' 'data/feature_deletion/bm_levine/seurat/embedding' 9 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/bm_levine/liger/metrics.csv' 'data/feature_deletion/bm_levine/sxc/orig' 'data/feature_deletion/bm_levine/liger/embedding' 9 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/bm_levine/fstmnn/metrics.csv' 'data/feature_deletion/bm_levine/sxc/orig' 'data/feature_deletion/bm_levine/fstmnn/embedding' 9 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/bm_levine/scan/metrics.csv' 'data/feature_deletion/bm_levine/sxc/orig' 'data/feature_deletion/bm_levine/scan/embedding' 9 &
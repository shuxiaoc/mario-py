/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/sxc_5k_20k-v3/metrics.csv' 'data/feature_deletion/murine/sxc_5k_20k-v3/orig' 'data/feature_deletion/murine/sxc_5k_20k-v3/embedding' 17 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/seurat/metrics.csv' 'data/feature_deletion/murine/sxc_5k_20k-v3/orig' 'data/feature_deletion/murine/seurat/embedding' 17 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/liger/metrics.csv' 'data/feature_deletion/murine/sxc_5k_20k-v3/orig' 'data/feature_deletion/murine/liger/embedding' 17 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/fstmnn/metrics.csv' 'data/feature_deletion/murine/sxc_5k_20k-v3/orig' 'data/feature_deletion/murine/fstmnn/embedding' 17 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/scan/metrics.csv' 'data/feature_deletion/murine/sxc_5k_20k-v3/orig' 'data/feature_deletion/murine/scan/embedding' 17 &

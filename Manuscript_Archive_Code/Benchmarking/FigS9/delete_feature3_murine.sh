/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/sxc/metrics.csv' 'data/feature_deletion/murine/sxc/orig' 'data/feature_deletion/murine/sxc/embedding' 16 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/seurat/metrics.csv' 'data/feature_deletion/murine/sxc/orig' 'data/feature_deletion/murine/seurat/embedding' 16 &
wait
/usr/bin/Rscript calculate_metricsVer2.R 'data/feature_deletion/murine/liger/metrics.csv' 'data/feature_deletion/murine/liger/orig' 'data/feature_deletion/murine/liger/embedding' 16 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/fstmnn/metrics.csv' 'data/feature_deletion/murine/sxc/orig' 'data/feature_deletion/murine/fstmnn/embedding' 16 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/murine/scan/metrics.csv' 'data/feature_deletion/murine/sxc/orig' 'data/feature_deletion/murine/scan/embedding' 16 &

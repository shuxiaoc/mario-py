/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/felix/sxc/metrics.csv' 'data/feature_deletion/felix/sxc/orig' 'data/feature_deletion/felix/sxc/embedding' 9 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/felix/seurat/metrics.csv' 'data/feature_deletion/felix/sxc/orig' 'data/feature_deletion/felix/seurat/embedding' 9 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/felix/liger/metrics.csv' 'data/feature_deletion/felix/sxc/orig2' 'data/feature_deletion/felix/liger/embedding' 9 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/felix/fstmnn/metrics.csv' 'data/feature_deletion/felix/sxc/orig' 'data/feature_deletion/felix/fstmnn/embedding' 9 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/felix/scan/metrics.csv' 'data/feature_deletion/felix/sxc/orig' 'data/feature_deletion/felix/scan/embedding' 9 &

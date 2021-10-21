/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab/sxc/metrics.csv' 'data/feature_deletion/human_crab/sxc/orig' 'data/feature_deletion/human_crab/sxc/embedding' 13 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab/seurat/metrics.csv' 'data/feature_deletion/human_crab/sxc/orig' 'data/feature_deletion/human_crab/seurat/embedding' 13 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab/liger/metrics.csv' 'data/feature_deletion/human_crab/sxc/orig' 'data/feature_deletion/human_crab/liger/embedding' 13 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab/fstmnn/metrics.csv' 'data/feature_deletion/human_crab/sxc/orig' 'data/feature_deletion/human_crab/fstmnn/embedding' 13 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab/scan/metrics.csv' 'data/feature_deletion/human_crab/sxc/orig' 'data/feature_deletion/human_crab/scan/embedding' 13 &

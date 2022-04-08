/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab_il4/sxc/metrics.csv' 'data/feature_deletion/human_crab_il4/sxc/orig' 'data/feature_deletion/human_crab_il4/sxc/embedding' 13 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab_il4/seurat/metrics.csv' 'data/feature_deletion/human_crab_il4/sxc/orig' 'data/feature_deletion/human_crab_il4/seurat/embedding' 13 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab_il4/liger/metrics.csv' 'data/feature_deletion/human_crab_il4/sxc/orig' 'data/feature_deletion/human_crab_il4/liger/embedding' 13 &
wait
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab_il4/fstmnn/metrics.csv' 'data/feature_deletion/human_crab_il4/sxc/orig' 'data/feature_deletion/human_crab_il4/fstmnn/embedding' 13 &
/usr/bin/Rscript calculate_metrics.R 'data/feature_deletion/human_crab_il4/scan/metrics.csv' 'data/feature_deletion/human_crab_il4/sxc/orig' 'data/feature_deletion/human_crab_il4/scan/embedding' 13 &

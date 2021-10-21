## script to make plots realted to benchmark

################################################## 
# currently this is for producting feature delete
##################################################

sxc_path = "../Benchmarking/data/feature_deletion/human_crab/sxc/metrics.csv"
seurat_path = "../Benchmarking/data/feature_deletion/human_crab/seurat/metrics.csv"
liger_path = "../Benchmarking/data/feature_deletion/human_crab/liger/metrics.csv"
fstmnn_path = "../Benchmarking/data/feature_deletion/human_crab/fstmnn/metrics.csv"
scan_path = "../Benchmarking/data/feature_deletion/human_crab/scan/metrics.csv"

# read metrics
a = read.csv(sxc_path)
colnames(a)[1] = "n_deleted"
colnames(a)[7] = "X1v1_acc"
b = read.csv(seurat_path)
c = read.csv(liger_path)
d = read.csv(fstmnn_path)
e = read.csv(scan_path)
colnames(e)[1] = "n_deleted"

# add method column
a$method = "sxc"
b$method = "seurat"
c$method = "liger"
d$method = "fstmnn"
e$method = "scan"

# columns to use
use = c("n_deleted","X1v1_acc","prop_remain","method",
        "sam_x","sam_y","slt_mix","slt_clust","slt_f1",
        "ari_mix","ari_clust","ari_f1","lisi_mix","lisi_clust","avg_mix")

############ match related
match = rbind(a[,use], b[,use])
match = rbind(match, c[,use])
match = rbind(match, d[,use])
match = rbind(match, e[,use])

# accu
p = ggplot(data=match, aes(x=n_deleted, y=X1v1_acc, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("accu") + theme(text = element_text(size=20)) +ylim(c(0.7,1))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/1.svg", plot=p, width=6, height=3)
#ggsave()
# prop remain
p = ggplot(data=match, aes(x=n_deleted, y=prop_remain, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic()+ ggtitle("remain")+ theme(text = element_text(size=20))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/2.svg", plot=p, width=6, height=3)

############ metrics related
match = match %>% mutate(sam_avg = (sam_x + sam_y)/2 ) 

## avg_sam
###
p = ggplot(data=match, aes(x=n_deleted, y=sam_avg, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("sam_avg") + theme(text = element_text(size=20)) +ylim(c(0.3,0.9))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/3.svg", plot=p, width=6, height=3)

## slt_mix
p = ggplot(data=match, aes(x=n_deleted, y=slt_mix, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("slt_mix") + theme(text = element_text(size=20)) +ylim(c(0.3,0.6))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/4.svg", plot=p, width=6, height=3)

## slt_clust
p = ggplot(data=match, aes(x=n_deleted, y=slt_clust, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("slt_clust") + theme(text = element_text(size=20)) + ylim(c(0.3,0.9))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/5.svg", plot=p, width=6, height=3)

## slt_f1
#####
p = ggplot(data=match, aes(x=n_deleted, y=slt_f1, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("slt_f1") + theme(text = element_text(size=20)) + ylim(c(0.45,0.55))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/6.svg", plot=p, width=6, height=3)
## ari_mix
p = ggplot(data=match, aes(x=n_deleted, y=ari_mix, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("ari_mix") + theme(text = element_text(size=20))+ ylim(c(0.4,0.6))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/7.svg", plot=p, width=6, height=3)
## ari_mix
p = ggplot(data=match, aes(x=n_deleted, y=ari_clust, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("ari_clust") + theme(text = element_text(size=20))+ ylim(c(0.6,1))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/8.svg", plot=p, width=6, height=3)
## ari_f1
####
p = ggplot(data=match, aes(x=n_deleted, y=ari_f1, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("ari_f1") + theme(text = element_text(size=20))+ ylim(c(0.5,0.65))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/9.svg", plot=p, width=6, height=3)
## lisi_mix
p = ggplot(data=match, aes(x=n_deleted, y=lisi_mix, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("lisi_mix") + theme(text = element_text(size=20))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/10.svg", plot=p, width=6, height=3)
## lisi_clust
p = ggplot(data=match, aes(x=n_deleted, y=lisi_clust, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("lisi_clust") + theme(text = element_text(size=20))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/11.svg", plot=p, width=6, height=3)
## avg_mix
#####
p = ggplot(data=match, aes(x=n_deleted, y=avg_mix, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("avg_mix") + theme(text = element_text(size=20)) + ylim(c(0.075,0.14))
ggsave(file="../Benchmarking/data/feature_deletion/human_crab/figures/12.svg", plot=p, width=6, height=3)

#################################
# now start the add noise 
#################################

sxc_path = "../Benchmarking/data/noise_addition/humancrab/sxc/metrics.csv"
seurat_path = "../Benchmarking/data/noise_addition/humancrab/seurat/metrics.csv"
liger_path = "../Benchmarking/data/noise_addition/humancrab/liger/metrics.csv"
fstmnn_path = "../Benchmarking/data/noise_addition/humancrab/fstmnn/metrics.csv"
scan_path = "../Benchmarking/data/noise_addition/humancrab/scan/metrics.csv"
# read metrics
a = read.csv(sxc_path)
colnames(a)[1] = "n_deleted"
colnames(a)[7] = "X1v1_acc"
## if avg_pval <0.05, stop match
a$X1v1_acc[c(6:9)] = NA
a$prop_remain[c(6:9)] = NA
b = read.csv(seurat_path)
colnames(b)[1] = "n_deleted"
c = read.csv(liger_path)
colnames(c)[1] = "n_deleted"
d = read.csv(fstmnn_path)
colnames(d)[1] = "n_deleted"
e = read.csv(scan_path)
colnames(e)[1] = "n_deleted"
# add method column
a$method = "sxc"
a$noise = c(0,0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
b$method = "seurat"
b$noise = c(0,0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
c$method = "liger"
c$noise = c(0,0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
d$method = "fstmnn"
d$noise = c(0,0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
e$method = "scan"
e$noise = c(0,0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
# columns to use
use = c("n_deleted","X1v1_acc","prop_remain","method","noise")
# match
match = rbind(a[,use], b[,use])
match = rbind(match, c[,use])
match = rbind(match, d[,use])
match = rbind(match, e[,use])
# accu
p = ggplot(data=match, aes(x=noise, y=X1v1_acc, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic() + ggtitle("accu") + theme(text = element_text(size=20)) +ylim(c(0.5,1))
ggsave(file="../Benchmarking/data/noise_addition/humancrab/humancrab_figs/15.svg", plot=p, width=6, height=3)
# prop remain
p = ggplot(data=match, aes(x=noise, y=prop_remain, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method)) + theme_classic()+ ggtitle("remain")+ theme(text = element_text(size=20))+ylim(c(0.6,1))
ggsave(file="../Benchmarking/data/noise_addition/humancrab/humancrab_figs/16.svg", plot=p, width=6, height=3)




#################################
# now start cell type deletion
#################################

sxc_path0.2 = "../Benchmarking/data/delete_cell_types/human_crab/sxc/metrics.csv"
seurat_path = "../Benchmarking/data/delete_cell_types/human_crab/seurat/metrics.csv"
fstmnn_path = "../Benchmarking/data/delete_cell_types/human_crab/fstmnn/metrics.csv"
scan_path = "../Benchmarking/data/delete_cell_types/human_crab/scan/metrics.csv"
# read metrics
a = read.csv(sxc_path0.2)
colnames(a)[1] = "n_deleted"
colnames(a)[7] = "X1v1_acc"
colnames(a)[9] = "celltype_error"
b = read.csv(seurat_path)
#c = read.csv(liger_path)
d = read.csv(fstmnn_path)
e = read.csv(scan_path)
colnames(e)[1] = "n_deleted"
# add method column
a$method = "sxc"
a2$method = "sxc2"
b$method = "seurat"
#c$method = "liger"
d$method = "fstmnn"
e$method = "scan"
# columns to use
use = c("n_deleted","X1v1_acc","prop_remain","method","celltype_error")
# match
match = rbind(a[,use], b[,use])
match = rbind(match, d[,use])
match = rbind(match, e[,use])
match[1,"celltype_error"] = "NA"
match = subset(match,match$n_deleted!=0)
match$n_deleted = as.factor(match$n_deleted)
match$celltype_error = as.numeric(match$celltype_error)
match$celltype_error = sqrt(1/(match$celltype_error))
#### figure making start

p <- ggplot(data=match, aes(x=n_deleted, y=celltype_error, fill=method)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  theme_classic()+ ggtitle("All")+ theme(text = element_text(size=20)) +ylab("error_avoidance_score")
ggsave(file="../Benchmarking/data/delete_cell_types/crab_figs/15.svg", plot=p, width=6, height=3)



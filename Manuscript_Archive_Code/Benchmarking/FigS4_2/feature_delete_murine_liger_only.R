path0 = "/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/feature_deletion/murine/sxc_5k_20k-v3/"
path = "/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/feature_deletion/murine/"
##### liger start
library(rliger)
## read in 
#path0 = "/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0616/data/feature_deletion/bmcite_levine/"
#path = "/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0616/data/feature_deletion/"
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
feature_2_delete = c("CD1632","F480","Ly6G","Ly6C","CD2135","CD11b","CD11c","CD5","CD19","IgM","IgD","CD27","TCR","NKp46","B220","CD3","CD4")
path_out = paste0(path,"liger/")
feature_2_delete = c("Full",feature_2_delete)
## container
coluse = c("n_deleted","1v1_acc","transfer_acc","prop_remain")
container = matrix(NA, nrow = length(feature_2_delete), ncol = length(coluse))
container = as.data.frame(container)
colnames(container) = coluse
## label
x_label=x[,"label"]
y_label=y[,"label"]
x [,"label"] <- NULL
y [,"label"] <- NULL
##
idx = 0
for (item in feature_2_delete){
  #print(item)
  idx = idx + 1
  if (item != "Full"){ # drop the current feature
    x [,item] <- NULL
    y [,item] <- NULL
  }
  ###### the matching and integration process
  x_value=x
  y_value=y
  #x_value_scale=scale(x_value) # scaled
  #y_value_scale=scale(y_value) # scaled
  x_value_scale=x_value # not scaled
  y_value_scale=y_value # not scaled
  x_value_scale=as.data.frame(x_value_scale)
  rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
  y_value_scale=as.data.frame(y_value_scale)
  rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
  # create seurat object datax
  ligerobj=createLiger( list(x = t(x_value_scale), y = t(y_value_scale) ) )
  ###Start integration, feature will drop so change this
  features=intersect(colnames(x),colnames(y))
  # default preprocessing
  ligerobj <- rliger::normalize(ligerobj)
  #ifnb_liger <- selectGenes(ifnb_liger, var.thresh = 0, alpha.thresh=1)
  ligerobj@var.genes=features #  just use all
  ligerobj <- scaleNotCenter(ligerobj)
  ligerobj <- optimizeALS(ligerobj, k = (length(features)-1))
  ligerobj <- quantile_norm(ligerobj)
  # remain porp
  container[idx,"prop_remain"] = NA
  # 1v1 match
  container[idx,"1v1_acc"] = NA
  # knn label trans
  container[idx,"transfer_acc"] = NA
  # embedding
  embedding = ligerobj@H.norm # pca scores
  container[idx,"n_deleted"] = idx -1 # just a useless column
  ###### end of the matching and integration process for current feature drop
  name_1 = paste0("embedding_x",idx-1,".csv")
  name_2 = paste0("embedding_y",idx-1,".csv")
  write.csv(embedding[c(1:5000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  write.csv(embedding[c(5001:24980),], paste0(path_out,name_2), row.names=FALSE) # need to decide output pca cell
  # some errors causing cell missing after reduction, possibily due to bad integration
}
write.csv(container,"/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/feature_deletion/murine/liger/metrics.csv")
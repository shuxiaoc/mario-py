## script for feature deletion after sxc python script
## computes seurat liger fstmnn

#### seurat
## read in 
path0 = "../Benchmarking/data/feature_deletion/human_crab/sxc/"
path = "../Benchmarking/data/feature_deletion/human_crab/"
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
feature_2_delete = c('BDCA3',
        'CCR7',
        'CD11b',
        'CD11c',
        'CD123',
        'CD14',
        'CD16',
        'CD161',
        'CD19',
        'CD1c',
        'CD235a',
        'CD3',
        'CD33')
path_out = paste0(path,"seurat/")
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
  print(item)
  idx = idx + 1
  if (item != "Full"){ # drop the current feature
    x [,item] <- NULL
    y [,item] <- NULL
  }
  ###### the matching and integration process
  library(Seurat)
  x_value=x
  y_value=y
  x_value_scale=scale(x_value) # scaled
  y_value_scale=scale(y_value) # scaled
  x_value_scale=as.data.frame(x_value_scale)
  rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
  y_value_scale=as.data.frame(y_value_scale)
  rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
  # create seurat object datax
  x_obj=CreateSeuratObject(counts=t(x_value_scale),assay="x")
  SetAssayData(object = x_obj, slot = "data", new.data = t(x_value_scale), assay="x")
  x_obj = ScaleData(x_obj)
  # add suerat object datay
  y_obj=CreateSeuratObject(counts=t(y_value_scale),assay="y")
  SetAssayData(object = y_obj, slot = "data", new.data = t(y_value_scale), assay="y")
  y_obj = ScaleData(y_obj)
  list_modality=list(x_obj,y_obj)
  ###Start integration, feature will drop so change this
  features=intersect(colnames(x),colnames(y))
  x_obj <- RunPCA(x_obj, npcs = (length(features)-1), verbose = FALSE, features = features)
  y_obj <- RunPCA(y_obj, npcs = (length(features)-1), verbose = FALSE, features = features)
  #### anchors for 1v1 match
  pre.anchors <- FindTransferAnchors(reference = y_obj, query = x_obj, 
                                     dims = 1:(length(features)-1), features = features, npcs = (length(features)-1))
  predictions <- TransferData(anchorset = pre.anchors, refdata = y_label, 
                              dims = 1:(length(features)-1))
  #### anchors for integration pca
  Int.anchors <- FindIntegrationAnchors(object.list = list_modality,
                                        dims = 1:(length(features)-1), anchor.features =features, k.filter = 10)
  xy_int <- IntegrateData(anchorset = Int.anchors, dims = 1:(length(features)-1), k.weight = 10)
  ### get embedding
  DefaultAssay(xy_int) <- "integrated"
  xy_int <- ScaleData(xy_int, verbose = FALSE)
  xy_int <- RunPCA(xy_int, npcs = (length(features)-1), verbose = FALSE) # need to change this to 
  matchmm=as.data.frame(pre.anchors@anchors)
  colnames(matchmm)=c("dY","dX","score")
  matchmm$dXlabel=x_label[matchmm$dX]
  matchmm$dYlabel=y_label[matchmm$dY]
  # remain porp
  rem_prop = length(unique(matchmm$dX)) / dim(x)[1] # 50%
  container[idx,"prop_remain"] = rem_prop
  # 1v1 match
  ovoM = sum(as.character(matchmm$dXlabel)==as.character(matchmm$dYlabel))/dim(matchmm)[1] # 0.9
  container[idx,"1v1_acc"] = ovoM
  # knn label trans
  knnM = sum(predictions$predicted.id==x_label)/dim(x)[1] # 0.94
  container[idx,"transfer_acc"] = knnM
  # embedding
  embedding = xy_int@reductions$pca@cell.embeddings # pca scores
  container[idx,"n_deleted"] = idx -1 # just a useless column
  ###### end of the matching and integration process for current feature drop
  name_1 = paste0("embedding_x",idx-1,".csv")
  name_2 = paste0("embedding_y",idx-1,".csv")
  write.csv(embedding[c(1:8000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  write.csv(embedding[c(8001:40000),], paste0(path_out,name_2), row.names=FALSE) # need to decide output pca cell
}
write.csv(container,"../Benchmarking/data/feature_deletion/human_crab/seurat/metrics.csv")

##### liger start
library(rliger)
## read in 
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
feature_2_delete = c('BDCA3',
        'CCR7',
        'CD11b',
        'CD11c',
        'CD123',
        'CD14',
        'CD16',
        'CD161',
        'CD19',
        'CD1c',
        'CD235a',
        'CD3',
        'CD33')
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
  write.csv(embedding[c(1:8000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  write.csv(embedding[c(8001:40000),], paste0(path_out,name_2), row.names=FALSE) # need to decide output pca cell
}
write.csv(container,"../Benchmarking/data/feature_deletion/human_crab/liger/metrics.csv")
###### fstmnn

library(batchelor)
## read in 
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
feature_2_delete = c('BDCA3',
        'CCR7',
        'CD11b',
        'CD11c',
        'CD123',
        'CD14',
        'CD16',
        'CD161',
        'CD19',
        'CD1c',
        'CD235a',
        'CD3',
        'CD33')
path_out = paste0(path,"fstmnn/")
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
  print(item)
  idx = idx + 1
  if (item != "Full"){ # drop the current feature
    x [,item] <- NULL
    y [,item] <- NULL
  }
  ###### the matching and integration process
  x_value=x
  y_value=y
  x_value_scale=scale(x_value) # scaled
  y_value_scale=scale(y_value) # scaled
  #x_value_scale=x_value # not scaled
  #y_value_scale=y_value # not scaled
  x_value_scale=as.data.frame(x_value_scale)
  rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
  y_value_scale=as.data.frame(y_value_scale)
  rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
  # fstmnn integration
  features=intersect(colnames(x),colnames(y))
  fmnn <- fastMNN(t(x_value_scale[features]), t(y_value_scale[features]),d=(length(features) -1 ))
  # fstmnn matching
  mnnMatch=findMutualNN(x_value_scale[features],y_value_scale[features], k1=20)
  mnnMatch=as.data.frame(mnnMatch)
  colnames(mnnMatch)=c("dx","dy")
  mnnMatch$dXlabel=x_label[mnnMatch$dx]
  mnnMatch$dYlabel=y_label[mnnMatch$dy]
  #print(dim(mnnMatch))
  # remain porp
  rem_prop = length(unique(mnnMatch$dx)) / dim(x)[1]
  container[idx,"prop_remain"] = rem_prop
  # 1v1 match
  ovoM = sum(as.character(mnnMatch$dXlabel)==as.character(mnnMatch$dYlabel))/dim(mnnMatch)[1]
  container[idx,"1v1_acc"] = ovoM
  # knn label trans
  container[idx,"transfer_acc"] = NA
  # embedding
  embedding = fmnn@int_colData$reducedDims@listData$corrected
  container[idx,"n_deleted"] = idx -1 # just a useless column
  ###### end of the matching and integration process for current feature drop
  name_1 = paste0("embedding_x",idx-1,".csv")
  name_2 = paste0("embedding_y",idx-1,".csv")
  write.csv(embedding[c(1:8000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  write.csv(embedding[c(8001:40000),], paste0(path_out,name_2), row.names=FALSE) # need to decide output pca cell
}
write.csv(container,"../Benchmarking/data/feature_deletion/human_crab/fstmnn/metrics.csv")

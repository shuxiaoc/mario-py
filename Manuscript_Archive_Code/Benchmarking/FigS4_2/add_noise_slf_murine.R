# script to create low quality data
# embedding not yet saved

########### seurat start
## read in 
path0 = "../Benchmarking/data/noise_addition/murine/sxc/"
path = "../Benchmarking/data/noise_addition/murine/"
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
path_out = paste0(path,"seurat/")
## container
noise_level = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
noise_level = c( 0 ,noise_level)
coluse = c("noise_idx","1v1_acc","transfer_acc","prop_remain")
container = matrix(NA, nrow = length(noise_level), ncol = length(coluse))
container = as.data.frame(container)
colnames(container) = coluse
## label
x_label=x[,"label"]
y_label=y[,"label"]
x [,"label"] <- NULL
y [,"label"] <- NULL
##
idx = 0
for (item in noise_level){
  print(item)
  idx = idx + 1
  ###### the matching and integration process
  library(Seurat)
  x_value=x
  y_value=y
  x_value_scale=scale(x_value) # scaled
  y_value_scale=scale(y_value) # scaled
  # dj make some noise
  noise_value_x = matrix( rnorm(dim(x_value)[1]*dim(x_value)[2],mean=0,sd=item), dim(x_value)[1], dim(x_value)[2]) 
  noise_value_y = matrix( rnorm(dim(y_value)[1]*dim(y_value)[2],mean=0,sd=item), dim(y_value)[1], dim(y_value)[2]) 
  # add noise
  x_value_scale = x_value_scale + noise_value_x
  y_value_scale = y_value_scale + noise_value_y
  # start
  x_value_scale=as.data.frame(x_value_scale)
  rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
  y_value_scale=as.data.frame(y_value_scale)
  rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
  # create seurat object datax
  x_obj=CreateSeuratObject(counts=t(x_value_scale),assay="x")
  SetAssayData(object = x_obj, slot = "data", new.data = t(x_value_scale), assay="x")
  # add suerat object datay
  y_obj=CreateSeuratObject(counts=t(y_value_scale),assay="y")
  SetAssayData(object = y_obj, slot = "data", new.data = t(y_value_scale), assay="y")
  list_modality=list(x_obj,y_obj)
  ###Start integration, feature will drop so change this
  features=intersect(colnames(x),colnames(y))
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
  ###### end of the matching and integration process for current feature drop
  name_1 = paste0("embedding_x",idx-1,".csv")
  name_2 = paste0("embedding_y",idx-1,".csv")
  #write.csv(embedding[c(1:7000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  #write.csv(embedding[c(7001:17000),], paste0(path_out,name_2), row.names=FALSE) # need to decide output pca cell
}
write.csv(container,"../Benchmarking/data/noise_addition/murine/seurat/metrics.csv")

############# fstmnn start
## read in 
#path = "/home/bkzhu/SNE-multi/figure_rcode/sup_simulation/code_0713/data/noise_addition/"
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
path_out = paste0(path,"fstmnn/")
## container
noise_level = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
noise_level = c( 0 ,noise_level)
coluse = c("noise_idx","1v1_acc","transfer_acc","prop_remain")
container = matrix(NA, nrow = length(noise_level), ncol = length(coluse))
container = as.data.frame(container)
colnames(container) = coluse
## label
x_label=x[,"label"]
y_label=y[,"label"]
x [,"label"] <- NULL
y [,"label"] <- NULL
##
idx = 0
for (item in noise_level){
  print(item)
  idx = idx + 1
  ###### the matching and integration process
  library(batchelor)
  x_value=x
  y_value=y
  x_value_scale=scale(x_value) # scaled
  y_value_scale=scale(y_value) # scaled
  # dj make some noise
  noise_value_x = matrix( rnorm(dim(x_value)[1]*dim(x_value)[2],mean=0,sd=item), dim(x_value)[1], dim(x_value)[2]) 
  noise_value_y = matrix( rnorm(dim(y_value)[1]*dim(y_value)[2],mean=0,sd=item), dim(y_value)[1], dim(y_value)[2]) 
  # add noise
  x_value_scale = x_value_scale + noise_value_x
  y_value_scale = y_value_scale + noise_value_y
  # start
  x_value_scale=as.data.frame(x_value_scale)
  rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
  y_value_scale=as.data.frame(y_value_scale)
  rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
  # fastmnn
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
  ovoM = sum(as.character(mnnMatch$dXlabel)== as.character(mnnMatch$dYlabel))/dim(mnnMatch)[1]
  container[idx,"1v1_acc"] = ovoM
  # knn label trans
  container[idx,"transfer_acc"] = NA
  # embedding
  embedding = fmnn@int_colData$reducedDims@listData$corrected
  container[idx,"n_deleted"] = idx -1 # just a useless column
  ###### end of the matching and integration process for current feature drop
  name_1 = paste0("embedding_x",idx-1,".csv")
  name_2 = paste0("embedding_y",idx-1,".csv")
  #write.csv(embedding[c(1:7000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  #write.csv(embedding[c(7001:17000),], paste0(path_out,name_2), row.names=FALSE) # need to decide output pca cell
}
write.csv(container,"../Benchmarking/data/noise_addition/murine/fstmnn/metrics.csv")

### liger
## read in 
x = read.csv(paste0(path0,"orig_x.csv"))
y = read.csv(paste0(path0,"orig_y.csv"))
## some parameters
path_out = paste0(path,"liger/")
## container
noise_level = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
noise_level = c( 0 ,noise_level)
coluse = c("noise_idx","1v1_acc","transfer_acc","prop_remain")
container = matrix(NA, nrow = length(noise_level), ncol = length(coluse))
container = as.data.frame(container)
colnames(container) = coluse
## label
x_label=x[,"label"]
y_label=y[,"label"]
x [,"label"] <- NULL
y [,"label"] <- NULL
##
idx = 0
for (item in noise_level){
  print(item)
  idx = idx + 1
  ###### the matching and integration process
  library(rliger)
  x_value=x
  y_value=y
  x_value_scale=scale(x_value) # scaled
  y_value_scale=scale(y_value) # scaled
  # dj make some noise
  noise_value_x = matrix( rnorm(dim(x_value)[1]*dim(x_value)[2],mean=0,sd=item), dim(x_value)[1], dim(x_value)[2]) 
  noise_value_y = matrix( rnorm(dim(y_value)[1]*dim(y_value)[2],mean=0,sd=item), dim(y_value)[1], dim(y_value)[2]) 
  # add noise
  x_value_scale = x_value_scale + noise_value_x
  y_value_scale = y_value_scale + noise_value_y
  #x_value_scale=x_value # not scaled
  #y_value_scale=y_value # not scaled
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
  #write.csv(embedding[c(1:7000),], paste0(path_out,name_1), row.names=FALSE) # need to decide output pca cell
  #write.csv(embedding[c(7001:17000),], paste0(path_out,name_2), row.names=FALSE)
}
write.csv(container,"../Benchmarking/data/noise_addition/murine/liger/metrics.csv")

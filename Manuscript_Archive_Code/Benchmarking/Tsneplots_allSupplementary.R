## bokai zhu
## 0830
## script related to produce the tsne plots presented in FiguresS2.2;S2.4;S3.1;S3.3;S4.2



############################### figureS2.2 ##############################################

# read in the MARIO files, all other method used same cells as MARIO have produced for fair comparison.
# (though MARIO have helped filtered out bad cell pairs for other methods already)

cite_cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/bmcite_cca10.csv")
lv32_cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/lv32_cca10.csv")

cite_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/bmcite_mario_matched.csv")
lv32_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/lv32_mario_matched.csv")

## first produce the preintegration tsne:
cite_value_scale=scale(cite_matched[,c(3:27)]) # only the protein features
lv_value_scale=scale(lv32_matched[,c(3:34)]) # only the protein features
shared=intersect(colnames(cite_value_scale),colnames(lv_value_scale))
cite_lv=rbind(cite_value_scale[,shared],lv_value_scale[,shared])
set.seed(42)
all_tsne_scale_10=Rtsne(cite_lv, check_duplicates = FALSE, num_threads = 10) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("cite",26822), rep("cytof",26822)))# lazy input of modality
p = ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with all shared pre")

## produce tsne for mario result
cite_cca=cite_cca[,-1] # first column is rowname
lv_cca=lv32_cca[,-1]
all_cca_scale_10=scale(rbind(cite_cca,lv_cca))
set.seed(42)
all_tsne_sxc10=Rtsne(all_cca_scale_10, check_duplicates = FALSE, num_threads = 10) # run with cca
all_tsne_sxc10_data <- data.frame(x = all_tsne_sxc10$Y[,1], y = all_tsne_sxc10$Y[,2])
label=as.factor(c(rep("cite",26822), rep("cytof",26822)))
p = ggplot(all_tsne_sxc10_data)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with mario")


## then produce seurat tsne
library(Seurat)
all=rbind(cite_value_scale[,shared],lv_value_scale[,shared]) # use scaled each individual dataset, performs better for seurat
all_t=t(all) #formatting for seurat
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
seuratobj=CreateSeuratObject(counts=all_t)
SetAssayData(object = seuratobj, slot = "data", new.data = all_t)
SetAssayData(object = seuratobj, slot = "scale.data", new.data = all_t)
groups=c(rep("cite",26822),rep("seurat",26822))
seuratobj <- AddMetaData(object = seuratobj, metadata = groups, col.name = "group")
seuratobj <- SplitObject(seuratobj, split.by = "group")
features=shared
seuratobj.anchors <- FindIntegrationAnchors(object.list = seuratobj, dims = 1:11,max.features =11, anchor.features =shared, k.filter = 10)
seuratobj.combined <- IntegrateData(anchorset = seuratobj.anchors ,dims = 1:11, k.weight = 10) # note NA will be produced in some cells
temp = seuratobj.combined@assays$integrated@data
temp[is.na(temp)] <- 0 # na value change to 0 for visualization

DefaultAssay(seuratobj.combined) <- "integrated"
seuratobj.combined=SetAssayData(object = seuratobj.combined, slot = "data", new.data =temp , assay="integrated") # re-input data
seuratobj.combined <- ScaleData(seuratobj.combined, verbose = FALSE) # scale the seurat corrected
seuratobj.combined <- RunPCA(seuratobj.combined, npcs = 30, verbose = FALSE)
seuratpca=seuratobj.combined@reductions$pca@cell.embeddings
set.seed(42)
all_tsne_seurat=Rtsne(seuratpca, check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_seurat_data <- data.frame(x = all_tsne_seurat$Y[,1], y = all_tsne_seurat$Y[,2])
label=as.factor(c(rep("cite",26822),rep("seurat",26822)))
p = ggplot(all_tsne_seurat_data)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with all shared seurat")

## produce liger tsne
library(rliger)
# should not use scaled for liger, use liger function to scale
all=rbind(cite_matched[,shared],lv32_matched[,shared]) # input non-scaled shared
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
all_t1=all_t[,c(1:26822)] # cite
all_t2=all_t[,c(26823:53644)] # cytof
ligerobj <- createLiger(list(cite =all_t1, lv32 = all_t2))
ligerobj <- rliger::normalize(ligerobj)
ligerobj@var.genes=shared #  just use all
ligerobj <- scaleNotCenter(ligerobj)
ligerobj <- optimizeALS(ligerobj, k=11)
ligerobj <- quantile_norm(ligerobj)
ligerpca=ligerobj@H.norm
set.seed(42)
all_tsne_liger=Rtsne(ligerpca[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_liger_data_all <- data.frame(x = all_tsne_liger$Y[,1], y = all_tsne_liger$Y[,2])
label=as.factor(c(rep("cite",26822),rep("seurat",26822)))
p = ggplot(all_tsne_liger_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("lr with 10cca") +theme_classic()
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/fig3/tsne_figs/liger_integration-update.pdf", plot=p, width=11, height=10)

## produce fastmnn tsne plot
#  for fastmnn scale for them
all=rbind(cite_value_scale[,shared],lv_value_scale[,shared]) # use scaled each individual dataset, performs better for fstmnn
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
all_t1=all_t[,c(1:26822)] # cite
all_t2=all_t[,c(26823:53644)] # cytof
## fstmnn
fmnnobj <- fastMNN(all_t1,all_t2)
fmmn_reduce = fmnnobj@int_colData$reducedDims@listData$corrected
set.seed(42)
all_tsne_fst_all=Rtsne(fmmn_reduce[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_fst_data_all <- data.frame(x = all_tsne_fst_all$Y[,1], y = all_tsne_fst_all$Y[,2])
label=as.factor(c(rep("cite",26822),rep("seurat",26822)))
p = ggplot(all_tsne_fst_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("lr with 10cca") + theme_classic()
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/fig3/tsne_figs/fstmnn_integration-update.pdf", plot=p, width=11, height=10)

## produce scanorama tsne plot
# scanorama run parameters please refer to the python script
# our r script just read in the results
cite_scan = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/bmcite_scan.csv", header = FALSE)
lv_scan = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/lv32_scan.csv", header = FALSE)
scanpy_pca=rbind(cite_scan,lv_scan)
set.seed(42)
all_tsne_scan=Rtsne(scanpy_pca[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_data_scan <- data.frame(x = all_tsne_scan$Y[,1], y = all_tsne_scan$Y[,2])
label=as.factor(c(rep("cite",26822), rep("cytof",26822)))
p = ggplot(all_tsne_data_scan)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with scannorama 10pca")
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/fig3/tsne_figs/scan_integration-update.pdf", plot=p, width=11, height=10)




############################### figureS2.4 ##############################################




## read in files from MARIO
tenx_cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/10xgenomic_cca10.csv")
felix_cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/felixetal_cca10.csv")

tenx_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/pbmc10x_mario_matched.csv")
felix_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/pbmcCytof_mario_matched.csv")

## MARIO tsne
all_cca_10=rbind(tenx_cca[,-1],felix_cca[,-1])
all_cca_10=all_cca_10[,c(1:10)]
all_cca_scale_10=scale(all_cca_10)
set.seed(42)
all_tsne_scale_10=Rtsne(all_cca_scale_10, check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("cite",4606), rep("cytof",4606)))
p = ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 15,stroke = 1, alpha =0.001), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")

## pre integration tsne plot
tenx_value_scale=scale(tenx_matched[,c(3:31)])
felix_value_scale=scale(felix_matched[,c(3:32)])
shared=intersect(colnames(tenx_value_scale),colnames(felix_value_scale))
tenx_felix=rbind(tenx_value_scale[,shared],felix_value_scale[,shared])
set.seed(42)
all_tsne_scale_10=Rtsne(tenx_felix, check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("cite",4606), rep("cytof",4606)))
p = ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 15,stroke = 1, alpha =0.001), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with all shared pre") + ylim(c(-50,50))

## seurat tsne plot
x_value_scale=tenx_value_scale # scaled
y_value_scale=felix_value_scale # scaled
x_value_scale=as.data.frame(x_value_scale)
rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
y_value_scale=as.data.frame(y_value_scale)
rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
# create seurat object datax
x_obj=CreateSeuratObject(counts=t(x_value_scale),assay="x")
SetAssayData(object = x_obj, slot = "data", new.data = t(x_value_scale), assay="x")
x_obj <- ScaleData(x_obj)
# add suerat object datay
y_obj=CreateSeuratObject(counts=t(y_value_scale),assay="y")
SetAssayData(object = y_obj, slot = "data", new.data = t(y_value_scale), assay="y")
y_obj <- ScaleData(y_obj)
list_modality=list(x_obj,y_obj)
features=intersect(colnames(x_value_scale),colnames(x_value_scale))
# anchors for integration pca
Int.anchors <- FindIntegrationAnchors(object.list = list_modality,
                                      dims = 1:10, anchor.features =features, k.filter = 10)
xy_int <- IntegrateData(anchorset = Int.anchors, dims = 1:10, k.weight = 10)
# in case NA produced, although in this case there were no NA produced
temp = xy_int@assays$integrated@data
temp[is.na(temp)] <- 0
# new obj in case NA issue
temp_obj = CreateSeuratObject(counts=as.matrix(temp),assay="temp")
SetAssayData(object = temp_obj, slot = "data", new.data = as.matrix(temp), assay="temp")
temp_obj <- ScaleData(temp_obj)
temp_obj <- RunPCA(temp_obj, npcs = 10, verbose = FALSE, features = features) # need to change this to 
embedding = temp_obj@reductions$pca@cell.embeddings
set.seed(41)
seurat_tsne_10=Rtsne(embedding, check_duplicates = FALSE, num_threads = 10) # run with cca
seurat_tsne_data_10 <- data.frame(x = seurat_tsne_10$Y[,1], y = seurat_tsne_10$Y[,2])
label=as.factor(c(rep("codex",4606), rep("cite",4606)))
p = ggplot(seurat_tsne_data_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 15,stroke = 1, alpha =0.001), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("seurat integration")


## rliger tsne plots
tenx_value = as.data.frame(tenx_matched[,c(3:31)]) # no outside scaling for liger, better performance
felix_value = as.data.frame(felix_matched[,c(3:32)])
rownames(tenx_value) = tenx_matched$X
rownames(felix_value) = felix_matched$X
x_value_scale = tenx_value
y_value_scale = felix_value
ligerobj=createLiger( list(x = t(x_value_scale[features]), y = t(y_value_scale[features]) ) )
features=intersect(colnames(x_value_scale),colnames(y_value_scale))
# default preprocessing
ligerobj <- rliger::normalize(ligerobj)
ligerobj@var.genes=features #  just use all
ligerobj <- scaleNotCenter(ligerobj)
ligerobj <- optimizeALS(ligerobj, k = 10)
ligerobj <- quantile_norm(ligerobj)
embedding = ligerobj@H.norm # pca scores

set.seed(41)
liger_tsne_10=Rtsne(embedding[,c(1:10)], check_duplicates = FALSE, num_threads = 10) # run with cca
liger_tsne_data_10 <- data.frame(x = liger_tsne_10$Y[,1], y = liger_tsne_10$Y[,2])
label=as.factor(c(rep("codex",4606), rep("cite",4606)))
p = ggplot(liger_tsne_data_10)  + 
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) + 
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("liger integration")

## fastmnn tsne plots
x_value_scale=tenx_value_scale # scaled
y_value_scale=felix_value_scale # scaled
fmnn <- fastMNN(t(x_value_scale[,features]), t(y_value_scale[,features]),d=20)
embedding = fmnn@int_colData$reducedDims@listData$corrected
set.seed(41)
fstmnn_tsne_10=Rtsne(embedding[,c(1:10)], check_duplicates = FALSE, num_threads = 10) # run with cca
fstmnn_tsne_data_10 <- data.frame(x = fstmnn_tsne_10$Y[,1], y = fstmnn_tsne_10$Y[,2])
label=as.factor(c(rep("codex",4606), rep("cite",4606)))
p = ggplot(fstmnn_tsne_data_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 15,stroke = 1, alpha =0.001), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("fstmnn integration")

## produce scanorama tsne plot
# scanorama run parameters please refer to the python script
# our r script just read in the results
cite_scan = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/10x_scan.csv", header = FALSE)
lv_scan = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/felix_scan.csv", header = FALSE)
scanpy_pca=rbind(cite_scan,lv_scan)
set.seed(42)
all_tsne_scan=Rtsne(scanpy_pca[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_data_scan <- data.frame(x = all_tsne_scan$Y[,1], y = all_tsne_scan$Y[,2])
label=as.factor(c(rep("cite",4606), rep("cytof",4606)))
p = ggplot(all_tsne_data_scan)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with scannorama 10pca")




############################### figureS3.1 ##############################################



# this from gca triology, which these cells are matched to all three datasets
wcct=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/wcctHuman_mario_matched_sub20k.csv")
zh=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/XhumanIFNG_mario_matched_sub20k.csv")
zr=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/XRhesusIFNG_mario_matched_sub20k.csv")
zc=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/XCynoIFNG_mario_matched_sub20k.csv")

trio_cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/gcca_sub20k.csv")

## tsne plot for MARIO
set.seed(42)
mario_tsne_scale_10=Rtsne(trio_cca[,-1], check_duplicates = FALSE, num_threads = 10) # run with cca
mario_tsne_data_scale_10 <- data.frame(x = mario_tsne_scale_10$Y[,1], y = mario_tsne_scale_10$Y[,2])
label=wzhrc$dataset
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(mario_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("MARIO integration with 10cca") + scale_color_manual(values =  clr )

## tsne plot for seurat
share=intersect(colnames(wcct),colnames(zh))
share=intersect(share,colnames(zr))
share=intersect(share,colnames(zc))
share=share[c(3:41)]
# use scale
all=rbind(scale(wcct[,share]),scale(zh[,share])) # scale each individual dataset
all=rbind(all,scale(zr[,share]))
all=rbind(all,scale(zc[,share]))
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
trio=CreateSeuratObject(counts=all_t,assay="RNA")
SetAssayData(object = trio, slot = "data", new.data = all_t, assay="RNA")
SetAssayData(object = trio, slot = "scale.data", new.data = all_t, assay="RNA")
groups=c(rep("wcct",20000),rep("zh",20000),rep("zr",20000),rep("zc",20000))
trio <- AddMetaData(object = trio, metadata = groups, col.name = "group")
obj.list <- SplitObject(trio, split.by = "group")
features=share
trio.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
trio.combined <- IntegrateData(anchorset = trio.anchors)
trio.combined <- ScaleData(trio.combined, verbose = FALSE)
trio.combined <- RunPCA(trio.combined, npcs = 30, verbose = FALSE)
trio_seuratpca=trio.combined@reductions$pca@cell.embeddings
set.seed(42)
all_tsne_trioseurat_all=Rtsne(trio_seuratpca, check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_trioseurat_data_all <- data.frame(x = all_tsne_trioseurat_all$Y[,1], y = all_tsne_trioseurat_all$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(all_tsne_trioseurat_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("seurat trio with 10cca") + scale_color_manual(values =  clr)+theme_classic()


## tsne plot for liger
# use not scaled
all=rbind(wcct[,share],zh[,share]) # no scale with liger use in package scaling
all=rbind(all,zr[,share])
all=rbind(all,zc[,share])
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
all_t1=all_t[,c(1:20000)]
all_t2=all_t[,c(20001:40000)]
all_t3=all_t[,c(40001:60000)]
all_t4=all_t[,c(60001:80000)]
trio_liger <- createLiger(list(wcct =all_t1, zh = all_t2, zr = all_t3, zc = all_t4 ))
trio_liger <- rliger::normalize(trio_liger)
trio_liger@var.genes=share #  just use all
trio_liger <- scaleNotCenter(trio_liger)
trio_liger <- optimizeALS(trio_liger, k=30)
trio_liger <- quantile_norm(trio_liger)
a=trio_liger@H.norm
set.seed(42)
all_tsne_triolr_all=Rtsne(a[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_triolr_data_all <- data.frame(x = all_tsne_triolr_all$Y[,1], y = all_tsne_triolr_all$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(all_tsne_triolr_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("lr trio with 10cca") + scale_color_manual(values =  clr)+theme_classic()


## tsne plot for fastmnn
all=rbind(scale(wcct[,share]),scale(zh[,share])) # scale each individual dataset
all=rbind(all,scale(zr[,share]))
all=rbind(all,scale(zc[,share]))
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
all_t1=all_t[,c(1:20000)]
all_t2=all_t[,c(20001:40000)]
all_t3=all_t[,c(40001:60000)]
all_t4=all_t[,c(60001:80000)]
fmnn_trio <- fastMNN(all_t1,all_t2,all_t3,all_t4)
fmmn_reduce_trio = fmnn_trio@int_colData$reducedDims@listData$corrected

## tsne plot for scanorama
sc1=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/wcct_scanpca.csv", header = FALSE)
sc2=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/zh_scanpca.csv", header = FALSE)
sc3=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/zr_scanpca.csv", header = FALSE)
sc4=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/zc_scanpca.csv", header = FALSE)
scall=rbind(sc1,sc2)
scall=rbind(scall,sc3)
scall=rbind(scall,sc4)
set.seed(42)
all_tsne_trioscan_all=Rtsne(scall[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_trioscan_data_all <- data.frame(x = all_tsne_trioscan_all$Y[,1], y = all_tsne_trioscan_all$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(all_tsne_trioscan_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("scanorama trio with 10cca") + scale_color_manual(values =  clr)+theme_classic()


## tsne plot for pre-integration dataset
wcct1_value_scale=as.data.frame(scale(wcct[,share])) # original data but scaled per dataset, as usual
zh_value_scale=as.data.frame(scale(zh[,share]))
zr_value_scale=as.data.frame(scale(zr[,share]))
zc_value_scale=as.data.frame(scale(zc[,share]))
wcct1_value_scale$dataset="wcct"
zh_value_scale$dataset="zh"
zr_value_scale$dataset="zr"
zc_value_scale$dataset="zc"
wzhrc=rbind(wcct1_value_scale,zh_value_scale,zr_value_scale,zc_value_scale)
set.seed(42)
pre_tsne_scale_10=Rtsne(wzhrc, check_duplicates = FALSE, num_threads = 10) # run with cca
pre_tsne_data_scale_10 <- data.frame(x = pre_tsne_scale_10$Y[,1], y = pre_tsne_scale_10$Y[,2])
label=wzhrc$dataset

clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(pre_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("pre integration") + scale_color_manual(values =  clr )



############################### figureS3.3 ##############################################


# this part of code should be near identical for the code related to figure S3.1, as the data format should be same
# this from gca triology, which these cells are matched to all three datasets
wcct=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/wcctHuman2_mario_matched_sub20k.csv")
zh=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/XhumanIL4_mario_matched_sub20k.csv")
zr=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/XRhesusIL4_mario_matched_sub20k.csv")
zc=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/XCynoIL4_mario_matched_sub20k.csv")

trio_cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/gcca_sub20k.csv")

## tsne plot for MARIO
set.seed(42)
mario_tsne_scale_10=Rtsne(trio_cca[,-1], check_duplicates = FALSE, num_threads = 10) # run with cca
mario_tsne_data_scale_10 <- data.frame(x = mario_tsne_scale_10$Y[,1], y = mario_tsne_scale_10$Y[,2])
label=wzhrc$dataset
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(mario_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("MARIO integration with 10cca") + scale_color_manual(values =  clr )

## tsne plot for seurat
share=intersect(colnames(wcct),colnames(zh))
share=intersect(share,colnames(zr))
share=intersect(share,colnames(zc))
share=share[c(3:41)]
# use scale
all=rbind(scale(wcct[,share]),scale(zh[,share])) # scale each individual dataset
all=rbind(all,scale(zr[,share]))
all=rbind(all,scale(zc[,share]))
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
trio=CreateSeuratObject(counts=all_t,assay="RNA")
SetAssayData(object = trio, slot = "data", new.data = all_t, assay="RNA")
SetAssayData(object = trio, slot = "scale.data", new.data = all_t, assay="RNA")
groups=c(rep("wcct",20000),rep("zh",20000),rep("zr",20000),rep("zc",20000))
trio <- AddMetaData(object = trio, metadata = groups, col.name = "group")
obj.list <- SplitObject(trio, split.by = "group")
features=share
trio.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
trio.combined <- IntegrateData(anchorset = trio.anchors)
trio.combined <- ScaleData(trio.combined, verbose = FALSE)
trio.combined <- RunPCA(trio.combined, npcs = 30, verbose = FALSE)
trio_seuratpca=trio.combined@reductions$pca@cell.embeddings
set.seed(42)
all_tsne_trioseurat_all=Rtsne(trio_seuratpca, check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_trioseurat_data_all <- data.frame(x = all_tsne_trioseurat_all$Y[,1], y = all_tsne_trioseurat_all$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(all_tsne_trioseurat_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("seurat trio with 10cca") + scale_color_manual(values =  clr)+theme_classic()


## tsne plot for liger
# use not scaled
all=rbind(wcct[,share],zh[,share]) # no scale with liger use in package scaling
all=rbind(all,zr[,share])
all=rbind(all,zc[,share])
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
all_t1=all_t[,c(1:20000)]
all_t2=all_t[,c(20001:40000)]
all_t3=all_t[,c(40001:60000)]
all_t4=all_t[,c(60001:80000)]
trio_liger <- createLiger(list(wcct =all_t1, zh = all_t2, zr = all_t3, zc = all_t4 ))
trio_liger <- rliger::normalize(trio_liger)
trio_liger@var.genes=share #  just use all
trio_liger <- scaleNotCenter(trio_liger)
trio_liger <- optimizeALS(trio_liger, k=30)
trio_liger <- quantile_norm(trio_liger)
a=trio_liger@H.norm
set.seed(42)
all_tsne_triolr_all=Rtsne(a[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_triolr_data_all <- data.frame(x = all_tsne_triolr_all$Y[,1], y = all_tsne_triolr_all$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(all_tsne_triolr_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("lr trio with 10cca") + scale_color_manual(values =  clr)+theme_classic()


## tsne plot for fastmnn
all=rbind(scale(wcct[,share]),scale(zh[,share])) # scale each individual dataset
all=rbind(all,scale(zr[,share]))
all=rbind(all,scale(zc[,share]))
all_t=t(all)
colnames(all_t)=as.character(c(1:dim(all_t)[2]))
all_t1=all_t[,c(1:20000)]
all_t2=all_t[,c(20001:40000)]
all_t3=all_t[,c(40001:60000)]
all_t4=all_t[,c(60001:80000)]
fmnn_trio <- fastMNN(all_t1,all_t2,all_t3,all_t4)
fmmn_reduce_trio = fmnn_trio@int_colData$reducedDims@listData$corrected

## tsne plot for scanorama
sc1=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/wcct_il4_scanpca.csv", header = FALSE)
sc2=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/zh_il4_scanpca.csv", header = FALSE)
sc3=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/zr_il4_scanpca.csv", header = FALSE)
sc4=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/zc_il4_scanpca.csv", header = FALSE)
scall=rbind(sc1,sc2)
scall=rbind(scall,sc3)
scall=rbind(scall,sc4)
set.seed(42)
all_tsne_trioscan_all=Rtsne(scall[,c(1:10)], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_trioscan_data_all <- data.frame(x = all_tsne_trioscan_all$Y[,1], y = all_tsne_trioscan_all$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(all_tsne_trioscan_data_all)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("scanorama trio with 10cca") + scale_color_manual(values =  clr)+theme_classic()


## tsne plot for pre-integration dataset
wcct1_value_scale=as.data.frame(scale(wcct[,share])) # original data but scaled per dataset, as usual
zh_value_scale=as.data.frame(scale(zh[,share]))
zr_value_scale=as.data.frame(scale(zr[,share]))
zc_value_scale=as.data.frame(scale(zc[,share]))
wcct1_value_scale$dataset="wcct"
zh_value_scale$dataset="zh"
zr_value_scale$dataset="zr"
zc_value_scale$dataset="zc"
wzhrc=rbind(wcct1_value_scale,zh_value_scale,zr_value_scale,zc_value_scale)
set.seed(42)
pre_tsne_scale_10=Rtsne(wzhrc, check_duplicates = FALSE, num_threads = 10) # run with cca
pre_tsne_data_scale_10 <- data.frame(x = pre_tsne_scale_10$Y[,1], y = pre_tsne_scale_10$Y[,2])
label=wzhrc$dataset

clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(pre_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("pre integration") + scale_color_manual(values =  clr )





############################### figureS4.2 ##############################################




## read in the MARIO result use same cells for other methods

codex_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCodex_mario_matched.csv")
cite_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCiteseq_mario_matched.csv")

codex_cca = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/cca_codex.csv")
cite_cca = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/cca_cite.csv")

#set.seed(41)
#randidx=sample(dim(codex_input)[1]/2,8000, replace = FALSE)
# make sure randix is reproducible:
#randidx = read.csv("/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/sup_related/randidx_from_r.csv")
#randidx = randidx$rand
codex_input_rand = codex_matched[randidx,]
cite_input_rand = cite_matched[randidx,]
shared=intersect(colnames(codex_input_rand), colnames(cite_input_rand))

## tsne plot for pre-integration
codex_input_rand_value=codex_input_rand[,shared[c(2:24)]]
codex_input_rand_value=scale(codex_input_rand_value)
cite_input_rand_value=cite_input_rand[,shared[c(2:24)]]
cite_input_rand_value=scale(cite_input_rand_value)
rand_all=rbind(codex_input_rand_value,cite_input_rand_value)

set.seed(41)
pre_tsne_scale_10=Rtsne(rand_all, check_duplicates = FALSE, num_threads = 10) # run with cca
pre_tsne_data_scale_10 <- data.frame(x = pre_tsne_scale_10$Y[,1], y = pre_tsne_scale_10$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000)))
p = ggplot(pre_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("pre integration") 

## tsne plot for MARIO
all_cca_sub = rbind(codex_cca[randidx,],cite_cca[randidx,])
set.seed(41)
all_tsne_sub=Rtsne(all_cca_sub[,-1], check_duplicates = FALSE, num_threads=10) # run with cca
all_tsne_data_sub<- data.frame(x = all_tsne_sub$Y[,1], y = all_tsne_sub$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000)))
p = ggplot(all_tsne_data_sub)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 1.2), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("MARIO integration with 10cca")

## tsne plot for Seurat
x_value_scale=codex_input_rand_value # scaled
y_value_scale=cite_input_rand_value # scaled
x_value_scale=as.data.frame(x_value_scale)
rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
y_value_scale=as.data.frame(y_value_scale)
rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
# create seurat object datax
x_obj=CreateSeuratObject(counts=t(x_value_scale),assay="x")
SetAssayData(object = x_obj, slot = "data", new.data = t(x_value_scale), assay="x")
x_obj <- ScaleData(x_obj)
# add suerat object datay
y_obj=CreateSeuratObject(counts=t(y_value_scale),assay="y")
SetAssayData(object = y_obj, slot = "data", new.data = t(y_value_scale), assay="y")
y_obj <- ScaleData(y_obj)
list_modality=list(x_obj,y_obj)
features=intersect(colnames(x_value_scale),colnames(x_value_scale))
## test run pca
y_obj = RunPCA(y_obj, features = features, npcs = 10)
x_obj = RunPCA(x_obj, features = features, npcs = 10)
#### anchors for integration pca
Int.anchors <- FindIntegrationAnchors(object.list = list_modality,
                                      dims = 1:10, anchor.features =features, k.filter = 10)
xy_int <- IntegrateData(anchorset = Int.anchors, dims = 1:10, k.weight = 10)
######WARNING########## in the integration process NaN is produced in the integrated data slot
temp = xy_int@assays$integrated@data
temp[is.na(temp)] <- 0
# new obj
temp_obj = CreateSeuratObject(counts=as.matrix(temp),assay="temp")
SetAssayData(object = temp_obj, slot = "data", new.data = as.matrix(temp), assay="temp")
temp_obj <- ScaleData(temp_obj)
temp_obj <- RunPCA(temp_obj, npcs = 10, verbose = FALSE, features = features) # need to change this to 
embedding = temp_obj@reductions$pca@cell.embeddings
set.seed(41)
seurat_tsne_10=Rtsne(embedding, check_duplicates = FALSE, num_threads = 10) # run with cca
seurat_tsne_data_10 <- data.frame(x = seurat_tsne_10$Y[,1], y = seurat_tsne_10$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000)))
p = ggplot(seurat_tsne_data_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("seurat integration")

## tsne plots for liger
codex_input_rand = codex_matched[randidx,]
cite_input_rand = cite_matched[randidx,]
shared=intersect(colnames(codex_input_rand), colnames(cite_input_rand))
codex_input_rand_value=codex_input_rand[,shared[c(2:24)]]
#codex_input_rand_value=scale(codex_input_rand_value)
cite_input_rand_value=cite_input_rand[,shared[c(2:24)]]
#cite_input_rand_value=scale(cite_input_rand_value)
x_value_scale = as.data.frame(codex_input_rand_value)
y_value_scale = as.data.frame(cite_input_rand_value)
rownames(x_value_scale)=paste0(as.character(rownames(x_value_scale)),"_x")
rownames(y_value_scale)=paste0(as.character(rownames(y_value_scale)),"_y")
ligerobj=createLiger( list(x = t(x_value_scale), y = t(y_value_scale) ) )
###Start integration, feature will drop so change this
features=intersect(colnames(x_value_scale),colnames(y_value_scale))
# default preprocessing
ligerobj <- rliger::normalize(ligerobj)
#ifnb_liger <- selectGenes(ifnb_liger, var.thresh = 0, alpha.thresh=1)
ligerobj@var.genes=features #  just use all
ligerobj <- scaleNotCenter(ligerobj)
ligerobj <- optimizeALS(ligerobj, k = 20)
ligerobj <- quantile_norm(ligerobj)
embedding = ligerobj@H.norm # pca scores
set.seed(41)
liger_tsne_10=Rtsne(embedding[,c(1:10)], check_duplicates = FALSE, num_threads = 10) # run with cca
liger_tsne_data_10 <- data.frame(x = liger_tsne_10$Y[,1], y = liger_tsne_10$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000)))
p = ggplot(liger_tsne_data_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("liger integration")


## tsne plots for fstmnn
fmnn <- fastMNN(t(scale(x_value_scale[features])), t(scale(y_value_scale[features])),d=20) # scale for better results
embedding = fmnn@int_colData$reducedDims@listData$corrected
set.seed(41)
fstmnn_tsne_10=Rtsne(embedding[,c(1:10)], check_duplicates = FALSE, num_threads = 10) # run with cca
fstmnn_tsne_data_10 <- data.frame(x = fstmnn_tsne_10$Y[,1], y = fstmnn_tsne_10$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000)))
p = ggplot(fstmnn_tsne_data_10)  + 
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("fstmnn integration")

## tsne plots for scanorama
scan_x = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/murine_codex_scan.csv")
scan_y = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/temp/murine_cite_scan.csv")
embedding = rbind(scan_x, scan_y)
set.seed(41)
scan_tsne_10=Rtsne(embedding[,-1], check_duplicates = FALSE, num_threads = 10) # run with cca
scan_tsne_data_10 <- data.frame(x = scan_tsne_10$Y[,1], y = scan_tsne_10$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000)))
p = ggplot(scan_tsne_data_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("scan integration")
### Bokai Zhu
### 0823-2021
##### script related to the production of Figure S3.2 pre and post Mario analysis
##### (integrative analysis on wcct human h1n1 -> Xspecies human rhesus cyno IL4 sitmulated data)




#################################### part 1 #######################################

######### part 1 is preprocessing and clustering/annotation of the raw cytof cells
######### files large and clustering is slow, does not recommend re-running this part of code
######### generally this part is not need for replication of the figures, only for record

# raw file is large, original fcs file from flowrepo + gating info from Zach+Dave+Han, share per request

#### pre-processing and clustering of the Xspecies human IL4 stimulated cytof dataset
basal=fread("/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/CrossSpecies_human_IL-4_annotated.csv")# raw file
basal_3d=subset(basal,basal$Donor %in% c("7826","7718","2810"))
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di",
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di",
           "Ir191Di","Ir193Di","Donor","Species","Condition","Filename","DNA1","DNA2")
basal_3d=as.data.frame(basal_3d)
basal_3d_clean=basal_3d[,!(colnames(basal_3d) %in% drops)]
gate_cols=basal_3d_clean[,c(41:55)]
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
basal_3d_clean$cluster.orig=annotation
basal_3d_clean=basal_3d_clean[,-c(41:55)]
#arcshin transform by cofactor 5
cofactor = 5
basal_3d_clean4Value=basal_3d_clean[,c(2:40)]
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
basal_3d_clean[,c(2:40)]=basal_3d_clean4Value_arch
# remove whole blood irrelevant cell types
basal_3d_clean_noNeu=basal_3d_clean
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Erythrocytes")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Platelets")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="CD4+CD8+") 
# random downsample
rows <- sample(120000) # sample 120k cells
basal_3d_clean_noNeu_small <- basal_3d_clean_noNeu[rows, ]
basal_3d_clean_noNeu_small1=basal_3d_clean_noNeu_small[c(1:60000),]
basal_3d_clean_noNeu_small2=basal_3d_clean_noNeu_small[c(60001:120000),]
### make seurat object and cluster / annotation
library(Seurat)
tm=as.data.frame(t(basal_3d_clean_noNeu_small1[,c(2:40)]))
data6 <- CreateSeuratObject(counts =tm)
data6 <- SetAssayData(object = data6, slot = "scale.data", new.data = as.matrix(tm))
data6 <- SetAssayData(object = data6, slot = "data", new.data = as.matrix(tm))
data6 <- RunPCA(data6, features = rownames(data6), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
DimPlot(data6, reduction = "pca_cytof")
ElbowPlot(data6, ndims = 50, reduction = "pca_cytof")
data6 <- RunTSNE(data6, dims = 1:10, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data6 <- FindNeighbors(data6, features = rownames(data6), dims = NULL)
cy.data <- GetAssayData(data6, slot = "data")
cy.dist <- dist(t(cy.data))
data6[["cy_snn"]] <- FindNeighbors(cy.dist)$snn
data6 <- FindClusters(data6, resolution = 0.2, graph.name = "cy_snn")
## add manual annotation
new.cluster.ids <- c("Neutrophil","Neutrophil","CD4 T",
                     "Neutrophil","MC","dirt","CD8 T","NK",
                     "Neutrophil","dirt","B","Neutrophil","Basophil","dirt")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)
basal_3d_clean_noNeu_small1$cluster.sr=Idents(data6)

### repeat the upper part code for the second dataframe clustering
basal_3d_clean_noNeu_smallall=rbind(basal_3d_clean_noNeu_small1,basal_3d_clean_noNeu_small2)
# remove unannotated cell types
basal_3d_clean_noNeu_smallall_nodirt=subset(basal_3d_clean_noNeu_smallall,basal_3d_clean_noNeu_smallall$cluster.sr!="dirt")
#write.csv(basal_3d_clean_noNeu_smallall_nodirt,"/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/zac_IL2_Hu_108538_forSup.csv")


#### pre-processing and clustering of the Xspecies rhesus IL4 stimulated cytof dataset
basal=fread("/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/CrossSpcies_NHP_IL-4_annotated.csv")# file with all NHP
basal_3d=subset(basal,basal$Donor %in% c("D00522","D06022","D06122")) # three rhesus donors
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di",
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di",
           "Ir191Di","Ir193Di","Donor","Species","Condition","Filename","DNA1","DNA2")
basal_3d=as.data.frame(basal_3d)
basal_3d_clean=basal_3d[,!(colnames(basal_3d) %in% drops)]
gate_cols=basal_3d_clean[,c(41:55)] # get gating
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
basal_3d_clean$cluster.orig=annotation
basal_3d_clean=basal_3d_clean[,-c(41:55)]
# before saving, both arcshin transform by cofactor 5
cofactor = 5
basal_3d_clean4Value=basal_3d_clean[,c(2:40)]
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
basal_3d_clean[,c(2:40)]=basal_3d_clean4Value_arch
# delete unused cell types in whole blood cell samples
basal_3d_clean_noNeu=basal_3d_clean
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Erythrocytes")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Platelets")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="CD4+CD8+") 
# random donwsample
rows <- sample(120000) # sample 120k cells
basal_3d_clean_noNeu_small <- basal_3d_clean_noNeu[rows, ]
#apperently 100k cluster at once is too much for seurat, lets do it in a linear fashion
basal_3d_clean_noNeu_small1=basal_3d_clean_noNeu_small[c(1:60000),]
basal_3d_clean_noNeu_small2=basal_3d_clean_noNeu_small[c(60001:120000),]
### make seurat object and cluster / annotation
tm=as.data.frame(t(basal_3d_clean_noNeu_small1[,c(2:40)]))
data6 <- CreateSeuratObject(counts =tm)
data6 <- SetAssayData(object = data6, slot = "scale.data", new.data = as.matrix(tm))
data6 <- SetAssayData(object = data6, slot = "data", new.data = as.matrix(tm))
data6 <- RunPCA(data6, features = rownames(data6), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
DimPlot(data6, reduction = "pca_cytof")
ElbowPlot(data6, ndims = 50, reduction = "pca_cytof")
data6 <- RunTSNE(data6, dims = 1:10, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data6 <- FindNeighbors(data6, features = rownames(data6), dims = NULL)
cy.data <- GetAssayData(data6, slot = "data")
cy.dist <- dist(t(cy.data))
data6[["cy_snn"]] <- FindNeighbors(cy.dist)$snn
data6 <- FindClusters(data6, resolution = 0.2, graph.name = "cy_snn")
new.cluster.ids <- c("Neutrophil","CD4 T","NK","CD8 T",
                     "B","dirt","Neutrophil","MC","CD8 T",
                     "Neutrophil","NK","dirt","CD4 T","Basophil")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)
basal_3d_clean_noNeu_small1$cluster.sr=Idents(data6) # add clustering annotation to the dataset

### repeat upper part for df2 clustering

# save output remove unannotated cells
basal_3d_clean_noNeu_small2$cluster.sr=Idents(data6)
basal_3d_clean_noNeu_smallall=rbind(basal_3d_clean_noNeu_small1,basal_3d_clean_noNeu_small2)
basal_3d_clean_noNeu_smallall_nodirt=subset(basal_3d_clean_noNeu_smallall,basal_3d_clean_noNeu_smallall$cluster.sr!="dirt")
#write.csv(basal_3d_clean_noNeu_smallall_nodirt,"/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/zac_IL2_rhesMac_110328_forSup.csv")

#### pre-processing and clustering of the Xspecies Cynolo IL4 stimulated cytof dataset

#basal=fread("/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/CrossSpcies_NHP_IL-4_annotated.csv") #already loaded previously
basal_3d=subset(basal,basal$Donor %in% c("D07282","D07292","D07322")) # three cyno donor 
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di",
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di",
           "Ir191Di","Ir193Di","Donor","Species","Condition","Filename","DNA1","DNA2")
basal_3d=as.data.frame(basal_3d)
basal_3d_clean=basal_3d[,!(colnames(basal_3d) %in% drops)]
gate_cols=basal_3d_clean[,c(41:55)]
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
basal_3d_clean$cluster.orig=annotation
basal_3d_clean=basal_3d_clean[,-c(41:55)]
# arcshin transform by cofactor 5
cofactor = 5
basal_3d_clean4Value=basal_3d_clean[,c(2:40)]
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
basal_3d_clean[,c(2:40)]=basal_3d_clean4Value_arch
basal_3d_clean_noNeu=basal_3d_clean
# remove irrelevant cell types in whole blood samples
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Erythrocytes")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Platelets")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="CD4+CD8+") 
rows <- sample(120000) # sample 120k cells
basal_3d_clean_noNeu_small <- basal_3d_clean_noNeu[rows, ]
# split 2 fast cluster
basal_3d_clean_noNeu_small1=basal_3d_clean_noNeu_small[c(1:60000),]
basal_3d_clean_noNeu_small2=basal_3d_clean_noNeu_small[c(60001:120000),]
### make seurat object and cluster / annotation
tm=as.data.frame(t(basal_3d_clean_noNeu_small1[,c(2:40)]))
data6 <- CreateSeuratObject(counts =tm)
data6 <- SetAssayData(object = data6, slot = "scale.data", new.data = as.matrix(tm))
data6 <- SetAssayData(object = data6, slot = "data", new.data = as.matrix(tm))
data6 <- RunPCA(data6, features = rownames(data6), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
DimPlot(data6, reduction = "pca_cytof")
ElbowPlot(data6, ndims = 50, reduction = "pca_cytof")
data6 <- RunTSNE(data6, dims = 1:10, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data6 <- FindNeighbors(data6, features = rownames(data6), dims = NULL)
cy.data <- GetAssayData(data6, slot = "data")
cy.dist <- dist(t(cy.data))
data6[["cy_snn"]] <- FindNeighbors(cy.dist)$snn
data6 <- FindClusters(data6, resolution = 0.2, graph.name = "cy_snn")
# cluster manual annotation
new.cluster.ids <- c("CD4 T","Neutrophil","dirt","Neutrophil",
                     "CD8 T","dirt","cMC","B","NK","NK","dirt","dirt","Basophil","dirt")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)
basal_3d_clean_noNeu_small1$cluster.sr=Idents(data6)

###repeat upper part for df2

# remove unannotated cells
basal_3d_clean_noNeu_smallall=rbind(basal_3d_clean_noNeu_small1,basal_3d_clean_noNeu_small2)
basal_3d_clean_noNeu_smallall_nodirt=subset(basal_3d_clean_noNeu_smallall,basal_3d_clean_noNeu_smallall$cluster.sr!="dirt")
#write.csv(basal_3d_clean_noNeu_smallall_nodirt,"/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/zac_IL2_cyno_90302__forSup.csv")

##############################################
########   finished prepping data  ###########
########   mario in python script  ###########
##############################################

#############################  part 2 ##################################

### this part is about the production of figures in figS3.2 for the manuscript

wcct=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/wcctHuman2_mario_matched_sub20k.csv") # wcct human
zh=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/XhumanIL4_mario_matched_sub20k.csv") # xspecies human
zr=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/XRhesusIL4_mario_matched_sub20k.csv") # xspecies rhesus
zc=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/XCynoIL4_mario_matched_sub20k.csv") # xspecies cyno
# general cca scores from MARIO
cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IL4/gcca_sub20k.csv", header = TRUE)

# start producing the tsne plots
library(Rtsne)
all_cca_scale_10=cca[,-1]
set.seed(42)
all_tsne_scale_10=Rtsne(all_cca_scale_10, check_duplicates = FALSE, num_threads = 10) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000))) # lazy modality
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
# tsne plot with modality color
p = ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct-zac-integration with 10cca") + scale_color_manual(values =  clr )
# tsne plot with only wcct
label=as.character(wcct$cluster.sr)
temp = all_tsne_data_scale_10[c(1:20000),]
temp = temp[label != "Basophils",] # clean out basophil, cell population to little, ignore
label = label[label != "Basophils"]
p = ggplot(temp)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.5), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") + 
  theme_classic() + scale_color_manual(values =clr )
# tsne plot with only xspecies human il4
label=as.character(zh$cluster.sr)
temp = all_tsne_data_scale_10[c(20001:40000),]
temp = temp[label != "Basophils",] # clean out basophil, cell population to little, ignore
label = label[label != "Basophils"]
p = ggplot(temp)  + 
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.5), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") + 
  theme_classic() + scale_color_manual(values =clr)
# tsne plot with only xspecies rhesus il4
label=as.character(zr$cluster.sr)
temp = all_tsne_data_scale_10[c(40001:60000),]
temp = temp[label != "Basophils",] # clean out basophil, cell population to little, ignore
label = label[label != "Basophils"]
p = ggplot(temp)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.5), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") +
  theme_classic() + scale_color_manual(values =clr )
# tsne plot with only xspecies cyno il4
label=as.character(zc$cluster.sr)
temp = all_tsne_data_scale_10[c(60001:80000),]
temp = temp[label != "Basophils",] # clean out basophil, cell population to little, ignore
label = label[label != "Basophils"]
p = ggplot(temp)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.5), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") +
  theme_classic() + scale_color_manual(values =clr )

### produce tsne plot pre-integration
wcct1_value_scale=as.data.frame(scale(wcct[,c(3:43)]))
zh_value_scale=as.data.frame(scale(zh[,c(3:41)]))
zr_value_scale=as.data.frame(scale(zr[,c(3:41)]))
zc_value_scale=as.data.frame(scale(zc[,c(3:41)]))
wcct1_value_scale$dataset="wcct"
zh_value_scale$dataset="zh"
zr_value_scale$dataset="zr"
zc_value_scale$dataset="zc"
share1=intersect(colnames(wcct1_value_scale),colnames(zh_value_scale))
share2=intersect(share1, colnames(zr_value_scale))
share3=intersect(share2, colnames(zc_value_scale))
# use shared features to produce
wzhrc=rbind(wcct1_value_scale[,share3],zh_value_scale[,share3],zr_value_scale[,share3],zc_value_scale[,share3])
set.seed(42)
pre_tsne_scale_10=Rtsne(wzhrc[,c(1:39)], check_duplicates = FALSE, num_threads = 10) # run with cca
pre_tsne_data_scale_10 <- data.frame(x = pre_tsne_scale_10$Y[,1], y = pre_tsne_scale_10$Y[,2])
label=wzhrc$dataset
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000))) #lazy modality
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
# tsne plot pre-integration
p = ggplot(pre_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct-zac-integration with 10cca") +
  scale_color_manual(values =  clr ) +theme_classic()


## now start to produce the feature plots
# first create seurat objects
# wcct
library(Seurat)
wcct1_value=wcct[,c(3:43)]
rownames(wcct1_value)=wcct$index
wcct1obj=CreateSeuratObject(counts=t(wcct1_value),assay="cytof")
SetAssayData(object = wcct1obj, slot = "data", new.data = t(wcct1_value), assay="cytof")
cca_embed=as.data.frame(all_tsne_data_scale_10[c(1:20000),])
colnames(cca_embed) <- paste0("ccatsne_", 1:2)
rownames(cca_embed)=wcct$index
wcct1obj[["ccatsne"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "ccatsne_", assay = "cytof")
wcct1obj@meta.data$anno.term=as.character(wcct$cluster.sr)
# xsp human il4
zh_value=zh[,c(3:41)]
zh_id_nonDup=make.names(as.character(zh$index), unique = TRUE)
rownames(zh_value)=zh_id_nonDup
zhobj=CreateSeuratObject(counts=t(zh_value),assay="cytof")
SetAssayData(object = zhobj, slot = "data", new.data = t(zh_value), assay="cytof")
cca_embed=as.data.frame(all_tsne_data_scale_10[c(20001:40000),])
colnames(cca_embed) <- paste0("ccatsne_", 1:2)
rownames(cca_embed)=zh_id_nonDup
zhobj[["ccatsne"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "ccatsne_", assay = "cytof")
zhobj@meta.data$anno.term=as.character(zh$cluster.sr)
# xsp rhesus il4
zr_value=zr[,c(3:41)]
zr_id_nonDup=make.names(as.character(zr$index), unique = TRUE)
rownames(zr_value)=zr_id_nonDup
zrobj=CreateSeuratObject(counts=t(zr_value),assay="cytof")
SetAssayData(object = zrobj, slot = "data", new.data = t(zr_value), assay="cytof")
cca_embed=as.data.frame(all_tsne_data_scale_10[c(40001:60000),])
colnames(cca_embed) <- paste0("ccatsne_", 1:2)
rownames(cca_embed)=zr_id_nonDup
zrobj[["ccatsne"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "ccatsne_", assay = "cytof")
zrobj@meta.data$anno.term=as.character(zr$cluster.sr)
# xsp cyno il4
zc_value=zc[,c(3:41)]
zc_id_nonDup=make.names(as.character(zc$index), unique = TRUE)
rownames(zc_value)=zc_id_nonDup
zcobj=CreateSeuratObject(counts=t(zc_value),assay="cytof")
SetAssayData(object = zcobj, slot = "data", new.data = t(zc_value), assay="cytof")
cca_embed=as.data.frame(all_tsne_data_scale_10[c(60001:80000),])
colnames(cca_embed) <- paste0("ccatsne_", 1:2)
rownames(cca_embed)=zc_id_nonDup
zcobj[["ccatsne"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "ccatsne_", assay = "cytof")
zcobj@meta.data$anno.term=as.character(zc$cluster.sr)

## start producing feature plots

DefaultAssay(wcct1obj) <- 'cytof'
p=FeaturePlot(wcct1obj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/wcct_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(wcct1obj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/wcct_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(wcct1obj) <- 'cytof'
p=FeaturePlot(wcct1obj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/wcct_response3_refg.png", plot=p, width=7, height=7)


DefaultAssay(zhobj) <- 'cytof'
p=FeaturePlot(zhobj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zh_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(zhobj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zh_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(zhobj) <- 'cytof'
p=FeaturePlot(zhobj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zh_response3_refg.png", plot=p, width=7, height=7)


DefaultAssay(zrobj) <- 'cytof'
p=FeaturePlot(zrobj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zr_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(zrobj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zr_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(zrobj) <- 'cytof'
p=FeaturePlot(zrobj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zr_response3_refg.png", plot=p, width=7, height=7)


DefaultAssay(zcobj) <- 'cytof'
p=FeaturePlot(zcobj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zc_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(zcobj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zc_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(zcobj) <- 'cytof'
p=FeaturePlot(zcobj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/sup/sup_plots/white_feature_version/zc_response3_refg.png", plot=p, width=7, height=7)

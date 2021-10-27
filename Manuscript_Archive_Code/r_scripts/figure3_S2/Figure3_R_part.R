### Bokai Zhu
### 0823-2021
##### script related to the production of Figure 2 pre and post Mario analysis
##### (integrative analysis on WCCT h1n1 virus / vaccine chanllenge + Xspecies Zachary cytof project)

#################################### part 1 #######################################

######### part 1 is preprocessing and clustering/annotation of the raw cytof cells
######### files large and clustering is slow, does not recommend re-running this part of code
######### generally this part is not need for replication of the figures, only for record

library(data.table) # file is large fast read
# raw file is large, original fcs file from flowrepo + gating info from Zach+Dave+Han, share per request
#### pre-processing and clustering of the wcct h1n1 cytof dataset
wcctg_full=fread("/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/WCCTG14-01_annotated.csv") # fast read
wcctg_3d=subset(wcctg_full,wcctg_full$ID %in% c(101,107,108)) # >900k cells
wcctg_3don_4d=subset(wcctg_3d,wcctg_3d$DAY==4) # 300k cell left at day 4
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di", # drop the unused channels
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di","Ir191Di",
           "Ir193Di","DAY","Filename","ID")
wcctg_3don_4d=as.data.frame(wcctg_3don_4d)
wcctg_3don_4d_clean=wcctg_3don_4d[,!(colnames(wcctg_3don_4d) %in% drops)]
# now clean up the annotation
gate_cols=wcctg_3don_4d_clean[,c(43:65)] # gating info
gate_cols <- gate_cols[c(1:19,21,20,22:23)] # reorder the annotation for cd8 because we want to find the lowest annotation
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
wcctg_3don_4d_clean$cluster.orig=annotation
wcctg_3don_4d_clean=wcctg_3don_4d_clean[,-c(43:65)] # clean out
# before saving, arcshin transform by cofactor 5
cofactor = 5
basal_3d_clean4Value=wcctg_3don_4d_clean[,c(2:42)]
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
wcctg_3don_4d_clean[,c(2:42)]=basal_3d_clean4Value_arch
# random downsample to 120k cells:
rows <- sample(120000) 
wcctg_3don_4d_clean_small <- wcctg_3don_4d_clean[rows, ]
# full size too large for clustering so sepearte into two and repeat the code
wcctg_3don_4d_clean_small1=wcctg_3don_4d_clean_small[c(1:60000),]
wcctg_3don_4d_clean_small2=wcctg_3don_4d_clean_small[c(60001:120000),]
#### start the seurat object and clustering
tm=as.data.frame(t(wcctg_3don_4d_clean_small2[,c(2:42)])) # only input the antibody signals
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
data6 <- FindClusters(data6, resolution = 0.1, graph.name = "cy_snn")
data6 <- RunTSNE(data6, dims = 1:25, method = "FIt-SNE", reduction = "pca_cytof",check_duplicates = FALSE)
DimPlot(data6, label = TRUE) + NoLegend()
# subcluster
data6 <- FindSubCluster(data6, cluster = 3, graph.name = "cy_snn",resolution=0.05)
data6_store=Idents(data6)
Idents(data6) <- "sub.cluster"
DimPlot(data6, label = TRUE) + NoLegend()
# human annotation
new.cluster.ids <- c("CD4 T","cMC","Neutrophil","MC","CD8 T","NK","dirt","dirt","MC","B","dirt","B","Basophil","dirt")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)
wcctg_3don_4d_clean_small2$cluster.sr=Idents(data6)
wcctg_3don_4d_clean_all=rbind(wcctg_3don_4d_clean_small1,wcctg_3don_4d_clean_small2) # combine the two sep files agin
wcctg_3don_4d_clean_all_term=subset(wcctg_3don_4d_clean_all, wcctg_3don_4d_clean_all$cluster.sr !="dirt") # remove unannotated cells
#write.csv(wcctg_3don_4d_clean_all_term,"/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/wcct4d_term_srcluster120k.csv")

#### pre-processing and clustering of the wcct h1n1 cytof dataset end


#### pre-processing and clustering of the cytof human dataset
# Xcross species IFNG human fcs + gating information, file share per request
basal=fread("/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/cross species human IFNg_annotated.csv")
# 22 milion cells, only chose 3 donors: 7826,7718,2810
basal_3d=subset(basal,basal$Donor %in% c("7826","7718","2810")) #600k cells
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di",
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di","Ir191Di",
           "Ir193Di","Donor","Species","Condition","Filename","DNA1","DNA2")
basal_3d=as.data.frame(basal_3d)
basal_3d_clean=basal_3d[,!(colnames(basal_3d) %in% drops)]
gate_cols=basal_3d_clean[,c(41:55)]
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
basal_3d_clean$cluster.orig=annotation
basal_3d_clean=basal_3d_clean[,-c(41:55)]
# before saving, both arcshin transform by cofactor 5
cofactor = 5
basal_3d_clean4Value=basal_3d_clean[,c(2:40)]
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
basal_3d_clean[,c(2:40)]=basal_3d_clean4Value_arch
# since this dataset is whole cell not pbmc, remove not used cell types for downstream using
basal_3d_clean_noNeu=basal_3d_clean
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Erythrocytes")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Platelets")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="CD4+CD8+") 
# random downsample to 120k cells
rows <- sample(120000)
basal_3d_clean_noNeu_small <- basal_3d_clean_noNeu[rows, ]
#same problem, smaller dataframe for faster clustering
basal_3d_clean_noNeu_small1=basal_3d_clean_noNeu_small[c(1:60000),]
basal_3d_clean_noNeu_small2=basal_3d_clean_noNeu_small[c(60001:120000),]

#### start the seurat object and clustering and annotation

tm=as.data.frame(t(basal_3d_clean_noNeu_small1[,c(2:40)])) # only input protein information
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
data6 <- FindClusters(data6, resolution = 0.1, graph.name = "cy_snn")
# some sub cell type clustering
data6 <- FindSubCluster(data6, cluster = 2, graph.name = "cy_snn",resolution=0.08)
Idents(data6) <- "sub.cluster"
# manual annotation
new.cluster.ids <- c("Neutrophil","CD4 T","dirt","MC","CD8 T","MC","Neutro_like","dirt","NK","Neutro_like2","MC","B","Basophil")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)

### then need to rerun the upper part for the second dataframe

# save the files with removed cells without annotations
basal_3d_clean_noNeu_small2$cluster.sr=Idents(data6)
basal_3d_clean_noNeu_smallall_nodirt=subset(basal_3d_clean_noNeu_smallall,basal_3d_clean_noNeu_smallall$cluster.sr!="dirt")
basal_3d_clean_noNeu_smallall_nodirt$cluster.sr[basal_3d_clean_noNeu_smallall_nodirt$cluster.sr=="Neutro_like"]="Neutrophil"
basal_3d_clean_noNeu_smallall_nodirt$cluster.sr[basal_3d_clean_noNeu_smallall_nodirt$cluster.sr=="Neutro_like2"]="Neutrophil"
#write.csv(basal_3d_clean_noNeu_smallall_nodirt,"/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/zac_ifgn_term_srcluster120k.csv")

#### pre-processing and clustering of the cytof NHP dataset
# Xcross species IFNG NHP fcs + gating information, file share per request
# this part is for rhesus macaque
basal=fread("/home/bkzhu/SNE-multi/figure_rcode/siv/zac_recluster/CrossSpcies_NHP_IFNb_annotated.csv")# nhp ifnb
basal_3d=subset(basal,basal$Donor %in% c("D00522","D06022","D06122")) # 600k, three rhesus macaque
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di",
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di","Ir191Di",
           "Ir193Di","Donor","Species","Condition","Filename","DNA1","DNA2")
basal_3d=as.data.frame(basal_3d)
basal_3d_clean=basal_3d[,!(colnames(basal_3d) %in% drops)]
gate_cols=basal_3d_clean[,c(41:55)] # get gating info
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
basal_3d_clean$cluster.orig=annotation
basal_3d_clean=basal_3d_clean[,-c(41:55)]
# before saving, both arcshin transform by cofactor 5
cofactor = 5
basal_3d_clean4Value=basal_3d_clean[,c(2:40)] # antibody values
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
basal_3d_clean[,c(2:40)]=basal_3d_clean4Value_arch
# since this is whole blood cell data remove irrelavent cells
basal_3d_clean_noNeu=basal_3d_clean
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Erythrocytes")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Platelets")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="CD4+CD8+") 
# random downsample dataset
rows <- sample(120000) # sample 120k cells
basal_3d_clean_noNeu_small <- basal_3d_clean_noNeu[rows, ]
#seperate data for fast clustering
basal_3d_clean_noNeu_small1=basal_3d_clean_noNeu_small[c(1:60000),]
basal_3d_clean_noNeu_small2=basal_3d_clean_noNeu_small[c(60001:120000),]

#### start seurat object and clustering annotation
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
data6 <- FindClusters(data6, resolution = 0.1, graph.name = "cy_snn")
# sub cluster cell type
data6.c3 <- subset(data6, idents = 7)
data6.c3 <- FindNeighbors(data6.c3, features = rownames(data6.c3), dims = NULL)
cy.data3 <- GetAssayData(data6.c3, slot = "data")
cy.dist3 <- dist(t(cy.data3))
data6.c3[["cy_snn"]] <- FindNeighbors(cy.dist3)$snn
data6.c3 <- FindClusters(data6.c3, resolution = 0.1, graph.name = "cy_snn")
data6$sub_cluster <- as.character(Idents(data6))
data6$sub_cluster[Cells(data6.c3)] <- paste("c7_",Idents(data6.c3))
# sub cluster cell type
data6.c3 <- subset(data6, idents = 8)
data6.c3 <- FindNeighbors(data6.c3, features = rownames(data6.c3), dims = NULL)
cy.data3 <- GetAssayData(data6.c3, slot = "data")
cy.dist3 <- dist(t(cy.data3))
data6.c3[["cy_snn"]] <- FindNeighbors(cy.dist3)$snn
data6.c3 <- FindClusters(data6.c3, resolution = 0.05, graph.name = "cy_snn")
data6_store=Idents(data6)
# manual annotation
Idents(data6) <- "sub_cluster"
new.cluster.ids <- c("B","Neutrophil","CD8 T","CD4 T","NK","CD8 T","MC","Neutrophil","dirt","MC","CD4 T","MC","NK")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)

#### rerun this part of code for clustering the second dataframe

basal_3d_clean_noNeu_small2$cluster.sr=Idents(data6)
basal_3d_clean_noNeu_smallall=rbind(basal_3d_clean_noNeu_small1,basal_3d_clean_noNeu_small2)
basal_3d_clean_noNeu_smallall_nodirt=subset(basal_3d_clean_noNeu_smallall,basal_3d_clean_noNeu_smallall$cluster.sr!="dirt")
table(basal_3d_clean_noNeu_smallall_nodirt$cluster.sr)
dim(basal_3d_clean_noNeu_smallall_nodirt) # 110k cells
basal_3d_clean_noNeu_smallall_nodirt$cluster.sr[basal_3d_clean_noNeu_smallall_nodirt$cluster.sr=="Neutro_like"]="Neutrophil"
basal_3d_clean_noNeu_smallall_nodirt$cluster.sr[basal_3d_clean_noNeu_smallall_nodirt$cluster.sr=="Neutro_like2"]="Neutrophil"
#write.csv(basal_3d_clean_noNeu_smallall_nodirt,"/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/zac_Rheus_ifgn_term_srcluster120k.csv")


#### pre-processing and clustering of the cytof NHP dataset
# Xcross species IFNG NHP fcs + gating information, file share per request
# this part is for cynomolgus
basal=fread("/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/CrossSpcies_NHP_IFNg_annotated.csv")# nhp ifng
basal_3d=subset(basal,basal$Donor %in% c("D07282","D07292","D07322")) # 769k, three cynomolgus monkey
drops <- c("Time","Event_length","Pd102Di","Pd104Di","Pd105Di",
           "Pd106Di","Pd108Di","Pd110Di","Ba138Di","Sm150Di","Ir191Di",
           "Ir193Di","Donor","Species","Condition","Filename","DNA1","DNA2")
basal_3d=as.data.frame(basal_3d)
basal_3d_clean=basal_3d[,!(colnames(basal_3d) %in% drops)]
gate_cols=basal_3d_clean[,c(41:55)] # get gating information
annotation=names(gate_cols)[max.col(gate_cols == 1)] # finding the max colm
basal_3d_clean$cluster.orig=annotation
basal_3d_clean=basal_3d_clean[,-c(41:55)]
# archsin transformation with cofactor 5
cofactor = 5
basal_3d_clean4Value=basal_3d_clean[,c(2:40)]
basal_3d_clean4Value_arch = asinh(basal_3d_clean4Value/cofactor)
basal_3d_clean[,c(2:40)]=basal_3d_clean4Value_arch
# whole blood sample, remove irrelevant cell types based on gating
basal_3d_clean_noNeu=basal_3d_clean
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Erythrocytes")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="Platelets")
basal_3d_clean_noNeu=subset(basal_3d_clean_noNeu, basal_3d_clean_noNeu$cluster.orig!="CD4+CD8+") 
# random downsample
rows <- sample(120000) # sample 120k cells
basal_3d_clean_noNeu_small <- basal_3d_clean_noNeu[rows, ]
# sep for fast clustering
basal_3d_clean_noNeu_small1=basal_3d_clean_noNeu_small[c(1:60000),]
basal_3d_clean_noNeu_small2=basal_3d_clean_noNeu_small[c(60001:120000),]
### start seurat object and clustering/ annotation
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
data6 <- FindClusters(data6, resolution = 0.1, graph.name = "cy_snn")
#sub cluster
data6.c3 <- subset(data6, idents = 7)
data6.c3 <- FindNeighbors(data6.c3, features = rownames(data6.c3), dims = NULL)
cy.data3 <- GetAssayData(data6.c3, slot = "data")
cy.dist3 <- dist(t(cy.data3))
data6.c3[["cy_snn"]] <- FindNeighbors(cy.dist3)$snn
data6.c3 <- FindClusters(data6.c3, resolution = 0.07, graph.name = "cy_snn")
# Generate a new column called sub_cluster in the metadata
data6$sub_cluster <- as.character(Idents(data6))
# Change the information of cells containing sub-cluster information
data6$sub_cluster[Cells(data6.c3)] <- paste("c7_",Idents(data6.c3))
data6_store=Idents(data6)
Idents(data6) <- "sub_cluster"
# add manual annotation
new.cluster.ids <- c("B","dirt","Neutrophil","CD4 T","NK"," CD8 T","Neutrophil","MC","dirt","MC","Neutrophil","dirt")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)

#### repeat the upper code for clustering of the second dataframe

# remove unanottated cell type and change annotation name
basal_3d_clean_noNeu_small2$cluster.sr=Idents(data6)
basal_3d_clean_noNeu_smallall=rbind(basal_3d_clean_noNeu_small1,basal_3d_clean_noNeu_small2)
basal_3d_clean_noNeu_smallall_nodirt=subset(basal_3d_clean_noNeu_smallall,basal_3d_clean_noNeu_smallall$cluster.sr!="dirt")
basal_3d_clean_noNeu_smallall_nodirt$cluster.sr[basal_3d_clean_noNeu_smallall_nodirt$cluster.sr=="Neutro_like"]="Neutrophil"
basal_3d_clean_noNeu_smallall_nodirt$cluster.sr[basal_3d_clean_noNeu_smallall_nodirt$cluster.sr=="Neutro_like2"]="Neutrophil"
#write.csv(basal_3d_clean_noNeu_smallall_nodirt,"/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/zac_Cyno_ifgn_term_srcluster120k.csv")



##############################################
########   finished prepping data  ###########
########   mario in python script  ###########
##############################################




#############################  part 2 ##################################

### this part is about the production of figures in fig3 for the manuscript

wcct=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/wcctHuman_mario_matched_sub20k.csv") # wcct human
zh=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/XhumanIFNG_mario_matched_sub20k.csv") # xspecies human
zr=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/XRhesusIFNG_mario_matched_sub20k.csv") # xspecies rhesus
zc=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/XCynoIFNG_mario_matched_sub20k.csv") # xspecies cyno
# general cca scores from MARIO
cca=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/Cytof-Xspecies-IFNG/gcca_sub20k.csv", header = TRUE)

## start production of the tsne plots
all_cca_scale_10=cca
set.seed(42)
all_tsne_scale_10=Rtsne(all_cca_scale_10[,-1], check_duplicates = FALSE, num_threads = 10) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000))) # lazy modality
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
# tsne with all four datasets
p = ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct-zac-integration with 10cca") +
  scale_color_manual(values =  clr )
# tsne with only wcct
label=as.character(wcct$cluster.sr)
p = ggplot(all_tsne_data_scale_10[c(1:20000),])  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0.1, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") + scale_color_brewer(palette = 'Set2' )
# tsne with only xspecies ifng human
label=as.character(zh$cluster.sr)
ggplot(all_tsne_data_scale_10[c(20001:40000),])  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") + scale_color_brewer(palette = 'Set2' )
# tsne with only xspecies ifng rhesus
label=as.character(zr$cluster.sr)
clr=c("#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
ggplot(all_tsne_data_scale_10[c(40001:60000),])  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") + scale_color_manual(values =clr )
# tsne with only xspecies ifng cyno
label=as.character(zc$cluster.sr)
ggplot(all_tsne_data_scale_10[c(60001:80000),])  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct- with 10cca") + scale_color_manual(values =clr )

## also need to produce the pre-integration plots
wcct1_value_scale=as.data.frame(scale(wcct[,c(3:43)] ))
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
wzhrc=rbind(wcct1_value_scale[,share3],zh_value_scale[,share3],zr_value_scale[,share3],zc_value_scale[,share3])
set.seed(42)
pre_tsne_scale_10=Rtsne(wzhrc[,c(1:39)], check_duplicates = FALSE, num_threads = 10) # run with cca
pre_tsne_data_scale_10 <- data.frame(x = pre_tsne_scale_10$Y[,1], y = pre_tsne_scale_10$Y[,2])
label=wzhrc$dataset
label=as.factor(c(rep("wcct-human",20000), rep("human-ifgn",20000), rep("rhesus-ifgn",20000), rep("cyno-ifgn",20000)))
clr=c("#f72585","#560bad","#4361ee","#9bf6ff")
p = ggplot(pre_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("wcct-zac-integration with 10cca") + scale_color_manual(values =  clr )

### start producing the feature plots

## first create the seurat objects
# wcct
wcct1_value=wcct[,c(3:43)] #input only the values
rownames(wcct1_value)=wcct$index
wcct1obj=CreateSeuratObject(counts=t(wcct1_value),assay="cytof")
SetAssayData(object = wcct1obj, slot = "data", new.data = t(wcct1_value), assay="cytof")
cca_embed=as.data.frame(all_tsne_data_scale_10[c(1:20000),]) # input the cca scores
colnames(cca_embed) <- paste0("ccatsne_", 1:2)
rownames(cca_embed)=wcct$index
wcct1obj[["ccatsne"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "ccatsne_", assay = "cytof")
wcct1obj@meta.data$anno.term=as.character(wcct$cluster.sr)
# xspecies ifng human
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
# xspecies ifng rhesus
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
# xspecies ifng zac
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

## now all seurat objects here, produce the feature plots
DefaultAssay(wcct1obj) <- 'cytof'
p=FeaturePlot(wcct1obj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/wcct_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(wcct1obj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/wcct_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(wcct1obj) <- 'cytof'
p=FeaturePlot(wcct1obj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/wcct_response3_refg.png", plot=p, width=7, height=7)


DefaultAssay(zhobj) <- 'cytof'
p=FeaturePlot(zhobj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zh_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(zhobj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zh_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(zhobj) <- 'cytof'
p=FeaturePlot(zhobj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zh_response3_refg.png", plot=p, width=7, height=7)


DefaultAssay(zrobj) <- 'cytof'
p=FeaturePlot(zrobj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zr_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(zrobj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zr_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(zrobj) <- 'cytof'
p=FeaturePlot(zrobj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zr_response3_refg.png", plot=p, width=7, height=7)


DefaultAssay(zcobj) <- 'cytof'
p=FeaturePlot(zcobj, features =  c("Ki67"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zc_response1_refg.png", plot=p, width=7, height=7)
p=FeaturePlot(zcobj, features =  c("STAT1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#6e1423"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zc_response2_refg.png", plot=p, width=7, height=7)
DefaultAssay(zcobj) <- 'cytof'
p=FeaturePlot(zcobj, features =  c("p38"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1, reduction = "ccatsne", pt.size = 1.2,cols=c("white","#641220"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/cytof_cytof/triology/tsne_plots/white_feature/zc_response3_refg.png", plot=p, width=7, height=7)

### now start to produce the corresponding boxplots

## normalize the values per dataset/ per column
wcct_scale = as.data.frame(scale(wcct[,c(3:43)]))
wcct_scale$cluster.sr = wcct$cluster.sr
wcct_scale$dataset = "wcct-human"
zh_scale = as.data.frame(scale(zh[,c(3:41)]))
zh_scale$cluster.sr = zh$cluster.sr
zh_scale$dataset = "IFNG-human"
zr_scale = as.data.frame(scale(zr[,c(3:41)]))
zr_scale$cluster.sr = zr$cluster.sr
zr_scale$dataset = "IFNG-rhesus"
zc_scale = as.data.frame(scale(zc[,c(3:41)]))
zc_scale$cluster.sr = zc$cluster.sr
zc_scale$dataset = "IFNG-crab-eating"
col_use = c("Ki67","STAT1","p38","cluster.sr","dataset")
all = rbind(wcct_scale[,col_use],zh_scale[,col_use])
all = rbind(all,zr_scale[,col_use])
all = rbind(all,zc_scale[,col_use])
# box plots 1
temp = subset(all, all$cluster.sr == "CD4 T")
temp$dataset <- factor(temp$dataset, levels = c("wcct-human", "IFNG-human", "IFNG-rhesus","IFNG-crab-eating"))
p <- ggplot(temp, aes(x=dataset, y=Ki67, fill = dataset)) + 
  geom_violin(lwd = 0.1) +
  scale_y_continuous(limits = quantile(temp$Ki67, c(0.05, 0.95))) + theme_classic() +
  theme(text = element_text(size=20))+
  geom_boxplot(outlier.shape = NA,lwd=0.1, width = 0.1, fill = "lightgrey")
# box plots 2
temp = subset(all, all$cluster.sr == "MC")
temp$dataset <- factor(temp$dataset, levels = c("wcct-human", "IFNG-human", "IFNG-rhesus","IFNG-crab-eating"))
p <- ggplot(temp, aes(x=dataset, y=STAT1, fill = dataset)) + 
  geom_violin(lwd = 0.1) +
  scale_y_continuous(limits = quantile(temp$Ki67, c(0.05, 0.95))) + theme_classic() +
  theme(text = element_text(size=20))+
  geom_boxplot(outlier.shape = NA,lwd=0.1, width = 0.1, fill = "lightgrey")
# box plot 3
temp = subset(all, all$cluster.sr == "MC")
temp$dataset <- factor(temp$dataset, levels = c("wcct-human", "IFNG-human", "IFNG-rhesus","IFNG-crab-eating"))
p <- ggplot(temp, aes(x=dataset, y=p38, fill = dataset)) + 
  geom_violin(lwd = 0.1) +
  scale_y_continuous(limits = quantile(temp$Ki67, c(0.05, 0.95))) + theme_classic() +
  theme(text = element_text(size=20))+
  geom_boxplot(outlier.shape = NA,lwd=0.1, width = 0.1, fill = "lightgrey")
# box plot 4
temp = subset(all, all$cluster.sr == "CD4 T")
temp$dataset <- factor(temp$dataset, levels = c("wcct-human", "IFNG-human", "IFNG-rhesus","IFNG-crab-eating"))
p <- ggplot(temp, aes(x=dataset, y=p38, fill = dataset)) + 
  geom_violin(lwd = 0.1) +
  scale_y_continuous(limits = quantile(temp$Ki67, c(0.05, 0.95))) + theme_classic() +
  theme(text = element_text(size=20))+
  geom_boxplot(outlier.shape = NA,lwd=0.1, width = 0.1, fill = "lightgrey")


### now start to produce distance distribution across datasets,  based on cca scores
# retrieve cca for each datasets
wcctcca=cca[c(1:20000),-1]
zhcca=cca[c(20001:40000),-1]
zrcca=cca[c(40001:60000),-1]
zccca=cca[c(60001:80000),-1]
# calculate the euclidean distance between pairs of cells
x=sqrt(apply((wcctcca-zhcca)^2,1,sum))
y=sqrt(apply((wcctcca-zrcca)^2,1,sum))
z=sqrt(apply((wcctcca-zccca)^2,1,sum))
# also calculate the distance between random cell pairs from human (wcct human vs wcct human)
randidnx=sample(20000, replace = FALSE)
wcctcca_random=wcctcca[randidnx,]
x_r=sqrt(apply((wcctcca-wcctcca_random)^2,1,sum))
# distance and ploting
df=data.frame(dist=c(x,y,z), data=c(rep("Human",20000),rep("Rhesus",20000),rep("Cyno",20000)))
p2=ggplot(df, aes(dist, fill = data, colour = data)) +geom_density(alpha = 0.3, size=0.5)+theme_classic() + xlim(0,4)
## with random
df=data.frame(dist=c(x,y,z,x_r), data=c(rep("Human",20000),rep("Rhesus",20000),rep("Cyno",20000),rep("Rand",20000)))
p2=ggplot(df, aes(dist, fill = data, colour = data)) +geom_density(alpha = 0.3, size=0.5)+theme_classic() + xlim(0,9)

### now start to producting the mathcing accuracy plots (bar plots)
# accuracy include wcct - human ; wcct - rhesus; wcct - cyno
wcct_zh=data.frame(label1=wcct$cluster.sr,label2=zh$cluster.sr)
wcct_zh$agree[wcct_zh$label1 == wcct_zh$label2]=1
wcct_zh$agree[wcct_zh$label1 != wcct_zh$label2]=0
wcct_zr=data.frame(label1=wcct$cluster.sr,label2=zr$cluster.sr) # we need to remove basophils
wcct_zr$label1=as.character(wcct_zr$label1)
wcct_zr=subset(wcct_zr,wcct_zr$label1!='Basophil')
wcct_zr$agree[wcct_zr$label1 == wcct_zr$label2]=1
wcct_zr$agree[wcct_zr$label1 != wcct_zr$label2]=0
wcct_zc=data.frame(label1=wcct$cluster.sr,label2=zc$cluster.sr) # we need to remove basophils
wcct_zc$label1=as.character(wcct_zc$label1)
wcct_zc=subset(wcct_zc,wcct_zc$label1!='Basophil')
wcct_zc$agree[wcct_zc$label1 == wcct_zc$label2]=1
wcct_zc$agree[wcct_zc$label1 != wcct_zc$label2]=0
# we also want to now how stable the accuracy could be so we subsample and calculate the accuracy

# loop:
library(caret)
cm <- confusionMatrix(factor(wcct_zc$label1),factor(wcct_zc$label2) , dnn = c("wcct","zh"))
list_ovacc=c()
n=5000
container=as.data.frame(cm$byClass)
container$batch=NA
container$type=NA
ipt=wcct_zh
loop=5
idx=0
for (i in c(1:loop)){
  idx=idx+1
  randidxx=sample(dim(ipt)[1], n)
  ipt_sub=ipt[randidxx,]
  cm=confusionMatrix(factor(ipt_sub$label1),factor(ipt_sub$label2) , dnn = c("lab1","lab2"))
  cm$byClass=as.data.frame(cm$byClass)
  cm$byClass$type=rownames(cm$byClass)
  cm$byClass$batch=idx
  container=rbind(container,cm$byClass)
  list_ovacc=c(list_ovacc,cm$overall[1])
}
container$match="wcct_zh"
container_zh=container
###### wcct-rhesus
cm <- confusionMatrix(factor(wcct_zc$label1),factor(wcct_zc$label2) , dnn = c("wcct","zh"))
container=as.data.frame(cm$byClass)
container$batch=NA
container$type=NA
ipt=wcct_zr
loop=5
idx=0
for (i in c(1:loop)){
  idx=idx+1
  randidxx=sample(dim(ipt)[1], n)
  ipt_sub=ipt[randidxx,]
  cm=confusionMatrix(factor(ipt_sub$label1),factor(ipt_sub$label2) , dnn = c("lab1","lab2"))
  cm$byClass=as.data.frame(cm$byClass)
  cm$byClass$type=rownames(cm$byClass)
  cm$byClass$batch=idx
  container=rbind(container,cm$byClass)
  list_ovacc=c(list_ovacc,cm$overall[1])
}
container$match="wcct_zr"
container_zr=container
######### wcct - cyno
cm <- confusionMatrix(factor(wcct_zc$label1),factor(wcct_zc$label2) , dnn = c("wcct","zh"))
container=as.data.frame(cm$byClass)
container$batch=NA
container$type=NA
ipt=wcct_zc
loop=5
idx=0
for (i in c(1:loop)){
  idx=idx+1
  randidxx=sample(dim(ipt)[1], n)
  ipt_sub=ipt[randidxx,]
  cm=confusionMatrix(factor(ipt_sub$label1),factor(ipt_sub$label2) , dnn = c("lab1","lab2"))
  cm$byClass=as.data.frame(cm$byClass)
  cm$byClass$type=rownames(cm$byClass)
  cm$byClass$batch=idx
  container=rbind(container,cm$byClass)
  list_ovacc=c(list_ovacc,cm$overall[1])
}
container$match="wcct_zc"
container_zc=container
all=rbind(container_zh[-c(1:6),],container_zr[-c(1:6),],container_zc[-c(1:6),])
## get all the accuracy quick summary
all2=all[,c(11:14)] # only the numbers
all2=subset(all2,all2$type!="Class: Basophil") # dont need to look at basophil since not in NHP
colnames(all2)=c("BalancedAccuracy","batch","type","match")
all2$match <- factor(all2$match, levels = c("wcct_zh", "wcct_zr", "wcct_zc"))
### small custome function for sumarization of the data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
## summary end
clr=c("#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
all3 <- data_summary(all2, varname="BalancedAccuracy", 
                     groupnames=c("match", "type"))
# now this produce the accuracy barplots in figure 3
p <- ggplot(data=all3, aes(x=match, y=BalancedAccuracy, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=(BalancedAccuracy-sd), ymax=(BalancedAccuracy+sd)) , width=.2,
                position=position_dodge(.9), size =0.1) +
  theme_minimal()+  scale_fill_manual(values = clr ) + coord_cartesian(ylim=c(0.5,1))



### Bokai Zhu
### 0822-2021
##### script related to the production of Figure S1.3 pre and post Mario analysis
##### (integrative analysis on 10x-5k pbmc and cytof pbmc felix et al cell report)

############################# part 1 ##################################

# part 1 is for preping the 10x and felix data before MARIO


###### prepping of 10x genomic pbmce download file location:
#Download location is https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3?
#from 10x genomics 5k pbmc health control
inputdata.10x.d1 <- Read10X_h5("/home/bkzhu/SNE-multi/cite-seq/5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5")
ab_counts <- inputdata.10x.d1$`Antibody Capture`[1:29,] # the last three abs are just control so disregard
# also save out the RNA info for downstream analysis
gene_counts <- inputdata.10x.d1$`Gene Expression` # the last three abs are just control
cellid = colnames(gene_counts)
gene_counts = as.data.frame(t(as.matrix(gene_counts)))
#write.csv(gene_counts, "/home/bkzhu/SNE-multi/10x-pbmc5k-cite-RNA.csv",quote = FALSE)

# Create Seurat object, for standard processing and clustering for the 10x data
data3 <- CreateSeuratObject(counts = ab_counts)
data3 <- NormalizeData(data3, assay = "RNA", normalization.method = "CLR")
data3 <- ScaleData(data3, assay = "RNA")
data3 <- RunPCA(data3, features = rownames(data3), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
                verbose = FALSE)
DimPlot(data3, reduction = "pca_adt")
ElbowPlot(data3, ndims = 50, reduction = "pca_adt")
# Now, we rerun tSNE using the PCA only on ADT (protein) levels.
data3 <- RunTSNE(data3, dims = 1:9, reduction = "pca_adt", reduction.key = "adtTSNE_", reduction.name = "tsne_adt")
data3 <- FindNeighbors(data3, features = rownames(data3), dims = NULL)
adt.data <- GetAssayData(data3, slot = "data")
adt.dist <- dist(t(adt.data))
data3[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
data3 <- FindClusters(data3, resolution = 0.2, graph.name = "adt_snn")
data3 <- RunTSNE(data3, dims = 1:25, method = "FIt-SNE", reduction = "pca_adt")
DimPlot(data3, label = TRUE) + NoLegend()
# check diff expression patterns for annotation
data3.small <- subset(data3, downsample = 1000)
# Find protein markers for all clusters, and draw a heatmap
adt.markers <- FindAllMarkers(data3.small, assay = "RNA", only.pos = TRUE)
DoHeatmap(data3.small, features = unique(adt.markers$gene), assay = "RNA", angle = 90) + NoLegend()
new.cluster.ids <- c("CD14 Mono","naive CD4 T","mem CD4 T","CD8 T","NK","B","DC","Other","Other")
names(new.cluster.ids) <- levels(data3)
data3 <- RenameIdents(data3, new.cluster.ids) # annotation of 10x pbmc created
counts=t(as.data.frame(as.matrix(data3@assays$RNA@counts)))
cluster.info=as.character(Idents(data3))
tenXdata=cbind(counts,cluster.info) # add the new annotation information to 10x dataset
colnames(tenXdata) = gsub(pattern = "-TotalSeqB", replacement = "", x = colnames(tenXdata)) # clean up feature names
# let call them dirt
genomic10x_sub = subset(tenXdata, tenXdata$cluster.info != "Other") # remove cell annotated as other
#write.csv(genomic10x_sub, "/home/bkzhu/SNE-multi/10x-pbmc5k-cite-clean.csv")
# this is the 10x 5k pbmc data used in matchting

#### start preping the cytof pbmc data (felix et al cell report)
library(flowCore)
library(FlowSOM)
fcs="/home/bkzhu/SNE-multi/cytof/cell_rerport/HD06_run1.fcs" # directly retrieved from the cell report paper flowrepository
fcsout = read.flowSet(files = fcs, transformation = F, truncate_max_range = F)
colnames(fcsout)=c("Time",markernames(fcsout))
expr_notrans2 = fsApply(fcsout, exprs)
norm_channel = "DNA1"
cutoff = 0.1
expr_notrans_filt2 = expr_notrans2[expr_notrans2[, norm_channel] > cutoff,drop=F] # remove cells with too little DNA not cells
df = expr_notrans_filt2
df = df/(median(df[,norm_channel])) # norm with DNA1 signal
df=as.data.frame(df)
df2=df[ , -which(names(df) %in% c('Time', 'Event_length', 'DNA1','DNA2',"empty","dead"))] # remove not used channles
df2=as.matrix(df2)
df2_raw=df2*(median(expr_notrans_filt2[,norm_channel])) 
cofactor = 5 # archsin transformation with cofactor 5
expr_notrans_filt_norm_median2 = asinh(df2/cofactor)
rng = colQuantiles(expr_notrans_filt_norm_median2, probs = c(0.005,0.995)) # remove low and high values
norm_rng = colQuantiles(expr_notrans_filt_norm_median2, probs = c(0, 1))
expr01_2 = t((t(expr_notrans_filt_norm_median2)-rng[,1])/(rng[,2]-rng[,1]))
expr01_2[expr01_2<0] = 0
expr01_2[expr01_2>1] = 1
raw_cytof=as.data.frame(t(df2_raw))
archsin_cytof=as.data.frame(t(expr01_2))
## random downsample the dataset to 50k cells
set.seed(1234)
sub_sample=sample(colnames(raw_cytof),50000)
raw_cytof_downsample=raw_cytof[,sub_sample]
archsin_cytof_downsample=archsin_cytof[,sub_sample]
# Create Seurat object and perform clustering and cell type annotation
data4 <- CreateSeuratObject(counts =raw_cytof_downsample)
data4 <- SetAssayData(object = data4, slot = "scale.data", new.data = as.matrix(archsin_cytof_downsample))
data4 <- SetAssayData(object = data4, slot = "data", new.data = as.matrix(archsin_cytof_downsample))
data4 <- RunPCA(data4, features = rownames(data4), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
DimPlot(data4, reduction = "pca_cytof")
ElbowPlot(data4, ndims = 50, reduction = "pca_cytof")
data4 <- RunTSNE(data4, dims = 1:7, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data4 <- FindNeighbors(data4, features = rownames(data4), dims = NULL)
cy.data <- GetAssayData(data4, slot = "data")
cy.dist <- dist(t(cy.data))
data4[["cy_snn"]] <- FindNeighbors(cy.dist)$snn
data4 <- FindClusters(data4, resolution = 0.2, graph.name = "cy_snn")
data4 <- RunTSNE(data4, dims = 1:25, method = "FIt-SNE", reduction = "pca_cytof",check_duplicates = FALSE)
DimPlot(data4, label = TRUE) + NoLegend()
# check diff expression pattern for manual annotation
data4.small <- subset(data4, downsample = 1000)
# Find protein markers for all clusters, and draw a heatmap
cy.markers <- FindAllMarkers(data4.small, assay = "RNA", only.pos = TRUE)
DoHeatmap(data4.small, features = unique(cy.markers$gene), assay = "RNA", angle = 90) + NoLegend()
new.cluster.ids <- c("naive CD4 T","B","CD8 T","mem CD4 T","Other","CD14 Mono", # manual annotation
                     "CD8 T","Other","DC","Other","Other","NK","DC","Other","Other","Other")
names(new.cluster.ids) <- levels(data4)
data4 <- RenameIdents(data4, new.cluster.ids)
felixdata=t(as.data.frame(as.matrix(raw_cytof_downsample)))
cluster.info=as.character(Idents(data4))
felix=cbind(felixdata,cluster.info)# add the clustering information to felix data
felix_sub = subset(felix,felix$cluster.info != "Other") # remove un-annotated cells
#write.csv(felix_sub, "/home/bkzhu/SNE-multi/sean-pbmc50k-cytof-clean.csv") # this is the input for mario


##############################################
########   finished prepping data  ###########
########   mario in python script  ###########
##############################################


########################  part 2 #############################

### part 2 include code producing figures2.3 figrures (tsne, feature plots, and mathcting accuracy)
# two modality protein files row matched by mario
tenx = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/pbmc10x_mario_matched.csv")
felix = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/pbmcCytof_mario_matched.csv")
# the reduced cca scores produced by mario
all_cca_scale_10 = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/pbmc/10xgenomic_felixetal_cca10.csv")
library(Rtsne)
set.seed(42)
all_tsne_scale_10=Rtsne(all_cca_scale_10[,-1], check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("cite",4606), rep("cytof",4606))) # lazy input of modality information
# produce the tsne plots based on madality
p = ggplot(all_tsne_data_scale_10)  + 
  geom_point(aes(x=x, y=y, color=label,size = 15,stroke = 1, alpha =0.001), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")
label=c(as.character(tenx$cluster.info), as.character(felix$cluster.info)) # annotation from citeseq and cytof clustering info
# produce tsne plot based on two dataset annotation
p=ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 15,stroke = 1, alpha =0.001), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")


## no start to produce the featureplots with citeseq protein/rna and cytof protein
# start to retrieve the rna counts for the matched 10x cells
rna = gene_counts
cite_id = as.character(tenx$X) # barcode
rna_matched = rna[cite_id,] # find barcode in rna
rna_matched = rna_matched[,-1] # delete barcode
rownames(rna_matched) = cite_id
## make seurat compatiable version, citeseq stack cytof cell
rna_matched_dum = rna_matched
rownames(rna_matched_dum) = rownames(felix)
rna_matched_comp = rbind(rna_matched,rna_matched_dum)
## cite
rownames(tenx) = tenx$X
tenx_dum = tenx
rownames(tenx_dum) = rownames(felix)
tenx_comp = rbind(tenx,tenx_dum)
## cytof
felix_dum = felix
rownames(felix_dum) = rownames(tenx)
felix_comp = rbind(felix_dum,felix)
library(Seurat)
# create seurat object
pbmc_obj=CreateSeuratObject(counts=t(rna_matched_comp),assay="RNA")
SetAssayData(object = pbmc_obj, slot = "data", new.data = t(rna_matched_comp), assay="RNA")
pbmc_obj <- NormalizeData(object = pbmc_obj, assay = "RNA")
pbmc_obj[["CITE"]]=CreateAssayObject(counts = t(tenx_comp[,c(3:31)])) # only input the protein values
SetAssayData(object = pbmc_obj, slot = "data", new.data = t(tenx_comp[,c(3:31)]), assay="CITE")
pbmc_obj <- ScaleData(pbmc_obj, assay = "CITE")
pbmc_obj[["CYTOF"]]=CreateAssayObject(counts = t(felix_comp[,c(3:32)])) # only input the protein values
SetAssayData(object = pbmc_obj, slot = "data", new.data = t(felix_comp[,c(3:32)]), assay="CYTOF")
pbmc_obj <- ScaleData(pbmc_obj, assay = "CYTOF")
cca_embed=all_tsne_data_scale_10
# load the reduced mario inetgrated cca scores to seurat object
colnames(cca_embed) <- paste0("int_", 1:2)
rownames(cca_embed)=rownames(felix_comp)
pbmc_obj[["int"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "int_", assay = "RNA")
pbmc_obj@meta.data$tech=c(rep("cite",4606),rep("cyto",4606)) # lazy input of modality
DimPlot(pbmc_obj, reduction = "int", pt.size = 0.5,group.by = "tech", label=TRUE)
# featureplot 1
#pbmc cite proteins
DefaultAssay(pbmc_obj) <- 'CITE'
p=FeaturePlot(pbmc_obj,
              features =  c("CD8a","CD4","CD19","CD14","CD45RO","CD45RA","CD27"),
              min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "int",
              pt.size = 1.5,cols=c("lightgrey","#05668D"))
#pbmc cytof proteins
DefaultAssay(pbmc_obj) <- 'CYTOF'
p=FeaturePlot(pbmc_obj,
              features =  c("CD8a","CD4","CD19","CD14","CD45RO","CD45RA","CD27"),
              min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "int",
              pt.size = 1.5,cols=c("lightgrey","#05668D"))
#pbmc citeseq rna
DefaultAssay(pbmc_obj) <- 'RNA'
p=FeaturePlot(pbmc_obj,
              features =  c("CD8A","IL7R","MS4A1","LYZ","S100A4","CCR7","GNLY"),
              min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "int",
              pt.size = 1.5,cols=c("lightgrey","#05668D"))


#### now start to make confusion matrix with balanced accuracy
library(caret)
library(dplyr)
temp = data.frame(tenx$cluster.info , felix$cluster.info) # input matched cells annotations
cm <- confusionMatrix(factor(temp$felix.cluster.info), factor(temp$tenx.cluster.info) , dnn = c( "Cytof","Cite-Seq"))
table <- data.frame(cm$table)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Cite.Seq == table$Cytof, "good", "bad")) %>%
  group_by(Cite.Seq) %>%
  mutate(prop = Freq/sum(Freq))
# produce the confustion matrix
p = ggplot(data = plotTable, mapping = aes(x = Cite.Seq, y = Cytof, fill = goodbad, alpha = prop)) +
  geom_tile() +
  theme_classic() +
  xlim(rev(levels(table$Cite.Seq)))+ theme(text = element_text(size=20))
# the balanced accuracy is retrieved from the cm object








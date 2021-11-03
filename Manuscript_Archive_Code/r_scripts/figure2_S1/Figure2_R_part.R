### Bokai Zhu
### 2021-0821
### script related to the production parts of Figure 2 (tsne plots and featureplots)
### and production of FigureS1.1 

##################### part 1 ########################

### script related to data preping for bmcite and levine 32
## bone marrow from bmcite, seurat package
library(SeuratData)
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
bm.cite=bm@assays$ADT@counts
rownames(bm.cite)[25]="HLA-DR" # change name matching
bm.cluster2=bm@meta.data$celltype.l2 # use l2 annotation from dataset
bm.cluster2=as.character(bm.cluster2)
# new col: merge clusters to lower level, for comparing with levine32
bm.cluster2[bm.cluster2 =="pDC"] <- "pDCs" # good
bm.cluster2[bm.cluster2 =="Prog_RBC"] <- "HSPCs" # good
bm.cluster2[bm.cluster2 =="gdT"] <- "toss" # need to toss
bm.cluster2[bm.cluster2 =="CD4 Naive"] <- "CD4 T" # good
bm.cluster2[bm.cluster2 =="CD4 Memory"] <- "CD4 T" # good
bm.cluster2[bm.cluster2 =="CD14 Mono"] <- "monocyte" # good
bm.cluster2[bm.cluster2 =="Naive B"] <- "B" # good
bm.cluster2[bm.cluster2 =="CD8 Naive"] <- "CD8 T" # good
bm.cluster2[bm.cluster2 =="CD8 Effector_2"] <- "CD8 T" # good
bm.cluster2[bm.cluster2 =="GMP"] <- "HSPCs"
bm.cluster2[bm.cluster2 =="CD8 Effector_1"] <- "CD8 T" # good
bm.cluster2[bm.cluster2 =="CD16 Mono"] <- "monocyte" # good
bm.cluster2[bm.cluster2 =="CD8 Memory_1"] <- "CD8 T" # good
bm.cluster2[bm.cluster2 =="MAIT"] <- "toss" # toss
bm.cluster2[bm.cluster2 =="Memory B"] <- "B" # good
bm.cluster2[bm.cluster2 =="cDC2"] <- "toss" # toss
bm.cluster2[bm.cluster2 =="CD56 bright NK"] <- "NK" # good
bm.cluster2[bm.cluster2 =="Treg"] <- "toss" # toss
bm.cluster2[bm.cluster2 =="Prog_B 2"] <- "B" # good
bm.cluster2[bm.cluster2 =="Prog_Mk"] <- "HSPCs" # good
bm.cluster2[bm.cluster2 =="CD8 Memory_2"] <- "CD8 T" # good
bm.cluster2[bm.cluster2 =="Plasmablast"] <- "B" # good
bm.cluster2[bm.cluster2 =="HSC"] <- "HSPCs" # good
bm.cluster2[bm.cluster2 =="LMPP"] <- "HSPCs" 
bm.cluster2[bm.cluster2 =="Prog_DC"] <- "HSPCs" # toss
bm.cluster2[bm.cluster2 =="Prog_B 1"] <- "B" # good
value=t(bm@assays$ADT@data) # seurat transformed
value=as.data.frame(as.matrix(value))
value$cluster.orig=unlist(bm@meta.data$celltype.l2)# org form data
value$cluster.term=bm.cluster2# org form data
value=subset(value, value$cluster.term!="toss") # get rid of irrelevant cell types
write.csv(raw, "/home/bkzhu/SNE-multi/figure_rcode/fig3/bmcite_29000.csv")



###### bone marrow cytof data from levine32 (levine et al 2015 cell)
library(flowCore)
fcs2="/home/bkzhu/SNE-multi/cytof/Levine_32dim.fcs" # archsine transformed from paper
fcsout2 = read.flowSet(files = fcs2, transformation = F, truncate_max_range = F) # read
expr_notrans3 = as.data.frame(fsApply(fcsout2, exprs))
expr_notrans3 = expr_notrans3[, -which(colnames(expr_notrans3) %in% c("Time","Cell_length","DNA1","DNA2","Viability",
                                                                      "file_number","event_number","individual")  )] # remove unused cols
gated=subset(expr_notrans3, !is.na(expr_notrans3$label))
from=c("Basophils","CD16-_NK_cells","CD16+_NK_cells","CD34+CD38+CD123-_HSPCs",
       "CD34+CD38+CD123+_HSPCs","CD34+CD38lo_HSCs","CD4_T_cells","CD8_T_cells",
       "Mature_B_cells","Monocytes","pDCs","Plasma_B_cells","Pre_B_cells","Pro_B_cells")
idx=0
for (item in from){
  idx=idx+1
  gated$label[gated$label==idx]=item
}
colnames(gated)[33]="cluster.orig"
### new col: merged clusters for comparison with bmcite
origclust=gated$cluster.orig
levine.cluster2=as.character(origclust)
levine.cluster2[levine.cluster2 =="CD16-_NK_cells"] <- "NK" # good
levine.cluster2[levine.cluster2 =="CD16+_NK_cells"] <- "NK" # good
levine.cluster2[levine.cluster2 =="CD34+CD38+CD123-_HSPCs"] <- "HSPCs" # good
levine.cluster2[levine.cluster2 =="CD34+CD38+CD123+_HSPCs"] <- "HSPCs" # good
levine.cluster2[levine.cluster2 =="CD34+CD38lo_HSCs"] <- "HSPCs" # good
levine.cluster2[levine.cluster2 =="CD4_T_cells"] <- "CD4 T" # good
levine.cluster2[levine.cluster2 =="CD8_T_cells"] <- "CD8 T" # good
levine.cluster2[levine.cluster2 =="Mature_B_cells"] <- "B" # good
levine.cluster2[levine.cluster2 =="Monocytes"] <- "monocyte" # good
levine.cluster2[levine.cluster2 =="Plasma_B_cells"] <- "B" # good
levine.cluster2[levine.cluster2 =="Pre_B_cells"] <- "B" # good
levine.cluster2[levine.cluster2 =="Pro_B_cells"] <- "B" # good
gated$cluster.term=levine.cluster2
gated=subset(gated, gated$cluster.term!="Basophils") # exclude Basophils
write.csv(gated, "/home/bkzhu/SNE-multi/figure_rcode/fig3/levine32_102977.csv")


#####################################
#### MARIO performed in python ######
#####################################


##################### part 2 #######################

### this part is the post integration analysis related to fig2 and figs2.1
library(Rtsne)
# read in files
# integrated cca scores
all_cca_scale_10 = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/bmcite_lv32_cca10.csv")
# matched files with original expression/ annotations etc
sxc_cite=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/bmcite_mario_matched.csv") # bmcite
sxc_lv=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/bonemarrow/lv32_mario_matched.csv") # levine32
########## start producing tsne plots
set.seed(42)
all_tsne_scale_10=Rtsne(all_cca_scale_10[,-1], check_duplicates = FALSE, num_threads = 10) # run with cca
all_tsne_data_scale_10 <- data.frame(x = all_tsne_scale_10$Y[,1], y = all_tsne_scale_10$Y[,2])
label=as.factor(c(rep("cite",26822), rep("cytof",26822))) # lazy input of modality information
# produce tsne colored by modality
p=ggplot(all_tsne_data_scale_10)  + 
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")
# color based on independent annotations from their own
label=c(as.character(sxc_cite$cluster.term), as.character(sxc_lv$cluster.term))
p=ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")
# color based on levine 32 annotations
label=c(as.character(sxc_lv$cluster.orig), as.character(sxc_lv$cluster.orig))
color=c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C" ,"#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928","#989898")
p=ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")+ scale_color_manual(values =color )
# color based on bmcite annotations
label=c(as.character(sxc_cite$cluster.orig), as.character(sxc_cite$cluster.orig))
#color=c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C" ,"#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928","#989898")
p=ggplot(all_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca")#+ scale_color_manual(values =color )
## also produce pre-integration tsne plot
# just need to use original data and perform tsne
cite_value=sxc_cite[,c(3:27)] # only get the values
lv_value=sxc_lv[,c(3:34)] # only get the values
cite_value_scale=scale(cite_value) # norm the values
lv_value_scale=scale(lv_value)
shared=intersect(colnames(cite_value_scale),colnames(lv_value_scale)) # only used shared features
cite_lv=rbind(cite_value_scale[,shared],lv_value_scale[,shared])
set.seed(42)
all_tsne_scale_10_pre=Rtsne(cite_lv, check_duplicates = FALSE, num_threads = 5) # run with cca
all_tsne_data_scale_10_pre <- data.frame(x = all_tsne_scale_10_pre$Y[,1], y = all_tsne_scale_10_pre$Y[,2])
label=as.factor(c(rep("cite",26822), rep("cytof",26822))) #lazy input
p = ggplot(all_tsne_data_scale_10_pre)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with all shared pre")

## tsne production finisehd


### cluster based on the integrated cca scores and color by cca clustering
# now we want to do clustering directly on cca values
# make seurat object
library(Seurat)
bm_id=sxc_cite$Unnamed..0 # id for bmcite
lv_id=sxc_lv$Unnamed..0 # id for lv32
lv_id_nonDup=make.names(as.character(lv_id), unique = TRUE) # make sure name no duplication
bm_rna=bm@assays$RNA@data
bm_rna_match=as.matrix(bm_rna[,as.character(bm_id)]) # get the matched cells' rna information
bm_rna_match=as.data.frame(t(bm_rna_match))
bm_rna_match_dumy=bm_rna_match
rownames(bm_rna_match_dumy)=lv_id_nonDup # for cytof , remake the row names (id) so the names do not duplicate
inte_rna=rbind(bm_rna_match,bm_rna_match_dumy)
# get citeseq pro
bm_cite=bm@assays$ADT@data
bm_cite_match=as.matrix(bm_cite[,as.character(bm_id)])
bm_cite_match=as.data.frame(t(bm_cite_match))
bm_cite_match_dumy=bm_cite_match
rownames(bm_cite_match_dumy)=as.character(lv_id_nonDup) # for y dataset
inte_proCITE=rbind(bm_cite_match,bm_cite_match_dumy)
# get cytof pro
lev32_value=sxc_lv[,c(2:33)]
lev32_match_dumy=lev32_value
rownames(lev32_value)=as.character(lv_id_nonDup)
rownames(lev32_match_dumy)=as.character(bm_id) # for x dataset
inte_proCYTOF=rbind(lev32_match_dumy,lev32_value)
# seurat object
inte_obj=CreateSeuratObject(counts=t(inte_rna),assay="RNA")
SetAssayData(object = inte_obj, slot = "data", new.data = t(inte_rna), assay="RNA")
inte_obj[["CITE"]]=CreateAssayObject(counts = t(inte_proCITE))
SetAssayData(object = inte_obj, slot = "data", new.data = t(inte_proCITE), assay="CITE")
inte_obj <- ScaleData(inte_obj, assay = "CITE")
inte_obj[["CYTOF"]]=CreateAssayObject(counts = t(inte_proCYTOF))
SetAssayData(object = inte_obj, slot = "data", new.data = t(inte_proCYTOF), assay="CYTOF")
inte_obj <- ScaleData(inte_obj, assay = "CYTOF")
# put pre computed tsne corrdinates into the seurat object
cca_embed=all_tsne_data_scale_10
colnames(cca_embed) <- paste0("int_", 1:2)
rownames(cca_embed)=rownames(inte_proCYTOF)
inte_obj[["int"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "int_", assay = "RNA")
inte_obj@meta.data$anno.bm=c(as.character(sxc_cite$cluster.orig),as.character(sxc_cite$cluster.orig))
inte_obj@meta.data$anno.lv=c(as.character(sxc_lv$cluster.orig),as.character(sxc_lv$cluster.orig))
inte_obj@meta.data$anno.term=c(as.character(sxc_cite$cluster.term),as.character(sxc_cite$cluster.term))
inte_obj@meta.data$tech=c(rep("cite",26822),rep("cyto",26822))# input the modality information
# input the cca scores into the seurat object
cca_embed=all_cca_scale_10[,-1]
colnames(cca_embed) <- paste0("cca_", 1:10)
rownames(cca_embed)=rownames(inte_proCYTOF)
inte_obj[["cca"]] <- CreateDimReducObject(embeddings = as.matrix(cca_embed), key = "cca_", assay = "RNA")
# perform clustering based on cca scores
inte_obj <- FindNeighbors(inte_obj, features = rownames(inte_obj),reduction='cca')
inte_obj <- FindClusters(inte_obj, reduction = "cca", resolution = 0.3, dims.use = 1:10)
inte_obj.c10 <- subset(inte_obj, idents = 10)
inte_obj.c10 <- FindNeighbors(inte_obj.c10, features = rownames(inte_obj.c10),reduction='cca')
inte_obj.c10 <- FindClusters(inte_obj.c10, reduction = "cca", resolution = 0.1, dims.use = 1:10)
inte_obj$sub_cluster <- as.character(Idents(inte_obj))
inte_obj$sub_cluster[Cells(inte_obj.c10)] <- paste("c10_",Idents(inte_obj.c10))
inte_obj_store=Idents(inte_obj)
Idents(inte_obj)=inte_obj$sub_cluster
inte_obj.c9 <- subset(inte_obj, idents = "9")
inte_obj.c9 <- FindNeighbors(inte_obj.c9, features = rownames(inte_obj.c9),reduction='cca')
inte_obj.c9 <- FindClusters(inte_obj.c9, reduction = "cca", resolution = 0.15, dims.use = 1:10)
inte_obj$sub_cluster <- as.character(Idents(inte_obj))
inte_obj$sub_cluster[Cells(inte_obj.c9)] <- paste("c9_",Idents(inte_obj.c9))
inte_obj_store2=Idents(inte_obj)
Idents(inte_obj)=inte_obj$sub_cluster
# produce tsne plots and color based on cca clustering results
label=Idents(inte_obj)
color=c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C" ,
        "#FDBF6F" ,"#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99", "#B15928",
        "#989898","#FA26A0","#3D5B59")
p=ggplot(all_tsne_data_scale_10)  + 
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0, alpha =0.001), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("integration with 10cca") + scale_color_manual(values =color )
## prodece feature maps for protein/RNA expression
# produce citeseq protein feature map
DefaultAssay(inte_obj) <- 'CITE'
p=FeaturePlot(inte_obj, features =  c("CD8a","CD4","CD19","CD14","CD45RO","CD28","CD56"),
              min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "int",
              pt.size = 1.5,cols=c("lightgrey","#05668D"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/fig3/tsne_figs/citepro_pred-test.svg", plot=p, width=49, height=7)
# produce cytof protein expression
DefaultAssay(inte_obj) <- 'CYTOF'
p=FeaturePlot(inte_obj, features =  c("CD8","CD4","CD19","CD14","CD11b","CD20","CD33"),
              min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "int",
              pt.size = 1.5,cols=c("lightgrey","#05668D"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/fig3/tsne_figs/cytopro_pred.svg", plot=p, width=49, height=7)
# produce rna expression from citeseq
DefaultAssay(inte_obj) <- 'RNA'
p=FeaturePlot(inte_obj, features =  c("CD8A","IL7R","MS4A1","LYZ","S100A4","CCR7","GNLY"),
              min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "int",
              pt.size = 1.5,cols=c("lightgrey","#05668D"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/fig3/tsne_figs/citeRNA_pred.svg", plot=p, width=49, height=7)

##### skip the feature plot production code for figS2.1, same exact code but with differnt features

## start produce the violin plots for figs2.1 RNA counts

# quick normalization function as Seurat, used here only for visualization
norm_log=function(df, scale_fac=10000){
  dfd=df/rowSums(df)
  dfd=dfd*scale_fac
  dfdl=log1p(dfd)
  return(dfdl)
}
bm_rna_match_dumy_lognorm=norm_log(bm_rna_match_dumy) # squash furthur the counts for visualization
bm_rna_match_dumy_lognorm$cluster.term=as.character(sxc_cite$cluster.term) # add cell type information
cyto_rna_ms4a1=bm_rna_match_dumy_lognorm[,c("MS4A1","LYZ","cluster.term")] #
cyto_rna_ms4a1_melt=melt(cyto_rna_ms4a1)
p = ggplot(cyto_rna_ms4a1_melt, aes(x=cluster.term, y=value)) + # rna violin plot1
  geom_violin(aes(fill = factor(variable)),lwd=0.1) + theme_bw()
# the second rna violin plot
cyto_rna_ms4a1=bm_rna_match_dumy_lognorm[,c("IL7R","GNLY","cluster.term")]
cyto_rna_ms4a1_melt=melt(cyto_rna_ms4a1)
p = ggplot(cyto_rna_ms4a1_melt, aes(x=cluster.term, y=value)) + # rna violin plot2
  geom_violin(aes(fill = factor(variable)))+ theme_bw()


## start producting the confusion matrix with the balanced accuracy
bm_full=sxc_cite # get the matched file for cell annotation matching accuracy
lv_full=sxc_lv # matched file for levine 32
library(caret)
cm <- confusionMatrix(factor(lv_full$cluster.term),factor(bm_full$cluster.term) , dnn = c("Cytof", "Cite-Seq"))
table <- data.frame(cm$table)

plotTable <- table %>%
  mutate(goodbad = ifelse(table$Cite.Seq == table$Cytof, "good", "bad")) %>%
  group_by(Cytof) %>%
  mutate(prop = Freq/sum(Freq))
p = ggplot(data = plotTable, mapping = aes(x = Cite.Seq, y = Cytof, fill = goodbad, alpha = prop)) + # the confusion plot
  geom_tile() +
  #geom_text(aes(label = Freq), vjust = .5, alpha = 0.5) +
  theme_classic() +
  xlim(rev(levels(table$Cite.Seq))) + theme(text = element_text(size=20))
# and the balanced accuracy was calculated by package caret, extracted from the cm object


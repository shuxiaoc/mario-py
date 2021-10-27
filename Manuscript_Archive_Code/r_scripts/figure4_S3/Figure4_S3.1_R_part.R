### Bokai Zhu
### 0824-2021
##### script related to the production of Figure 4 and FigureS3.1 pre and post Mario analysis
##### (integrative analysis on codex1 murine cell paper data + citeseq murine spleen from totalvi paper)

#################################### part 1 #######################################

######### part 1 is preprocessing and clustering/annotation of the codex / citeseq cells

## part 1.1 is about codex data clean up and clustering
## data retrieved directly from Yury Goltsev et al, tiff files of the balbc1 murine codex imaging

## data was first re-segmented via Msemer whole cell segmentation as described in the material and methods
## Then lateral spill-over and aggregation of B220 signal was corrected via REDSEA

## then this script start with reading the single cell extracted, complemented and aggreation removed fcs file of codex balbc 1 data.

library(flowCore) # for fcs reading
## this file is a fcs file with single cell extracted, complemented and aggreation removed fcs file of codex balbc 1 data
fcs1=read.flowSet("/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/yh_red1_noah/dataRedSeaScaleSizeFCS_all63_cCmp_b220fix_s7.fcs")
expr1 = as.data.frame(fsApply(fcs1, exprs))
expr1 = subset(expr1,expr1$nucl>80) # filter cell with small nucl value removing non-cell segmentations
expr1 = as.matrix(expr1)
rng = colQuantiles(expr1[, 3:(dim(expr1)[2]-1)], probs = c(0.005,0.995)) # cap value between 0.5% - > 99.5%
expr1_sc = t((t(expr1[, 3:(dim(expr1)[2]-1)])-rng[,1])/(rng[,2]-rng[,1])) # norm
expr1_sc[expr1_sc<0] = 0
expr1_sc[expr1_sc>1] = 1 # capped values
# cluster the codex cells based on these markers
cluster_channel=c("CD45","Ly6C","TCR","Ly6G","CD19",
                  "CD169","CD3","CD8a","F480","CD11c",
                  "CD27","CD31","CD4","IgM","B220","ERTR7",
                  "MHCII","CD35","CD2135","NKp46","CD1632",
                  "CD90","CD5","CD79b","IgD","CD11b","CD106")
lineage_markers = expr1_sc[,cluster_channel]
tm=as.data.frame(t(lineage_markers))
data6 <- CreateSeuratObject(counts =tm)
data6 <- SetAssayData(object = data6, slot = "scale.data", new.data = as.matrix(tm))
data6 <- SetAssayData(object = data6, slot = "data", new.data = as.matrix(tm))
data6 <- RunPCA(data6, features = rownames(data6), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
DimPlot(data6, reduction = "pca_cytof")
ElbowPlot(data6, ndims = 20, reduction = "pca_cytof")
data6 <- RunTSNE(data6, dims = 1:15, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data6 <- FindNeighbors(data6, features = rownames(data6), dims = NULL)
cy.data <- GetAssayData(data6, slot = "data")
cy.dist <- dist(t(cy.data))
data6[["cy_snn"]] <- FindNeighbors(cy.dist)$snn
data6 <- FindClusters(data6, resolution = 1.5, graph.name = "cy_snn")
# manual annotation
new.cluster.ids <- c("CD4 T","CD8 T","B","B","CD4 T",
                     "B","Mph","empty","DC","gc B",
                     "Mph","B","Mph","fibro","mz Mph",
                     "B","NK","endo","mix","empty",
                     "dirt","empty","mono","fibro",
                     "naive B","neutrophil","blood vessel",
                     "empty","other T","empty")
names(new.cluster.ids) <- levels(data6)
idstore=Idents(data6)
data6 <- RenameIdents(data6, new.cluster.ids)
# clean up the annotations
fdf=as.data.frame(expr1)
fdf$cluster.sr=as.character(Idents(data6))
fdf$cluster.term=fdf$cluster.sr
ann=as.character(fdf$cluster.term)
ann[ann=="endo"]<-"Other"
ann[ann=="gc B"]<-"B"
ann[ann=="fibro"]<-"Other"
ann[ann=="neutrophil"]<-"Neutrophils"
ann[ann=="naive B"]<-"B"
ann[ann=="mz Mph"]<-"Mph"
ann[ann=="blood vessel"]<-"Other"
ann[ann=="other T"]<-"Other"
fdf$cluster.term=ann
# due to factors related to spatial location and antibody staining
# DC/ Macrophage need to be reclustered individually due to strong b cell signal containmination:
codex_mphr1=subset(fdf,fdf$cluster.term=="Mph")
codex_DC1=subset(fdf,fdf$cluster.term=="DC")
codex_mhdcr1=rbind(codex_mphr1,codex_DC1)
# cluster only based on dc/mph related channels
cluster_channelr2=c("CD45","CD19","CD169","CD8a","F480","CD11c",
                    "CD4","IgM","B220","MHCII","CD35","CD2135",
                    "CD1632","CD106","IgD")
lineage_markersr2 = codex_mhdcr1[,cluster_channelr2]
lineage_markersr2=scale(lineage_markersr2)
tm2=as.data.frame(t(lineage_markersr2))
data7 <- CreateSeuratObject(counts =tm2)
data7 <- SetAssayData(object = data7, slot = "scale.data", new.data = as.matrix(tm2))
data7 <- SetAssayData(object = data7, slot = "data", new.data = as.matrix(tm2))
data7 <- RunPCA(data7, features = rownames(data7), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
data7 <- RunTSNE(data7, dims = 1:10, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data7 <- FindNeighbors(data7, features = rownames(data7), dims = NULL)
cy.data7 <- GetAssayData(data7, slot = "data")
cy.dist7 <- dist(t(cy.data7))
data7[["cy_snn"]] <- FindNeighbors(cy.dist7)$snn
data7 <- FindClusters(data7, resolution = 1, graph.name = "cy_snn")
new.cluster.ids <- c("mix","Mph","DC","B mph mix","dirt","Other","Mph","B mph mix","B","B","DC","Mph","empty","dirt","B")
names(new.cluster.ids) <- levels(data7)
idstore=Idents(data7)
data7 <- RenameIdents(data7, new.cluster.ids)
lineage_markersr2=as.data.frame(lineage_markersr2)
lineage_markersr2$cluster.sr=as.character(Idents(data7))
# furthur cleanup the naming of the clusters
fdf2=fdf
fdf2$cluster.sr[fdf2$cluster.term=="Mph"]=lineage_markersr2$cluster.sr[1:11067]
fdf2$cluster.sr[fdf2$cluster.term=="DC"]=lineage_markersr2$cluster.sr[11068:14626]
ann=as.character(fdf2$cluster.sr)
ann[ann=="mz Mph"]<-"Mph"
ann[ann=="B mph mix"]<-"mix"
ann[ann=="mix"]<-"mix"
ann[ann=="Endothelial"]<-"Other"
ann[ann=="gc B"]<-"B"
ann[ann=="fibroblast"]<-"Other"
ann[ann=="Neutrophil"]<-"Neutrophils"
ann[ann=="Naive B"]<-"B"
ann[ann=="mz Mph"]<-"Mph"
ann[ann=="endo"]<-"Other"
ann[ann=="gc B"]<-"B"
ann[ann=="fibro"]<-"Other"
ann[ann=="neutrophil"]<-"Neutrophils"
ann[ann=="naive B"]<-"B"
ann[ann=="mz Mph"]<-"Mph"
ann[ann=="blood vessel"]<-"Other"
ann[ann=="other T"]<-"Other"
fdf2$cluster.term=ann
# caution: strong aggregation of NK marker in fov34, remove them
fdf2$cluster.term[fdf2$cluster.term=="NK" & fdf2$PointNum==34]="dirt"
# save out all the cells for future ploting need
#write.csv(fdf2,"/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/forMatch_noahCodexComp_b220Fix/c")
# disregard the dirt cells for downstream analysis
fdf2_clean=fdf2
fdf2_clean=subset(fdf2_clean,fdf2_clean$cluster.term!="dirt")
fdf2_clean=subset(fdf2_clean,fdf2_clean$cluster.term!="empty")
fdf2_clean=subset(fdf2_clean,fdf2_clean$cluster.term!="mix")
fdf2_clean=subset(fdf2_clean,fdf2_clean$cluster.term!="Other")
#write.csv(fdf2_clean,"/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/forMatch_noahCodexComp_b220Fix/clean_cells.csv")


### next we prep the murine spleen citeseq dataset 
## For citeseq data mouse spleen health comes from https://www.nature.com/articles/s41592-020-01050-x#data-availability
## github location: https://github.com/YosefLab/totalVI_reproducibility
## data retrieving via python skipped very straightforward

cite_murine_meta=read.csv("/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/cite_data/Mouse206-meta.csv") # meta data retrieved from github
cite_murine_pro=read.csv("/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/cite_data/Mouse206-pro.csv") # protein data retreived from github
cite_murine_meta$cellID=cite_murine_pro$cellID
# only use the cell from spleen
cite_murine_meta=subset(cite_murine_meta,cite_murine_meta$sample=="Spleen")
cite_murine_pro2=cite_murine_pro[cite_murine_pro$cellID %in% cite_murine_meta$cellID,]
cite_murine_rna2=cite_murine_rna[cite_murine_rna$cellID %in% cite_murine_meta$cellID,]
cite_murine_pro2$cluster.orig=cite_murine_meta$cell_type
# dont need to recluster cell
# just directly the annotated information from totalvi paper
#clean cite murine spleen annotation names
cite_murine_meta$cell_type_clean=NA
cite_murine_ann=as.character(cite_murine_meta$cell_type)
cite_murine_ann[cite_murine_ann =="Activated CD4 T"] <- "CD4 T"
cite_murine_ann[cite_murine_ann =="B1 B"] <- "B"
cite_murine_ann[cite_murine_ann =="B-CD4 T cell doublets"] <- "Mixed"
cite_murine_ann[cite_murine_ann =="B-CD8 T cell doublets"] <- "Mixed"
cite_murine_ann[cite_murine_ann =="B doublets"] <- "B"
cite_murine_ann[cite_murine_ann =="B-macrophage doublets"] <- "Mixed"
cite_murine_ann[cite_murine_ann =="CD122+ CD8 T"] <- "CD8 T"
cite_murine_ann[cite_murine_ann =="CD4 T"] <- "CD4 T"
cite_murine_ann[cite_murine_ann =="CD8 T"] <- "CD8 T"
cite_murine_ann[cite_murine_ann =="cDC1s"] <- "DC"
cite_murine_ann[cite_murine_ann =="cDC2s"] <- "DC"
cite_murine_ann[cite_murine_ann =="Cycling B/T cells"] <- "Cycling"
cite_murine_ann[cite_murine_ann =="Erythrocytes"] <- "Erythrocytes"
cite_murine_ann[cite_murine_ann =="GD T"] <- "otherT"
cite_murine_ann[cite_murine_ann =="ICOS-high Tregs"] <- "otherT"
cite_murine_ann[cite_murine_ann =="Ifit3-high B"] <- "B"
cite_murine_ann[cite_murine_ann =="Ifit3-high CD4 T"] <- "CD4 T"
cite_murine_ann[cite_murine_ann =="Ifit3-high CD8 T"] <- "CD8 T"
cite_murine_ann[cite_murine_ann =="Low quality B cells"] <- "B"
cite_murine_ann[cite_murine_ann =="Low quality T cells"] <- "otherT"
cite_murine_ann[cite_murine_ann =="Ly6-high mono"] <- "Mono"
cite_murine_ann[cite_murine_ann =="Ly6-low mono"] <- "Mono"
cite_murine_ann[cite_murine_ann =="Mature B"] <- "B"
cite_murine_ann[cite_murine_ann =="Migratory DCs"] <- "DC"
cite_murine_ann[cite_murine_ann =="MZ/Marco-high macrophages"] <- "Mph"
cite_murine_ann[cite_murine_ann =="Neutrophils"] <- "Neutrophils"
cite_murine_ann[cite_murine_ann =="NK"] <- "NK"
cite_murine_ann[cite_murine_ann =="MZ B"] <- "B"
cite_murine_ann[cite_murine_ann =="NKT"] <- "NKT"
cite_murine_ann[cite_murine_ann =="pDCs"] <- "DC"
cite_murine_ann[cite_murine_ann =="Plasma B"] <- "B"
cite_murine_ann[cite_murine_ann =="Red-pulp macrophages"] <- "Mph"
cite_murine_ann[cite_murine_ann =="T doublets"] <- "Mixed"
cite_murine_ann[cite_murine_ann =="Transitional B"] <- "B"
cite_murine_ann[cite_murine_ann =="Tregs"] <- "otherT"
cite_murine_pro2$cluster.info=cite_murine_ann
# change citeseq murine feature names clean up
cite_proName=gsub(x = colnames(cite_murine_pro2), pattern = "_.*", replacement = "")
cite_proName2=gsub(x = cite_proName, pattern = "/.*", replacement = "") 
colnames(cite_murine_pro2)=cite_proName2
# remove irrelevant cell types
cite_murine_pro2=cite_murine_pro2[cite_murine_pro2$cluster.info!="Mixed",]
ite_murine_pro2=cite_murine_pro2[cite_murine_pro2$cluster.info!="Erythrocytes",]
cite_murine_pro2=cite_murine_pro2[cite_murine_pro2$cluster.info!="NKT",]
cite_murine_pro2=cite_murine_pro2[cite_murine_pro2$cluster.info!="otherT",]
cite_murine_pro2=cite_murine_pro2[cite_murine_pro2$cluster.info!="Cycling",]
# duplicate the dataset to 2x to maintain the cell size for matching (rare cell types)


##############################################
########   finished prepping data  ###########
########   mario in python script  ###########
##############################################


#############################  part 2 ##################################


### this part is about the production of figures in fig4 for the manuscript

## read in the cca scores and matched files
cca_codex=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/cca_codex.csv", header = TRUE)
cca_cite=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/cca_cite.csv", header = TRUE)
codex_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCodex_mario_matched.csv")
cite_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCiteseq_mario_matched.csv")

## first we produce the tsne plots
set.seed(42)
# subsample 5k cells for visulization, otherwise unfair for citeseq
all_cca = rbind(cca_codex,cca_cite)
randidx=sample(dim(all_cca)[1]/2,8000, replace = FALSE)
all_cca_codex=all_cca[c(1:36673),-1] # lazy
all_cca_cite=all_cca[c(36674:73346),-1] # lazy
all_cca_codex_sub=all_cca_codex[randidx,] # random 5k
all_cca_cite_sub=all_cca_cite[randidx,]
all_cca_sub=rbind(all_cca_codex_sub,all_cca_cite_sub)
# tsne plot production
library(Rtsne)
all_tsne_sub=Rtsne(all_cca_sub[,c(1:10)], check_duplicates = FALSE, num_threads=10) # run with cca
all_tsne_data_sub<- data.frame(x = all_tsne_sub$Y[,1], y = all_tsne_sub$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000))) # lazy modality
# tsne plot with modality information
p = ggplot(all_tsne_data_sub)  +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 1.2), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("integration with 10cca")
# plot tsne plot with individual dataset annotation
cite_matched$cluster.info=as.character(cite_matched$cluster.info)
label=c(as.character(codex_matched$cluster.term[randidx]),as.character(cite_matched$cluster.info[randidx]))
# tsne plot with annotations
p = ggplot(all_tsne_data_sub) +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 0.8, alpha =1), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("integration with 10cca")
## start producing tsne plot pre-integration
codex_input_rand = codex_matched[randidx,] # get the same random cells
cite_input_rand = cite_matched[randidx,]
shared=intersect(colnames(codex_input_rand), colnames(cite_input_rand))
codex_input_rand_value=codex_input_rand[,shared[c(2:24)]] # only use shared features
codex_input_rand_value=scale(codex_input_rand_value)
cite_input_rand_value=cite_input_rand[,shared[c(2:24)]]
cite_input_rand_value=scale(cite_input_rand_value)
rand_all=rbind(codex_input_rand_value,cite_input_rand_value)
set.seed(42)
pre_tsne_scale_10=Rtsne(rand_all, check_duplicates = FALSE, num_threads = 10) # run with cca
pre_tsne_data_scale_10 <- data.frame(x = pre_tsne_scale_10$Y[,1], y = pre_tsne_scale_10$Y[,2])
label=as.factor(c(rep("codex",8000), rep("cite",8000))) # lazy
p = ggplot(pre_tsne_data_scale_10)  +
  geom_point(aes(x=x, y=y,color=label,size = 10,stroke = 0.8, alpha =1), cex = 0.3) +
  labs( x = "tsne 1", y = "tsne 2") + theme_bw() + ggtitle("pre integration") 


## next produce the codex annotation figures

## first job is to plot the bablc 1 spleen and color based on annotations from codex dataset alone
## this if figure 4 annotation
library(ggplot2)

## annotation plot start
codex_file=read.csv("../data/murine_spleen/all_clusters.csv") # file that contains annotation for all cells (include non-annotated)
for (i in c(1:63)){ # 63 fovs
  # read in the label map
  pathlabel=paste0("../data/murine_spleen/msemer-seg/point",i,"/nuclCD45_1.6_mpp0.8/point",i,"_feature_0.tif",sep='')
  labelmap=readTIFF(pathlabel)
  labelmap_melt=melt(labelmap)
  # read segementation from mesemer
  segpath=paste0("../data/murine_spleen/msemer-seg/point",i,"/nuclCD45_1.6_mpp0.8/point",i,"_outline.tif",sep='')
  segmap=readTIFF(segpath)
  segmap_melt=melt(segmap)
  segmap_melt$value = 1 - segmap_melt$value # reverse 0 1
  # read annotation and change the level
  labelmap_melt$value=labelmap_melt$value*segmap_melt$value # make boundary 0
  # joint the labels
  codex_filesub=subset(codex_file,codex_file$PointNum==i)
  temp2=left_join(labelmap_melt,codex_filesub, by=c("value"="cellLabelInImage"))
  temp2$cluster.term=as.character(temp2$cluster.term)
  temp2$cluster.term[is.na(temp2$cluster.term)] ="not cell"
  temp2$cluster.term <- factor(temp2$cluster.term,
                               levels = c("not cell", "dirt", "empty",
                                          "Other","mix","B","CD4 T","CD8 T",
                                          "DC","Mono","Mph","Neutrophils","NK"))
  # change colors corresponding to cell types
  clr_pl=c("black","black","grey","grey","grey","#771122","#77AADD","#114477","#44AA77","#CC99BB","#DDDD77","#117744","#774411")
  names(clr_pl) <- c("not cell", "dirt", "empty","Other","mix","B","CD4 T","CD8 T","DC","mono","Mph","Neutrophils","NK")
  colScale <- scale_fill_manual(name = "grp",values = clr_pl)
  # change direction for y for ploting
  temp2$Var1=abs(temp2$Var1-dim(labelmap)[1])
  # make ploot
  temp_plot = ggplot(temp2, aes(x = Var2, y = Var1, fill = cluster.term)) +
    geom_tile() +
    colScale +
    labs(x=NULL, y=NULL, title=NULL) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() +
    coord_fixed() +
    labs(x=NULL, y=NULL, title=NULL) +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position='none',
          panel.background = element_blank(),
          panel.grid = element_blank(),
          title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(rep(-1.25,4),'null'),
          plot.margin = unit(rep(-1.25,4),'lines'),
          axis.ticks.length = unit(0,'cm'),
          axis.ticks.margin = unit(0,'cm'),
          panel.border = element_rect(colour = 'black', fill=NA, size=1))
  # save the images
  core_path="/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/yh_red1_noah/codex_gd/" # output path please ignore
  tiff(file = paste0(core_path,"allcodexcells_cluster_terrm", i, ".tiff"), width = dim(labelmap)[2], height = dim(labelmap)[1])
  print(temp_plot)
  dev.off()
}
## annotation plot ended

## next color the same plot with the knn label transfer annotation from the citeseq dataset

## start label transfer plotting
clr_pl=c("black","grey","#771122","#77AADD","#114477","#44AA77","#CC99BB","#DDDD77","#117744","#774411")
names(clr_pl) <- c("not cell","not matched","B","CD4 T","CD8 T","DC","Mono","Mph","Neutrophils","NK")
codex_input=read.csv("./data/murine_spleen/MurineCodex_mario_labeltrans.csv") # file with label transfer annotations
codex_file=read.csv("../data/murine_spleen/all_clusters.csv") # file that contains annotation for all cells (include non-annotated)
#####
for (i in c(1:63)){
  # read cell labels
  pathlabel=paste0("../data/murine_spleen/msemer-seg/point",i,"/nuclCD45_1.6_mpp0.8/point",i,"_feature_0.tif",sep='')
  labelmap=readTIFF(pathlabel)
  labelmap_melt=melt(labelmap)
  # read segementation 
  segpath=paste0("../data/murine_spleen/msemer-seg/point",i,"/nuclCD45_1.6_mpp0.8/point",i,"_outline.tif",sep='')
  segmap=readTIFF(segpath)
  segmap_melt=melt(segmap)
  segmap_melt$value = 1 - segmap_melt$value # reverse 0 1
  # read annotation and change the level
  labelmap_melt$value=labelmap_melt$value*segmap_melt$value # make boundary 0
  # joint the labels
  codex_filesub=subset(codex_file,codex_file$PointNum==i)
  temp2=left_join(labelmap_melt,codex_filesub, by=c("value"="cellLabelInImage"))
  temp2$cluster.term=as.character(temp2$cluster.term)
  temp2$cluster.term[is.na(temp2$cluster.term)] ="not cell"
  temp2$cluster.term[temp2$cluster.term=="dirt"] ="not cell"
  # add matching info post integration with knn
  codex_filesubMatch=subset(codex_input,codex_input$PointNum==i)
  codex_filesubMatch=codex_filesubMatch[,c("cellLabelInImage","cluster.term","cluster.sr","cite_cluster2")] # get citeseq labels
  colnames(codex_filesubMatch)=c("cellLabelInImage","match.term","match.sr","knn.type") # rename to avoid duplication
  temp3=left_join(temp2,codex_filesubMatch, by=c("value"="cellLabelInImage"))
  temp3$knn.type=as.character(temp3$knn.type)
  temp3$knn.type[temp3$cluster.term=='not cell'] ="not cell"
  temp3$knn.type[is.na(temp3$knn.type)] ="not matched"
  # level annotation
  temp3$knn.type <- factor(temp3$knn.type,
                           levels =  c("not cell","not matched","B","CD4 T","CD8 T","DC","Mono","Mph","Neutrophils","NK"))
  colScale <- scale_fill_manual(name = "grp",values = clr_pl)
  # change direction for y
  temp3$Var1=abs(temp3$Var1-dim(labelmap)[1])
  # make ploot
  temp_plot = ggplot(temp3, aes(x = Var2, y = Var1, fill = knn.type)) +
    geom_tile() +
    colScale +
    labs(x=NULL, y=NULL, title=NULL) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() +
    coord_fixed() +
    labs(x=NULL, y=NULL, title=NULL) +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position='none',
          panel.background = element_blank(),
          panel.grid = element_blank(),
          title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(rep(-1.25,4),'null'),
          plot.margin = unit(rep(-1.25,4),'lines'),
          axis.ticks.length = unit(0,'cm'),
          axis.ticks.margin = unit(0,'cm'),
          panel.border = element_rect(colour = 'black', fill=NA, size=1))
  # save the images
  core_path="/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/yh_red1_noah/labeltrans-randinput/"
  tiff(file = paste0(core_path,"allcodexcells_cluster_term", i, ".tiff"), width = dim(labelmap)[2], height = dim(labelmap)[1])
  print(temp_plot)
  dev.off()
}
## label transfer plotting ended

## next is to produce the rna plotting
## only give the example to produce one RNA rest are all the same
library(foreach)
library(doParallel)
registerDoParallel(10) # fast production with multiple cores
#### input
codex_input=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCodex_mario_matched.csv") # matched codex cells
codex_file=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/all_clusters.csv") # file with all cell annotation
cite_input=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCiteseq_mario_matched.csv") # matched citeseq cells
#### also get the corresponding rna information from the citeseq cells
cite_rna_all=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/Mouse206-rna.csv") # rna counts retrieved from totalvi paper
cite_rna_all_matched=subset(cite_rna_all, cite_rna_all$cellID %in% cite_input$cellID) # rna counts for matched citeseq cells
#### simplyfiy
codex_file=codex_file[,c("cellLabelInImage","PointNum","cluster.sr","cluster.term")]
### prep citeseq info  #######
rna_target="Il7r"
rna=c("cellID","Il7r","Ms4a1","Cxcr5","Myc")
protein=c("cellID","cluster.info")
rna_name=rna[c(2:length(rna))]# not choosing cellID
cite_rna_all_matched[,2:13554]=cite_rna_all_matched[,2:13554]/rowSums(cite_rna_all_matched[,2:13554])# relative counts
cite_rna_all_matched[,2:13554]=cite_rna_all_matched[,2:13554]*10000 # scale factor is 1e4, actually does not matter
# quantile the target
target_count=cite_rna_all_matched[,rna_target]
target_count0=target_count[target_count>0] # quantile of > 0
target_cap=quantile(target_count0,0.8) # cap at 0.8
cite_rna_all_matched[cite_rna_all_matched[,rna_target]>target_cap,rna_target]=target_cap
cite_input_pro_rna=left_join(cite_input[,protein],cite_rna_all_matched[,c("cellID",rna_target)], by=c("cellID"="cellID"))

### for loop start for rna level on spleen
for (j in c(1:7)){
  result = foreach( i = 1:9, .combine=rbind) %dopar% {
    #for (i in c(1:9)){
    ptnum=(j-1)*9+i
    print(ptnum)
    pathlabel=paste0("../data/murine_spleen/msemer-seg/point",ptnum,"/nuclCD45_1.6_mpp0.8/point",ptnum,"_feature_0.tif",sep='')
    labelmap=readTIFF(pathlabel)
    labelmap_melt=melt(labelmap)
    # read segementation 
    segpath=paste0("../data/murine_spleen/msemer-seg/point",ptnum,"/nuclCD45_1.6_mpp0.8/point",ptnum,"_outline.tif",sep='')
    segmap=readTIFF(segpath)
    segmap_melt=melt(segmap)
    segmap_melt$value = 1 - segmap_melt$value # reverse 0 1
    # read annotation and change the level
    labelmap_melt$value=labelmap_melt$value*segmap_melt$value # make boundary 0
    # joint the labels
    codex_filesub=subset(codex_file,codex_file$PointNum==ptnum)
    temp2=left_join(labelmap_melt,codex_filesub, by=c("value"="cellLabelInImage"))
    temp2$cluster.term=as.character(temp2$cluster.term)
    temp2$cluster.term[is.na(temp2$cluster.term)] ="not cell"
    # add matching info, adding from codex basic info
    codex_filesubMatch=subset(codex_input,codex_input$PointNum==ptnum)
    codex_filesubMatch$rowid=rownames(codex_filesubMatch)
    codex_filesubMatch=codex_filesubMatch[,c("cellLabelInImage","cluster.term","cluster.sr","rowid")]
    colnames(codex_filesubMatch)=c("cellLabelInImage","match.term","match.sr","rowid")
    temp3=left_join(temp2,codex_filesubMatch, by=c("value"="cellLabelInImage"))
    temp3$match.term=as.character(temp3$match.term)
    temp3$match.term[temp3$cluster.term=='not cell'] ="not cell"
    temp3$match.term[is.na(temp3$match.term)] ="not matched"
    # join citeseq data with rna and protein depends
    cite_input_2=cite_input_pro_rna
    cite_input_2$rowid=rownames(cite_input_2)
    temp4=left_join(temp3,cite_input_2, by=c("rowid"="rowid"))
    # change direction for y
    temp4$Var1=abs(temp4$Var1-dim(labelmap)[1])
    # add tile pixels
    temp4$Var2=temp4$Var2+(7-j)*1344
    temp4$Var1=temp4$Var1+(9-i)*1008
    temp4
  }
  # save each column
  ##### make plot ########
  temp_plot = ggplot(result, aes(x = Var2, y = Var1, fill = Il7r)) + # plotting Il7r now
    geom_tile() +
    #colScale +
    scale_fill_gradient(low = "black", high = "#ECF87F", na.value = "black")+
    labs(x=NULL, y=NULL, title=NULL) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() +
    coord_fixed() +
    labs(x=NULL, y=NULL, title=NULL) +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position='none',
          panel.background = element_blank(),
          panel.grid = element_blank(),
          title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(rep(-1.25,4),'null'),
          plot.margin = unit(rep(-1.25,4),'lines'),
          axis.ticks.length = unit(0,'cm'),
          axis.ticks.margin = unit(0,'cm'),
          panel.border = element_rect(colour = 'black', fill=NA, size=1))
  # save the images
  core_path="/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/yh_red1_noah/RCcap_rna/Il7r/"
  tiff(file = paste0(core_path,"Il7r_", j, ".tiff"), width = dim(labelmap)[2], height = dim(labelmap)[1]*9)
  print(temp_plot)
  dev.off()
}
### for loop end for rna level on spleen


### next part start producing the b cell subtyping RNA and spatial

## first step is to focus on b cell and do reclustering based on b cell markers

# read in again the matched cells 
codex_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCodex_mario_matched.csv")
cite_matched=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/MurineCiteseq_mario_matched.csv")
# get the corresponding rna information
# get rna info
cite_rna_all=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/Mouse206-rna.csv") # same rna file from totalvi previous
cite_rna_all_matched=subset(cite_rna_all, cite_rna_all$cellID %in% cite_matched$cellID) # rna for matched citeseq cells
cite_input_pro_rna=left_join(cite_matched,cite_rna_all_matched, by=c("cellID"="cellID")) #  citeseq simplified info
cite_input_rna=cite_input_pro_rna[,-c(1:216)] # file only contains rna
# one thing to notice, matching did not use cd45 but we want to use it during clustering
# so we need to retrieved it from the original file as MARIO input
codex_orig = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/murine_spleen/clean_cells_neutroflip.csv")
codex_matched_allchannel=left_join(codex_matched,codex_orig, by =c("Unnamed..0.1"="X"))
#### b cell subtypeing with codex proteins
codex_b=subset(codex_matched_allchannel,codex_matched_allchannel$cluster.term.x=="B") # subset matched b cell
# b cell related markers
bmarkers=c("CD19.x","CD1632.x","IgD.x","CD27.x","CD79b.x","IgM.x","B220.x","MHCII.x","CD35.x","CD2135.x","CD45")
codex_b_matched_pro=codex_b[,bmarkers]
# get codex b cell matched citeseq rna
nums <- unlist(lapply(cite_input_rna, is.numeric)) 
nums_cite_input_rna=cite_input_rna[,nums]
nums_cite_input_rna_ann=cbind(nums_cite_input_rna,ct=codex_matched$cluster.term)
codex_b_matched_rna=subset(nums_cite_input_rna_ann,nums_cite_input_rna_ann$ct == "B")
codex_b_matched_rna$ct <-NULL
# then make the b cells into a seurat object
rownames(codex_b_matched_pro)=as.character(rownames(codex_b_matched_pro))
rownames(codex_b_matched_rna)=as.character(rownames(codex_b_matched_pro))
### seurat object
# input matched rna
codexB_obj=CreateSeuratObject(counts=t(codex_b_matched_rna),assay="RNA") # this is raw counts
codexB_obj <- NormalizeData(object = codexB_obj, assay="RNA")
codexB_obj <- ScaleData(codexB_obj,assay="RNA")
# input codex with b cell markers
codexB_obj[["codex"]]=CreateAssayObject(counts = t(codex_b_matched_pro))
SetAssayData(object = codexB_obj, slot = "data", new.data = t(codex_b_matched_pro), assay="codex")
codexB_obj=SetAssayData(object = codexB_obj, slot = "scale.data", new.data = t(codex_b_matched_pro), assay="codex")
DefaultAssay(codexB_obj) <- "codex"
codexB_obj <- RunPCA(codexB_obj, features = rownames(codexB_obj), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                     verbose = FALSE)
DimPlot(codexB_obj, reduction = "pca_cytof")
ElbowPlot(codexB_obj, ndims = 10, reduction = "pca_cytof")
## start clustering based on codex b cell markers
codexB_obj <- RunTSNE(codexB_obj, dims = 1:5, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
codexB_obj <- FindNeighbors(codexB_obj, features = rownames(codexB_obj), dims = NULL)
cy.codexB_obj <- GetAssayData(codexB_obj, slot = "data")
cy.codexB_obj <- dist(t(cy.codexB_obj))
codexB_obj[["cy_snn"]] <- FindNeighbors(cy.codexB_obj)$snn
codexB_obj <- FindClusters(codexB_obj, resolution = 0.2, graph.name = "cy_snn") # 0.2
DimPlot(codexB_obj, label = TRUE) + NoLegend()
## manual annotation of the b cell subtypes
new.cluster.ids <- c("Naive","Transitional","Mature","GC","dirt")
names(new.cluster.ids) <- levels(codexB_obj)
idstore=Idents(codexB_obj)
codexB_obj <- RenameIdents(codexB_obj, new.cluster.ids)

## produce the ridge plot presented in figure S4.1
p = RidgePlot(codexB_obj, features = bmarkers, ncol = 5)

### produce the heatmap presented in figure 4 for b cell subtype
codexB_obj.small <- subset(codexB_obj, downsample = 5000)
cy.markers <- FindAllMarkers(codexB_obj.small, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.25)
# dont need to show the DE genes for cluster 4 (dirt)
all_bgenes=cy.markers$gene
all_bgenes_sub=all_bgenes[1:156] # genes after 156 is DE for dirt
p = DoHeatmap(codexB_obj.small, features = all_bgenes_sub, assay = "RNA", angle = 90)

## then produce the spatial balbc1 spleen colored with sub types of b cells
# input needed for plotting
codex_b_plot = codex_b[,c("cellLabelInImage","PointNum")]
codex_b_plot$sub = Idents(codexB_obj)
colnames(codex_b_plot) = c("cellLabelInImage","PointNum","sub")
library(tiff)
clr_pl=c("black","#424242","#ff537f","#94bf00" ,"#43ff4c","#54daff","#424242")
names(clr_pl) <- c("not cell","not matched","Naive","Transitional","Mature","GC","dirt")
codex_input = codex_b_plot
codex_file=read.csv("../data/murine_spleen/all_clusters.csv") 
## b cell subtype spatial plot start
for (i in c(1:63)){ # intotal 63 fovs
  # read cell labels
  pathlabel=paste0("../data/murine_spleen/msemer-seg/point",i,"/nuclCD45_1.6_mpp0.8/point",i,"_feature_0.tif",sep='')
  labelmap=readTIFF(pathlabel)
  labelmap_melt=melt(labelmap)
  # read segementation 
  segpath=paste0("../data/murine_spleen/msemer-seg/point",i,"/nuclCD45_1.6_mpp0.8/point",i,"_outline.tif",sep='')
  segmap=readTIFF(segpath)
  segmap_melt=melt(segmap)
  segmap_melt$value = 1 - segmap_melt$value # reverse 0 1
  # read annotation and change the level
  labelmap_melt$value=labelmap_melt$value*segmap_melt$value # make boundary 0
  # joint the labels
  codex_filesub=subset(codex_file,codex_file$PointNum==i)
  temp2=left_join(labelmap_melt,codex_filesub, by=c("value"="cellLabelInImage"))
  temp2$cluster.term=as.character(temp2$cluster.term)
  temp2$cluster.term[is.na(temp2$cluster.term)] ="not cell"
  temp2$cluster.term[temp2$cluster.term=="dirt"] ="not cell"
  codex_filesubMatch=subset(codex_input,codex_input$PointNum==i) # merging with the files with b cell subtype information
  temp3=left_join(temp2,codex_filesubMatch, by=c("value"="cellLabelInImage"))
  temp3$sub=as.character(temp3$sub) # adding the b cell subtypes
  temp3$sub[temp3$cluster.term=='not cell'] ="not cell"
  temp3$sub[is.na(temp3$sub)] ="not matched"
  # level annotation
  temp3$sub <- factor(temp3$sub,
                      levels =  c("not cell","not matched","Naive","Transitional","Mature","GC","dirt"))
  colScale <- scale_fill_manual(name = "grp",values = clr_pl)
  # change direction for y
  temp3$Var1=abs(temp3$Var1-dim(labelmap)[1])
  # make ploot
  temp_plot = ggplot(temp3, aes(x = Var2, y = Var1, fill = sub)) +
    geom_tile() +
    colScale +
    labs(x=NULL, y=NULL, title=NULL) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() +
    coord_fixed() +
    labs(x=NULL, y=NULL, title=NULL) +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position='none',
          panel.background = element_blank(),
          panel.grid = element_blank(),
          title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(rep(-1.25,4),'null'),
          plot.margin = unit(rep(-1.25,4),'lines'),
          axis.ticks.length = unit(0,'cm'),
          axis.ticks.margin = unit(0,'cm'),
          panel.border = element_rect(colour = 'black', fill=NA, size=1))
  # save the images
  core_path="/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/yh_red1_noah/bsub2/"
  tiff(file = paste0(core_path,"bsub_image", i, ".tiff"), width = dim(labelmap)[2], height = dim(labelmap)[1])
  print(temp_plot)
  dev.off()
}
## b cell subtype spatial plot end

### now start to product the TSNE based on codex with matched RNA overlaying onto it
### this figure is presented in figure S4.1

## first perform tsne only on codex proteins
set.seed(42)
all_tsne_cdxaln=Rtsne(codex_matched[,c(3:27)], check_duplicates = FALSE, num_threads=10) # tsne on codex markers
all_tsne_cdxaln_df<- data.frame(x = all_tsne_cdxaln$Y[,1], y = all_tsne_cdxaln$Y[,2])
label=as.character(codex_matched$cluster.term)
p = ggplot(all_tsne_cdxaln_df) +
  geom_point(aes(x=x, y=y, color=label,size = 30,stroke = 1.2), cex = 1) +
  labs( x = "tsne 1", y = "tsne 2") + theme_classic() + ggtitle("tsne with codex protein alone")
## create seurat object with codex protein generated tsne coordinations + RNA information
# seurat object
rownames(cite_input_rna)=as.character(rownames(cite_input_rna))
rownames(codex_matched)=as.character(rownames(codex_matched))
codex_obj=CreateSeuratObject(counts=t(cite_input_rna),assay="RNA") # this is raw counts, as previously prepared
SetAssayData(object = codex_obj, slot = "data", new.data = t(cite_input_rna), assay="RNA")
codex_obj <- NormalizeData(object = codex_obj, assay="RNA", normalization.method="CLR")
# codex protein input, not actually used in the plotting
codex_obj[["codex"]]=CreateAssayObject(counts = t(codex_matched[,c(3:27)]))
SetAssayData(object = codex_obj, slot = "data", new.data = t(codex_matched[,c(3:27)]), assay="codex")
# input the tsne coordinates into the seurat object
tsne_embed=all_tsne_cdxaln_df
colnames(tsne_embed) <- paste0("tsne_", 1:2)
rownames(tsne_embed)=rownames(codex_matched)
codex_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne_embed), key = "tsne_", assay = "RNA")
codex_obj@meta.data$anno.cdx=as.character(codex_matched$cluster.term)
DimPlot(codex_obj, reduction = "tsne", pt.size = 0.5,group.by = "anno.cdx", label=TRUE)
# rna feature plots row1
DefaultAssay(codex_obj) = "RNA"
p = FeaturePlot(codex_obj, features =  c("Cd8a","Cd4","Il7r","Ms4a1","Flt3","Sepp1","Klrd1"),
                min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "tsne",
                pt.size = 0.3,cols=c("lightgrey","#05668D"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/sup_related/plots/codex_feature_update1.png", plot=p, width=49, height=7)
# protein features from codex
DefaultAssay(codex_obj) = "codex"
p = FeaturePlot(codex_obj, features =  c("CD8a","CD4","CD3","B220","CD11c","F480","NKp46"),
                min.cutoff = "q05", max.cutoff = "q95", ncol = 7, reduction = "tsne",
                pt.size = 0.3,cols=c("lightgrey","#05668D"))
#ggsave(file="/home/bkzhu/SNE-multi/figure_rcode/mouse_cite/sup_related/plots/codex_feature_update2.png", plot=p, width=49, height=7)

### then the last part of the figure S4.1 is the confustion plot with accuracy information

library(caret)
cm <- confusionMatrix(factor(cite_matched$cluster.info),factor(codex_matched$cluster.term) , dnn = c("citeseq", "codex"))
table <- data.frame(cm$table)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$codex == table$citeseq, "good", "bad")) %>%
  group_by(codex) %>%
  mutate(prop = Freq/sum(Freq))
p = ggplot(data = plotTable, mapping = aes(x = codex, y = cite, fill = goodbad, alpha = prop)) +
  geom_tile() +
  theme_classic() +
  xlim(rev(levels(table$codex))) + theme(text = element_text(size=20))
## and the balanced accuracy can be retreived from the cm caret object


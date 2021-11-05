### Bokai Zhu
### 2021-0902
### script related to the production parts of Figure 5 
### COVID lung CODEX imaging and COVID BALF CITE-seq data




##################### part 1 ########################

## to be noticed @0902, the covid Lung CODEX imaging is not public yet
## need to update when avaliable publically

#### first we prep the COVID lung CODEX imaging data, where we input the csv file with ready-to-go CODEX information
############## these files are not avalible to public ##################
codextma1=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/codex/TMA1_complete.csv")
codextma2=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/codex/TMA2_complete.csv")
codextma3=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/codex/TMA3_complete.csv")
codextma4=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/codex/TMA4_complete.csv")
codexann=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/codex/TMA_codex_clusters_annotable_final.csv", sep = "\t")
# the tma contains non-covid samples too, need to use only covid samples
codextma1=subset(codextma1,codextma1$region %in% c(14,16,5,7,15)) # select based on tma core position
codextma2=subset(codextma2,codextma2$region %in% c(4,5,8,9,10,13,14,15,16,18,21,26,3,7,12,20,19,17,23,24,25,22,30,29,28))
codextma3=subset(codextma3,codextma3$region %in% c(1,2,4,6,7,8,9,11,12,13,14,17,18,19,20,22,25,27,3,16,21,24,26,28,29))
codextma4=subset(codextma4,codextma4$region %in% c(6,7,8,11,12,16,17,19,21,22,24,26,27,29,3,5,10,14,18,23,28))
# note the 'region' column is only unique per tma, so create a new column with unique values
library(plyr)
codextma1$region_new <- mapvalues(codextma1$region, 
                                  from=c(14,16,5,7,15), to=c(1,2,47:49))
codextma2$region_new <- mapvalues(codextma2$region, 
                                  from=c(4,5,8,9,10,13,14,15,16,18,21,26,3,7,12,20,19,17,23,24,25,22,30,29,28), to=c(3:14,50:62))
codextma3$region_new <- mapvalues(codextma3$region, 
                                  from=c(1,2,4,6,7,8,9,11,12,13,14,17,18,19,20,22,25,27,3,16,21,24,26,28,29), to=c(15:32,63:69))
codextma4$region_new <- mapvalues(codextma4$region, 
                                  from=c(6,7,8,11,12,16,17,19,21,22,24,26,27,29,3,5,10,14,18,23,28), to=c(33:46,70:76))
# add cell type annotations from Yury
codextma1$annotation <- codexann$final.annotation[match(codextma1$ClusterID, codexann$ClusterID)]
codextma2$annotation <- codexann$final.annotation[match(codextma2$ClusterID, codexann$ClusterID)]
codextma3$annotation <- codexann$final.annotation[match(codextma3$ClusterID, codexann$ClusterID)]
codextma4$annotation <- codexann$final.annotation[match(codextma4$ClusterID, codexann$ClusterID)]
codextma1$tma=1
codextma2$tma=2
codextma3$tma=3
codextma4$tma=4
codextmaall=rbind(codextma1,codextma2)
codextmaall=rbind(codextmaall,codextma3)
codextmaall=rbind(codextmaall,codextma4)
# clean up naming system
label=codextmaall$annotation
label=gsub(".*CD4 T.*", "CD4+ T-cell", label, perl=TRUE) #  bin some annootations
label=gsub(".*CD8 T.*", "CD8+ T-cell", label, perl=TRUE)
label=gsub(".*mphs.*", "Macrophage", label, perl=TRUE)
label=gsub(".*gdT cells*", "gd T-cell", label, perl=TRUE)
label=gsub(".*neutrophils*", "Neutrophil", label, perl=TRUE)
label=gsub(".*plasma cells*", "Plasma cell", label, perl=TRUE)
label=gsub(".*mphs.*", "mphs", label, perl=TRUE)
label=gsub(".*B cells*", "B cell", label, perl=TRUE)
codextmaall$cluster.term=label
# here we are only interested in using the macrophages in the dataset
codextmaall_mph=subset(codextmaall,codextmaall$cluster.term == "Macrophage")
#write.csv(codextmaall,"/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/CODEX_covid19_all.csv")
#write.csv(codextmaall_mph,"/home/bkzhu/SNE-multi/figure_rcode/covid/codex/tma1234_covidclus1_forMatch_macrophage.csv")


###### now we start to prep the COVID balf cite-seq dataset
## data retrieved from the COVID19 cell atlas under VIB/Gent hospital
## reading of the h5ad file was in python using package anndata but since only one line we are skipping it for simplicity.

citebelg_clean=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/covidbel-protein.csv") # read from h5ad, not shared due to size
# first we bin some columns in the citeseq dataset since same protein used multiple clones/ change protein names
citebelg_clean$CD138=citebelg_clean$CD138.A0055+citebelg_clean$CD138.A0831
citebelg_clean$CD138.A0055 <- NULL
citebelg_clean$CD138.A0831 <- NULL
citebelg_clean$ICOS=citebelg_clean$CD278
citebelg_clean$CD278 <- NULL
citebelg_clean$CD56=citebelg_clean$CD56.A0047+citebelg_clean$CD56.A0084
citebelg_clean$CD56.A0047 <- NULL
citebelg_clean$CD56.A0084 <- NULL
citebelg_clean$CD3=citebelg_clean$CD3.A0034+citebelg_clean$CD3.A0049
citebelg_clean$CD3.A0034 <- NULL
citebelg_clean$CD3.A0049 <- NULL
citebelg_clean$CD4=citebelg_clean$CD4.A0045+citebelg_clean$CD4.A0072+citebelg_clean$CD4.A0922
citebelg_clean$CD4.A0045 <- NULL
citebelg_clean$CD4.A0072 <- NULL
citebelg_clean$CD4.A0922 <- NULL
citebelg_clean$CD45=citebelg_clean$CD45.A0048+citebelg_clean$CD45.A0391
citebelg_clean$CD45.A0048 <- NULL
citebelg_clean$CD45.A0391 <- NULL
citebelg_clean$CD66a=citebelg_clean$CD66a.c.e
citebelg_clean$CD66a.c.e <- NULL
citebelg_clean$CD38=citebelg_clean$CD38.A0389+citebelg_clean$CD38.A0410
citebelg_clean$CD38.A0389 <- NULL
citebelg_clean$CD38.A0410 <- NULL
citebelg_clean$PD1=citebelg_clean$CD279
citebelg_clean$CD279 <- NULL
citebelg_clean$PD.L1=citebelg_clean$CD274
citebelg_clean$CD274 <- NULL
citebelg_clean$MMR=citebelg_clean$CD206
citebelg_clean$CD206 <- NULL

# second during analysis we found some unclear annotation between neutrophils and macrophages in the citeseq dataset
# thus we decided to do our own clustering and annotations on this dataset
# meanwhile there are some batch effects in the data, so we clustered per patient

## only cluster based on the protein features
citeblg_protein=citeblg_protein[,c(2:293)]
citeblg_protein_value=citeblg_protein[,c(2:279, 282:292)]
rownames(citeblg_protein_value)=citeblg_protein$X



##### cluster for each patient, start with patient cov002
temp=citeblg_protein
COV002=dplyr::filter(temp, grepl('COV002', X)) # 18660 cells
COV002_value=COV002
COV002_value$celltype<-NULL
COV002_value$cluster.term<-NULL
tm=as.data.frame(t(COV002))
data6 <- CreateSeuratObject(counts =tm)
data6 <- SetAssayData(object = data6, slot = "scale.data", new.data = as.matrix(tm))
data6 <- SetAssayData(object = data6, slot = "data", new.data = as.matrix(tm))
data6 <- RunPCA(data6, features = rownames(data6), reduction.name = "pca_cytof", reduction.key = "pca_cytof_", 
                verbose = FALSE)
DimPlot(data6, reduction = "pca_cytof")
ElbowPlot(data6, ndims = 50, reduction = "pca_cytof")
data6 <- RunTSNE(data6, dims = 1:15, reduction = "pca_cytof", reduction.key = "cyTSNE_", reduction.name = "tsne_cy",check_duplicates = FALSE)
data6 <- FindNeighbors(data6, features = rownames(data6), dims = NULL)
cy.data <- GetAssayData(data6, slot = "data")
cy.dist <- dist(t(cy.data))
data6[["cy_snn"]] <- FindNeighbors(cy.dist)$snn
data6 <- FindClusters(data6, resolution = 1, graph.name = "cy_snn")
new.cluster.ids <- c("mph","cd8","mph","cd4","cd4","plasma","plasma","NK","gd","gd","pdc","plasma","other","cdc")
names(new.cluster.ids) <- levels(data6)
data6 <- RenameIdents(data6, new.cluster.ids)
##### repeat for cov013/015/024/034/036/037 patients, and save them out


# after seurat clustering, get the protein expression matrix with annotation column, bin these patients back together
covp1=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV002_anno-sr.csv")
covp2=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV013_anno-sr.csv")
covp3=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV034_anno-sr.csv")
covp4=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV024_anno-sr.csv")
covp5=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV015_anno-sr.csv")
covp6=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV036_anno-sr.csv")
covp7=read.csv("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/COV037_anno-sr.csv")
allcovp=rbind(covp1,covp2)
allcovp=rbind(allcovp,covp3)
allcovp=rbind(allcovp,covp4)
allcovp=rbind(allcovp,covp5)
allcovp=rbind(allcovp,covp6)
allcovp=rbind(allcovp,covp7)
allcovp_macrophage=subset(allcovp,allcovp$cluster.sr == "mph") # 16k macrophage

## to now, both codex macrophages and covid macrophages are all prepared


#####################################
#### MARIO performed in python ######
#####################################



##################### part 2 #######################


### this part is the post integration analysis related to fig5 and figS5

## first we cluster the codex macropahges based on matched RNA expression levels (C1q abc expression levels)
library(dplyr)
library(data.table)
codex=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/lung_mph_mario_matched.csv")
cite=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/balf_mph_mario_matched.csv")
# link rna info
# rna counts extracted from the vib ghent balf h5ad file
cite_rna=as.data.frame(fread("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/covidbel-rna.csv"))
cite_rna_protein=left_join(cite,cite_rna, by=c("X"="V1")) # join the dataset with rna levels, by barcode
cite_rna_only=cite_rna_protein[,-c(1:296)] # remove the citseq proteins

complement_data = cite_rna_only[,c("C1QA","C1QB","C1QC")]
## cluster base on c1q abc expression level
library(stats)
library(factoextra)
res <- hcut(complement_data, k = 2, stand = TRUE) # only need to cluster two out
comp_clust = res$cluster
codex$c1qgroup = comp_clust
## cluster finished
cell_type_label = comp_clust
df = complement_data
mat_markers = sapply(df, function(df) (df-mean(df))/sd(df)) # z norm
cell_types = unique(comp_clust)
cell_type_mat = matrix(0,nrow=length(cell_types),ncol=ncol(complement_data))
for (idx in c(1:length(cell_types))){
  print(idx)
  row_match = cell_type_label %in% cell_types[idx]
  print(length(which(row_match)))
  cell_type_mat[idx,] = colMeans(mat_markers[row_match,])
}
# Heat map of znormed RNA expression level
library(RColorBrewer)
library(gplots)
colors = c(seq(-1.282,1.282,length=22)) # cap low and high values
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))
# heatmap of c1q abc expression (c1q high and c1q low)
hm2_call = heatmap.2(cell_type_mat,col=my_palette,breaks=colors,density.info="none",trace="none",
                     Rowv=F,Colv=T,dendrogram="col",symm=F,labRow=cell_types,labCol=colnames(mat_markers),
                     margins=c(2*(dim(cell_type_mat)[2]/dim(cell_type_mat)[1]),3),scale="none",cexRow=1,
                     cexCol=1,rowsep=c(0:20),colsep=c(0:33),sepcolor="black",sepwidth=c(0.0001,0.0001))

## next heat map of codex protein for c1q high and low macrophages
cell_type_label = comp_clust
use = c("MMR","CD16","MMP.9","HLA.DR","Ki.67","CD163","CD11b","PD.L1","CD68","CD45RA","CD36") # only use mph relevant markers
df = codex[,use]
mat_markers = codex[,use]
mat_markers = scale(mat_markers)
cell_types = unique(comp_clust)
cell_type_mat = matrix(0,nrow=length(cell_types),ncol=ncol(mat_markers))
for (idx in c(1:length(cell_types))){
  print(idx)
  row_match = cell_type_label %in% cell_types[idx]
  print(length(which(row_match)))
  cell_type_mat[idx,] = colMeans(mat_markers[row_match,])
}
# Heat map of znormed codex proteins
library(RColorBrewer)
library(gplots)
colors = c(seq(-0.5,0.5,length=22)) # cap low and high values
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))
hm2_call = heatmap.2(cell_type_mat,col=my_palette,breaks=colors,density.info="none",
                     trace="none",Rowv=F,Colv=T,dendrogram="col",symm=F,labRow=cell_types,
                     labCol=colnames(mat_markers),margins=c(2*(dim(cell_type_mat)[2]/dim(cell_type_mat)[1]),3),
                     scale="none",cexRow=1,cexCol=1,rowsep=c(0:20),colsep=c(0:33),sepcolor="black",sepwidth=c(0.0001,0.0001))

## heatmap of citeseeq protein for c1q high and low macrophages
use = c("MMR","CD16","HLA.DR","CD163","CD11b","CD68","CD45RA","CD36","CD169","CD209","Siglec8","CD25") # mph relevant antibodies
df = cite[,use]
mat_markers = as.matrix(df)
mat_markers = scale(mat_markers)
cell_types = unique(comp_clust)
cell_type_mat = matrix(0,nrow=length(cell_types),ncol=ncol(mat_markers))
for (idx in c(1:length(cell_types))){
  print(idx)
  row_match = cell_type_label %in% cell_types[idx]
  print(length(which(row_match)))
  cell_type_mat[idx,] = colMeans(mat_markers[row_match,])
}
# Heat map of znomred citeseq protein values
library(RColorBrewer)
library(gplots)
colors = c(seq(-0.5,0.5,length=22)) # cap high and low values
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))
hm2_call = heatmap.2(cell_type_mat,col=my_palette,breaks=colors,density.info="none",
                     trace="none",Rowv=F,Colv=T,dendrogram="col",symm=F,labRow=cell_types,
                     labCol=colnames(mat_markers),margins=c(2*(dim(cell_type_mat)[2]/dim(cell_type_mat)[1]),3),
                     scale="none",cexRow=1,cexCol=1,rowsep=c(0:20),colsep=c(0:33),sepcolor="black",sepwidth=c(0.0001,0.0001))

## now analysis per patient

## first rank c1q high percentage of macrophage per patient

patient_region = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/patient_region.csv")
codex_patient = left_join(codex, patient_region, by = c("tma", "region")) # add matched file with patient info
codex_patient = as_tibble(codex_patient)
# calculate the percentage of c1q high mph among high and low mph
perc_cores = codex_patient %>%
  group_by(tma,region) %>%
  dplyr::count(c1qgroup) %>%
  dplyr::mutate(countT= sum(n)) %>%
  dplyr::mutate(perc= n / countT) 

perc_cores_patient = left_join(perc_cores, patient_region, by = c("tma", "region")) # need this for plotting
perc_cores_patient = subset(perc_cores_patient, perc_cores_patient$c1qgroup == 2) # only show high percentage on plot
p <- ggplot(perc_cores_patient, aes(x=reorder(patient, - perc), y=perc)) + 
  geom_boxplot(lwd = 0.8) + theme_classic()
p + ylab("perc macro high C1Q complex") + ggtitle("All Patient rank by macrophage complement percentage")

## use the same c1q mph percentage per patient rank and plot out neutrophil percentage
# original codex file will all the cells:
codex_allCelltypes = read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/CODEX_covid19_all.csv")

neutro_perc = codex_allCelltypes %>% dplyr::filter(!is.na(cluster.term)) %>% # dont need NA cells
  group_by(tma,region) %>% 
  dplyr::summarise(neutro_count= sum(cluster.term == "Neutrophil") , n = n()) %>%
  mutate( neu_perc = neutro_count / n) %>%
  left_join(patient_region , by = c("tma", "region"))
# use the same rank previously
neutro_perc$patient <- factor(neutro_perc$patient,                          
                                levels = c("b-18","b-9","b-24", "b-21", "b-10",
                                           "b-12","b-1","b-7","b-6","b-3","b-22","b-5",
                                           "b-17","b-20",
                                           "b-11","b-25","b-2","b-4","b-28","b-14","b-23","b-13","b-16"))
p <- ggplot(neutro_perc, aes(x=patient, y=neu_perc)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.8) + theme_classic()
p + ylab("perc neutrophils of all cells in core") + ggtitle("All Patient neutrophil infiltration use complement rank")

## directly compare per core neutrophil percentage vs c1q high macrophage percentage
comp_neutro = left_join (neutro_perc, perc_cores_patient, by = c("tma", "region"))
p = ggplot(comp_neutro, aes(x = log(neu_perc), y = log(perc))) +
  stat_summary(fun.data=mean_cl_normal) + theme_classic() +
  geom_smooth(method='lm', formula= y~x) + xlab("neutrophil percent") + ylab("complement macro percent") +
  geom_point(aes(color = patient.x))


## production of pseudo plot with c1q high and low macrophage positions:
# need codex all cells with c1q clustering info
common= colnames(codextmaall)[2:70]# the columns to use to match
## NOTE!!!: please do not use csv saved and read file here
## saved and read files will lose precise number info to match
## thus will lose some cells, bad code
## this is not a safe code to run
codex_allCelltypes_update=left_join(codextmaall, codex, by = common) # codex all cells with c1q clustering info 0,1,NA
coluse = c("x","y","tma.x","region","c1qgroup","cluster.term.x")
df_simp = codex_allCelltypes_update[,coluse]
df_simp$c1qgroup = df_simp$c1qgroup +1
df_simp[is.na(df_simp$c1qgroup),"c1qgroup"] <- 0 # three labels 0 not predicted, 1 low, 2 high
### example: plot pseudo tma3 region 18
temp = subset(df_simp, df_simp$tma ==3 & df_simp$region == 8 & df_simp$cluster.term.x =="Macrophage" & df_simp$c1qgroup!=0)
temp$c1qgroup = as.factor(temp$c1qgroup)
temp2  = subset(df_simp, df_simp$tma ==3 & df_simp$region == 8 &df_simp$c1qgroup!=1 &df_simp$c1qgroup!=2)
p = ggplot() +
  geom_point(data = temp, aes(x, y ,colour=as.factor(c1qgroup)), size =5, alpha = 0.5, stroke = 0.001) +
  geom_point(data = temp2, aes(x, y, colour=as.factor(c1qgroup)), size =2, alpha = 0.15, stroke = 0.001)+
  theme_classic() + scale_colour_manual(values = c("grey","#76B947","#741AAC"))
### this code repeat for other cores

## now produce the radius anchor plot

############################## anchor plot function, from sizun et al immunity paper ####################################

get_distance_matrix <- function(df_input, id_col, x_col, y_col, celltype_col, anchor_col, region_col, keep_anchor = FALSE, split_celltype = TRUE, columns_to_use,
                                distance_threshold = 320, bins = 10, num_cores = NULL, um_pixel = 100/320){
  # df_input = dataframe with marker expressions, celltypes, etc
  # id_col = column name that corresponds to each unique cell ID
  # x_col = column name that corresponds to x axis location of cells
  # y_col = column name that corresponds to y axis location of cells
  # celltype_col = includes option to group by celltype, if null will not group by the cell types, otherwise is column name that corresponds to celltype annotations
  # anchor_col = column name that has T/F boolean annotations for cells to use as anchor columns
  # distance_treshold = distance to consider
  # bins = number of bins to split distance and assign cells into
  # num_cores = number of cores to use during parallelization, default to using everything but 1
  start_time <- Sys.time()
  cat(paste0("Beginning distance matrix calculations: ", start_time, "\n"))
  
  # Area distance calculations
  r_max = distance_threshold * um_pixel
  
  bin_area = 3.14/bins*(r_max)^2
  bin_conversion = data.frame(bin_num = seq(1, bins, 1)) %>%
    # Get bin radii
    mutate(bin_radius_min = (bin_num-1)*r_max/bins,
           bin_radius_max = (bin_num)*r_max/bins) %>%
    # Calculate area of each bin
    mutate(bin_area = 3.14*(bin_radius_max^2 - bin_radius_min^2))
  
  # Get unique regions in the dataset
  unique_regions = unique(df_input[[region_col]])
  
  # Set up parallelization
  if(is.null(num_cores)){
    num_cores = parallel::detectCores() - 1
  }
  library(doSNOW)
  library(foreach)
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(unique_regions), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Calculate distance matrix
  df_distance_matrix <- foreach(fov = 1:length(unique_regions), .packages = c("tidyverse"), .combine = "rbind", .options.snow = opts)%dopar%{
    # Subset larger dataframe by region
    df_point_subset <- df_input %>%
      dplyr::filter(!!as.symbol(region_col) == unique_regions[fov])
    
    subset_anchor_cells = which(df_point_subset[[anchor_col]] == TRUE)
    
    # Check to make sure anchor cells are present
    if (length(subset_anchor_cells) > 0){ # anchor cells are present
      # Calculate df_anchor
      df_anchor = do.call(rbind, lapply(1:length(subset_anchor_cells), function(i){
        new_x = df_point_subset[[x_col]] - as.double(df_point_subset[subset_anchor_cells[i],][x_col])
        new_y = df_point_subset[[y_col]] - as.double(df_point_subset[subset_anchor_cells[i],][y_col])
        # Calculate radius around cells
        dist_from_anchor = spatstat.geom::crossdist(new_x, new_y, 0, 0, squared = F)
        
        # Threshold anchor cells to keep
        if(keep_anchor == FALSE){
          # Don't include anchor cells in downstream analysis
          cells_to_keep = which(dist_from_anchor <= distance_threshold & dist_from_anchor != 0)
        } else {
          # Include anchor cell
          cells_to_keep = which(dist_from_anchor <= distance_threshold)
        }
        kept_cells = sort(dist_from_anchor[cells_to_keep],
                          decreasing = FALSE,
                          index.return = TRUE)
        cell_id = df_point_subset[subset_anchor_cells[i],][[id_col]] # Anchor cell ID
        
        if(length(kept_cells$x) > 0){
          # make sure there are kept cells
          df_anchor_marker = data.frame(anchor_cell_id = cell_id,
                                        anchor_cell_type = df_point_subset[which(df_point_subset[[id_col]] == cell_id),][celltype_col],
                                        kept_cell_id = df_point_subset[cells_to_keep[kept_cells$ix], 1],
                                        kept_cell_distance = dist_from_anchor[cells_to_keep[kept_cells$ix]],
                                        row.names = NULL)
          colnames(df_anchor_marker) <- c("anchor_cell_id", "anchor_cell_type", "kept_cell_id", "kept_cell_distance")
          df_anchor_marker <- df_anchor_marker %>%
            # Join kept cell marker expressions
            left_join(df_point_subset, by = c("kept_cell_id" = id_col)) %>%
            # Convert from pixel to um
            mutate(kept_cell_distance = kept_cell_distance * um_pixel)
          
          # Append bin information to distances
          df_anchor_marker$bin = sapply(1:nrow(df_anchor_marker), function(x){
            dist = df_anchor_marker$kept_cell_distance[x]
            bin_num = sapply(dist, function(dist){
              ind = with(bin_conversion, which(dist >= bin_radius_min & dist < bin_radius_max))
              bin = bin_conversion$bin_num[ind]
              return(bin)
            })
          })
          df_anchor_marker$bin <- as.numeric(df_anchor_marker$bin)
          # assign cells at r_max to highest bin
          df_anchor_marker$bin[is.na(df_anchor_marker$bin)] <- max(df_anchor_marker$bin, na.rm = TRUE)
          
          
          # calculate marker summaries
          if(split_celltype == TRUE){
            marker_summary <- df_anchor_marker %>%
              gather(key = "marker", value = "expr", colnames(df_input[,columns_to_use])) %>%
              # Get total number of cells present in bin and the total expression in bin
              group_by(anchor_cell_id, anchor_cell_type, bin, !!as.symbol(celltype_col), !!as.symbol(region_col), marker) %>%
              summarize(cell_counts = n(),
                        sum_expr = sum(expr, na.rm = TRUE)) %>%
              ungroup() %>%
              # Append bin areas
              left_join(bin_conversion %>% 
                          select(bin_num, bin_area),  by = c("bin" = "bin_num")) %>%
              mutate(normalized = (sum_expr/cell_counts)/bin_area) %>%
              select(-cell_counts, -sum_expr, -bin_area) %>%
              spread(key = marker, value = normalized, fill = 0)
          }
          if(split_celltype == FALSE){
            marker_summary = df_anchor_marker %>%
              gather(key = "marker", value = "expr", colnames(df_input[,columns_to_use])) %>%
              group_by(anchor_cell_id, anchor_cell_type, bin, !!as.symbol(region_col), marker) %>%
              summarize(cell_counts = n(), 
                        sum_expr = sum(expr, na.rm = TRUE),
                        mean_expr = mean(expr, na.rm = TRUE)) %>%
              ungroup() %>%
              # reorganize output table
              gather(key = "metric", value = "value", cell_counts:mean_expr) %>%
              spread(key = "marker", value = "value") 
          }
          
        } else {
          marker_summary = data.frame()
        }
        return(marker_summary)
      }))
    } else { # no anchor cells present
      df_anchor = data.frame()
    }
    return(df_anchor)
  }
  
  close(pb)
  stopCluster(cl)
  # End information
  end_time <- Sys.time()
  cat(paste0("Completed distance matrix calculations: ", end_time, "\n"))
  cat(paste0("Elapsed time: ", end_time - start_time, "\n"))
  
  return(df_distance_matrix)
}


################################################################################### function finished

## start the anchor production
library(dplyr)
library(tidyverse)
outputdir = "temp" # not used
id_col = "rowid"
x_col = "x" # x pos
y_col = "y" # y pos
celltype_col = "cluster.term.x" # annotation
region_col = "region_new" # unique region specifing which core
split_celltype = FALSE
keep_anchor = FALSE
distance_threshold = 264 # codex 20x 264 pixels == 100um
um_pixel = 100/distance_threshold
bins = 6 # radiu 100um seperate by 6 parts
num_cores = 10

## some addtional radius functions

r_max = distance_threshold * um_pixel
bin_area = 3.14/bins*(r_max)^2
bin_conversion = data.frame(bin_num = seq(1, bins, 1)) %>%
  # Get bin radii
  mutate(bin_radius_min = (bin_num-1)*r_max/bins,
         bin_radius_max = (bin_num)*r_max/bins) %>%
  # Calculate area of each bin
  mutate(bin_area = 3.14*(bin_radius_max^2 - bin_radius_min^2))
##

## This analysis need one hot style cell type clustering result (annotation) in the matrix
library(data.table)
library(mltools)
codex_allCelltypes_update2 = subset(codex_allCelltypes_update, !is.na(codex_allCelltypes_update$cluster.term.x)) # exclude NA cell types
clusterm.x = codex_allCelltypes_update2$cluster.term.x
codex_allCelltypes_update_onehot = codex_allCelltypes_update2 %>% mutate(value = 1)  %>% spread(cluster.term.x , value,  fill = 0 ) 
codex_allCelltypes_update_onehot$cluster.term.x = clusterm.x # resume the column disregared during one hot
codex_allCelltypes_update_onehot = as.tibble(codex_allCelltypes_update_onehot)
use = c("region_new","c1qgroup","x","y","tma.x","cluster.term.x","B cell",
        "CD38pos CD45pos","CD4+ T-cell","CD8+ T-cell","dirt",
        "gd T-cell","lymphatic vessels","Macrophage","myeloid cells",
        "Neutrophil","NK","Plasma cell ","pneumocytes","smooth muscle",
        "stromal cells","unknown","vasculature cells","VISTA myeloid cells")
codex_allCelltypes_update_onehot = codex_allCelltypes_update_onehot[,use]
# start the process
df_master_clean <- codex_allCelltypes_update_onehot %>% arrange(region_new) # codex all cells with c1q clustering info
df_master_clean$rowid = as.numeric(rownames(df_master_clean))
df_master_clean$comp_clust  = NA
df_master_clean$comp_clust[df_master_clean$c1qgroup == 2] = TRUE
df_master_clean$comp_clust[df_master_clean$c1qgroup == 1] = FALSE # turn c1q annotation into NA/TRUE/FALSE
#########
write.csv(df_master_clean, "/home/bkzhu/SNE-multi/figure_rcode/covid/test.csv") # save for repeat use


###
anchor_col = "comp_clust"
df_master_clean <- read_csv("/home/bkzhu/SNE-multi/figure_rcode/covid/test.csv")
df_input <- df_master_clean
df_input = subset(df_input,!(df_input$comp_clust == FALSE & df_input$cluster.term.x == "Macrophage")) # get rid of low and NA macrophage
columns_to_use = colnames(df_input)[8:25] # using cell type columns one hot columns
df_input2 = as.tibble(df_input)
## actual function
H = get_distance_matrix(df_input = df_input2,
                        id_col = id_col, x_col = x_col, y_col = y_col, celltype_col = celltype_col, 
                        columns_to_use = columns_to_use,
                        anchor_col = anchor_col, region_col = region_col, split_celltype = split_celltype,
                        keep_anchor = keep_anchor, distance_threshold = distance_threshold, bins = bins, 
                        num_cores = num_cores,
                        um_pixel = um_pixel)
# scale values by cell numbers
all_target_columns = colnames(H)[6:23] # only take the cell type info
cell_counts = rowSums(H[,all_target_columns])
H[,all_target_columns] = H[,all_target_columns] / cell_counts

######
## now calculate the anchor plot for c1q low macrophages
# reverse the anchors this time
df_master_clean <- read_csv("/home/bkzhu/SNE-multi/figure_rcode/covid/test.csv")
df_master_clean_reverse <- df_master_clean
df_master_clean_reverse[,anchor_col]=!df_master_clean_reverse[,anchor_col]
df_master_clean_reverse[,anchor_col][is.na(df_master_clean_reverse[,anchor_col])] <- FALSE
df_input <- df_master_clean_reverse
df_input = subset(df_input,!(df_input$comp_clust == FALSE & df_input$cluster.term.x == "Macrophage")) # get rid of high and NA macrophage

columns_to_use = colnames(df_input)[8:25] # using cell type columns one hot columns
df_input = as.tibble(df_input)

L = get_distance_matrix(df_input = df_input,
                        id_col = id_col, x_col = x_col, y_col = y_col, celltype_col = celltype_col, 
                        columns_to_use = columns_to_use,
                        anchor_col = anchor_col, region_col = region_col, split_celltype = split_celltype,
                        keep_anchor = keep_anchor, distance_threshold = distance_threshold, bins = bins, 
                        num_cores = num_cores,
                        um_pixel = um_pixel)
cell_counts = rowSums(L[,columns_to_use])
L[,all_target_columns] = L[,all_target_columns] / cell_counts

###### not anchor targets input
# Get celltype column names

## note do not load plyr here since code not specified function via dplyr::
celltype_columns = 8:25 # same only look at cell type information
celltype_columns = colnames(df_input[,celltype_columns])
celltypes_summaryH <- H %>%
  # select only the celltype columns
  select(anchor_cell_id, anchor_cell_type, bin, region_new, metric, all_of(celltype_columns)) %>%
  # Get the distribution counts of each celltype in a bin
  dplyr::filter(metric == "sum_expr") %>%
  gather(key = "celltype", value = "cell_count", celltype_columns) %>%
  # Append bin data
  left_join(bin_conversion %>%
              select(bin_num, bin_area), by = c("bin" = "bin_num")) %>%
  mutate(cell_density = cell_count/1) %>% # 1 was binarea for cell density version, changed 0617
  select(-cell_count, -bin_area) %>%
  # Get statistics for each celltype and bin
  group_by(anchor_cell_type, bin , celltype) %>%
  dplyr::summarize(avg = mean(cell_density, na.rm = TRUE),
            std_dev = sd(cell_density, na.rm = TRUE),
            count = n()) %>%
  ungroup() %>%
  # Calculate standard error and confidence intervals
  mutate(se = std_dev/sqrt(count),
         lower_ci = avg - qt(1 - ((1 - 0.95)/2), count - 1) * se,
         upper_ci = avg + qt(1 - ((1 - 0.95)/2), count - 1) * se) %>%
  ungroup()

celltypes_summaryL <- L %>%
  # select only the celltype columns
  select(anchor_cell_id, anchor_cell_type, bin, region_new, metric, all_of(celltype_columns)) %>%
  # Get the distribution counts of each celltype in a bin
  dplyr::filter(metric == "sum_expr") %>%
  gather(key = "celltype", value = "cell_count", celltype_columns) %>%
  # Append bin data
  left_join(bin_conversion %>%
              select(bin_num, bin_area), by = c("bin" = "bin_num")) %>%
  mutate(cell_density = cell_count/1) %>% # 1 was binarea for cell density version, changed 0617
  select(-cell_count, -bin_area) %>%
  # Get statistics for each celltype and bin
  group_by(anchor_cell_type, bin , celltype) %>%
  dplyr::summarize(avg = mean(cell_density, na.rm = TRUE),
            std_dev = sd(cell_density, na.rm = TRUE),
            count = n()) %>%
  ungroup() %>%
  # Calculate standard error and confidence intervals
  mutate(se = std_dev/sqrt(count),
         lower_ci = avg - qt(1 - ((1 - 0.95)/2), count - 1) * se,
         upper_ci = avg + qt(1 - ((1 - 0.95)/2), count - 1) * se) %>%
  ungroup()

#######
celltypes_summaryH$group=anchor_col
celltypes_summaryL$group=paste0("No_",anchor_col)
HL_sum=rbind(celltypes_summaryH,celltypes_summaryL)
###### ploting
HL_sum_sub = subset(HL_sum,HL_sum$celltype == "Macrophage")

## ggplot production
p = ggplot(data = HL_sum_sub,
       mapping = aes(x = factor(bin), color = group, group = group)) +
  geom_line(mapping = aes(y = avg)) +
  geom_ribbon(mapping = aes( y = avg, ymin = lower_ci, ymax = upper_ci, fill = group),
              alpha = 0.2, linetype = "dashed") +
  theme_bw() + ggtitle("complement subcluster macrophages to macrophages")



############### cell cell interaction analysis ############

### used same function presented in sizun jiang et al immunity:
source("/home/bkzhu/SNE-multi/figure_rcode/covid/cell_cell_interaction/get_triangulation_distances.R")
source("/home/bkzhu/SNE-multi/figure_rcode/covid/cell_cell_interaction/iterate_triangulation_distance.R")
filedir = "/home/bkzhu/SNE-multi/figure_rcode/covid/cell_cell_interaction/"
library(dplyr)
library(tidyverse)
# read in data frame containing all the cells in tma123 covid patients, with c1q info
df_master_clean <- read_csv("/home/bkzhu/SNE-multi/figure_rcode/covid/test.csv")
# some pre processing, delete dirt and unknown cell types
df_master_clean = subset(df_master_clean,df_master_clean$cluster.term.x !="unknown")
df_master_clean = subset(df_master_clean,df_master_clean$cluster.term.x !="dirt")
## disregard vague annotations in this analysis
df_master_clean = subset(df_master_clean,df_master_clean$cluster.term.x !="CD38pos CD45pos")
df_master_clean$cluster.term.x = as.character(df_master_clean$cluster.term.x)
df_master_clean$cluster.term.x[df_master_clean$cluster.term.x == "VISTA myeloid cells"] = "myeloid cells" # change name
coluse = c("comp_clust", "region_new","x","y","cluster.term.x")
df_master_clean = df_master_clean[,coluse]
# change the c1q labels
df_master_clean$comp_clust <- as.integer(df_master_clean$comp_clust)
df_master_clean$comp_clust = df_master_clean$comp_clust + 1
df_master_clean[is.na(df_master_clean)] <- 0 # three labels 0 not predicted, 1 low, 2 high
############### first low
df_master_clean_low = df_master_clean
df_master_clean_low = df_master_clean_low %>% mutate(region_new_id = paste0(region_new, "_low"))
df_master_clean_low = subset(df_master_clean_low, df_master_clean_low$comp_clust !=2) # delete mphage c1q high
df_master_clean_low = subset(df_master_clean_low,
                             !(df_master_clean_low$cluster.term.x == "Macrophage" & df_master_clean_low$comp_clust==0)) # remove not predicted mph
df_master_clean_low$meta = "low"
############### high
df_master_clean_high = df_master_clean
df_master_clean_high = df_master_clean_high %>% mutate(region_new_id = paste0(region_new, "_high"))
df_master_clean_high = subset(df_master_clean_high, df_master_clean_high$comp_clust !=1) # delete mphage c1q low
df_master_clean_high = subset(df_master_clean_high,
                              !(df_master_clean_high$cluster.term.x == "Macrophage" & df_master_clean_high$comp_clust==0)) # remove not predicted mph
df_master_clean_high$meta = "high"
############### high low combine
df_master_clean_HL =  rbind(df_master_clean_high,df_master_clean_low)
df_master_clean_HL = df_master_clean_HL %>% rowid_to_column("data_index")
patient_metadata = df_master_clean_HL[,c("region_new_id","meta")]
##### observed trangulation distance
library(parallel)
fname = "triangulation_observed_c1qHL_temp.csv"
if(file.exists(paste0(filedir, fname))){
  print(paste0("reading ", fname))
  tri_dist_celltype <- read_csv(paste0(filedir, fname))
} else {
  tri_dist_celltype <- get_triangulation_distances(df_input = df_master_clean_HL,
                                                   id = "data_index",
                                                   x_pos = "x",
                                                   y_pos = "y",
                                                   cell_type = "cluster.term.x",
                                                   region = "region_new_id")
  write_csv(tri_dist_celltype,
            paste0(filedir, "triangulation_observed_c1qHL_temp.csv"))
}
tri_dist_celltype$group = gsub(".*_", "", tri_dist_celltype$region_new_id)
##### calculate iteration trangulation distance
fname = "triangulation_iteration_c1qHL_temp.csv"
if(file.exists(paste0(filedir, fname))){
  print(paste0("reading ", fname))
  tri_dist_iteration_celltype <- read_csv(paste0(filedir, fname))
} else {
  tri_dist_iteration_celltype <- iterate_triangulation_distance(df_input = df_master_clean_HL,
                                                                id = "data_index",
                                                                x_pos = "x",
                                                                y_pos = "y",
                                                                cell_type = "cluster.term.x",
                                                                region = "region_new_id")
  write_csv(tri_dist_celltype,
            paste0(filedir, "triangulation_iteration_c1qHL_temp.csv"))
}

tri_dist_iteration_celltype$group = gsub(".*_", "", tri_dist_iteration_celltype$region_new_id)

#### start the plotting process
#Set distance threshold for observed cell-cell interactions
distance_threshold = 264 # corresponds to 100um
# Reformat observed dataset
ggdf_observed_celltype <- tri_dist_celltype %>%
  dplyr::filter(distance <= distance_threshold) %>%
  # Calculate observed mean distance and list values
  group_by(celltype1, celltype2, group) %>%
  summarize(observed = list(distance),
            observed_mean = mean(unlist(observed), na.rm = TRUE)) %>%
  ungroup()
##### Reformat exepcted dataset
ggdf_expected_celltype <- tri_dist_iteration_celltype %>%
  gather(key = "iteration", value = "mean_distance", -c(celltype1, celltype2, region_new_id))
ggdf_expected_celltype <- ggdf_expected_celltype %>% 
  mutate( group = gsub(".*_", "", region_new_id) ) %>%
  # Calculate expected mean distance and list values
  group_by(celltype1, celltype2, group)
ggdf_expected_celltype$mean_distance = as.numeric(ggdf_expected_celltype$mean_distance) # NA exist but is taken care of 
ggdf_expected_celltype <- ggdf_expected_celltype %>%  
  summarize(expected = list(mean_distance),
            expected_mean = mean(mean_distance, na.rm = TRUE)) %>%
  ungroup() 
# Calculate pvalues
ggdf_pvals_celltypes <- ggdf_expected_celltype %>%
  left_join(ggdf_observed_celltype,
            by = c("celltype1", "celltype2", "group")) %>%
  # Calculate wilcoxon test between observed and expected distances
  group_by(celltype1, celltype2, group) %>%
  mutate(pvalue = wilcox.test(unlist(expected), unlist(observed), exact = FALSE)$p.value) %>%
  ungroup() %>%
  select(-observed, -expected)
# Calculate log fold enrichment
ggdf <- ggdf_pvals_celltypes %>%
  # Calculate log fold enrichment
  mutate(logfold_group = log2(observed_mean/expected_mean),
         interaction = paste0(celltype1, " --> ", celltype2)) %>%
  # Filter out pvalues less than 0.05
  dplyr::filter(pvalue < 0.05)
ggdf2 = subset(ggdf, ggdf$celltype1 == "Macrophage" | ggdf$celltype2 == "Macrophage")
# Create plots
outdir <- paste0(outputdir, "Dumbbell Plots/")
dir.create(outdir, showWarnings = FALSE)
# Assign cell type colors
ggdf3 <- ggdf2 %>%
  select(interaction, group, logfold_group) %>%
  spread(key = group, value = logfold_group) %>%
  mutate(difference = (high - low)) %>%
  dplyr::filter(!is.na(difference)) %>%
  arrange(low)
ord_diff = ggdf3$interaction
# Assign interaction order
ggdf2$interaction <- factor(ggdf2$interaction,
                            levels = ord_diff)
# Assign group factor
ggdf2$group <- factor(ggdf2$group,
                      levels = c("high", "low"))
### make dumb bell plots
p = ggplot(data = ggdf2 %>%
         dplyr::filter(!is.na(interaction))) +
  geom_vline(mapping = aes(xintercept = 0), linetype = "dashed") +
  geom_line(mapping = aes(x = logfold_group, y = interaction),
            na.rm = TRUE) +
  geom_point(mapping = aes(x = logfold_group, y = interaction, fill = group), 
             size = 4, shape = 22, stroke = 0.5, na.rm = TRUE) +
  scale_fill_manual(values=c("#FF4500", "#145DA0")) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


########## gene expression related analysis ############
# read again , but not previously already loaded can skip
codex=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/lung_mph_mario_matched.csv") 
cite=read.csv("/home/bkzhu/SNE-multi/figure_rcode/figure_related_code_submit/data/COVID-19/balf_mph_mario_matched.csv")
# link rna info
# large rna file ignore reading again
#cite_rna=as.data.frame(fread("/home/bkzhu/SNE-multi/figure_rcode/covid/cite/covidbel-rna.csv"))
codex_regall=codex
cite_regall=cite
cite_regall_rna_protein=left_join(cite_regall,cite_rna, by=c("X"="V1"))
cite_regall_rna_only=cite_regall_rna_protein[,-c(1:296)]
cite_rna_protein=left_join(cite,cite_rna, by=c("X"="V1"))
# add row name required for seurat
rows = paste0(as.character(c(1:dim(cite_regall_rna_only)[1])),"_")
mph_obj=CreateSeuratObject(counts=t(cite_regall_rna_only),assay="RNA")
SetAssayData(object = mph_obj, slot = "data", new.data = t(cite_regall_rna_only), assay="RNA")
Idents(mph_obj) = as.character(comp_clust) # just the c1q high and low clustering result as group variable
###### this will just take ~2hr
de.rna_full_more = FindMarkers(mph_obj, ident.1 = "1", ident.2 = "2", min.cells.group = 1, 
                               min.cells.feature = 1,
                               min.pct = 0,
                               logfc.threshold = 0,
                               only.pos = FALSE)
##### then we based on the DE genes, hand pick some interesting genes
comp12_temp_simp_values =  (cite_regall_rna_only[,c("C1QA","C1QB","C1QC","CXCL8", "CCL2", "IL1B",
                                      "S100A12", "S100A8", "CCL3", "CCL4", "SPP1",
                                      "NAMPT", "CCL7","F13A1", "PLTP", "FOLR2",
                                      "RNASE1", "LGMN", "TMEM176B", "MS4A6A", "CCL18","TLR4","TREM2")])
# range norm the counts
comp12_temp_simp_values = apply(comp12_temp_simp_values, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
comp12_temp_simp_values = as.data.frame(comp12_temp_simp_values)
comp12_temp_simp_values$type = comp_clust # group based on high and low clustering result
comp12_temp_simp_melt_znorm = melt(comp12_temp_simp_values, id.vars = c("type"),
                                   measure.vars = c("C1QA","C1QB","C1QC","CXCL8", "CCL2", "IL1B",
                                                    "S100A12", "S100A8", "CCL3", "CCL4", "SPP1",
                                                    "NAMPT", "CCL7","F13A1", "PLTP", "FOLR2",
                                                    "RNASE1", "LGMN", "TMEM176B", "MS4A6A", "CCL18","TLR4","TREM2"))

## patched violin plots
library(patchwork)
#polar_genes = "CCL2"
polar_genes = c("C1QA","C1QB","C1QC","CXCL8", "CCL2", "IL1B",
                "S100A8", "CCL3", "SPP1",
                "NAMPT", "CCL7","F13A1", "PLTP", "FOLR2",
                "RNASE1", "LGMN", "TMEM176B", "MS4A6A", "CCL18","TLR4","TREM2")
plot_list = list()
for (item in polar_genes){
  temp = subset(comp12_temp_simp_melt_znorm, comp12_temp_simp_melt_znorm$variable == item)
  temp$type = as.character(temp$type)
  p = ggplot(temp, aes(x=type, y=value, fill = type)) +
    geom_violin(scale = 'width', lwd = 0.2) +
    #geom_boxplot(width=0.05, outlier.shape = NA, color = "white") +
    coord_flip()+
    theme_void()+
    theme(legend.position = "none")+
    scale_fill_manual(values=c("#76B947","#741AAC"))
  plot_list[[item]] = p
}
wrap_plots(plot_list)+ plot_layout(nrow =3)


####### go term analysis #######

#  high expression in c1q high genes
#  high expression in c1q low genes
de.rna_full_more_minus = subset(de.rna_full_more, -log10(de.rna_full_more$p_val) > 2 & de.rna_full_more$avg_logFC < -0.5)
de.rna_full_more_plus = subset(de.rna_full_more, -log10(de.rna_full_more$p_val) > 2 & de.rna_full_more$avg_logFC > 0.5)
# save out these two list and use panther go term analysis
# http://geneontology.org/

# read in the panther results
c1qH_go = read.table("/home/bkzhu/SNE-multi/figure_rcode/covid/go_files/c1qhigh.txt", sep = "\t", header = TRUE)
c1qL_go = read.table("/home/bkzhu/SNE-multi/figure_rcode/covid/go_files/c1qlow.txt", sep = "\t", header = TRUE)
# filter based on fold change >10
c1qH_go = subset(c1qH_go,c1qH_go$upload_1..fold.Enrichment.>10)
c1qL_go = subset(c1qL_go,c1qL_go$upload_1..fold.Enrichment.>10)
# remove go term tail
c1qH_go$GO.biological.process.complete = gsub(' \\(.*\\)', '',c1qH_go$GO.biological.process.complete)
c1qL_go$GO.biological.process.complete = gsub(' \\(.*\\)', '',c1qL_go$GO.biological.process.complete)
# arrange by fold enrichment
c1qH_go$GO.biological.process.complete <- factor(c1qH_go$GO.biological.process.complete, 
                                                 levels = c1qH_go$GO.biological.process.complete[order(c1qH_go$upload_1..fold.Enrichment.)])

c1qL_go$GO.biological.process.complete <- factor(c1qL_go$GO.biological.process.complete, 
                                                 levels = c1qL_go$GO.biological.process.complete[order(c1qL_go$upload_1..fold.Enrichment.)])


## go term plots
p = ggplot(data = c1qH_go[c(1:30,], aes(x = upload_1..fold.Enrichment., y = GO.biological.process.complete, 
                               color = log10(upload_1..FDR.), size = upload_1..85.)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO for C1qH")

p2 = ggplot(data = c1qL_go[c(1:30,], aes(x = upload_1..fold.Enrichment., y = GO.biological.process.complete, 
                                color = log10(upload_1..FDR.), size = upload_1..106.)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO for C1qL")


####### ISG analysis ######


### list of ~ 50 isg related to covid response is retrieved from cell journal
# 'The interferon landscape along the respiratory tract impacts the severity of COVID-19'
polar_genes = c("MAX","ST3GAL4","GSDMD","NAPA", "STAT1", "SPATS2L",
                "B4GALT5", "ELF1", "RAB27A", "IFIT1", "IFIT5",
                "IFIT3", "DDX60","STAT2", "BST2",
                "UBD", "SUSD3", "REC8", "ETV6","CNP","USP18",
                "ARNTL","CLEC4D","GBP3","NRN1","TAGAP","PIK3AP1",
                "MYD88","ISG15","MLKL","RSAD2","CCND3","SPATA13",
                "ISG20","ZBP1","IFITM3","RGS22","FGD2","FZD5","CASP7",
                "APOL4","ERLIN1","MSR1","HSPA8","LY6E",
                "GNB4","TRIM21","RAB39A",
                "ANGPTL4","FNDC4","DNAJC6","APOL2","IL11","IL4I1")
temp12 = cite_regall_rna_only[,polar_genes]
temp12 = scale(temp12) # znorm before plotting
cell_type_label = comp_clust # add grouping information
cell_type_matdf = as.data.frame(cell_type_mat)
heatmap(as.matrix(cell_type_matdf))
## make heatmap plot
colors = c(seq(-0.3,0.3,length=16))
my_palette <- rev(colorRampPalette(brewer.pal(6,"YlOrRd"))(n = 15))
hm2_call = heatmap.2(cell_type_mat,col=my_palette,breaks=colors,density.info="none",
                     trace="none",Rowv=F,Colv=T,dendrogram="col",symm=F,labRow=cell_types,
                     labCol=colnames(mat_markers),margins=c(0.5*(dim(cell_type_mat)[2]/dim(cell_type_mat)[1]),2),
                     scale="none",cexRow=1,cexCol=1,rowsep=c(0:20),colsep=c(0:33),sepcolor=NA,sepwidth=c(0.0001,0.0001))
                           
                           

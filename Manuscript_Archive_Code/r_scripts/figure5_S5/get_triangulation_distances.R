# Description:
# For each cell, get the deldir cells and interactions in a dataframe per FOV

# Input parameters
# df_input: dataframe containing unique cell id, X position, Y position, cell type annotation, and region FOV
# ID: unique cell id (there should be no duplicates in the dataframe)
# X position: X location
# Y position: Y location
# cell_type: manually annotated cell type/comparison groups
# region: imaging region/FOV

get_triangulation_distances <- function(df_input, id, x_pos, y_pos, cell_type, region, num_cores = NULL){
  
  # get unique regions
  unique_regions = unique(df_input[[region]])
  
  # Set up parallelization
  library(doSNOW)
  library(foreach)
  if(is.null(num_cores)){
    num_cores <- detectCores() - 1 # get available cores, minus one
  }
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(unique_regions), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Select only necessary columns
  df_input <- df_input %>%
    select(!!as.symbol(id), !!as.symbol(x_pos), !!as.symbol(y_pos), !!as.symbol(cell_type), !!as.symbol(region))
  
  # Calculate triangulation distances
  triangulation_distances <- foreach(i = 1:length(unique_regions), .packages = c("deldir", "tidyverse"), .combine = "rbind", .options.snow = opts)%dopar%{
    # SUBSET DATASET
    subset <- df_input %>%
      dplyr::filter(!!as.symbol(region) == unique_regions[i]) %>%
      mutate(uniqueID = paste0(!!as.symbol(id), "-",
                               !!as.symbol(x_pos), "-",
                               !!as.symbol(y_pos)),
             XYcellID = paste0(!!as.symbol(x_pos),"_", !!as.symbol(y_pos)))
    
    # COMPUTE RDELAUN DISTANCES
    vtress = deldir(subset[[x_pos]], subset[[y_pos]])
    rdelaun_result = vtress$delsgs
    rdelaun_result = rdelaun_result %>%
      mutate(cell1ID = paste0(x1, "_", y1),
             cell2ID = paste0(x2, "_", y2))
    
    # ANNOTATE RDELAUN  RESULTS
    rdelaun_result <- left_join(rdelaun_result, subset, by = c("cell1ID" = "XYcellID")) %>%
      rename(celltype1 = !!as.symbol(cell_type)) %>%
      select(-!!as.symbol(x_pos), -!!as.symbol(y_pos), -!!as.symbol(region), -uniqueID)
    rdelaun_result <- left_join(rdelaun_result, subset, by = c("cell2ID" = "XYcellID")) %>%
      rename(celltype2 = !!as.symbol(cell_type)) %>%
      select(x1, y1, celltype1, !!as.symbol(paste0(id, ".x")),
             x2, y2, celltype2, !!as.symbol(paste0(id, ".y")),
             !!as.symbol(region)) %>%
      mutate(distance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
      select(!!as.symbol(region), 
             !!as.symbol(paste0(id, ".x")), celltype1, x1, y1,
             !!as.symbol(paste0(id, ".y")), celltype2, x2, y2,
             distance)
    colnames(rdelaun_result) <- c(region, 
                                  "celltype1_index", "celltype1", "celltype1_X", "celltype1_Y",
                                  "celltype2_index", "celltype2", "celltype2_X", "celltype2_Y", 
                                  "distance")
    
    rm(subset); rm(vtress)
    return(rdelaun_result)
  }
  close(pb)
  stopCluster(cl)
  
  return(triangulation_distances)
}

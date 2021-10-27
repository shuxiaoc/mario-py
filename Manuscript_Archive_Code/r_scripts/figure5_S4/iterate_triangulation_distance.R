# Description:
# Random permutation version of triangulation distances to evaluate how the observed distances compare to an expected normal distribution.

# Input parameters
# df_input: dataframe containing unique cell id, X position, Y position, cell type annotation, and region FOV
# ID: unique cell id (there should be no duplicates in the dataframe)
# X position: X location
# Y position: Y location
# cell_type annotated cell type/comparison groups
# region: imaging region/FOV

# Default number of iterations = 1000

iterate_triangulation_distance <- function(df_input, num_iterations = 1000, id, x_pos, y_pos, cell_type, region, num_cores = NULL){
  # Get unique regions
  unique_regions = unique(df_input[[region]])
  
  # Setup parallelization
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
  
  # Subset dataframe with only necessary columns
  df_input <- df_input %>%
    select(!!as.symbol(id), !!as.symbol(x_pos), !!as.symbol(y_pos), !!as.symbol(cell_type), !!as.symbol(region))
  
  # START CALCULATING CELL DISTANCES  
  output_results <- foreach(i = 1:length(unique_regions), .packages = c("deldir", "tidyverse"), .combine = "rbind", .options.snow = opts)%dopar%{
    y <- list()
    for (permutations in 1:num_iterations){
      subset <- df_input %>%
        dplyr::filter(!!as.symbol(region) == unique_regions[i]) %>%
        mutate(uniqueID = paste0(!!as.symbol(id), "-",
                                 !!as.symbol(x_pos), "-",
                                 !!as.symbol(y_pos)),
               XYcellID = paste0(!!as.symbol(x_pos),"_", !!as.symbol(y_pos)))
      # Randomize cell type annotations
      shuffled_annotations <- data.frame(subset[[cell_type]])
      set.seed(permutations + 1234) # change seed with every permutation
      rows <- sample(nrow(shuffled_annotations)) # get new row order
      shuffled_annotations <- data.frame(shuffled_annotations[rows,])
      colnames(shuffled_annotations) <- c("random_annotations")
      
      shuffled_subset <- cbind(subset, shuffled_annotations) %>%
        select(-!!as.symbol(cell_type))
      
      rm(subset); rm(rows); rm(shuffled_annotations)
      
      # COMPUTE RDELAUN DISTANCES
      vtress = deldir(shuffled_subset[[x_pos]], shuffled_subset[[y_pos]])
      rdelaun_result = vtress$delsgs
      rdelaun_result = rdelaun_result %>%
        mutate(cell1ID = paste0(x1, "_", y1),
               cell2ID = paste0(x2, "_", y2))
      
      # ANNOTATE RDELAUN  RESULTS
      rdelaun_result <- left_join(rdelaun_result, shuffled_subset, by = c("cell1ID" = "XYcellID")) %>%
        rename(celltype1 = random_annotations) %>%
        select(-!!as.symbol(x_pos), -!!as.symbol(y_pos), -!!as.symbol(region), -!!as.symbol(id))
      rdelaun_result <- left_join(rdelaun_result, shuffled_subset, by = c("cell2ID" = "XYcellID")) %>%
        rename(celltype2 = random_annotations) %>%
        select(x1, y1, x2, y2, celltype1, celltype2, !!as.symbol(region))
      
      rm(shuffled_subset); rm(vtress)
      
      # CALCULATE DISTANCE BETWEEN THE RDELAUN RESULTS
      rdelaun_result <- rdelaun_result %>%
        mutate(distance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
        group_by(celltype1, celltype2) %>%
        summarize(mean_distance = mean(distance)) %>%
        mutate(permutation_num = permutations,
               temp_name = unique_regions[i])
      colnames(rdelaun_result) <- c("celltype1", "celltype2", "mean_distance", "permutation_num", region)
      
      y <- rbind(y, rdelaun_result)
    }
    y <- y %>%
      spread(key = permutation_num, value = mean_distance)
    return(y)
  }
  close(pb)
  stopCluster(cl)
  return(output_results)
}
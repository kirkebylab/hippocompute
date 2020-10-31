# Load packages ----
library(dplyr)
library(tidyr)



load_all_reference_data <- function() {
  # return results in list
  dsets <- list()
  dsets[[ "h9_v0" ]] <- read.csv("data/h9_v0.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "rc17_v0" ]] <- read.csv("data/rc17_v0.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  return(dsets)
}


process_reference_data <- function(data_raw, hk_genes) {
    # load h9 reference values
    
    # extract h9 housekeeping genes and reference ct values
    hkg_list <- list()
    for (g in hk_genes) {
        # TODO: investigate any side-effect, i.e. is the right hk gene chosen first?
        v <- data_raw[data_raw$gene == g,]$ct[1] # get first ct value as ref
        hkg_list[[ g ]] <- v # store ct value
    }
    
    # compute delta.ct for each housekeeping gene
    delta_list = list()
    for (g in names(hkg_list)) {
        data <- data.frame(data_raw) # copy dataframe
        v <- hkg_list[[ g ]] # get housekeeping ct value
        
        data <- data %>%
            mutate(ct=ct - v) %>%
            group_by(gene) %>%
            summarise(avg.delta.ct=mean(ct))
        
        delta_list[[ g ]] <- data
    }
    
    return(delta_list)
}


process_plate <- function(in_file, col_names, row_names, reference_data, house_keeping_genes) {
  # in_file will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  if (is.null(in_file))
    return(NULL)
  
  # TODO: get header
  # header <- read.table(in_file$datapath, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
  
  tab <- read.csv(in_file$datapath, header=TRUE, sep='\t', skip=1)
  tab <- dplyr::mutate(tab, Pos_x = substr(Pos, 2, 4), Pos_y = substr(Pos, 1,1))
  tab <- tidyr::pivot_wider(tab[,c("Pos_y", "Pos_x", "Cp")], names_from=Pos_x, values_from=Cp)
  tab <- as.data.frame(tab)
  rownames(tab) <- tab[,"Pos_y"]
  tab <- tab[-1] # drop Pos_y col
  
  # TODO: Make separate function to contain logic
  col_labels <- unlist(strsplit(col_names, split="\n|\t| ")) # unlist cast to vector
  row_labels <- unlist(strsplit(row_names, split="\n|\t| ")) # unlsit cast to vector
  
  # stop if too many labels
  if (length(row_labels) > dim(tab)[1]) {
    stop()
  }
  if (length(row_labels) > dim(tab)[2]) {
    stop()
  }
  
  # fill up with "empty" labels, if necessary
  if (length(row_labels) < dim(tab)[1]) {
    diff <- dim(tab)[1] - length(row_labels)
    row_labels <- c(row_labels, replicate(diff, "empty"))
  }
  if (length(col_labels) < dim(tab)[2]) {
    diff <- dim(tab)[2] - length(col_labels)
    row_labels <- c(col_labels, replicate(diff, "empty"))
  }
  
  # Replace 0's with 35
  # Replace NA's with 35
  # make notification, see: https://gallery.shinyapps.io/116-notifications/
  # TODO: move notification?
  mask_zeros <- (tab == 0.)
  mask_nas <- is.na(tab)
  # mask_differences is computed later.
  
  if (any(mask_zeros, na.rm=TRUE)) {
    notification_warning_zeros <<- showNotification(
      "WARNING: 0's detected in raw values. Replacing with 35.",
      duration = 10, 
      closeButton = FALSE,
      type = "warning"
    )
  }
  if (any(mask_nas)) {
    notification_warning_nas <<- showNotification(
      "WARNING: NA's detected in raw values. Replacing with 35.",
      duration = 10,
      closeButton = FALSE,
      type = "warning"
    )
  }
  
  tab[mask_zeros | mask_nas] <- 35.
  
  # save raw ct values for qc
  tab_raw_ct <- tab
  rownames(tab_raw_ct) <- paste0(rownames(tab), '_', row_labels)
  colnames(tab_raw_ct) <- paste0(colnames(tab), '_', col_labels)
  
  # infer which axis contains replicates labels for group by
  if (any(duplicated(col_labels)) & !any(duplicated(row_labels))) {
    # if replicates are in columns
    # transpose tab, and add replicates as new col: "sample"
    # R can only groupby based on column values - not column names
    
    tab <- as.data.frame(t(tab))
    colnames(tab) <- row_labels # since tab is transposed, row_labs are now col_labs
    tab <- cbind(tab, sample=col_labels)
  } else if (any(duplicated(row_labels)) & !any(duplicated(col_labels))) {
    colnames(tab) <- col_labels
    tab <- cbind(tab, sample=row_labels)
  } else {
    # throw error
    stop("Could not infer which axis contains replicate samples or primers. Please check that only one axis contains replicate labels.")
  }
  
  # compute standard deviation per replicate
  mask_std <- tab %>%
   dplyr::group_by(sample) %>%
   dplyr::transmute(dplyr::across(.fns=sd, na.rm=TRUE)) %>%
   as.data.frame()

  mask_std <- mask_std[,-1] # drop sample col
  mask_std <- (mask_std > 2.) # TODO: change hard coded threshold. Ask AK.
  
  if (any(mask_std)) { # TODO: place this better
   notification_warning_std <<- showNotification(
     "WARNING: Some replicates have a high standard deviation.",
     duration = 10,
     closeButton = FALSE,
     type = "warning"
   )
  }
  
  
  # compute average Ct intensity per replicate
  sample_order <- unique(tab$sample) # preserve input row order
  tab$sample <- factor(tab$sample, levels = sample_order)
  tab <- tab %>%
    group_by(sample) %>% 
    summarise(across(.fns=mean, na.rm=TRUE)) %>%
    as.data.frame()
  rownames(tab) <- tab[,"sample"]
  tab <- tab[,-1] # drop sample col
  
  # check that primers are in cols and samples in rows
  if (any(rownames(tab) %in% reference_data[[ 1 ]]$gene)) {
    tab <- as.data.frame(t(tab)) # transpose table
  }
  
  # compute deltas
  # for each house-keeping primer in data: compute delta for all samples
  # i.e. delta = sample1_primer1 - sample1_ACTB
  house_keeping_genes <- names(reference_data) # get names from list of df's

  power_acc <- NULL
  for (g in house_keeping_genes) {
    delta_sample <- data.frame(tab, check.names=FALSE) # check.names prevents R changing '-' to '.'
    v <- delta_sample[,g] # get housekeeping gene ct vector
    delta_sample <- delta_sample - v # subtract houskeeping across data
    
    # calculate power
    # for each house-keeping primer in data: compute power in relation to h9 reference
    # i.e. power = 2 ^ -(sample1_primer1avg_actbdelta - h9_primer1avg_actbdelta)
    # N.B. no reference results in NA
    delta_h9 <- reference_data[[ g ]]
    delta_h9 <- data.frame(delta_h9[delta_h9$gene %in% colnames(delta_sample),])
    rownames(delta_h9) <- delta_h9$gene
    delta_h9 <- delta_h9[colnames(delta_sample),] # get h9 avg ct vals in correct order
    diff <- sweep(delta_sample, 2, delta_h9$avg.delta.ct, '-') # subtract h9 avg.delta.ct from each row
    power <- 2 ** -diff
    
    if (is.null(power_acc)) {
      power_acc <- power
    } else {
      power_acc <- power_acc + power
    }
  }
  
  # calculate mean_power
  # i.e. mean_power = avg(sample1_primer1_actbpower, sample1_primer1_gapdhpower)
  mean_power <- power_acc / length(house_keeping_genes)
  
  # return results in list
  results <- list()
  results[[ "data_processed "]] <- mean_power
  results[[ "data_raw_ct" ]] <- tab_raw_ct
  results[[ "mask_zeros" ]] <- mask_zeros
  results[[ "mask_nas" ]] <- mask_nas
  results[[ "mask_std" ]] <- mask_std
  return(results)
}

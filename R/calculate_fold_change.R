
calculate_fold_change <- function(df, reference_data) {
  # it should cbind or rbind the plates, depending on inferred shape
  # then add col_labels and row_labels
  # add replicate_name (sample) column
  # do the computation as usual
  
  # ensure colnames are unique
  if (any(duplicated(colnames(df)))) {
    colnames(df) <- make.unique(colnames(df), sep='.')
  }
  
  # compute average Ct intensity per replicate
  df <- df %>%
    group_by(sample) %>% 
    summarise(across(.fns=mean, na.rm=TRUE), .groups="drop") %>%
    as.data.frame()
  rownames(df) <- df[,"sample"]
  df <- df[,-1] # drop sample col
  # check that primers are in cols and samples in rows
  if (any(rownames(df) %in% reference_data[[ 1 ]]$gene)) {
    df <- as.data.frame(t(df)) # transpose table
  }
  
  # compute deltas
  # for each house-keeping primer in data: compute delta for all samples
  # i.e. delta = sample1_primer1 - sample1_ACTB
  house_keeping_genes <- names(reference_data) # get names from list of df's
  
  power_acc <- NULL
  for (g in house_keeping_genes) {
    # check that hk gene is in colnames
    if (! (g %in% colnames(df))) {
      stop(paste("Housekeeping gene ", g, " not found. Check selected housekeeping genes."))
    }
    
    delta_sample <- data.frame(df, check.names=FALSE) # check.names prevents R changing '-' to '.'
    
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
  
  return(mean_power)
}


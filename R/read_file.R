
read_file <- function(name, file, col_labels, row_labels, replicates_in_cols=NULL) {
  # file input will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  
  # stop if wrong number of labels
  if (length(col_labels) != 24) {
    stop()
  }
  if (length(row_labels) != 16) {
    stop()
  }
  # infer which axis contains replicate labels, if not specified
  if (is.null(replicates_in_cols)) {
    replicates_in_cols <- any(duplicated(col_labels))
  }
  
  # read txt file, but skip header. TODO: get later for metadata
  df <- read.csv(file, header=TRUE, sep='\t', skip=1)
  
  # extract plate coordinates and pivot dataframe
  df <- dplyr::mutate(df, Pos_x = substr(Pos, 2, 4), Pos_y = substr(Pos, 1,1))
  df <- tidyr::pivot_wider(df[,c("Pos_y", "Pos_x", "Cp")], names_from=Pos_x, values_from=Cp)
  df <- as.data.frame(df)
  rownames(df) <- df[,"Pos_y"]
  df <- df[-1] # drop Pos_y col
  
  # Replace 0's with 35
  # Replace NA's with 35
  mask_zeros <- (df == 0.)
  mask_nas <- is.na(df)
  
  df[mask_zeros | mask_nas] <- 35.
  
  # save raw ct values for qc
  df_raw <- df
  rownames(df_raw) <- paste0(rownames(df), '_', row_labels)
  colnames(df_raw) <- paste0(colnames(df), '_', col_labels)
  
  
  if (replicates_in_cols) {
    # if replicates are in columns
    # transpose df, and add replicates as new col: "sample"
    # R can only groupby based on column values - not column names
    df <- as.data.frame(t(df))
    # colnames(df) <- row_labels # since df is transposed, row_labs are now col_labs
    df <- cbind(df, sample=col_labels) # add sample column with labels for group by
  } else {
    # colnames(df) <- col_labels # add col_labels
    df <- cbind(df, sample=row_labels)
  }
  
  # compute standard deviation per replicate
  mask_absdiff <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::transmute(dplyr::across(.fns=mean, na.rm=TRUE)) %>%
    as.data.frame()
  
  mask_absdiff <- mask_absdiff[,-1] # drop sample col
  mask_absdiff <- abs(mask_absdiff - df[,-dim(df)[2]])
  mask_absdiff <- mask_absdiff >= .5 # TODO: evaluate threshold
  
  # transpose to shape 24x16 if necessary
  if (replicates_in_cols) {
    mask_absdiff <- t(mask_absdiff)
  }
  
  # For 2 samples. a SD of .3 is equivalent to a difference from the mean of .5,
  # which works with 2 samples, but not 3
  # TODO: calculate threshold depending on number of replicates. Discuss with AK.
  mask_absdiff <- (mask_absdiff > .5)
  
  # return results in list
  results <- list()
  results[[ "name" ]] <- name
  results[[ "df" ]] <- df_raw
  results[[ "mask_zeros" ]] <- mask_zeros
  results[[ "mask_nas" ]] <- mask_nas
  results[[ "mask_absdiff" ]] <- mask_absdiff
  return(results)
}

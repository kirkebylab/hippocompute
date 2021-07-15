
process_plates <- function(files, col_names, row_names, reference_data, replicates_in_cols=NULL) {
  # files will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  print("process_plates()")
  print(files)
  print(col_names)
  print(row_names)
  print(head(reference_data))
  print(replicates_in_cols)
  
  # -- check input ------------------------------------------------------------
  if (is.null(files))
    return(NULL)
  
  col_labels <- unlist(strsplit(col_names, split="\\s+")) # unlist cast to vector
  row_labels <- unlist(strsplit(row_names, split="\\s+")) # unlsit cast to vector
  
  # infer which axis contains replicate labels, if not specified
  if (is.null(replicates_in_cols)) {
    replicates_in_cols <- any(duplicated(col_labels))
  }
  
  # check number of labels
  if ((length(col_labels) %% 24) != 0) {
    stop("Number of column labels is not a multiple of 24.")
  }
  if ((length(row_labels) %% 16) != 0) {
    stop("Number of column labels is not a multiple of 16.")
  }
  
  
  # -- read files and compute qc -----------------------------------------------
  
  # iterate over files and process them
  n_plates_x <- as.integer(length(col_labels) / 24)
  n_plates_y <- as.integer(length(row_labels) / 16)
  plates <- list()
  df_row <- list()
  for (i in 1:n_plates_y) {
    ix <- (i-1)*16
    plate_rows <- row_labels[(ix+1):(ix+16)]
    df_col <- list()
    for (j in 1:n_plates_x) {
      jx <- (j-1)*24
      plate_cols <- col_labels[(jx+1):(jx+24)]
      px <- i+j-1
      plates[[px]] <- read_file(files$name[px], files$datapath[px], plate_cols, plate_rows, replicates_in_cols)
      
      df <- plates[[px]]$df
      if (replicates_in_cols) {
        # if transposing, doesn't go here
        df <- as.data.frame(t(df))
        colnames(df) <- plate_rows
      } else {
        colnames(df) <- plate_cols
      }
      
      df_col[[j]] <- df
    }
    if (replicates_in_cols) {
      df_row[[i]] <- do.call(rbind, df_col)
    } else {
      df_row[[i]] <- do.call(cbind, df_col)
    }
    
  }
  
  # transpose df if needed
  if (replicates_in_cols) {
    df <- do.call(cbind, df_row)
    # for some reason the plate is built in the wrong direction
    # transpose df, and add replicates as new col: "sample"
    # R can only groupby based on column values - not column names
    df <- cbind(df, sample=col_labels) # add sample column with labels for group by
  } else {
    df <- do.call(rbind, df_row)
    df <- cbind(df, sample=row_labels)
  }
  
  # preserve input row order
  sample_order <- unique(df$sample)
  df$sample <- factor(df$sample, levels = sample_order)
  
  
  # -- compute fold change ----------------------------------------------------
  df_foldchange <- calculate_fold_change(df, reference_data)
  
  # return results in list
  results <- list()
  results[[ "fold_change" ]] <- df_foldchange
  results[[ "raw" ]] <- plates
  return(results)
}


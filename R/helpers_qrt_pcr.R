# Load packages ----
library(dplyr)
library(tidyr)


ceiling_to_multiple <- function(number, multiple) {
    # ceil number to nearest multiple
    return((floor(number / multiple) * multiple) + multiple)
}


load_all_reference_data <- function() {
  # return results in list
  dsets <- list()
  dsets[[ "h9_v0" ]] <- read.csv("data/h9_v0.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9_v1" ]] <- read.csv("data/h9_v1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
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
            summarise(avg.delta.ct=mean(ct), .groups="drop")
        
        delta_list[[ g ]] <- data
    }
    
    return(delta_list)
}


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
  mask_std <- df %>%
   dplyr::group_by(sample) %>%
   dplyr::transmute(dplyr::across(.fns=sd, na.rm=TRUE)) %>%
   as.data.frame()

  mask_std <- mask_std[,-1] # drop sample col
  
  # transpose to shape 24x16 if necessary
  if (replicates_in_cols) {
    mask_std <- t(mask_std)
  }
  
  # For 2 samples. a SD of .3 is equivalent to a difference from the mean of .5,
  # which works with 2 samples, but not 3
  # TODO: calculate threshold depending on number of replicates. Discuss with AK.
  mask_std <- (mask_std > .5)

  # return results in list
  results <- list()
  results[[ "name" ]] <- name
  results[[ "df" ]] <- df_raw
  results[[ "mask_zeros" ]] <- mask_zeros
  results[[ "mask_nas" ]] <- mask_nas
  results[[ "mask_std" ]] <- mask_std
  return(results)
}


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


process_plates <- function(files, col_names, row_names, reference_data, replicates_in_cols=NULL) {
  # files will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.

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
              df <- as.data.frame(t(df))
              colnames(df) <- plate_rows
          } else {
              colnames(df) <- plate_cols
          }

          df_col[[j]] <- df
      }
      df_row[[i]] <- do.call(cbind, df_col)
  }

  df <- do.call(rbind, df_row)

  # transpose df if needed
  if (replicates_in_cols) {
    # transpose df, and add replicates as new col: "sample"
    # R can only groupby based on column values - not column names
    df <- cbind(df, sample=col_labels) # add sample column with labels for group by
  } else {
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


make_datatables_ct <- function(plates) {
  # process a list of plates (16x24 raw ct values)
  # each plate has masks that are used to apply coloring to the datatable
  dt_list <- list()
  for (i in 1:length(plates)) {
      name_obj <- paste0("plate_", i)
      p <- plates[[i]]

      df <- p$df
      mask <- matrix(0L, nrow = dim(df)[1], ncol = dim(df)[2])
      mask[p$mask_std] <- 1
      mask[p$mask_zeros] <- 2
      mask[p$mask_nas] <- 3

      # set table cell color according to replacement, see:
      # https://stackoverflow.com/questions/50798941/r-shiny-rendertable-change-cell-colors
      # https://stackoverflow.com/questions/42569302/format-color-of-shiny-datatable-dt-according-to-values-in-a-different-dataset
      df <- datatable(
        data = cbind(df,mask),
        caption = p$name,
        extensions = c("Buttons"),
        options=list(
          dom = 'Bfrtip',
          pageLength = 25,
          columnDefs = list(list(visible=FALSE, targets=c((1+ncol(df)):(ncol(df)+ncol(mask))))),
          buttons = list('copy',
                        'csv',
                        list(extend = 'excel', filename="hippocompute_qpcr_ct", title = NULL),
                        'colvis')
        ),
        selection = "single"
      ) %>%
        formatStyle(
          1:ncol(df),
          valueColumns=(1+ncol(df)):(ncol(df)+ncol(mask)),
          backgroundColor=styleEqual(c(1,2,3), c("khaki", "lightsalmon", "lightcoral"))
        ) %>%
        formatRound(1:ncol(df), 2)

      dt_list[[name_obj]] <- df
  }
  return(dt_list)
}


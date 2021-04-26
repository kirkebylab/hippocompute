
make_datatables_ct <- function(plates) {
  # process a list of plates (16x24 raw ct values)
  # each plate has masks that are used to apply coloring to the datatable
  if (length(plates) <= 0) {
    stop(paste("No plates to make datatables from."))
  }
  
  dt_list <- list()
  for (i in 1:length(plates)) {
      name_obj <- paste0("plate_", i)
      p <- plates[[i]]

      df <- p$df
      mask <- matrix(0L, nrow = dim(df)[1], ncol = dim(df)[2])
      mask[p$mask_absdiff] <- 1
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


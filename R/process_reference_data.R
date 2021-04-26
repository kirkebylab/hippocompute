
process_reference_data <- function(data_raw, hk_genes) {
  # load h9 reference values
  # TODO: refactor this function - 2x loops not needed + update names c.f. article
  # extract h9 housekeeping genes and reference ct values
  hkg_list <- list()
  for (g in hk_genes) {
    v <- mean(data_raw[data_raw$gene == g,]$ct)
    hkg_list[[ g ]] <- v # store ct value
  }
  
  # compute delta.ct for each housekeeping gene
  delta_list = list()
  for (g in names(hkg_list)) {
    data <- data.frame(data_raw) # copy dataframe
    v <- hkg_list[[ g ]] # get housekeeping ct value
    
    data <- data %>%
      mutate(ct=ct - v) %>% # subtract avg of hk gene
      group_by(gene) %>%    # group by gene
      summarise(avg.delta.ct=mean(ct), .groups="drop") # compute mean
    
    delta_list[[ g ]] <- data
  }
  
  return(delta_list)
}


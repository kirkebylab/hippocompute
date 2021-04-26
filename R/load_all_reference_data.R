
load_all_reference_data <- function() {
  # return results in list
  dsets <- list()
  dsets[[ "h9_v0_2018" ]] <- read.csv("data/h9_v0_2018.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9_v0_2020" ]] <- read.csv("data/h9_v0_2020.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9_v1_2020" ]] <- read.csv("data/h9_v1_2020.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "rc17_v1_2020" ]] <- read.csv("data/rc17_v1_2020.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9-rc17_v1_2020" ]] <- rbind(dsets[[ "h9_v1_2020" ]], dsets[[ "rc17_v1_2020" ]])
  return(dsets)
}



load_all_reference_data <- function() {
  # return results in list
  dsets <- list()
  dsets[[ "h9_v0_2018" ]] <- read.csv("data/h9_v0_2018.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9_v0_2020" ]] <- read.csv("data/h9_v0_2020.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9_v1_2020" ]] <- read.csv("data/h9_v1_2020.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "rc17_v1_2020" ]] <- read.csv("data/rc17_v1_2020.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
  dsets[[ "h9-rc17_v1_2020" ]] <- rbind(dsets[[ "h9_v1_2020" ]], dsets[[ "rc17_v1_2020" ]])
  
  #additional datasets
  #250212 added the kolf2-1_2025 ref dataset (alrik)
  dsets[["kolf2-1_2025"]] <- read.csv("data/kolf2-1_2025.tsv", sep='\t', stringsAsFactors=FALSE)
  
  #250214 added the UKBI011_A2_2025 and UKBi011_A_184_E11 refs dataset (alrik)
  dsets[["UKBI011_A2_2025"]] <- read.csv("data/UKBI011_A2_2025.tsv", sep='\t', stringsAsFactors=FALSE)
  dsets[["UKBi011_A_184_E11_2025"]] <- read.csv("data/UKBi011_A_184_E11_2025.tsv", sep='\t', stringsAsFactors=FALSE)
  
  return(dsets)
}


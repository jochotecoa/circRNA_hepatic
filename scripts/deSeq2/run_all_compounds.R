args = commandArgs(trailingOnly=TRUE)

library(dplyr)

if (length(args)==0) {
  stop("At least one argument must be supplied (input folder).n", call.=FALSE)
} 

if (args[1] == 'circ' ) {
  input_folder = path_data = 'data/salmon_circBase/'
  script = 'scripts/deSeq2/deseq2_circ.R'
}

if (args[1] == 'gene') {
  input_folder = path_data = 'data/salmon_circBase/'
  script = 'scripts/deSeq2/deseq2_genes.R'
}

if (args[1] == 'mirna') {
  input_folder = path_data = 'data/mirge2_mirna/'
  script = 'scripts/deSeq2/deseq2_miR.R'
}


compounds = list.files(path = input_folder) %>% 
  gsub(pattern = '_total_quant.sf', replacement = '')

compounds = compounds[!grepl('Con|con', compounds)] %>% 
  gsub(pattern = '.miR.Counts.rds', replacement = '')


for (comp in compounds) {
  source(file = script, local = T, echo = T)
}

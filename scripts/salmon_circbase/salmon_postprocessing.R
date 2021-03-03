library(dplyr)
source('scripts/functions/functions_JOA.R')
setwd('/ngs-data/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase')
comps = list.dirs(recursive=F)

for (comp in comps) {
  setwd(comp)
  if (length(list.files(pattern='quant.sf', recursive = T)) == 0) {
    print(paste(comp, 'has no quantification files'))
    setwd('../')
    next()
  }
  if (length(list.files(pattern = 'total')) > 0) {
    print(paste(comp, 'already has a merged file'))
    setwd('../')
    next()
  }
  fil = list.dirs(full.names = F, recursive = F)[1]
  cmp = fil %>% strsplit('_') %>% .[[1]] %>% .[1]
  big_file = mergeFiles()
  big_file = big_file[!duplicated(as.list(big_file))]
  colnames(big_file) = colnames(big_file) %>% 
    gsub(pattern = '_quant/quant.sf', replacement = '')
  saveRDS(big_file, file = paste0(cmp, '_total_quant.sf'))
  setwd('../')
}

a = list.files(pattern='total_quant.sf', recursive=T)
for (b in a) {
  file.copy(b, '/share/script/hecatos/juantxo/circRNA_project/data/salmon_circBase/')
}


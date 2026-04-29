source('scripts/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble'))

maindir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/hepatic/'

output_dir = '/share/analysis/hecatos/juantxo/circRNA/salmon_circBase/'

dir.create(output_dir, recursive = T)

salmon_output_dir = '/share/analysis/hecatos/juantxo/mRNA/quant_salmon/Homo_sapiens.GRCh38.cdna.ncrna.circbase/hepatic/'

salmon_output_compounds = list.dirs(path = salmon_output_dir, full.names = T, recursive = F)

for (sal_out_comp in salmon_output_compounds) {
  
  quant_file = mergeFiles(path = sal_out_comp, 
                          files_patt = 'quant.sf', 
                          recursive = T, header = T)
  
  comp = sal_out_comp %>% 
    gsub(salmon_output_dir, '', x = .) %>% 
    gsub('/', '', x = .) %>% 
    tolower()
  
  colnames(quant_file) = colnames(quant_file) %>% 
    gsub(pattern = sal_out_comp, replacement = '') %>% 
    gsub(pattern = '_quant/quant.sf', replacement = '') %>% 
    gsub(pattern = '^.*\\/\\/', replacement = '')
  
  
  saveRDS(object = quant_file, 
          file = paste0(output_dir, comp, '.rds'))
  
}


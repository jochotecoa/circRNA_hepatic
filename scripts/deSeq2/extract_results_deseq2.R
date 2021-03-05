

DEmRNA_fdr = data.frame(ensembl_gene_id = character())

# contrast = list(resultsNames(dds)[2])
contrast = c('conditions', paste(comp, dose, sep = '_'), control)

method_dir = paste('output/data', method, sep = '/')
gene_level_dir = paste(method_dir, 'gene_level', sep = '/')

if (exists('genes')) {
  method_dir = gene_level_dir
}

norm_counts_dir = paste(method_dir, 'normalized_counts', sep = '/')
dr_dir = paste(method_dir, 'DESeqResults', sep = '/')


for (dir_str in c('method_dir', 'norm_counts_dir', 'dr_dir')) {
  dir_str = get(dir_str)
  if (!dir.exists(dir_str)) {
    dir.create(dir_str, recursive = T)
  }
}

norm_counts = DESeq2::counts(object = dds, normalized = T)

# for (contrast in contrast_list) {
if (filtering) {
  res <- results(dds, contrast = contrast)
} else {
  res <- results(dds, 
                 independentFiltering=F, cooksCutoff = F, 
                 contrast = contrast) 
}
DEmRNA_fdr_tmp = data.frame(ensembl_gene_id = rownames(res), 
                            res$padj < 0.05)
contrast_name = paste(contrast[2], 'vs', contrast[3], sep = '_')
colnames(DEmRNA_fdr_tmp)[2] = contrast_name
print(paste("There are", sum(DEmRNA_fdr_tmp[, 2], na.rm = T), 
            "gene(s) that have a padj < 0.05"))
# DEmRNA_fdr = merge.data.frame(DEmRNA_fdr, DEmRNA_fdr_tmp, 
#                               by = 'ensembl_gene_id', all = T)


norm_counts_filename = paste0(norm_counts_dir, '/', contrast_name, '.rds')
saveRDS(object = norm_counts, file = norm_counts_filename)

results_filename = paste0(dr_dir, '/', contrast_name, '.rds')
saveRDS(object = res, file = results_filename)
# }

# DECs = DEmRNA_fdr %>%
#   filter(grepl(':', ensembl_gene_id)) #%>%
# # filter(conditions_UNTR_vs_Therapeutic == T | conditions_UNTR_vs_Toxic == T)
# 
# saveRDS(object = DECs, file = paste0('output/data/', method, '/DECs/',
#                                      contrast_name, '_', comp,
#                                      '__controlling_for_', 'timepoint_',
#                                      'replicate'))

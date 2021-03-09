<<<<<<< HEAD
source('scripts/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'pbmcapply'))

# Extract gene data -------------------------------------------------------

compounds = c('5FU', 'AMI', 'CEL', 'DOC', 'Dox', 'Epi', 'MXT', 'PTX')

pb_2 <- progressBar(max = length(compounds))
i = 1

if (!file.exists('../circ_prot_all.rds')) {
  for (comp in compounds) {
    # source('scripts/postDESeq2_analysis/assign_direction_degs_all.R')
    if (!exists('gene_changes_ori')) {
      gene_changes_ori = readRDS('data/de_changes/gene_changes.rds')
    }
    
    gene_changes = gene_changes_ori %>% 
      dplyr::select(contains(!!comp))
    
    gene_changes = gene_changes[gene_changes[,1] == 0 ,, F]
    gene_changes = gene_changes[gene_changes[,1] == gene_changes[,2],, F]
    
    
    gene_changes = gene_changes[-ncol(gene_changes)]
    colnames(gene_changes) = paste0(comp, '_gene')
    
    # Combine gene data with gene-miRNA table ---------------------------------
    
    
    if (!exists('mirna_transcript_strict_gene')) {
      mirna_transcript_strict_gene = readRDS('data/mti/mrna_mir/mirna_transcript_strict_gene.rds')
    }
    
    mirna_genes_ori = mirna_transcript_strict_gene %>% 
      dplyr::select(miRBase_ID, ensembl_gene_id) %>% 
      unique() 
    
    
    mirna_genes = merge.data.frame(x = mirna_genes_ori, by = 'ensembl_gene_id', 
                                   y = rownames_to_column(gene_changes, 'ensembl_gene_id'))
    
    # Extract miRNA data ------------------------------------------------------
    
    
    # source('scripts/postDESeq2_analysis/assign_direction_demis_all.R')
    if (!exists('mirna_changes_filt_ori')) {
      mirna_changes_filt_ori = readRDS('data/de_changes/mirna_changes_filt.rds')
    }
    mirna_changes_filt = mirna_changes_filt_ori %>% 
      dplyr::select(contains(!!comp)) %>% 
      na.omit()
    
    if (comp == 'MXT') {
      mirna_changes_filt = mirna_changes_filt_ori %>% 
        dplyr::select(contains('Mitoxantrone'))
    }
    
    if (comp == 'PTX') {
      mirna_changes_filt = mirna_changes_filt_ori %>% 
        dplyr::select(contains('Paclitaxel'))
    }
    
    mirna_changes_filt = mirna_changes_filt[mirna_changes_filt[,1] == 0 ,, F]
    mirna_changes_filt = mirna_changes_filt[mirna_changes_filt[,1] == mirna_changes_filt[,2],, F]
    
    mirna_changes_filt = mirna_changes_filt[-ncol(mirna_changes_filt)]
    colnames(mirna_changes_filt) = paste0(comp, '_mirna')
    
    # Combine miRNA data with gene-miRNA data ---------------------------------
    
    mirna_genes = merge.data.frame(x = mirna_genes, by = 'miRBase_ID', 
                                   y = rownames_to_column(mirna_changes_filt, 'miRBase_ID'))
    
    # mirna_genes = mirna_genes[mirna_genes[,3] == -(mirna_genes[,4]),, F]
    
    
    # Combine miRNA and circRNA data ------------------------------------------
    if (!file.exists('../ciRmiR_unique_highsponging_split.rds')) {
      mir_circ_ori = readRDS('../ciRmiR_unique_highsponging.rds')
      
      mir_circ_ori = mir_circ_ori$Var1 %>% 
        as.character %>% 
        sapply(strsplit, ' ') %>% 
        as.data.frame %>% 
        t
    } else {
      mir_circ_ori = readRDS('../ciRmiR_unique_highsponging_split.rds') %>% 
        .[-3]
    }
    colnames(mir_circ_ori) = c('miRBase_ID', 'circBase_ID')
    
    # mir_circ = mir_circ_ori[mir_circ_ori$miRBase_ID %in% mirna_genes$miRBase_ID,, F]
    if (!exists('circrna_changes_ori')) {
      circrna_changes_ori = readRDS('data/de_changes/circrna_changes.rds')
    }
    
    circrna_changes = circrna_changes_ori %>% 
      dplyr::select(contains(!!comp))
    
    circrna_changes = circrna_changes[circrna_changes[,1] != 0 ,, F]
    circrna_changes = circrna_changes[circrna_changes[,1] == circrna_changes[,2],, F]
    
    circrna_changes = circrna_changes[-ncol(circrna_changes)]
    colnames(circrna_changes) = paste0(comp, '_circrna')
    
    
    
    
    # upinallrows = rowSums(x = circrna_changes, na.rm = T) == ncol(circrna_changes) 
    # downinallrows = rowSums(x = circrna_changes, na.rm = T) == -ncol(circrna_changes) 
    # 
    # circrna_changes = upinallrows %>% names() %>% .[upinallrows]
    
    
    mir_circ = circrna_changes %>% 
      rownames_to_column(var = 'circBase_ID') %>% 
      merge.data.frame(x = ., 
                       y = mir_circ_ori, 
                       by = 'circBase_ID')
    
    gene_mir_circ = merge.data.frame(mirna_genes, mir_circ, by = 'miRBase_ID')
    
    if (ncol(gene_mir_circ) == 7) {
      gene_mir_circ = gene_mir_circ[-7]
    }
    
    gene_mir_circ = gene_mir_circ[, c(2,1,5,3,4,6)]
    
    # gene_mir_circ = gene_mir_circ[gene_mir_circ[,4] == gene_mir_circ[,6],, F]
    
    
    
    # Combine proteomics with other omics -------------------------------------
    
    
    library(biomaRt)
    
    mart = openMart2018()
    
    genes = gene_mir_circ$ensembl_gene_id %>% unique()
    
    genes_proteins = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
                           filters = 'ensembl_gene_id', 
                           values = genes, 
                           mart = mart)
    
    # gmc_prot_ori = merge.data.frame(x = gene_mir_circ_all, y = genes_proteins, 
    #                                 by = 'ensembl_gene_id')
    
    protein_changes_ori = readRDS('data/de_changes/protein_changes.rds')
    
    protein_changes = protein_changes_ori %>%
      column_to_rownames() %>% 
      dplyr::select(contains(!!toupper(comp))) %>% 
      rownames_to_column('uniprot_gn')
    
    colnames(protein_changes) = gsub(pattern = 'direff_quant75_data/de_changes/proteomics/|.rds', 
                                     replacement = '', 
                                     x = colnames(protein_changes))
    
    protein_changes = column_to_rownames(protein_changes, 'uniprot_gn') %>% 
      naToZero()
    
    protein_changes = protein_changes[protein_changes[,1] != 0 ,, F]
    protein_changes = protein_changes[protein_changes[,1] == protein_changes[,2],, F]
    
    protein_changes = protein_changes[-ncol(protein_changes)]
    colnames(protein_changes) = paste0(comp, '_protein')
    
    
    
    genes_proteins_changes = protein_changes %>% 
      rownames_to_column(var = 'uniprot_gn') %>% 
      merge.data.frame(x = ., 
                       y = genes_proteins, 
                       by = 'uniprot_gn')
    
    
    gene_mir_circ_prot = merge.data.frame(gene_mir_circ, genes_proteins_changes, 
                                          by = 'ensembl_gene_id') 
    
    
    gene_mir_circ_prot = gene_mir_circ_prot[gene_mir_circ_prot[,6] == gene_mir_circ_prot[,8],] %>% 
      dplyr::select(ensembl_gene_id, miRBase_ID, circBase_ID, uniprot_gn,
                    everything())
    
    
    gene_mir_circ_prot$genemircirc = apply(gene_mir_circ_prot[, 1:4], 1, paste, collapse = '~')
    
    gene_mir_circ_prot = gene_mir_circ_prot[, -1:-4]
    
    if (comp == '5FU') {
      
      gene_mir_circ_prot_all = gene_mir_circ_prot
      saveRDS(object = gene_mir_circ_prot_all, file = '../circ_prot_all.rds')
      
      
    } else {
      gene_mir_circ_prot_all = merge.data.frame(x = gene_mir_circ_prot_all, y = gene_mir_circ_prot, 
                                                by = 'genemircirc', all = T)
      saveRDS(object = gene_mir_circ_prot_all, file = '../circ_prot_all.rds')
    }
    setTxtProgressBar(pb_2, i)
    if (i == length(compounds)) {
      close(pb_2)
    }
    i = i + 1
    
  }
  
  omic_ids = gene_mir_circ_prot_all[,1] %>% 
    strsplit('~') %>% 
    as.data.frame() %>% 
    t()
  
  gene_mir_circ_prot_all = cbind.data.frame(gene_mir_circ_prot_all, omic_ids)
  
  colnames(gene_mir_circ_prot_all)[(ncol(gene_mir_circ_prot_all)-3):ncol(gene_mir_circ_prot_all)] = 
    c('ensembl_gene_id', 'miRBase_ID', 'circBase_ID', 'uniprot_gn')
  
  gene_mir_circ_prot_all = remove_rownames(gene_mir_circ_prot_all)
  
  saveRDS(object = gene_mir_circ_prot_all, file = '../circ_prot_all.rds')
  
} else {
  gene_mir_circ_prot_all = readRDS(file = '../circ_prot_all.rds')
  
}


# Clean sponge data -------------------------------------------------------

circrna_dec_spo_exp = readRDS('data/de_changes/circrna_dec_spo_exp.rds')

# circrna_dec_spo_exp = subset.data.frame(circrna_dec_spo_exp, circrna_dec_spo_exp$median_TPM_de >= 1)

circrna_dec_spo_exp = circrna_dec_spo_exp %>% 
  dplyr::select(-matches('_vs_|rank'))


# Merge both data ---------------------------------------------------------

if (!exists('gene_mir_circ_prot_all')) {
  gene_mir_circ_prot_all = readRDS('../circ_prot_all.rds')
  
}

gene_mir_circ_prot_all = gene_mir_circ_prot_all[, !duplicated(colnames(gene_mir_circ_prot_all))]

# gene_mir_all = gene_mir_circ_prot_all %>% 
#   dplyr::select(ensembl_gene_id, miRBase_ID)

gene_mir_circ_prot_all$mircirc = paste(gene_mir_circ_prot_all$miRBase_ID, 
                                  gene_mir_circ_prot_all$circBase_ID, 
                                  sep = '~')
  # dplyr::select(-miRBase_ID, -ensembl_gene_id)
circrna_dec_spo_exp$mircirc = paste(circrna_dec_spo_exp$miRBase_ID, 
                                    circrna_dec_spo_exp$circBase_ID, 
                                    sep = '~')


gene_mir_circ_prot_all = gene_mir_circ_prot_all %>% 
  dplyr::select(-miRBase_ID, -circBase_ID)


# circrna_quantile = readRDS(file = 'data/quantile_DE_expression/circ_0.083.rds') %>% 
#   rownames_to_column('circBase_ID')
# gene_quantile = readRDS(file = 'data/quantile_expression/quant0.05_gene_norm_counts.rds') %>% 
#   rownames_to_column('ensembl_gene_id')
# gene_col = colnames(gene_quantile)[1]
# mirna_quantile = readRDS(file = 'data/quantile_expression/quant0.05_mirna_norm_counts.rds') %>% 
#   rownames_to_column('miRBase_ID')
# mirna_col = colnames(mirna_quantile)[1]
# protein_quantile = readRDS(file = 'data/quantile_expression/quant0.05_protein_norm_counts.rds') %>% 
#   rownames_to_column('uniprot_gn')
# protein_col = colnames(protein_quantile)[1]
circrna_iq = readRDS("data/interquantile_ratio_expression/iq_ratio_circrna_norm_counts.rds") %>%
  rownames_to_column('circBase_ID')
mirna_iq = readRDS("data/interquantile_ratio_expression/iq_ratio_mirna_norm_counts.rds") %>%
  rownames_to_column('miRBase_ID')
gene_iq = readRDS("data/interquantile_ratio_expression/iq_ratio_gene_norm_counts.rds") %>%
  rownames_to_column('ensembl_gene_id')
protein_iq = readRDS("data/interquantile_ratio_expression/iq_ratio_protein_norm_counts.rds") %>%
  rownames_to_column('uniprot_gn')



gmcp_spo = merge.data.frame(x = gene_mir_circ_prot_all, y = circrna_dec_spo_exp, 
                            by = 'mircirc') %>% 
  dplyr::select(-mircirc, -genemircirc)

gmcp_spo = gmcp_spo %>% 
  merge.data.frame(y = circrna_iq, by = 'circBase_ID', all.x = T) %>% 
  merge.data.frame(y = gene_iq, by = 'ensembl_gene_id', all.x = T) %>% 
  merge.data.frame(y = mirna_iq, by = 'miRBase_ID', all.x = T) %>% 
  merge.data.frame(y = protein_iq, by = 'uniprot_gn', all.x = T) 
  # merge.data.frame(y = circrna_iq, by = 'circBase_ID', all.x = T) 
  

gmcp_spo$quant_DE_circ = gmcp_spo %>% dplyr::select(contains('quant_DE')) %>% rowSums(na.rm = T)


# gmcp_spo$rank_expr_circ = rank(-gmcp_spo$median_TPM_de)
# gmcp_spo$rank_expr_gene = rank(-gmcp_spo$quant0.05_gene_norm_counts_all)
# gmcp_spo$rank_expr_mir = rank(-gmcp_spo$quant0.05_mirna_norm_counts_all)
# gmcp_spo$rank_expr_prot = rank(-gmcp_spo$quant0.05_protein_norm_counts_all)
# gmcp_spo$ranking = rank(90*gmcp_spo$rank_expr_circ + (10/3)*gmcp_spo$rank_expr_gene + (10/3)*gmcp_spo$rank_expr_mir + (10/3)*gmcp_spo$rank_expr_prot)

top_circs_df = gmcp_spo

top_circs_df$ID = apply(top_circs_df[, 1:4], 1, paste, collapse = '~')

compounds = c('5FU', 'AMI', 'CEL', 'DOC', 'Dox', 'Epi', 'MXT', 'PTX')

best_cases = NULL

for (comp in compounds) {

  col_iq_total_comp = paste0('iq_total_', comp)
  if (comp %in% c('Dox', 'Epi')) {
    control = 'DF2|Fluct'
  } else {
    control = 'DMSO_01|Con_DMSO|ConDMSO'
  }

  regex_comp = paste(comp, control, sep = '|')

  top_circs_comp = top_circs_df %>%
    select(matches(regex_comp)) %>%
    select(contains('iq_ratio'))

  stopifnot(ncol(top_circs_comp) == 8)
  top_circs_df[, col_iq_total_comp] = rowSums(top_circs_comp)

  comp_circ = paste0(comp, '_circrna')

  top_circs_df[, comp_circ] = naToZero(top_circs_df[, comp_circ])
  comp_de = top_circs_df[top_circs_df[, comp_circ] != 0, ]
  best_case = comp_de[, col_iq_total_comp] %>% which.min() %>% comp_de$ID[.]
  best_cases = c(best_cases, best_case)

  best_case_df = comp_de[, c('ID', col_iq_total_comp)]
  colnames(best_case_df) = c('ID', 'col_iq_total_comp')
  best_case_df$comp = comp
  if (comp == '5FU') {
    best_cases_df = best_case_df
  } else {
    best_cases_df = rbind.data.frame(best_cases_df, best_case_df)
  }

}

omic_ids = best_cases_df[,1] %>% 
  strsplit('~') %>% 
  as.data.frame() %>% 
  t()

best_cases_df = cbind.data.frame(best_cases_df, omic_ids)

top_circs_df %>% 
  droplevels() %>% 
  saveRDS(file = '../dput_circ_prot.rds')





n_mti = gmcp_spo$circBase_ID %>% 
  as.character %>% 
  table %>% 
  as.data.frame

colnames(n_mti) = c('circBase_ID', 'N_MTI')

gmcp_spo = merge.data.frame(x = gmcp_spo, 
                                       by = 'circBase_ID', 
                                       y = n_mti)

gmcp_spo$rank_n_mti = rank(gmcp_spo$N_MTI) %>% 
  as.numeric()

gmcp_spo[order(gmcp_spo$median_TPM_de, decreasing = T),] %>% head


gmcp_spo$ranking = gmcp_spo$rank_tpm_de + gmcp_spo$rank_n_mti

# 
# # Extract protein data ----------------------------------------------------
# 
# 
# library(biomaRt)
# 
# mart = openMart2018()
# 
# genes = gmc_spo$ensembl_gene_id %>% unique()
# 
# genes_proteins = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
#                        filters = 'ensembl_gene_id', 
#                        values = genes, 
#                        mart = mart)
# 
# # gmc_prot_ori = merge.data.frame(x = gene_mir_circ_prot_all, y = genes_proteins, 
# #                                 by = 'ensembl_gene_id')
# 
# protein_changes_ori = readRDS('data/de_changes/protein_changes.rds')
# 
# protein_changes = protein_changes_ori %>%
#   column_to_rownames() %>% 
#   rownames_to_column('uniprot_gn')
# 
# protein_changes_ori = gsub('direff_quant75_data/de_changes/proteomics/', 
#                        '', 
#                        colnames(protein_changes_ori))
# 
# gmc_prot = merge.data.frame(x = gmc_prot_ori, 
#                             y = protein_changes, 
#                             by = 'uniprot_gn')
# 
# a = gmc_prot[gmc_prot['5FU_The_vs_ConDMSO_gene'] == 1 & gmc_prot['direff_quant75_data/de_changes/proteomics/5FU_The.rds'] == 1,, F]
# 
# colnames(gmc_prot) = c('uniprot_gn', 'ensembl_gene_id', 'miRBase_ID', 'gene_change', 'mirna_change', 'circBase_ID', 'circrna_change', 'N_seeds', 'protein_change')
# 
# gmc_prot$protein_change = naToZero(gmc_prot$protein_change)
# 
# gmc_prot_de = gmc_prot %>% 
#   filter(protein_change != 0,
#          protein_change == circrna_change,
#          gene_change == circrna_change,
#          mirna_change != gene_change)
# 
# 
# 
# 
# gmc_spo_prot = gmc_spo %>% 
#   subset.data.frame(gmc_spo$ensembl_gene_id %in% gmc_prot_de$ensembl_gene_id)
# 
# gmc_spo_prot = gmc_spo_prot[!grepl('N_MTI', colnames(gmc_spo_prot))]
# 
# n_mti = gmc_spo_prot$circBase_ID %>% 
#   as.character %>% 
#   table %>% 
#   as.data.frame
# 
# colnames(n_mti) = c('circBase_ID', 'N_MTI')
# 
# gmc_spo_prot = merge.data.frame(x = gmc_spo_prot, by = 'circBase_ID', 
#                                        y = n_mti)
# 
# 
# # Which genes are related to these top circRNAs? --------------------------
# 
# # gmc_spo %>% 
# #   filter(miRBase_ID == 'hsa-miR-1248',
# #          circBase_ID == 'hsa_circ_0055985') %>% 
# # head()
# #          ensembl_gene_id == 'ENSG00000025800')


top_circs = gmcp_spo[order(gmcp_spo$median_TPM_de, decreasing = T), 'circBase_ID'] %>% 
  unique() %>% 
  .[1:10] %>% 
  droplevels() %>% 
  as.character()

top_mirs = asf %>% 
  dplyr::filter(grepl(top_circs, ))


top1 = gmcp_spo %>% 
  dplyr::filter(circBase_ID == !!top_circs[1]) %>% 
  dplyr::select(circBase_ID, miRBase_ID, ensembl_gene_id, uniprot_gn)



=======
source('scripts/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'tibble', 'pbmcapply'))

# Extract gene data -------------------------------------------------------

compounds = c('5FU', 'AMI', 'CEL', 'DOC', 'Dox', 'Epi', 'MXT', 'PTX')

pb_2 <- progressBar(max = length(compounds))
i = 1

if (!file.exists('../circ_prot_all.rds')) {
  for (comp in compounds) {
    # source('scripts/postDESeq2_analysis/assign_direction_degs_all.R')
    if (!exists('gene_changes_ori')) {
      gene_changes_ori = readRDS('data/de_changes/gene_changes.rds')
    }
    
    gene_changes = gene_changes_ori %>% 
      dplyr::select(contains(!!comp))
    
    gene_changes = gene_changes[gene_changes[,1] == 0 ,, F]
    gene_changes = gene_changes[gene_changes[,1] == gene_changes[,2],, F]
    
    
    gene_changes = gene_changes[-ncol(gene_changes)]
    colnames(gene_changes) = paste0(comp, '_gene')
    
    # Combine gene data with gene-miRNA table ---------------------------------
    
    
    if (!exists('mirna_transcript_strict_gene')) {
      mirna_transcript_strict_gene = readRDS('data/mti/mrna_mir/mirna_transcript_strict_gene.rds')
    }
    
    mirna_genes_ori = mirna_transcript_strict_gene %>% 
      dplyr::select(miRBase_ID, ensembl_gene_id) %>% 
      unique() 
    
    
    mirna_genes = merge.data.frame(x = mirna_genes_ori, by = 'ensembl_gene_id', 
                                   y = rownames_to_column(gene_changes, 'ensembl_gene_id'))
    
    # Extract miRNA data ------------------------------------------------------
    
    
    # source('scripts/postDESeq2_analysis/assign_direction_demis_all.R')
    if (!exists('mirna_changes_filt_ori')) {
      mirna_changes_filt_ori = readRDS('data/de_changes/mirna_changes_filt.rds')
    }
    mirna_changes_filt = mirna_changes_filt_ori %>% 
      dplyr::select(contains(!!comp)) %>% 
      na.omit()
    
    if (comp == 'MXT') {
      mirna_changes_filt = mirna_changes_filt_ori %>% 
        dplyr::select(contains('Mitoxantrone'))
    }
    
    if (comp == 'PTX') {
      mirna_changes_filt = mirna_changes_filt_ori %>% 
        dplyr::select(contains('Paclitaxel'))
    }
    
    mirna_changes_filt = mirna_changes_filt[mirna_changes_filt[,1] == 0 ,, F]
    mirna_changes_filt = mirna_changes_filt[mirna_changes_filt[,1] == mirna_changes_filt[,2],, F]
    
    mirna_changes_filt = mirna_changes_filt[-ncol(mirna_changes_filt)]
    colnames(mirna_changes_filt) = paste0(comp, '_mirna')
    
    # Combine miRNA data with gene-miRNA data ---------------------------------
    
    mirna_genes = merge.data.frame(x = mirna_genes, by = 'miRBase_ID', 
                                   y = rownames_to_column(mirna_changes_filt, 'miRBase_ID'))
    
    # mirna_genes = mirna_genes[mirna_genes[,3] == -(mirna_genes[,4]),, F]
    
    
    # Combine miRNA and circRNA data ------------------------------------------
    if (!file.exists('../ciRmiR_unique_highsponging_split.rds')) {
      mir_circ_ori = readRDS('../ciRmiR_unique_highsponging.rds')
      
      mir_circ_ori = mir_circ_ori$Var1 %>% 
        as.character %>% 
        sapply(strsplit, ' ') %>% 
        as.data.frame %>% 
        t
    } else {
      mir_circ_ori = readRDS('../ciRmiR_unique_highsponging_split.rds') %>% 
        .[-3]
    }
    colnames(mir_circ_ori) = c('miRBase_ID', 'circBase_ID')
    
    # mir_circ = mir_circ_ori[mir_circ_ori$miRBase_ID %in% mirna_genes$miRBase_ID,, F]
    if (!exists('circrna_changes_ori')) {
      circrna_changes_ori = readRDS('data/de_changes/circrna_changes.rds')
    }
    
    circrna_changes = circrna_changes_ori %>% 
      dplyr::select(contains(!!comp))
    
    circrna_changes = circrna_changes[circrna_changes[,1] != 0 ,, F]
    circrna_changes = circrna_changes[circrna_changes[,1] == circrna_changes[,2],, F]
    
    circrna_changes = circrna_changes[-ncol(circrna_changes)]
    colnames(circrna_changes) = paste0(comp, '_circrna')
    
    
    
    
    # upinallrows = rowSums(x = circrna_changes, na.rm = T) == ncol(circrna_changes) 
    # downinallrows = rowSums(x = circrna_changes, na.rm = T) == -ncol(circrna_changes) 
    # 
    # circrna_changes = upinallrows %>% names() %>% .[upinallrows]
    
    
    mir_circ = circrna_changes %>% 
      rownames_to_column(var = 'circBase_ID') %>% 
      merge.data.frame(x = ., 
                       y = mir_circ_ori, 
                       by = 'circBase_ID')
    
    gene_mir_circ = merge.data.frame(mirna_genes, mir_circ, by = 'miRBase_ID')
    
    if (ncol(gene_mir_circ) == 7) {
      gene_mir_circ = gene_mir_circ[-7]
    }
    
    gene_mir_circ = gene_mir_circ[, c(2,1,5,3,4,6)]
    
    # gene_mir_circ = gene_mir_circ[gene_mir_circ[,4] == gene_mir_circ[,6],, F]
    
    
    
    # Combine proteomics with other omics -------------------------------------
    
    
    library(biomaRt)
    
    mart = openMart2018()
    
    genes = gene_mir_circ$ensembl_gene_id %>% unique()
    
    genes_proteins = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
                           filters = 'ensembl_gene_id', 
                           values = genes, 
                           mart = mart)
    
    # gmc_prot_ori = merge.data.frame(x = gene_mir_circ_all, y = genes_proteins, 
    #                                 by = 'ensembl_gene_id')
    
    protein_changes_ori = readRDS('data/de_changes/protein_changes.rds')
    
    protein_changes = protein_changes_ori %>%
      column_to_rownames() %>% 
      dplyr::select(contains(!!toupper(comp))) %>% 
      rownames_to_column('uniprot_gn')
    
    colnames(protein_changes) = gsub(pattern = 'direff_quant75_data/de_changes/proteomics/|.rds', 
                                     replacement = '', 
                                     x = colnames(protein_changes))
    
    protein_changes = column_to_rownames(protein_changes, 'uniprot_gn') %>% 
      naToZero()
    
    protein_changes = protein_changes[protein_changes[,1] != 0 ,, F]
    protein_changes = protein_changes[protein_changes[,1] == protein_changes[,2],, F]
    
    protein_changes = protein_changes[-ncol(protein_changes)]
    colnames(protein_changes) = paste0(comp, '_protein')
    
    
    
    genes_proteins_changes = protein_changes %>% 
      rownames_to_column(var = 'uniprot_gn') %>% 
      merge.data.frame(x = ., 
                       y = genes_proteins, 
                       by = 'uniprot_gn')
    
    
    gene_mir_circ_prot = merge.data.frame(gene_mir_circ, genes_proteins_changes, 
                                          by = 'ensembl_gene_id') 
    
    
    gene_mir_circ_prot = gene_mir_circ_prot[gene_mir_circ_prot[,6] == gene_mir_circ_prot[,8],] %>% 
      dplyr::select(ensembl_gene_id, miRBase_ID, circBase_ID, uniprot_gn,
                    everything())
    
    
    gene_mir_circ_prot$genemircirc = apply(gene_mir_circ_prot[, 1:4], 1, paste, collapse = '~')
    
    gene_mir_circ_prot = gene_mir_circ_prot[, -1:-4]
    
    if (comp == '5FU') {
      
      gene_mir_circ_prot_all = gene_mir_circ_prot
      saveRDS(object = gene_mir_circ_prot_all, file = '../circ_prot_all.rds')
      
      
    } else {
      gene_mir_circ_prot_all = merge.data.frame(x = gene_mir_circ_prot_all, y = gene_mir_circ_prot, 
                                                by = 'genemircirc', all = T)
      saveRDS(object = gene_mir_circ_prot_all, file = '../circ_prot_all.rds')
    }
    setTxtProgressBar(pb_2, i)
    if (i == length(compounds)) {
      close(pb_2)
    }
    i = i + 1
    
  }
  
  omic_ids = gene_mir_circ_prot_all[,1] %>% 
    strsplit('~') %>% 
    as.data.frame() %>% 
    t()
  
  gene_mir_circ_prot_all = cbind.data.frame(gene_mir_circ_prot_all, omic_ids)
  
  colnames(gene_mir_circ_prot_all)[(ncol(gene_mir_circ_prot_all)-3):ncol(gene_mir_circ_prot_all)] = 
    c('ensembl_gene_id', 'miRBase_ID', 'circBase_ID', 'uniprot_gn')
  
  gene_mir_circ_prot_all = remove_rownames(gene_mir_circ_prot_all)
  
  saveRDS(object = gene_mir_circ_prot_all, file = '../circ_prot_all.rds')
  
} else {
  gene_mir_circ_prot_all = readRDS(file = '../circ_prot_all.rds')
  
}


# Clean sponge data -------------------------------------------------------

circrna_dec_spo_exp = readRDS('data/de_changes/circrna_dec_spo_exp.rds')

# circrna_dec_spo_exp = subset.data.frame(circrna_dec_spo_exp, circrna_dec_spo_exp$median_TPM_de >= 1)

circrna_dec_spo_exp = circrna_dec_spo_exp %>% 
  dplyr::select(-matches('_vs_|rank'))


# Merge both data ---------------------------------------------------------

if (!exists('gene_mir_circ_prot_all')) {
  gene_mir_circ_prot_all = readRDS('../circ_prot_all.rds')
  
}



# gene_mir_all = gene_mir_circ_prot_all %>% 
#   dplyr::select(ensembl_gene_id, miRBase_ID)

gene_mir_circ_prot_all$mircirc = paste(gene_mir_circ_prot_all$miRBase_ID, 
                                  gene_mir_circ_prot_all$circBase_ID, 
                                  sep = '~')
  # dplyr::select(-miRBase_ID, -ensembl_gene_id)
circrna_dec_spo_exp$mircirc = paste(circrna_dec_spo_exp$miRBase_ID, 
                                    circrna_dec_spo_exp$circBase_ID, 
                                    sep = '~')


gene_mir_circ_prot_all = gene_mir_circ_prot_all %>% 
  dplyr::select(-miRBase_ID, -circBase_ID)


circrna_quantile = readRDS(file = 'data/quantile_DE_expression/circ_0.083.rds') %>% 
  rownames_to_column('circBase_ID')
gene_quantile = readRDS(file = 'data/quantile_expression/quant0.05_gene_norm_counts.rds') %>% 
  rownames_to_column('ensembl_gene_id')
gene_col = colnames(gene_quantile)[1]
mirna_quantile = readRDS(file = 'data/quantile_expression/quant0.05_mirna_norm_counts.rds') %>% 
  rownames_to_column('miRBase_ID')
mirna_col = colnames(mirna_quantile)[1]
protein_quantile = readRDS(file = 'data/quantile_expression/quant0.05_protein_norm_counts.rds') %>% 
  rownames_to_column('uniprot_gn')
protein_col = colnames(protein_quantile)[1]
circrna_iq = readRDS("data/interquantile_ratio_expression/iq_ratio_circrna_norm_counts.rds") %>% 
  rownames_to_column('circBase_ID')



gmcp_spo = merge.data.frame(x = gene_mir_circ_prot_all, y = circrna_dec_spo_exp, 
                            by = 'mircirc') %>% 
  dplyr::select(-mircirc, -genemircirc)

gmcp_spo = gmcp_spo %>% 
  merge.data.frame(y = circrna_quantile, by = 'circBase_ID', all = T) %>% 
  merge.data.frame(y = gene_quantile, by = 'ensembl_gene_id', all = T) %>% 
  merge.data.frame(y = mirna_quantile, by = 'miRBase_ID', all = T) %>% 
  merge.data.frame(y = protein_quantile, by = 'uniprot_gn', all = T) %>% 
  merge.data.frame(y = circrna_iq, by = 'circBase_ID', all = T) 
  

gmcp_spo$quant_DE_circ = gmcp_spo %>% dplyr::select(contains('quant_DE')) %>% rowSums(na.rm = T)


gmcp_spo$rank_expr_circ = rank(-gmcp_spo$median_TPM_de)
gmcp_spo$rank_expr_gene = rank(-gmcp_spo$quant0.05_gene_norm_counts_all)
gmcp_spo$rank_expr_mir = rank(-gmcp_spo$quant0.05_mirna_norm_counts_all)
gmcp_spo$rank_expr_prot = rank(-gmcp_spo$quant0.05_protein_norm_counts_all)
gmcp_spo$ranking = rank(90*gmcp_spo$rank_expr_circ + (10/3)*gmcp_spo$rank_expr_gene + (10/3)*gmcp_spo$rank_expr_mir + (10/3)*gmcp_spo$rank_expr_prot)

top_circs_df = gmcp_spo #%>% 
  # dplyr::filter(quant_DE_circ > 0,
  #               quant0.05_gene_norm_counts_all > 0)#,
                # quant0.05_mirna_norm_counts_all > 0,
                # quant0.05_protein_norm_counts_all > 0) 

top_circs_df %>% dput

top_circs = top_circs_df$circBase_ID[order(top_circs_df$median_TPM_de)] %>% 
  unique() %>% 
  .[1:10] %>% 
  droplevels()

n_mti = gmcp_spo$circBase_ID %>% 
  as.character %>% 
  table %>% 
  as.data.frame

colnames(n_mti) = c('circBase_ID', 'N_MTI')

gmcp_spo = merge.data.frame(x = gmcp_spo, 
                                       by = 'circBase_ID', 
                                       y = n_mti)

gmcp_spo$rank_n_mti = rank(gmcp_spo$N_MTI) %>% 
  as.numeric()

gmcp_spo[order(gmcp_spo$median_TPM_de, decreasing = T),] %>% head


gmcp_spo$ranking = gmcp_spo$rank_tpm_de + gmcp_spo$rank_n_mti

# 
# # Extract protein data ----------------------------------------------------
# 
# 
# library(biomaRt)
# 
# mart = openMart2018()
# 
# genes = gmc_spo$ensembl_gene_id %>% unique()
# 
# genes_proteins = getBM(attributes = c('ensembl_gene_id', 'uniprot_gn'), 
#                        filters = 'ensembl_gene_id', 
#                        values = genes, 
#                        mart = mart)
# 
# # gmc_prot_ori = merge.data.frame(x = gene_mir_circ_prot_all, y = genes_proteins, 
# #                                 by = 'ensembl_gene_id')
# 
# protein_changes_ori = readRDS('data/de_changes/protein_changes.rds')
# 
# protein_changes = protein_changes_ori %>%
#   column_to_rownames() %>% 
#   rownames_to_column('uniprot_gn')
# 
# protein_changes_ori = gsub('direff_quant75_data/de_changes/proteomics/', 
#                        '', 
#                        colnames(protein_changes_ori))
# 
# gmc_prot = merge.data.frame(x = gmc_prot_ori, 
#                             y = protein_changes, 
#                             by = 'uniprot_gn')
# 
# a = gmc_prot[gmc_prot['5FU_The_vs_ConDMSO_gene'] == 1 & gmc_prot['direff_quant75_data/de_changes/proteomics/5FU_The.rds'] == 1,, F]
# 
# colnames(gmc_prot) = c('uniprot_gn', 'ensembl_gene_id', 'miRBase_ID', 'gene_change', 'mirna_change', 'circBase_ID', 'circrna_change', 'N_seeds', 'protein_change')
# 
# gmc_prot$protein_change = naToZero(gmc_prot$protein_change)
# 
# gmc_prot_de = gmc_prot %>% 
#   filter(protein_change != 0,
#          protein_change == circrna_change,
#          gene_change == circrna_change,
#          mirna_change != gene_change)
# 
# 
# 
# 
# gmc_spo_prot = gmc_spo %>% 
#   subset.data.frame(gmc_spo$ensembl_gene_id %in% gmc_prot_de$ensembl_gene_id)
# 
# gmc_spo_prot = gmc_spo_prot[!grepl('N_MTI', colnames(gmc_spo_prot))]
# 
# n_mti = gmc_spo_prot$circBase_ID %>% 
#   as.character %>% 
#   table %>% 
#   as.data.frame
# 
# colnames(n_mti) = c('circBase_ID', 'N_MTI')
# 
# gmc_spo_prot = merge.data.frame(x = gmc_spo_prot, by = 'circBase_ID', 
#                                        y = n_mti)
# 
# 
# # Which genes are related to these top circRNAs? --------------------------
# 
# # gmc_spo %>% 
# #   filter(miRBase_ID == 'hsa-miR-1248',
# #          circBase_ID == 'hsa_circ_0055985') %>% 
# # head()
# #          ensembl_gene_id == 'ENSG00000025800')


top_circs = gmcp_spo[order(gmcp_spo$median_TPM_de, decreasing = T), 'circBase_ID'] %>% 
  unique() %>% 
  .[1:10] %>% 
  droplevels() %>% 
  as.character()

top_mirs = asf %>% 
  dplyr::filter(grepl(top_circs, ))


top1 = gmcp_spo %>% 
  dplyr::filter(circBase_ID == !!top_circs[1]) %>% 
  dplyr::select(circBase_ID, miRBase_ID, ensembl_gene_id, uniprot_gn)



>>>>>>> 895cd48bf7209b7d4188cedf2a43f9e49e5ad781

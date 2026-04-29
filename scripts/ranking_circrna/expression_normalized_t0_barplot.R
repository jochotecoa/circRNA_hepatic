source('scripts/functions/functions_JOA.R')
forceLibrary(c('dplyr', 'UpSetR', 'tibble', 'pbmcapply'))
iq_ratio <- function(x, prba = F) {
  y = x #%>% 
  # dplyr::select(matches(paste0(col, collapse = '.*')))
  
  # colname = paste0('iq_ratio_', paste0(col, collapse = '_'), collapse = '_')
  # x[, colname] = NA
  ratio_vector = NULL
  if (prba) {  pb_514 <- progressBar(max = nrow(x))}
  for (rwo in rownames(x)) {
    expr_values = y[rwo,]
    expr_values = expr_values %>% unlist() %>% as.numeric()
    expr_values = naToZero(expr_values)
    quantiles = quantile(x = expr_values, probs = c(1/7, 0.5, 7/8))
    ratio = (quantiles[3]-quantiles[1])/quantiles[2]
    ratio = as.numeric(ratio)
    ratio_vector = c(ratio_vector, ratio)
    # x[rwo, colname] = xa
    if (prba) { 
      i = grep(rwo, rownames(x))
      setTxtProgressBar(pb_514, i)
      if (i == nrow(x)) {
        close(pb_514)
      }
    }
  }
  # x[is.na(x[, colname]), colname] = 1000
  ratio_vector[is.na(ratio_vector)] = 1000
  return(ratio_vector)
}

treatm_files = list.files('output/data/salmon_circBase/DESeqResults/', ) %>% 
  subset(x = ., grepl(x = ., 'T0'))
DEcircs_all = data.frame()


pb_2 = progressBar(max = length(treatm_files), style = 'ETA')

for (treatm_file in treatm_files) {
  
  res = 
    readRDS(
      paste0(
        'output/data/salmon_circBase/DESeqResults/', 
        treatm_file
      )
    )
  norm_counts_ori = 
    readRDS(
      paste0(
        'output/data/salmon_circBase/normalized_counts/', 
        treatm_file
      )
    )
  
  treatment = gsub(pattern =  '.rds', replacement = '', x = treatm_file)
  
  norm_counts = norm_counts_ori
  
  comp = colnames(norm_counts)[ncol(norm_counts)] %>%
    strsplit(split = '_') %>%
    sapply('[', 1) #
  
  
  first_treat_col = grepl(comp, colnames(norm_counts))
  # conCols = 1:ncol(cts_control)
  treatmCols = grepl(comp, colnames(norm_counts))
  conCols = !treatmCols
  
  dim_DEs = dim(res)
  
  DEs = res %>% 
    as.data.frame %>% 
    dplyr::filter(padj < 0.1)
  
  norm_counts = norm_counts[rownames(norm_counts) %in% rownames(DEs), ]
  all.equal(rownames(norm_counts), rownames(DEs)) %>% stopifnot()
  dim_DEs = rbind(dim_DEs, dim(DEs)) 
  
  DEs$baseMedian = apply(norm_counts, 1, median)
  
  DEs = DEs %>% 
    dplyr::filter(baseMedian > 0)
  
  norm_counts = norm_counts[rownames(norm_counts) %in% rownames(DEs), ]
  all.equal(rownames(norm_counts), rownames(DEs)) %>% stopifnot()
  dim_DEs = rbind(dim_DEs, dim(DEs)) 
  
  DEs$conFQ = apply(norm_counts[, conCols], 1, quantile, 0.25)
  
  DEs = DEs %>% 
    dplyr::filter(conFQ > 0)
  
  norm_counts = norm_counts[rownames(norm_counts) %in% rownames(DEs), ]
  all.equal(rownames(norm_counts), rownames(DEs)) %>% stopifnot()
  dim_DEs = rbind(dim_DEs, dim(DEs)) 
  
  DEs$treatmFQ = apply(norm_counts[, treatmCols], 1, quantile, 0.25)
  
  DEs = DEs %>% 
    dplyr::filter(treatmFQ > 0)
  
  norm_counts = norm_counts[rownames(norm_counts) %in% rownames(DEs), ]
  all.equal(rownames(norm_counts), rownames(DEs)) %>% stopifnot()
  dim_DEs = rbind(dim_DEs, dim(DEs)) 
  
  
  DEs$conTQ = apply(norm_counts[, conCols], 1, quantile, 0.75)
  DEs$treatmTQ = apply(norm_counts[, treatmCols], 1, quantile, 0.75)
  
  
  DEs = DEs %>% 
    dplyr::filter(treatmTQ < conFQ  | treatmFQ > conTQ,)
  
  norm_counts = norm_counts[rownames(norm_counts) %in% rownames(DEs), ]
  all.equal(rownames(norm_counts), rownames(DEs)) %>% stopifnot()
  dim_DEs = rbind(dim_DEs, dim(DEs)) 
  
  
  DEs$iq_ratio_con = iq_ratio(norm_counts[, conCols])
  
  DEs = DEs %>% 
    dplyr::filter(iq_ratio_con < 1)
  
  norm_counts = norm_counts[rownames(norm_counts) %in% rownames(DEs), ]
  all.equal(rownames(norm_counts), rownames(DEs)) %>% stopifnot()
  dim_DEs = rbind(dim_DEs, dim(DEs)) 
  
  
  DEcircs = DEs[grepl('hsa', rownames(DEs)), ]
  
  clean_rnames = DEcircs %>% rownames() %>% cleanCircNames()
  rownames(DEcircs) = clean_rnames
  stopifnot(all.equal(rownames(DEcircs), clean_rnames))
  
  DEcircs_all = rbind.data.frame(DEcircs_all, rownames_to_column(DEcircs))
  # DEcircs %>% .[sort(abs(.[, 'log2FoldChange']), decreasing = T), ] %>% head()

  
  
  if (exists('DEcircs_all_top')) {
    comp_dose = colnames(norm_counts)[ncol(norm_counts)] %>%
      strsplit(split = '_') %>%
      sapply('[', 1:2) %>% 
      paste(collapse='_')
    
    last_con_col = grep('VPA_Tox_000_3', colnames(norm_counts)) #
    first_treat_col = last_con_col + 1
    
    groups = c(
      rep('control', last_con_col),
      rep('treatment', ncol(norm_counts) - last_con_col)
    )
    
    control_colnames = norm_counts %>% colnames() %>% .[1:last_con_col]
    if (comp %in% c('Dox', 'Epi')) {
      compound_concols = grep(pattern = 'con_DF2', control_colnames)
      comp_concols = control_colnames[compound_concols]
    } else {
      compound_concols = grep(pattern = comp_dose, control_colnames)
      if (length(compound_concols) == 0) {
        compound_concols = grep(pattern = comp, control_colnames)
      }
      comp_concols = control_colnames[compound_concols]
    }
    reordered_cols = c(control_colnames[-compound_concols],
                       control_colnames[compound_concols],
                       colnames(norm_counts)[first_treat_col:ncol(norm_counts)])
    norm_counts = norm_counts[, reordered_cols]
    
    DEcircs = DEcircs[rownames(DEcircs) %in% DEcircs_all_top, ]
    
    for (circ in rownames(DEcircs)) {
      test = norm_counts[grep(circ, rownames(norm_counts)), ]
      
      
      colos = groups %>%
        gsub(pattern = 'treatment', replacement = 'pink') %>%
        gsub(pattern = 'control', replacement = 'gray')
      colos[colnames(norm_counts) %in% comp_concols] = 'lightgray'
      
      plot_dir = '/share/analysis/hecatos/juantxo/circRNA/plots'
      barplot_dir = paste0(plot_dir, '/barplot/hepatic')
      barplot_dir = barplot_dir %>% 
        list.files(pattern = circ, full.names = T)
      dir.create(barplot_dir, recursive = T)
      
      png(filename = paste0(barplot_dir, '/', treatment, '.png'), width = 1024, height = 768)
      barplot(test, las = 2, cex.names = 0.6, col = colos, main = circ)
      dev.off()
      
      boxplot_dir = paste0(plot_dir, '/boxplot/hepatic')
      boxplot_dir = boxplot_dir %>% 
        list.files(pattern = circ, full.names = T)
      
      dir.create(boxplot_dir, recursive = T)
      
      png(filename = paste0(boxplot_dir, '/', treatment, '.png'), width = 1024, height = 768)
      boxplot(test ~ groups)
      dev.off()
      
      
      
    }
    
  }
  
  
  
  setTxtProgressBar(pb_2, grep(treatm_file, treatm_files))
  if(grep(treatm_file, treatm_files) == length(treatm_files)){close(pb_2)}
  
  
}

DEcircs_all = DEcircs_all[order(DEcircs_all$baseMedian, decreasing = T), ]
DEcircs_all = DEcircs_all[!duplicated(DEcircs_all$rowname), ]
DEcircs_all_top = DEcircs_all$rowname[1:100]
DEcircs_all_top_rank = paste(1:100, DEcircs_all_top, sep = '_')

for (circ in DEcircs_all_top_rank) {
  plot_dir = '/share/analysis/hecatos/juantxo/circRNA/plots/'

  barplot_dir = paste0(plot_dir, '/barplot/hepatic/', circ)
  dir.create(barplot_dir, recursive = T)

  boxplot_dir = paste0(plot_dir, '/boxplot/hepatic/', circ)
  dir.create(boxplot_dir, recursive = T)

  # Sys.sleep(1)

}

norm_counts = data.frame()
  

for (treatm_file in treatm_files) {
  
  # res = 
  #   readRDS(
  #     paste0(
  #       'output/data/salmon_circBase/DESeqResults/', 
  #       treatm_file
  #     )
  #   )
  norm_counts_ori = 
    readRDS(
      paste0(
        'output/data/salmon_circBase/normalized_counts/', 
        treatm_file
      )
    )
  
  a = norm_counts_ori[, 1:45] %>% 
    as.data.frame() %>% 
    rownames_to_column %>% 
    dplyr::filter(grepl('hsa', rowname)) 
  b = apply(a[, -1], 1, median)
  
  print(a$rowname[which.max(b)])
  print(max(b))
  
  
}
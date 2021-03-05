# /share/script/hecatos/juantxo/circRNA_project/
# Libraries and functions -------------------------------------------------

source('scripts/functions/functions_JOA.R')
forceLibrary(c('biomaRt', "tximport", "dplyr", "DESeq2", "grid", "ggplot2", 
               "pheatmap", "BiocParallel"))
register(MulticoreParam(20))
getCtsMiRNA <- function(path, comp) {
  quant_file = list.files(path = path, pattern = comp, full.names = T) %>% 
      readRDS()

  quant_file_counts = quant_file  %>% 
    dplyr::select(contains('.fastq'))
  # Rename colnames
  colnames(quant_file_counts) = colnames(quant_file_counts) %>% 
    gsub(pattern = '_R1.fastq', replacement = '') %>% 
    gsub(pattern = '.fastq', replacement = '')
  # Take away samples with low sequencing depth
  quant_file_counts = quant_file_counts %>% 
    dplyr::select(!matches('lowseqdepth'))
  # DESeq2 needs round numbers
  cts = apply(quant_file_counts, 2, round) %>% as.data.frame()
  
  rownames(cts) = quant_file$miRNA
  
  if (any(grepl('miRNAtotal', rownames(cts)))) {
    total_row = grep('miRNAtotal', rownames(cts))
    cts = cts[-total_row, ]
  }
  
  return(cts)
}
getSampleNames <- function(fullSampleName) {
  tablenames = fullSampleName %>% colnames() %>% strsplit('_') %>% 
    as.data.frame %>% t  
  tmpnt = ncol(tablenames) - 1
  smpl = ncol(tablenames) 
  stringnames = paste0(tablenames[, tmpnt], '_', tablenames[, smpl])
  return(stringnames)
}
getTimepoints <- function(fullSampleName) {
  tablenames = fullSampleName %>% colnames() %>% strsplit('_') %>% 
    as.data.frame %>% t  
  tmpnt = ncol(tablenames) - 1
  stringnames = tablenames[, tmpnt]
  return(stringnames)
}
getReplicates <- function(fullSampleName) {
  tablenames = fullSampleName %>% colnames() %>% strsplit('_') %>% 
    as.data.frame %>% t  
  rplct = ncol(tablenames)
  stringnames = tablenames[, rplct]
  return(stringnames)
}
rmDuplicatedColumns <- function(df) {
  df = df[!duplicated(as.list(df))]
  return(df)
}

# Input data --------------------------------------------------------------

stopifnot("No compound specified"= exists('comp'))
stopifnot("No path for the input data specified"= exists('path_data'))

# Parameters

method = path_data %>% strsplit(split = '/') %>% .[[1]] %>% .[2]
setSizeFactorToOne = F # Default: FALSE
filtering = T # Default: TRUE

# Analysis ----------------------------------------------------------------


if (grepl('DAU|DOX|EPI|IDA|CON_0.1_DMSO', toupper(comp))) {
  control = 'Con_Flu_DMSO_Gal'
} else {
  control = 'Con_0.1_DMSO'
}

cts_control = getCtsMiRNA(path = path_data, comp = control) %>% 
  rmDuplicatedColumns() %>% 
  filterSamplesBySeqDepth()
cts_treatm = getCtsMiRNA(path = path_data, comp = comp) %>% 
  rmDuplicatedColumns() %>% 
  filterSamplesBySeqDepth()

if (nrow(cts_control) != nrow(cts_treatm)) {
  cts_control = cts_control[rownames(cts_control) %in% rownames(cts_treatm), ]
  cts_treatm = cts_treatm[rownames(cts_treatm) %in% rownames(cts_control), ]
  
  warning('Some miRNAs are found only in 1 of the compared compounds')
}


cts_the = cts_treatm %>% 
  dplyr::select(contains('The'))

cts_tox = cts_treatm %>% 
  dplyr::select(contains('Tox'))

for (cts_dose_cha in c('cts_The', 'cts_Tox')) {
  
  if (!grepl('con', tolower(comp))) {
    dose = cts_dose_cha %>% gsub(pattern = 'cts_', replacement = '')
    
    cts_dose = get(x = tolower(cts_dose_cha))
    
  } else {
    cts_dose = cts_treatm
  }
  
  if (ncol(cts_control) > ncol(cts_dose)) {
    
    cts_control = cts_control[, 1:ncol(cts_dose), F]
    
  } else if (ncol(cts_control) < ncol(cts_dose)) {
    
    cts_dose = cts_dose[, 1:ncol(cts_control), F]
    
  }
  
  stopifnot(ncol(cts_control) == ncol(cts_dose))
  stopifnot(nrow(cts_control) == nrow(cts_dose))
  
  cts = cbind.data.frame(cts_control, cts_dose)
  
  # Metadata ----------------------------------------------------------------
  
  source(file = 'scripts/deSeq2/metadata_deseq2.R', local = T, echo = T)
  
  # DESeq2 pipeline ---------------------------------------------------------
  
  source(file = 'scripts/deSeq2/deseq2.R', local = T, echo = T)
  
  # Extract results from a DESeq analysis -----------------------------------
  
  source(file = 'scripts/deSeq2/extract_results_deseq2.R', local = T, echo = T)
  
  if (!grepl('con', tolower(comp))) {
    break()
  }
  
}
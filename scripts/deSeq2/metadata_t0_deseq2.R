conditions = repls = NULL

if (!exists('dose')) {
  dose = ''
}

  repls = rep(control, ncol(cts_control))
  conditions = c(conditions, repls)
  repls = rep(paste0(comp, '_', dose), ncol(cts_dose))
  conditions = c(conditions, repls)
  
# timepoints = NULL
# 
# for (cts_table in c('cts_control', 'cts_the', 'cts_tox')) {
#   
#   timepoints = cts_table %>% get() %>% getTimepoints() %>% c(timepoints, .)
#   
# }
# 
# replicates = NULL
# 
# for (cts_table in c('cts_control', 'cts_the', 'cts_tox')) {
#   
#   replicates = cts_table %>% get() %>% getReplicates() %>% c(replicates, .)
#   
# }


coldata = data.frame(row.names = colnames(cts), 
                     conditions = conditions)
                     # timepoints = timepoints,
                     # replicates = replicates)

rownames(coldata) %in% colnames(cts) %>% all() %>% stopifnot()
all(rownames(coldata) == colnames(cts)) %>% stopifnot()

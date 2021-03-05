# DESeq2 pipeline ---------------------------------------------------------

# txi = tximport(files = list.files(pattern = 'quant.sf', recursive = T)[-1], type = "salmon", tx2gene = tx2gene)
dds = DESeqDataSetFromMatrix(countData = cts, 
                             colData = coldata, 
                             design = ~ conditions)

if (setSizeFactorToOne) {
  sizeFactors(dds) <- rep(1, ncol(dds)) # this is the important line
}

print(paste0('Running ', dose, ' dose of ', comp, ' against ', control))

dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

# Independent filtering can be turned off by setting independentFiltering to FALSE.
# cooksCutoff: count outliers


# perform the normalization of the read count 
norm_data <- counts(dds,normalized=TRUE)               

# norm_data_name <- merge(norm_data,rnaseq_id["gene_name"],by="row.names",all.x=TRUE)
# rownames(norm_data_name) <- norm_data_name$Row.names
# norm_data_name <- norm_data_name[,-1]

norm_read_count <- colSums(norm_data)
# barplot(norm_read_count, las=2, cex.names=0.6)
# barplot(log2(norm_read_count), las=2, cex.names=0.6)

## Set the pvalue to 0.05
# res <- results(dds, alpha=0.05, contrast=c("status", "Response" , "No_response"))
# res <- res[order(res$padj),]

# PCA ---------------------------------------------------------------------


# perform a PCA plot. You can change the condition assessed by the pca plot by changing the value between [] from 1 to 3 (corresponding to the column number of the metadata file)
# rld <- vst(dds)
# pca = plotPCA(rld, intgroup = names(colData(dds)))
# filename = paste0('output/plots/PCA_', comp, '.png')
# ggsave(filename = filename, device = "png", plot = pca)



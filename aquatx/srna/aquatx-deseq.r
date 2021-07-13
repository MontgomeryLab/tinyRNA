#!/usr/bin/env Rscript --vanilla

## The module for running DESeq2 for small RNA sequencing data

library(DESeq2)
library(lattice)

#### ---- Get the command line arguments ---- ####

args <- commandArgs(trailingOnly = TRUE)

# Throw an error if there are more args than expected
if (length(args) > 4){
  stop(sprintf("Too many arguments given. Only --input-files and --sample-names are accepted.
        Files and names must be comma-separated without spaces. The bad arguments given were 

        %s

        Did you add a space in a list?", args[5:length(args)]))
} else if (length(args) == 0){
  stop("No arguments given. The following arguments are accepted:
       
       --input-file <count_file>
              A text file containing a table of features x samples of the run to 
              process by DESeq2. The [...]feature_counts.csv output of aquatx-count is expected here.
      
       --outfile-prefix <outfile>
              Name of the output files to write. These will be created:
                  1. Normalized count table of all samples
                  2. Differential gene expression table per comparison

       --pca
             Use this option to produce principle component analysis plots
             using the DESeq2 library

       ")
} else if (length(args) < 4){
  stop(sprintf("Not enough arguments given. Only parsed %d arguments.", length(args)))
}

# Grab all the arg information
arg_pos <- 1
plot_pca = FALSE
while (arg_pos <= length(args)){
  
  if (args[arg_pos] == '--input-file'){
    count_file <- args[arg_pos + 1]
  } else if (args[arg_pos] == '--outfile-prefix'){
    out_pref <- args[arg_pos + 1]
  } else if (args[arg_pos] == '--pca') {
    plot_pca <- TRUE
  } else if (!(args[arg_pos - 1] %in% c('--input-file', '--outfile-prefix', '--pca'))){
    stop(sprintf("This argument %s%s is not accepted. Did you accidentally include a space in your file or
         sample names?", args[arg_pos - 1], args[arg_pos]))
  }
  arg_pos <- arg_pos + 1
}

#### ---- Set up the parameters ---- ####
counts <- read.csv(count_file,row.names = 1)
# The first column (Feature ID) is absorbed as the index
# The remaining two columns (Feature Name and Feature Class) need to be dropped
counts <- counts[3:length(counts)]
counts <- data.frame(sapply(counts, as.integer), row.names = rownames(counts))
# Create a data frame matching the file, name, condition
samples <- colnames(counts)
condition <- rep('none', length(samples))
for (i in 1:length(samples)){
  sample <- as.character(samples[i])
  condition[i] <- strsplit(sample, '_rep_')[[1]][1]
}

# Create the deseqdataset
sample_table <- data.frame(row.names=samples, condition=condition)
deseq_ds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_table, design = ~ condition)

#### ---- Run DESeq2 & write outputs ---- ####

# Create the DESeq Object
deseq_run <- DESeq(deseq_ds)

# Get the normalized counts
deseq_counts <- counts(deseq_run, normalized=TRUE)
write.csv(deseq_counts, paste(out_pref, "norm_counts.csv", sep="_"))

# Create & retrieve all possible comparisons
all_comparisons <- t(combn(unique(condition), 2))
for (i in 1:nrow(all_comparisons)){
  comparison <- all_comparisons[i,] 
  deseq_res <- results(deseq_run, c("condition", comparison[1], comparison[2]))
  deseq_res <- deseq_res[order(deseq_res$padj),]

  plt <- plotPCA(rlog(deseq_ds))
  trellis.device(device="pdf", file=paste(out_pref, "cond1", comparison[1], "cond2", comparison[2], "pca_plot.pdf", sep="_"))
  print(plt)
  dev.off()
  write.csv(deseq_res, paste(out_pref, "cond1", comparison[1], "cond2", comparison[2], "deseq_table.csv", sep="_"))
}


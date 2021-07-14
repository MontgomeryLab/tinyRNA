#!/usr/bin/env Rscript --vanilla

#### ---- Get the command line arguments ---- ####

args <- commandArgs(trailingOnly = TRUE)
usage <- "The following arguments are accepted:

       --input-file <count_file>
              A text file containing a table of features x samples of the run to
              process by DESeq2. The [...]feature_counts.csv output of aquatx-count is expected here.

       --outfile-prefix <outfile>
              Name of the output files to write. These will be created:
                  1. Normalized count table of all samples
                  2. Differential gene expression table per comparison
                  3. A PCA plot per comparison, if --pca is also provided.

       --pca
              Optional. This will produce principle component analysis plots
              using the DESeq2 library. Output files are PDF format.

       "

## Throw an error if an unexpected number of arguments is provided
if (length(args) > 5){
  stop(gettextf("Too many arguments given. %d arguments were parsed.
  Was there an unquoted space in your arguments?\n\n%s", length(args), usage))
} else if (length(args) == 0){
  stop(gettextf("No arguments given.\n\n%s", usage))
} else if (length(args) < 4){
  stop("Not enough arguments given. Only parsed %d arguments.\n\n%s", length(args), usage)
}

## Assign arg variables
count_file <- args[match('--input-file', args) + 1]
out_pref <- args[match('--outfile-prefix', args) + 1]
plot_pca <- '--pca' %in% args

## Make sure that argument parameters weren't skipped
if (out_pref %in% c('--input-file', '--pca')){
  stop(gettextf("Please be sure to include your outfile prefix.\n\n%s", usage))
} else if (count_file %in% c('--outfile-prefix', '--pca')){
  stop(gettextf("Please be sure to include your count file.\n\n%s", usage))
}

## Now that command line args have been validated, load libraries for performing DEG and writing PCA plots
library(DESeq2)
library(lattice)

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

  if (plot_pca){
    plt <- plotPCA(rlog(deseq_ds))
    trellis.device(device="pdf", file=paste(out_pref, "cond1", comparison[1], "cond2", comparison[2], "pca_plot.pdf", sep="_"))
    print(plt)
    dev.off()
  }

  write.csv(deseq_res, paste(out_pref, "cond1", comparison[1], "cond2", comparison[2], "deseq_table.csv", sep="_"))
}


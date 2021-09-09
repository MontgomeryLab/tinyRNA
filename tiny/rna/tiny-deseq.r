#!/usr/bin/env Rscript

#### ---- Validate the provided command line arguments ---- ####

args <- commandArgs(trailingOnly = TRUE)
usage <- "The following arguments are accepted:

       --input-file <count_file>
              A text file containing a table of features x samples of the run to
              process by DESeq2. The [...]feature_counts.csv output of tinyrna-count is expected here.

       --outfile-prefix <outfile>
              Name of the output files to write. These will be created:
                  1. Normalized count table of all samples
                  2. Differential gene expression table per comparison
                  3. A PCA plot per comparison, if --pca is also provided.

       --pca
              Optional. This will produce principle component analysis plots
              using the DESeq2 library. Output files are PDF format.

       --drop-zero
              Optional. Prior to performing analysis, this will drop all
              rows/features which have a zero count in all libraries."

## Returns the provided data as a dataframe with a corresponding classes column, regardless of order
df_with_classes <- function(data){
  return(
    transform(
      merge(classes, data.frame(data), by=0),
      row.names = Row.names,
      Row.names = NULL
    )
  )
}

## Throw an error if an unexpected number of arguments is provided
if (length(args) > 6){
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
drop_zero <- '--drop-zero' %in% args

## Make sure that required argument parameters weren't skipped
if (out_pref %in% c('--input-file', '--pca', '--drop-zero')){
  stop(gettextf("Please be sure to include your outfile prefix.\n\n%s", usage))
} else if (count_file %in% c('--outfile-prefix', '--pca', '--drop-zero')){
  stop(gettextf("Please be sure to include your count file.\n\n%s", usage))
}

## Now that command line args have been validated, load libraries for performing DEG and writing PCA plots
library(DESeq2)
library(lattice)
library(utils)


#### ---- Read in inputs ---- ####


## Resolve relative paths and perform tilde expansion for input file path
count_file <- normalizePath(count_file)

## DESeq2 prefers R-safe column names. We want to preserve the original col names to use in final outputs
orig_sample_names <- colnames(read.csv(count_file, check.names = FALSE, row.names = 1, nrows = 1))[-(1:2)]

## Read counts CSV with sanitized column names for handling. Subset classes.
counts <- read.csv(count_file, row.names = 1)
classes <- data.frame(counts[2])

## Subset counts table to drop Feature Name and Feature Class columns before integer sapply
counts <- data.frame(sapply(counts[-(1:2)], as.integer), row.names = rownames(counts))

## Remove rows containing all zeros if user requests
if (drop_zero){
  counts <- counts[rowSums(counts[]) > 0, ]
}


#### ---- Set up the parameters ---- ####


## Get samples and corresponding conditions as a named vector
## sampleCondition: "names" are R safe sample names, members are the original sample names
samples <- orig_sample_names
sampleCondition <- rep('none', length(samples))
for (i in seq_along(samples)){
  sample <- as.character(samples[i])
  condition <- strsplit(sample, '_rep_')[[1]][1]
  sampleCondition[i] <- condition
  names(sampleCondition)[i] <- make.names(condition)
}

## Create the deseqdataset
sample_table <- data.frame(row.names=make.names(samples), condition=factor(names(sampleCondition)))
deseq_ds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = sample_table, design = ~ condition)


#### ---- Run DESeq2 & write outputs ---- ####


## Create the DESeq Object
deseq_run <- DESeq2::DESeq(deseq_ds)

## Produce PCA plot if requested
if (plot_pca){
    plt <- plotPCA(DESeq2::rlog(deseq_ds))
    trellis.device(device="pdf", file=paste(out_pref, "pca_plot.pdf", sep="_"))
    print(plt)
}

## Get normalized counts and write them to CSV with original sample names in header
deseq_counts <- df_with_classes(counts(deseq_run, normalized=TRUE))
colnames(deseq_counts) <- c("Feature Class", orig_sample_names)
write.csv(deseq_counts, paste(out_pref, "norm_counts.csv", sep="_"))

## Create & retrieve all possible comparisons
all_comparisons <- t(combn(unique(names(sampleCondition)), 2))
for (i in seq_len(nrow(all_comparisons))){
  comparison <- all_comparisons[i,]

  deseq_res <- DESeq2::results(deseq_run, c("condition", comparison[1], comparison[2]))
  result_df <- df_with_classes(deseq_res[order(deseq_res$padj),])
  colnames(result_df)[1] <- "Feature Class"

  write.csv(
    result_df,
    paste(out_pref,
          # Resolve original condition names for use in output filename
          "cond1", sampleCondition[[comparison[1]]],
          "cond2", sampleCondition[[comparison[2]]],
          "deseq_table.csv", sep="_")
  )
}


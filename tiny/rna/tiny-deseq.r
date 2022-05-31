#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- "
Required arguments:

       --input-file <count_file>
              A text file containing a table of features x samples of the run to
              process by DESeq2. The [...]feature_counts.csv output of tinyrna-count is expected here.

       --outfile-prefix <outfile>
              Name of the output files to write. These will be created:
                  1. Normalized count table of all samples
                  2. Differential gene expression table per comparison
                  3. A PCA plot per comparison, if --pca is also provided.
                  
Optional arguments:

       --control <control_condition>
              If the control condition is specified, comparisons will only 
              be made between the control and experimental conditions.

       --pca
              This will produce principle component analysis plots using 
              the DESeq2 library. Output files are PDF format.

       --drop-zero
              Prior to performing analysis, this will drop all rows/features 
              which have a zero count in all samples."

# Increase max error string length by the length of the usage string
# 8170: https://github.com/wch/r-source/blob/tags/R-4-1-1/src/main/options.c#L626
options(warning.length = min(getOption("warning.length") + nchar(usage), 8170))

# Suppress conversion to scientific notation
options(scipen = 999)


#### ---- Validate commandline arguments ---- ####


## Throw an error if an unexpected number of arguments is provided
if (length(args) > 8){
  stop(gettextf("Too many arguments given. %d arguments were parsed.
  Was there an unquoted space in your arguments?\n\n%s", length(args), usage))
} else if (length(args) < 4){
  stop(gettextf("Not enough arguments given. Only parsed %d arguments.\n\n%s", length(args), usage))
}

## Assign arg variables
count_file <- args[match('--input-file', args) + 1]
out_pref <- args[match('--outfile-prefix', args) + 1]
control_grp <- args[match('--control', args) + 1]
has_control <- !(control_grp %in% NA)
plot_pca <- '--pca' %in% args
drop_zero <- '--drop-zero' %in% args

## Make sure that the user included parameters for arguments that require them
all_args <- append(args[startsWith(args, '--')], NA)
if (count_file %in% all_args[all_args != '--input-file']){
  stop(gettextf("Please be sure to include your count file.\n\n%s", usage))
} else if (out_pref %in% all_args[all_args != '--outfile-prefix']){
  stop(gettextf("Please be sure to include your outfile prefix.\n\n%s", usage))
} else if ('--control' %in% args && control_grp %in% all_args[all_args != '--control']){
  stop(gettextf("Please be sure to include your control condition name.\n\n%s", usage))
}


#### ---- Read in inputs ---- ####


## Resolve relative paths and perform tilde expansion for input file path
count_file <- normalizePath(count_file)

## DESeq2 prefers R-safe column names. We want to preserve original sample names to use in final outputs
orig_sample_names <- colnames(read.csv(count_file, check.names = FALSE, row.names = 1, nrows = 1))[0:-3]

## Read counts CSV with sanitized column names for handling. Subset classes and aliases.
counts <- read.csv(count_file)

## Set rownames to concatenated Feature ID and Tag, then drop these columns
rownames(counts) <- split(counts[, c('Feature.ID', 'Tag')], seq(nrow(counts)))
counts <- subset(counts, select = -c(Feature.ID, Tag))

## Copy feature metadata and drop these columns, leaving only integer columns
metadata <- data.frame("Feature Name" = counts[[1]], "Feature Class" = counts[[2]], row.names = rownames(counts), check.names = FALSE)
counts <- data.frame(sapply(counts[0:-2], as.integer), row.names = rownames(counts))

## Remove rows containing all zeros if user requests
if (drop_zero){
  counts <- counts[rowSums(counts[]) > 0, ]
}


#### ---- Set up the parameters ---- ####


## Get sample conditions as a named vector where names are R safe names
sampleConditions <- sapply(strsplit(as.character(orig_sample_names), '_rep_'), '[[', 1)
names(sampleConditions) <- make.names(sampleConditions)

## Ensure control group is present in samples, if specified
if (has_control && !(control_grp %in% sampleConditions)){
  stop("The control condition was not found among your samples.")
}

## Now that inputs have been validated and read, load libraries
library(DESeq2)
library(ggplot2)
library(utils)

## Returns a new dataframe with metadata columns prepended
df_with_metadata <- function(classless_df){
  return(data.frame(
    merge(metadata, data.frame(classless_df), by=0),
    row.names = "Row.names",
    check.names = FALSE
  ))
}

restore_multiindex <- function(base_df){
  base_df[] <- sapply(base_df, as.character)
  base_matrix <- as.matrix(base_df)

  # Parse row names back to "MultiIndex" comprised of Feature ID and Tag
  split_rn <- lapply(rownames(base_matrix), function(x) eval(parse(text=x)))
  split_mx <- do.call(rbind, split_rn)

  # Assign the MultiIndex columns and drop rownames
  multiidx_mx <- cbind(split_mx, base_matrix)
  colnames(multiidx_mx) <- c("Feature ID", "Tag", colnames(multiidx_mx))
  rownames(multiidx_mx) <- NULL

  return(multiidx_mx)
}

## Create the deseqdataset
sample_table <- data.frame(row.names=names(orig_sample_names), condition=factor(names(sampleConditions)))
deseq_ds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = sample_table, design = ~ condition)


#### ---- Run DESeq2 & write outputs ---- ####


## Create the DESeq Object
deseq_run <- DESeq2::DESeq(deseq_ds)

## Produce PCA plot if requested
if (plot_pca){
  vst_res <- tryCatch(
    DESeq2::vst(deseq_ds),
    error=function(e) {
      return(DESeq2::varianceStabilizingTransformation(deseq_ds))
  })

  plt <- DESeq2::plotPCA(vst_res) + ggplot2::theme(aspect.ratio = 1)
  ggplot2::ggsave(paste(out_pref, "pca_plot.pdf", sep="_"))
}

## Get normalized counts and write them to CSV with original sample names in header
deseq_counts <- df_with_metadata(DESeq2::counts(deseq_run, normalized=TRUE))
colnames(deseq_counts)[0:-2] <- orig_sample_names
write.csv(
  restore_multiindex(deseq_counts),
  paste(out_pref, "norm_counts.csv", sep="_"),
  row.names=FALSE,
  quote=1:4,
  na=""
)

if (has_control){
  # Comparison is the cartesian product of control and experimental conditions
  r_safe_control_name <- make.names(control_grp)
  all_comparisons <- as.matrix(expand.grid(
    r_safe_control_name,
    unique(names(sampleConditions[sampleConditions != r_safe_control_name]))
  ))
} else {
  # Comparison is all combinations of conditions
  all_comparisons <- t(combn(unique(names(sampleConditions)), 2))
}

write_dge_table <- function (dge_df, cond1, cond2){
  write.csv(
    restore_multiindex(dge_df),
    paste(out_pref,
          "cond1", cond1,
          "cond2", cond2,
          "deseq_table.csv", sep="_"),
    row.names=FALSE,
    quote=1:4,
    na=""
  )
}

## Perform condition comparisons
for (i in seq_len(nrow(all_comparisons))){
  comparison <- all_comparisons[i,]

  deseq_res <- DESeq2::results(deseq_run, c("condition", comparison[2], comparison[1]))
  result_df <- df_with_metadata(deseq_res[order(deseq_res$padj),])

  # Resolve original condition names for use in output filename
  cond1 <- sampleConditions[[comparison[1]]]
  cond2 <- sampleConditions[[comparison[2]]]
  write_dge_table(result_df, cond1, cond2)

  if (!has_control){
    # Save DGE table for reverse comparison
    flip_sign <- c("log2FoldChange", "stat")
    result_df[flip_sign] <- result_df[flip_sign] * -1
    write_dge_table(result_df, cond2, cond1)
  }
}

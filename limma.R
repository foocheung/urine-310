#' Perform analysis on = expression data
#' Usage:
#'   analyze_protein_expression.R [--input_file=<input_file>] [--output_file=<output_file>]
#'
#' Options:
#'   -h --help             Show this help message and exit.
#'   --input_file=<input_file>   Path to the input file (annotated_foo_S2_stats.txt by default).
#'   --output_file=<output_file> Path to the output file (results.txt by default).
#'

library(tidyverse)
library(limma)
library(fgsea)
library(docopt)

# Define command line arguments
doc <- "Perform analysis on gene expression data"
args <- docopt(doc)

# Read input file
input_file <- ifelse(is.null(args$'--input_file'), "./annotated_foo_S2_stats.txt", args$'--input_file')
list <- read_tsv(input_file)

# Load data and perform initial processing
a <- read_tsv("../Step1_Matrix_28-Mar-2023_10_53_54.txt")
a <- a[, names(a) %in% c("SampleID", "TimePoint", "Group", "Subject", dput(list$EntrezGeneSymbol))] %>%
  filter(Group %in% c("BB"))

b <- a %>%
  select(1:4) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  mutate("group" = paste0(TimePoint, "_", Group, "_", Subject)) %>%
  arrange(group)

b$TimePoint <- gsub("-", "_", b$TimePoint)

design <- model.matrix(~ 0 + TimePoint, data = b)

# Extract counts and perform correlation
counts <- a[7:length(a)]
rownames(counts) <- paste0(a$TimePoint, "_", a$Group, "_", a$Subject)

count2 <- counts %>%
  rownames_to_column %>%
  gather(var, value, -rowname) %>%
  spread(rowname, value)

rownames(count2) <- count2$var
count3 <- count2[2:length(count2)]
rownames(count3) <- count2$var
count3 <- count3 %>% mutate_if(is.character, as.numeric)
count3 <- log(count3)
rownames(count3) <- count2$var

corfit <- duplicateCorrelation(count3, design, block = b$Subject)

cont.matrix <- makeContrasts(
  delta1 = (TimePointPost_Op) - (TimePointPre_Op),
  levels = design
)

fit <- lmFit(count3, design = design, block = b$Subject, correlation = corfit$consensus)
contrast_fit <- contrasts.fit(fit, contrasts = cont.matrix)
ebays_fit <- eBayes(contrast_fit, robust = TRUE)

res <- as.data.frame(topTable(ebays_fit, coef = 1, number = Inf))

# Load pathways and perform GSEA
tm <- gmtPathways("./ReactomePathways.gmt")
ch <- res %>% arrange(t)
res2 <- cbind(rownames(ch), as.numeric(as.character(ch$t)))
res2 <- as.tibble(res2)
res2$V2 <- as.numeric(res2$V2)
ranks <- deframe(as.data.frame(res2))

fgseaRes <- fgsea(pathways = tm, stats = ranks, maxSize = 500)

fgseaRes_IBD <- fgseaRes %>% arrange(padj)

# Write results to output file
output_file <- ifelse(is.null(args$'--output_file'), "results.txt", args$'--output_file')
write_tsv(fgseaRes_IBD, output_file)

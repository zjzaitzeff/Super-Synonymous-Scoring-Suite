args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
#grabs output file
user_lib <- file.path(Sys.getenv("HOME"), "R", "libs")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

options(repos = c(CRAN = "https://cran.rstudio.com/"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("atSNP")

library(atSNP)
data(encode_library)
length(encode_motif)

# https://www.bioconductor.org/packages/3.0/data/annotation/html/BSgenome.Hsapiens.NCBI.GRCh38.html
install.packages("BiocManager")
BiocManager::install("atSNP")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

# loading SNPs (https://www.bioconductor.org/packages/release/bioc/vignettes/atSNP/inst/doc/atsnp-vignette.html#load-the-snp-data)
data(example)
# remove "chr" from chr column entries since the example was built for BSgenome.Hsapiens.UCSC.GRCh38, which is not available for this version of BiocManager
snp_tbl$chr <- sub("^chr", "", snp_tbl$chr)
snp_file <- file.path(output_dir, "test_snp_file.txt")
write.table(snp_tbl, file = snp_file,
            row.names = FALSE, quote = FALSE)
snp_info <- LoadSNPData(snp_file, genome.lib = "BSgenome.Hsapiens.NCBI.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)
ncol(snp_info$sequence) == nrow(snp_tbl)
snp_info$rsid.rm

write("\n", file = snp_file, append = TRUE)
#affinity scores
atsnp.scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 1)
write.table(atsnp.scores$motif.scores, file = snp_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE)

write("\n", file = snp_file, append = TRUE)
#p-values
atsnp.result <- ComputePValues(motif.lib = motif_library, snp.info = snpInfo, motif.scores = atsnp.scores$motif.scores, ncores = 1, testing.mc=TRUE)
write.table(atsnp.result, file = snp_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE)

#just the p-values
p_values_only <- head(atsnp.result[order(atsnp.result$pval_rank), c("snpid", "motif", "pval_ref", "pval_snp", "pval_rank")])
write("\n", file =snp_file, append = TRUE)
write.table(p_values_only, file = snp_file, sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE)

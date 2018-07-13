#!/usr/bin/env Rscript
library(dplyr)
library(methods)


args <- commandArgs(TRUE)

sample_id <- dir(file.path(args[1]))
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(args[1], id))
kal_dirs

s2c <- read.table(args[2], header = TRUE, stringsAsFactors=FALSE) 
s2c <- s2c[order(s2c$run_accession),]
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

s2c
# Extract gene or transcripts infromation from Ensemble data source
ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = ensembl)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Load the kallisto processed data into the object
so <- sleuth::sleuth_prep(s2c, target_mapping = t2g, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, aggregation_column = 'ens_gene')#

so <- sleuth::sleuth_fit(so, ~condition, 'full')
so <- sleuth::sleuth_fit(so, ~1, 'reduced')
so <- sleuth::sleuth_lrt(so, 'reduced', 'full')
full_results <- sleuth::sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(full_results, qval <= 0.05)
write.table(sleuth_significant, paste("sleuth_significant.csv", sep = "\t"))

# Generate the whole gene table
gene_table <- sleuth::kallisto_table(so)
write.table(gene_table, paste("gene_table.csv"), sep = "\t")
# Saving for shinny
save(so, file=paste("sleuth_object.so"))
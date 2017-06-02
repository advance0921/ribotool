
devtools::use_package("tidyverse")
devtools::use_package("devtools")
devtools::use_package("magrittr")
devtools::use_package("roxygen2")

# Import the codon count table ------------------------------------------------
import_codon_count_data <- function(codoncount_data, levels_sample) {
    # This file contains all the genes in E. coli. We only use genes that are coded for protein and exclude pseudogenes
    Genedata <- EcoliAnnotation %>% filter(Type == "aa") %>% filter(!grepl("pseudogene", Protein))
    # Load codon count
    data <- readr::read_delim(codoncount_data, delim = "\t") %>% semi_join(., Genedata, by = c(gene_name = "Gene")) %>% mutate(sample = factor(sample,
        levels = levels_sample))
    return(data)
}

# Check the total count number from each sample that are mapped to coding regions
# ------------------------------------------------ check_raw_data_seqnumber <- function(raw_data){ output <- raw_data %>%
# group_by(condition) %>% summarise(seqcounts = sum(codon_count_sum)) return(output) }


# This function takes codon count data and normalize it to RPM, RPKM and TPM -----------------------------------------
normalize_seqdata <- function(codoncount_data){

    output <- codoncount_data %>% group_by(sample) %>% mutate(RPM = codon_count_sum/sum(codon_count_sum)*10^6) %>%
        group_by(sample,gene_name) %>%
        mutate(RPKM = RPM/max(codon_index)*10^3, TPM = codon_count_sum/max(codon_index)*10^3) %>%
        group_by(sample) %>%
        mutate(TPM = TPM/sum(TPM)*10^6) %>%
        group_by(sample,gene_name) %>%
        mutate(RPM_normgene = RPM/sum(RPM),
               RPKM_normgene = RPKM/sum(RPKM),
               TPM_normgene = TPM/sum(TPM)) %>%
        select(-codon_count_sum,-contains("position"))
    return(output)
}


# This function sums up all the counts within one gene in each sample -----------------------------------------
#Nest gene and sample first. It will just sum up the data
summarize_seqdata <- function(codoncount_data){
    output_summary <- codoncount_data %>%
        summarise(gene_len = n(),RPM = sum(RPM),RPKM = sum(RPKM), TPM = sum(TPM))
    return(output_summary)
}

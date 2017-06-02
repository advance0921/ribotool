
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


# Check the total count number from each sample that are mapped to coding regions
# ------------------------------------------------ normalize_raw_data <- function(){ print('hello! World!') }


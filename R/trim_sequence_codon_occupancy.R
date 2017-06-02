
devtools::use_package("tidyverse")
devtools::use_package("devtools")
devtools::use_package("magrittr")
devtools::use_package("roxygen2")
devtools::use_package("cowplot")

# This function trim the end of sequence. It is useful because for codon occupancy analysis during elongation, we need to remove effect from initiation and termination  ----------------------------------------------------
# Nest each sample first before running this function
trim_sequence <- function(data,trim_start,trim_end,genelength_thresh,TPM_thresh){
    trimmed_data <- data %>% group_by(gene_name) %>%
        filter(sum(TPM)>=TPM_thresh) %>%
        filter(n()>=genelength_thresh) %>%
        filter(codon_index %in% seq(trim_start + 1 ,max(codon_index)-trim_end))
    return(trimmed_data)
}

# This function takes gene sequence, splits it into two halves and sum up the data  ----------------------------------------------------
summarise_two_halves <- function(input_data){
    # First we label which codon index belongs to first or second half of the gene
    labeled_data <- input_data %>% mutate( gene_position = ifelse(
        (codon_index <= floor((min(codon_index) + max(codon_index))/2)), "first" , "second" ))

    # Then we call up the summarize_seqdata from RiboTool to summarise it
    labeled_data_summary <- labeled_data %>% group_by(gene_position) %>% nest() %>%
        mutate(summary = map(data,~RiboTool::summarize_seqdata(.))) %>% select(-data) %>% unnest(summary)
    return(labeled_data_summary)
}


# This function calculate codon occupancy  ----------------------------------------------------




# This function plot the codon occupancy analysis and highlight the codon with its corresponding amino acid based on the significance cut-off ----------------------------------------------------






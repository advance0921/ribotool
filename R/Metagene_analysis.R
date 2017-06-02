
devtools::use_package("tidyverse")
devtools::use_package("devtools")
devtools::use_package("magrittr")
devtools::use_package("roxygen2")
devtools::use_package("cowplot")

# This function takes the data of count for each position, and position which user wants to calculate metagene (50 codons
# from start and stop codon, for example) and returns the ribosome distribution----------------------------------------------------
cal_metagene <- function(data, norm_position) {

    data <- data %>% group_by(sample) %>% mutate(RPM_renorm = RPM/sum(RPM),
                                                 RPKM_renorm = RPKM/sum(RPKM),
                                                 TPM_renorm = TPM/sum(TPM)) %>%
        group_by(gene_id, gene_name, sample) %>%
        mutate(interval = cut(codon_index, c(0, norm_position + 1, max(codon_index) -
        norm_position + 1, max(codon_index)), labels = c("start", "mid", "stop"))) %>%
        filter(interval %in% c("start", "stop")) %>%
        group_by(gene_id, gene_name, interval) %>%
        mutate(dist_from_tip = case_when(interval == "start" ~ as.double(codon_index)-2,
                                         interval == "stop" ~ as.double(max(codon_index) - codon_index))) %>%
        filter(dist_from_tip >= 0) %>% group_by(sample, interval) %>% nest() %>%
        mutate(ribo_distri = map(data, ~RiboTool::cal_ribosome_distribution(.))) %>% select(-data)
    return(data)
}

# This function can calculate the ribosome distribution based on sum or median---------------------------------------------
cal_ribosome_distribution <- function(data) {
    ribo_distri <- data %>% group_by(dist_from_tip) %>%
        summarise(RPM_renorm_sum = sum(RPM_renorm),
                  RPKM_renorm_sum = sum(RPKM_renorm),
                  TPM_renorm_sum = sum(TPM_renorm),
                  median_normgene = median(RPM_normgene,na.rm=TRUE),
                  mean_normgene = mean(RPM_normgene,na.rm=TRUE))
    return(ribo_distri)
}


# This function takes the metagene distribution and plot them based on whether it's start or stop part of the gene----------------------
plot_ribo_dist <- function(interval_position, data, plot_yvalue, palette, plot_position, plot_yscale,aspect_ratio) {
    switch(as.character(interval_position),
           start = {
               plot_style <- list(theme(legend.position = "top", aspect.ratio = aspect_ratio),
                                  scale_y_continuous(limits = plot_yscale, labels = scales::percent),
                                  scale_x_continuous(limits = c(0,plot_position + 1)),
                                  ylab("Percentage of ribosome total count"), xlab("Distance from start codon (AA)"))},
           stop = {
               plot_style <- list(theme(legend.position = "top", aspect.ratio = aspect_ratio),
                                  scale_x_reverse(limits = c(plot_position + 1,0)), scale_y_continuous(position = "right",
                                                                        limits = plot_yscale, labels = scales::percent),
                                  xlab("Distance from stop codon (AA)"))
           })
    ribo_dist_plot <- data %>% filter(dist_from_tip < plot_position) %>% ggplot(aes_string(x = "dist_from_tip", y = plot_yvalue, col = "sample")) + geom_line(size = 1) +
        scale_color_manual(values = palette) + plot_style
    return(ribo_dist_plot)
}

# This function simply plot the start and end together using the cowplot plot grid------------------
plot_start_end <- function(plot_data) {

    plot_start <- plot_data$ribodist_plot[[which(plot_data$interval == "start")]]
    plot_stop <- plot_data$ribodist_plot[[which(plot_data$interval == "stop")]]
    ribo_plot <- cowplot::plot_grid(plot_start, plot_stop, align = "h")
    return(ribo_plot)
}


# This function takes the summary data of riboseq and rnaseq and return translational efficiency ----------------------------------------------------

calc_TE <- function(cds,rna){
    te <- cds %>% inner_join(.,rna,by= c("gene_name","sample","gene_len")) %>% mutate(TransEff = RPKM.x/RPKM.y,TransEff_TPM = TPM.x/TPM.y) %>% select(-contains(".x"),-contains(".y")) %>% IDPmisc::NaRV.omit()
    return(te)
}

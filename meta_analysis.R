

rm(list=ls())

library(tidyverse)
library(data.table)
library(janitor)
library(knitr)
library(kableExtra)
library(ggpubr)
library(gridExtra)
library(here)
library(meta)
library(grid)
library(tableone)

slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter

################################################################################
# helper functions

# function to format meta::metagen output for data.frame
format_metagen <- function(x) {
  # x input is metagen object
  
  # get meta fixed and random estimates 
  out <- bind_rows(summary(x)$random, summary(x)$fixed) %>% select(1:6)
  colnames(out) <- c("estimate", "std.error", "conf.low", "conf.high", "statistic", "p.value")
  
  # get heterogeneity diagnoistics
  diag <- x[c("Q", "df.Q", "pval.Q", "tau2", "lower.tau2", "upper.tau2", 
              "I2", "lower.I2", "upper.I2", "H", "lower.H", "upper.H")] %>% data.frame %>% 
    mutate_all(as.numeric)
  diag <- data.frame(diag[1:3], diag[4], 
                     tau_ci=paste0("[", round(diag[5], 2), "-", round(diag[6], 2), "]"), 
                     diag[7], 
                     I2_ci=paste0("[", round(diag[8], 2), "-", round(diag[9], 2), "]"), 
                     diag[10], 
                     H_ci=paste0("[", round(diag[11], 2), "-", round(diag[12], 2), "]"))
  diag <- rbind(diag, c(rep(NA, length(diag))))
  
  # return output
  return(data.frame(dataset=c("meta_random", "meta_fixed"), out, diag))
}


################################################################################
# meta anlysis: vitd IV

# read in results
path <- "manuscript/code_to_upload/for_repo/"
vitd_results <- readRDS(file=here(file.path(path, "results/regression_vitd_wgrs_by_dataset.rds")))

# meta anlysisis for estiamtes adjusted for cov --------
results_vitd_adj <- vitd_results$results_vitd_adj

meta_vitd_adj <- results_vitd_adj %>% filter(grepl("v", term)) %>% mutate(effect=name) %>% 
  group_by(effect) %>% 
  group_modify(~metagen(TE=.x$estimate, seTE=.x$std.error, studlab=.x$dataset) %>% 
                 format_metagen) 
meta_vitd_adj


#  meta-analysis (fixed and random) for estiamtes stratified by sex --------

results_vitd_sex_adj <- vitd_results$results_vitd_sex_adj

meta_vitd_sex_adj <- results_vitd_sex_adj %>% filter(grepl("v", term)) %>% mutate(effect=name) %>% 
  group_by(effect, sex) %>% 
  group_modify(~metagen(TE=.x$estimate, seTE=.x$std.error, studlab=.x$dataset) %>% 
                 format_metagen) 

# save
saveRDS(list(meta_vitd_adj=meta_vitd_adj, 
             meta_vitd_sex_adj=meta_vitd_sex_adj, 
             results_vitd_adj=results_vitd_adj, 
             results_vitd_sex_adj=results_vitd_sex_adj),
        file=here(file.path(path, "results/meta_analysis_vitd_wgrs.rds")))



################################################################################
# meta anlysis: VDR IVs

# read in VDR only results
path <- "manuscript/code_to_upload/for_repo/"
results <- readRDS(file=here(file.path(path, "results/regression_vdr_by_dataset.rds")))
results_swed <- readRDS(file=here(file.path(path, "results/regression_vdr_by_dataset_sweden.rds")))

# check overlapping VDR SNPs in both datasets
vdr_snps_in_both <- intersect(results$results_vdr_adj$term, results_swed$results_vdr_adj$term)
length(vdr_snps_in_both)



# meta anlysisis for estiamtes adjusted for cov --------

# combine kpnc/gera, ukb, sweden results
results_vdr <- bind_rows(results$results_vdr_adj, results_swed$results_vdr_adj)


# meta-analysis (fixed and random)
meta_vdr_adj <- results_vdr %>% mutate(effect=term) %>% 
  group_by(effect) %>% 
  group_modify(~metagen(TE=.x$estimate, seTE=.x$std.error, studlab=.x$dataset) %>% 
                 format_metagen %>% slice(1)) %>% ungroup %>% 
  mutate(p_fdr=p.adjust(p.value, method="BH"), 
         p_bonf=p.adjust(p.value, method="bonf")) %>% 
  mutate(ci=paste0("[", round(conf.low,2), ", ", round(conf.high,2), "]")) %>% 
  select(effect, dataset, logOR=estimate, ci, p.value, p_fdr, p_bonf, everything())


# meta anlysisis for estiamtes stratified by sex --------

# combine results
results_vdr_sex_adj <- bind_rows(results$results_vdr_sex_adj, results_swed$results_vdr_sex_adj)

# meta-analysis (fixed and random)
meta_vdr_sex_adj <- results_vdr_sex_adj %>% mutate(effect=term) %>% 
  group_by(effect, sex) %>% 
  group_modify(~metagen(TE=.x$estimate, seTE=.x$std.error, studlab=.x$dataset) %>% 
                 format_metagen %>% slice(1)) %>% group_by(sex) %>% 
  ungroup %>% group_by(sex) %>% 
  mutate(p_fdr=p.adjust(p.value, method="BH"), 
         p_bonf=p.adjust(p.value, method="bonf")) %>% 
  mutate(ci=paste0("[", round(conf.low,2), ", ", round(conf.high,2), "]")) %>% 
  select(effect, sex, dataset, logOR=estimate, ci, p.value, p_fdr, p_bonf, everything())


# save --------
saveRDS(list(meta_vdr_adj=meta_vdr_adj, 
             meta_vdr_sex_adj=meta_vdr_sex_adj, 
             results_vdr=results_vdr, 
             results_vdr_sex_adj=results_vdr_sex_adj), 
        file=here(file.path(path, "results/meta_analysis_vdr.rds")))

################################################################################
# meta anlysis: VDR x vitd

# read in regression results
path <- "manuscript/code_to_upload/for_repo/"
results <- readRDS(file=here(file.path(path, "results/regression_vdr_vitd_by_dataset.rds")))
results_swed <- readRDS(file=here(file.path(path, "results/regression_vdr_vitd_by_dataset_sweden.rds")))

# meta anlysisis for estiamtes adjusted for cov --------

# combine kpnc/gera, ukb, sweden results
results_vdr_vitd <- bind_rows(results$results_vdr_vitd_adj, 
                              results_swed$results_vdr_vitd_adj) %>% 
  filter(term=="snp:vitd_value") %>% 
  filter(snp %in% vdr_snps_in_both) %>% 
  mutate(vitd=sub("_test_dose", "", vitd))

# meta analysis across all combinations of VDR and vitD IVs
# random and fixed
meta_vdr_vitd_adj <- results_vdr_vitd %>% 
  mutate(effect=paste0(vitd, "_", snp)) %>% 
  group_by(snp, vitd) %>% 
  group_modify(~metagen(TE=.x$estimate, seTE=.x$std.error, studlab=.x$dataset) %>% 
                 format_metagen %>% slice(1)) %>% ungroup %>% 
  mutate(flip=ifelse(grepl("flip", vitd), "flip", "no_flip")) %>% 
  group_by(flip) %>% 
  mutate(p_fdr=p.adjust(p.value, method="fdr"), 
         p_bonf=p.adjust(p.value, method="bonf")) %>% 
  mutate(ci=paste0("[", round(conf.low,2), ", ", round(conf.high,2), "]")) %>% 
  mutate(effect=paste0(vitd, "_", snp)) %>% 
  select(effect, snp, vitd, dataset, logOR=estimate, ci, p.value, p_fdr, p_bonf, everything())


# meta anlysisis for estiamtes stratified by sex --------

# combine kpnc/gera, ukb, sweden results
results_vdr_vitd_sex <- bind_rows(results$results_vdr_vitd_sex_adj, 
                                  results_swed$results_vdr_vitd_sex_adj) %>% 
  filter(snp %in% vdr_snps_in_both) %>% 
  filter(term=="snp:vitd_value") %>% 
  filter(snp %in% vdr_snps_in_both) %>% 
  mutate(vitd=sub("_test_dose", "", vitd)) %>% 
  mutate(sex=ifelse(grepl("GSA|Omni", dataset), sex-1, sex))

# meta analysis across all combinations of VDR and vitD IVs
# fixed and random

meta_vdr_vitd_sex_adj <- results_vdr_vitd_sex %>% 
  mutate(effect=paste0(vitd, "_", snp)) %>% 
  group_by(snp, vitd, sex) %>% 
  group_modify(~metagen(TE=.x$estimate, seTE=.x$std.error, studlab=.x$dataset) %>% 
                 format_metagen %>% slice(1)) %>% group_by(sex) %>% 
  ungroup %>% group_by(sex, vitd) %>% 
  mutate(effect=paste0(vitd, "_", snp)) %>% 
  mutate(p_fdr=p.adjust(p.value, method="BH"), 
         p_bonf=p.adjust(p.value, method="bonf")) %>% 
  mutate(ci=paste0("[", round(conf.low,2), ", ", round(conf.high,2), "]")) %>% 
  select(effect, sex, snp, vitd, dataset, logOR=estimate, ci, p.value, p_fdr, p_bonf, everything())

# save --------
saveRDS(list(meta_vdr_vitd_adj=meta_vdr_vitd_adj, 
             meta_vdr_vitd_sex_adj=meta_vdr_sex_adj, 
             results_vdr_vitd=results_vdr_vitd,
             results_vdr_vitd_sex=results_vdr_vitd_sex), 
        file=here(file.path(path, "results/meta_analysis_vdr_vitd.rds")))






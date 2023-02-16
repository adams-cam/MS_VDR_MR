

rm(list=ls())

library(tidyverse)
library(data.table)
library(janitor)
library(here)
library(broom)
library(glm2)
library(tableone)
library(kableExtra)
library(flextable)
library(officer)
library(conflicted)

select <- dplyr::select
filter <- dplyr::filter

# read in covarites and vdr and 25OHD IV data
path <- "manuscript/code_to_upload/for_repo/"
cov <- readRDS(file=here(file.path(path, "data/pheno.rds")))
vdr <- readRDS(file=here(file.path(path, "data/vdr_geno_iv.rds")))
vitd <- readRDS(file=here(file.path(path, "data/vitd_IVs.rds")))

# what dataset to do analysis on
# 
cur_dataset <- c("kpnc", "gera")


################################################################################
# make cov

# filter to dataset of interest
cov <- cov %>% filter(dataset==cur_dataset)

# make age category variable
cov <- cov %>% group_by(dataset) %>% 
  mutate(age_cat=cut(age, breaks=quantile(age, probs=c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm=T), 
                     right=F, labels=F))
with(cov, table(is.na(age), dataset))
with(cov, table(age_cat, dataset))


################################################################################
# vitd


# continuous vitd --------

# join data
dat <- vitd %>% inner_join(cov) %>% distinct() 

# get vars to iterate over
cols <- dat %>% select(matches("vD")) %>% colnames %>% grep("flip", ., invert=T, value=T)

# crude 
results_vitd_crude <- 
  dat %>% pivot_longer(cols=cols) %>% group_by(dataset, name) %>% 
  mutate(value=scale(value, scale=F)) %>% 
  group_modify(~glm(case ~ value, data=.x, family="binomial") %>% 
                 tidy(exponentiate=F, conf.int=T) %>% filter(grepl("v", term)))

# adjusted
vars <- c("case", "sex", "age_cat", "drb1_1501_bin", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
dat %>% select(vars) %>% is.na %>% colSums
results_vitd_adj <- 
  dat %>% pivot_longer(cols=cols) %>% mutate(dataset2=dataset) %>% 
  group_by(dataset, name) %>% 
  mutate(value=scale(value, scale=F)) %>% 
  group_modify(~{
    f <- "case ~ value + factor(age_cat) + sex + drb1_1501_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"
    glm(f, data=.x, family="binomial") %>% 
                 tidy(exponentiate=F, conf.int=T) %>% filter(grepl("v", term))
})


# adjusted stratified by sex
results_vitd_sex_adj <- 
  dat %>% pivot_longer(cols=cols) %>% mutate(dataset2=dataset) %>% 
  group_by(dataset, name, sex) %>% 
  mutate(value=scale(value, scale=F)) %>% 
  group_modify(~{
    f <- "case ~ value + factor(age_cat) + drb1_1501_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"
    glm(f, data=.x, family="binomial") %>% 
      tidy(exponentiate=F, conf.int=T) %>% filter(grepl("v", term))
  })

saveRDS(list(results_vitd_crude=results_vitd_crude, results_vitd_adj=results_vitd_adj, 
             results_vitd_sex_adj=results_vitd_sex_adj, dat=dat), 
        file=paste0("results/vitd_wgrs_regression_", cur_dataset, ".rds"))



################################################################################
# vdr

# make function for getting tidied output of glm only for paramters of interest
# much faster when needing profile confidence intervals from adjusted models
# this second version is only for the MS ~ VDR models

tidy2 <- function(x, pat) {
  #x <- m; pat<- "sex"
  ci <- suppressMessages(confint(x, parm=grep(pat, names(coef(x))), trace=F) %>% t() %>% data.frame %>% 
                           rename(conf.low=1, conf.high=2))
  tidy(x) %>% filter(grepl(pat, term)) %>% bind_cols(ci)
}


# merge data
dat <- vdr %>% inner_join(cov, by="FID") %>% distinct() %>% 
  mutate(dataset2=dataset)

# mean impute missing values within dataset (if genotype data)
vdr %>% is.na %>% colSums
dat <- dat %>% group_by(dataset) %>% 
  mutate(across(matches("rs"), ~replace_na(.x, mean(.x, na.rm=T)))) %>% ungroup

# get SNPs to iterate over
cols <- dat %>% select(contains("rs")) %>% colnames; length(cols)

# crude
f <- "case ~ snp"
results_vdr_crude <- lapply(cols, function(vdr_snp) {
  out <- dat %>% mutate(snp=unlist(dat[, vdr_snp])) %>% 
    group_by(dataset) %>% 
    group_modify(~glm(case ~ snp, data=.x, family="binomial") %>% 
                   tidy2(pat="snp") %>% filter(grepl("snp", term))) %>% 
    #mutate(conf.low=estimate+qnorm(0.025)*std.error, conf.high=qnorm(0.975)*std.error+estimate) %>% 
    mutate(term=vdr_snp)
  return(out)
})
results_vdr_crude <- results_vdr_crude %>% do.call(rbind, .)
results_vdr_crude %>% head

# adjusted
f <- "case ~ snp + sex + factor(age_cat) + drb1_1501_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"
results_vdr_adj <- lapply(cols, function(vdr_snp) {
  #vdr_snp=cols[1]
  out <- dat %>% mutate(snp=unlist(dat[, vdr_snp])) %>% 
    group_by(dataset) %>% 
    group_modify(~{
       glm(f, data=.x, family="binomial") %>% 
        tidy2(pat="snp") %>% filter(grepl("snp", term))}) %>% 
    mutate(term=vdr_snp)
    return(out)
})
results_vdr_adj <- results_vdr_adj %>% do.call(rbind, .)
results_vdr_adj %>% head

# stratfied by sex
f <- "case ~ snp + factor(age_cat) + drb1_1501_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"
results_vdr_sex_adj <- lapply(cols, function(vdr_snp) {
  #vdr_snp=cols[1]
  out <- dat %>% mutate(snp=unlist(dat[, vdr_snp])) %>% 
    group_by(dataset, sex) %>% 
    group_modify(~{glm(f, data=.x, family="binomial") %>% 
                   tidy2(pat="snp") %>% filter(grepl("snp", term))}) %>% 
    mutate(term=vdr_snp)
  return(out)
})
results_vdr_sex_adj <- results_vdr_sex_adj %>% do.call(rbind, .)


# save
saveRDS(list(results_vdr_adj=results_vdr_adj, results_vdr_crude=results_vdr_crude, 
             results_vdr_sex_adj=results_vdr_sex_adj, dat=dat), 
        file=paste0("results/vdr_regression", cur_dataset, ".rds"))

################################################################################
# vitd x vdr

# join data
dat <- vdr %>% inner_join(cov, by="FID") %>% distinct() %>% 
  inner_join(vitd)


# get values to iterate over
vitd_cols <- dat %>% select(matches("vD")) %>% colnames
vdr_cols <- dat %>% select(contains("rs")) %>% colnames

# make dat long against vitd vars
dat_long <- dat %>% pivot_longer(cols=all_of(vitd_cols)) %>% rename(vitd=name, vitd_value=value) 

# crude
length(vdr_cols); length(vitd_cols)
f <- "case ~ snp*vitd_value"
results_vdr_vitd_crude <- lapply(vdr_cols, function(vdr_snp) {
  #vdr_snp=cols[1];
  out <- dat_long %>% mutate(snp=unlist(dat_long[, vdr_snp])) %>% 
    group_by(dataset, vitd) %>% 
    group_modify(~{glm(f, data=.x, family="binomial") %>% tidy2(pat="snp|value") %>% 
    mutate(snp=vdr_snp)}) %>% select(dataset, vitd, snp, everything())
  return(out)
})
results_vdr_vitd_crude <- results_vdr_vitd_crude %>% do.call(rbind, .)
results_vdr_vitd_crude %>% head

# adj
results_vdr_vitd_adj <- lapply(vdr_cols, function(vdr_snp) {
  #vdr_snp=cols[1]
  out <- dat_long  %>% mutate(dataset2=dataset) %>% mutate(snp=unlist(dat_long[, vdr_snp])) %>% 
    group_by(dataset, vitd) %>% 
    group_modify(~{
      f <- "case ~ snp*vitd_value + factor(age_cat) + sex + drb1_1501_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"
      glm(f, data=.x, family="binomial") %>% 
        tidy2(pat="snp|value") %>% mutate(snp=vdr_snp)}) %>% 
    select(dataset, vitd, snp, everything())
  return(out)
})
results_vdr_vitd_adj <- results_vdr_vitd_adj %>% do.call(rbind, .)

# stratified by sex
results_vdr_vitd_sex_adj <- lapply(vdr_cols, function(vdr_snp) {
  out <- dat_long  %>% mutate(dataset2=dataset) %>% mutate(snp=unlist(dat_long[, vdr_snp])) %>% 
    group_by(dataset, vitd, sex) %>% 
    group_modify(~{
      f <- "case ~ snp*vitd_value + factor(age_cat) + drb1_1501_bin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"
      glm(f, data=.x, family="binomial") %>% 
        tidy2(pat="snp|value") %>% mutate(snp=vdr_snp)}) %>% 
    select(dataset, vitd, snp, everything())
  return(out)
})
results_vdr_vitd_sex_adj <- results_vdr_vitd_sex_adj %>% do.call(rbind, .)


# save
saveRDS(list(results_vdr_vitd_crude=results_vdr_vitd_crude, results_vdr_vitd_adj=results_vdr_vitd_adj, 
             results_vdr_vitd_sex_adj=results_vdr_vitd_sex_adj, dat=dat), 
        file=paste0("results/vdr_vitd_regression", cur_dataset, ".rds"))








###############################################################################
#                                                                             #
#                             Run MR analyses                                 #
#                                                                             #
###############################################################################

###############################################################################
#                                 Set up                                      #
###############################################################################

# Clear the work environment
rm(list = ls())

# Setting digits
options(digits = 10)

# Required libraries
## install.packages("remotes")
## remotes::install_github('MRCIEU/TwoSampleMR')
x <- c("dplyr", "purrr", "data.table", "TwoSampleMR", "stringr", "tidyr")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AUTOIMMUNE_DIR"), "working/")
data_dir <- file.path(home_dir, "data/")
gwas_dir <- file.path(data_dir, "gwas/")
exp_data_dir <- file.path(data_dir, "exposure_dat/")
outcome_data_dir <- file.path(data_dir, "/outcome/")
out_dir <- paste0(home_dir, "/results/")
setwd(home_dir)

# Load functions
source("scripts/functions/generate_outdat_with_proxies.R")

# define MHC region on chr6
# according to GRC: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
# defined as chr6:28,477,797-33,448,354 on build 37
mhc_start <- 28477797
mhc_end <- 33448354

###############################################################################
#                        Load exposures & outcomes                            #
###############################################################################

# Exposures
opengwas_list <- read.csv(paste0(gwas_dir, "/opengwas_list.csv"))
opengwas_list$opengwas_id[opengwas_list$phenotype_id == "ms"] <- "ms"
all_exposures <- c(opengwas_list$phenotype_id)

# Outcomes
primary_outcomes <- read.table(paste0(outcome_data_dir, "primary_outcomes.txt"))$V1
secondary_outcomes <- read.table(paste0(outcome_data_dir, "secondary_outcomes.txt"))$V1 
all_outcomes <- c(primary_outcomes, secondary_outcomes)
# exclude outcomes not available in current data release
all_outcomes <- all_outcomes[!(all_outcomes %in% c("conmal", "conmal_chd"))]
# continuous
continuous_outcomes <- c("zbw_all", "ga_all")
# full outcome names
outcome_clean <- data.frame(outcome = c(all_outcomes),
    outcome_clean = c("Miscarriage", "Stillbirth", "Hypertensive disorders of pregnancy",
    "Gestational diabetes mellitus", "Preterm birth", "Small for gestational age",
    "Large for gestational age", "Low Apgar score at 5 minutes", "NICU admission",
    "Preterm birth (spontaneous)",  "Gestational hypertension", "Preeclampsia",
    "Birthweight", "Low birthweight", "High birthweight", "Gestational age"))

###############################################################################
#                            Load Outcome datasets                            #
###############################################################################

# read in full file for maternal-fetal duos GWAS
full_duos <- fread(paste0(home_dir, "/data/duos_out_dat.txt")) %>%
    mutate(chr.outcome = chr, pos.outcome = pos,
            effect_allele.outcome = effect_allele,
            other_allele.outcome = other_allele, eaf.outcome = NA,
            ncase.outcome = NA, ncontrol.outcome = NA,
            outcome = Phenotype, mr_keep.outcome = TRUE,
            pval_origin.outcome = "reported", id.outcome = Phenotype,
            data_source.outcome = "textfile")  %>%
    # filter to outcomes of interest only
    filter(outcome %in% all_outcomes)

# read in full file for study level GWAS
full_stu <- fread(paste0(home_dir, "/data/stu_out_dat.txt")) %>%
        mutate(chr.outcome = chr, pos.outcome = pos,
            effect_allele.outcome = effect_allele,
            other_allele.outcome = other_allele, eaf.outcome = eaf,
            beta.outcome = beta, se.outcome = se, pval.outcome = pval,
            samplesize.outcome = samplesize, ncase.outcome = ncase,
            ncontrol.outcome = ncontrol, outcome = Phenotype,
            mr_keep.outcome = TRUE, pval_origin.outcome = "reported",
            id.outcome = Phenotype, data_source.outcome = "textfile") %>%
    # filter to outcomes of interest only
    filter(outcome %in% all_outcomes)

# list of included studies
studies <- unique(full_stu$study)

###############################################################################
#                              Run MR Analyses                                #
###############################################################################

# select methods
mr_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")

# dataframes to write results to
mhc_record <- data.frame()
mr_res <- data.frame()
study_loo_res <- data.frame()
q_stats <- data.frame()
egger_intercept <- data.frame()

# appending sample sizes to results as they are generated
# function for leave-one-out sample sizes, calculating sample without study_name
calculate_loo_sample_size <- function(study_name){
    stu_res_ivw %>%
    filter(study != study_name & outcome == outcome_name) %>%
    group_by(outcome) %>%
    summarise(
        sum_min_n = sum(min_n, na.rm = TRUE),
        sum_mean_n = sum(mean_n, na.rm = TRUE),
        sum_median_n = sum(median_n, na.rm = TRUE),
        sum_max_n = sum(max_n, na.rm = TRUE),
        sum_median_ncases = sum(median_ncases, na.rm = TRUE),
        sum_median_ncontrol = sum(median_ncontrol, na.rm = TRUE)) %>%
    mutate_if(is.numeric, round) %>%
    mutate(analysis = paste(study_name, "left out")) %>%
    select(-outcome)
}

# for each exposure
for(exposure_name in all_exposures){

    ### 0 - read in or filter datasets
    
    # exposure data
    exp_dat <- read.csv(paste0(exp_data_dir, exposure_name,
        "_pval_5e_08_clumped.csv"))
    # read in proxies
    proxies <- read.csv(paste0(exp_data_dir, exposure_name,
        "_proxies.csv"))
    all_snps <- c(exp_dat$SNP, proxies$rsid)

    # outcome datasets, retrieving all exposure SNPs & proxies
    # outcome data - meta-analysis
    ma_outdat <- read_outcome_data(
        filename = paste0(home_dir, "/data/autoimmune/outcome/ma_out_dat.txt"),
        snps = all_snps)
    # outcome data - single studies
    stu_outdat <- full_stu %>%
        # filter to only exposure SNPs
        filter(SNP %in% all_snps) %>%
        select(c(colnames(ma_outdat), "study")) %>%
        as.data.frame()
    # outcome data - maternal GWAS (adjusted for fetal genotype)
    duos_outdat <- full_duos %>%
        mutate(beta.outcome = beta_mat_donuts, se.outcome = se_mat_donuts,
            pval.outcome = p_mat_donuts, samplesize.outcome = n_mat) %>%
        filter(SNP %in% all_snps) %>%
        select(c(colnames(ma_outdat))) %>%
        as.data.frame()
    # outcome data - fetal GWAS (unadjusted)
    fetal_outdat <- full_duos %>%
        mutate(beta.outcome = beta_off, se.outcome = se_off,
            pval.outcome = p_off, samplesize.outcome = n_off) %>%
        filter(SNP %in% all_snps) %>%
        select(c(colnames(ma_outdat))) %>%
        as.data.frame()
    # outcome data - fetal GWAS (adjusted for maternal genotype)
    fetal_adj_outdat <- full_duos %>%
        mutate(beta.outcome = beta_off_donuts, se.outcome = se_off_donuts,
            pval.outcome = p_off_donuts, samplesize.outcome = n_off) %>%
        filter(SNP %in% all_snps) %>%
        select(c(colnames(ma_outdat))) %>%
        as.data.frame()

    ### run analysis on each outcome, i.e. each exposure-outcome pair
    for(outcome_name in all_outcomes){

        print(paste(exposure_name, outcome_name))

        ### 1 - run primary mr analysis using meta-analysis results
        ma_outdat_tmp <- generate_outdat_with_proxies(exp_dat, ma_outdat,
            outcome_name, proxies)
        # harmonise
        ma_dat <- harmonise_data(exp_dat, ma_outdat_tmp)
        # mr
        print("Running MR")
        mr_res_tmp <- mr(ma_dat, method_list = mr_methods) %>% 
         mutate(exposure = exposure_name, id.outcome = outcome_name,
            analysis = "Main", study = "ma")
        # sample size
        ma_sample_sizes <- ma_outdat_tmp %>%
            group_by(outcome) %>%
            summarise(min_n = min(samplesize.outcome),
            mean_n = mean(samplesize.outcome),
            median_n = median(samplesize.outcome),
            max_n = max(samplesize.outcome),
            median_ncases = median(ncase.outcome),
            median_ncontrol = median(ncontrol.outcome)) %>%
            mutate_if(is.numeric, round) %>% select(-outcome)
        mr_res_tmp <- cbind(mr_res_tmp, ma_sample_sizes)
        mr_res <- rbind(mr_res, mr_res_tmp)
        # cochrane's Q statistic
        q_stats_tmp <- mr_heterogeneity(ma_dat) %>% 
            mutate(exposure = exposure_name, outcome = outcome_name,
            analysis = "Main", study = "ma")
        q_stats <- rbind(q_stats, q_stats_tmp)
        # MR-egger intercept test
        egger_intercept_tmp <- mr_pleiotropy_test(ma_dat) %>% 
            mutate(exposure = exposure_name, outcome = outcome_name, study = "ma")
        egger_intercept <- rbind(egger_intercept, egger_intercept_tmp)

        ### 2 - leave one SNP out
        mr_res_loo_snp <- mr_leaveoneout(ma_dat, method = mr_ivw) %>%
         mutate(exposure = exposure_name, id.outcome = outcome_name,
            analysis = paste0("Leave one SNP out: ", SNP),
            method = "Inverse variance weighted",
            # re-format
            pval = p, nsnp = NA, study = "loo SNP") %>% 
         select(-c("SNP", "samplesize", "p"))
        mr_res_loo_snp <- cbind(mr_res_loo_snp, ma_sample_sizes)
        mr_res <- rbind(mr_res, mr_res_loo_snp)

        ### 3 - run mr excluding instruments in HLA / MHC region
        mhc_snps <- ma_dat$SNP[(ma_dat$chr.exposure == 6 &
            ma_dat$pos.exposure >= mhc_start & ma_dat$pos.exposure <= mhc_end)]

        if(length(mhc_snps) == 0){
            print("No exposure instruments in MHC region")
            
            # record that no SNPs were excluded
            mhc_record_tmp <- data.frame(exposure = exposure_name, outcome = outcome_name,
                MHC_instruments_excluded = NA)
            mhc_record <- rbind(mhc_record, mhc_record_tmp)

        } else {
            print(paste0(length(mhc_snps), " exposure SNPs in MHC region - running sensitivity MR exluding these SNPs."))
            # record which SNPs were excluded
            mhc_record_tmp <- data.frame(exposure = exposure_name, outcome = outcome_name,
                MHC_instruments_excluded = mhc_snps)
            mhc_record <- rbind(mhc_record, mhc_record_tmp)

            # run mr again without HLA
            ma_dat_no_mhc <- ma_dat[!(ma_dat$SNP %in% mhc_snps),]
            mr_res_tmp <- mr(ma_dat_no_mhc, method_list = mr_methods) %>% 
             mutate(exposure = exposure_name, id.outcome = outcome_name,
                 analysis = "Excluding MHC SNPs", study = "mhc")
            mr_res_tmp <- cbind(mr_res_tmp, ma_sample_sizes)
            mr_res <- rbind(mr_res, mr_res_tmp)
        }


        ### 4 a - run mr for single study outcome data
        stu_res <- data.frame()
        for(study_name in studies){
            # filter outcome data to study
            stu_outdat_all_tmp <- stu_outdat %>% 
                filter(study == study_name) %>% select(-study)
            # filter to outcome of interest, add proxies where available
            stu_outdat_tmp <- generate_outdat_with_proxies(exp_dat,
                stu_outdat_all_tmp, outcome_name, proxies)
            if(dim(stu_outdat_tmp)[1] == 0){
                print(paste("Outcome", outcome_name, "not available in", study_name))
            } else if(dim(stu_outdat_tmp)[1] > 0) {
            # harmonise
            stu_outdat_tmp$id.outcome <- stu_outdat_tmp$id.outcome[1]
            stu_dat <- harmonise_data(exp_dat, stu_outdat_tmp)
            # run MR analysis
            stu_res_tmp <- mr(stu_dat, method_list = c("mr_ivw")) %>% 
                mutate(analysis = paste(study_name, "only"), exposure = exposure_name,
                    outcome = outcome_name, study = study_name)
            # calulcate sample size
            stu_sample_sizes <- stu_outdat_tmp %>%
                group_by(outcome) %>%
                summarise(min_n = min(samplesize.outcome),
                mean_n = mean(samplesize.outcome),
                median_n = median(samplesize.outcome),
                max_n = max(samplesize.outcome),
                median_ncases = median(ncase.outcome),
                median_ncontrol = median(ncontrol.outcome)) %>% 
                mutate_if(is.numeric, round) %>% select(-outcome)
            stu_res_tmp <- cbind(stu_res_tmp, stu_sample_sizes)
            stu_res <- rbind(stu_res, stu_res_tmp)
            }
        }

        mr_res <- rbind(mr_res, stu_res)

        ### 4 b - meta-analysis leaving one study out
        stu_res_ivw <- filter(stu_res, method == "Inverse variance weighted")
        # format data for LOO
        study_loo_dat <- metafor::rma(yi = b, sei = se, data = stu_res_ivw,
          slab = study, method = "FE")
        # run LOO removing one study at a time
        print("Running leave one study out (LOO) MR analysis")
        study_loo_res_tmp <- metafor::leave1out(study_loo_dat)
        study_loo_res_tmp_formatted <- study_loo_res_tmp %>% 
            as.data.frame %>%
            tibble::rownames_to_column(., var = "study") %>%
            mutate(se = estimate/zval,
            # to remove numbers after study:
            study = str_remove(study, "[.][0-9]*"),
            analysis = paste(study, "left out"),
            exposure = exposure_name,
            outcome = outcome_name)
        # sample sizes
        loo_sample_sizes <- data.frame()
        for(study_name in unique(stu_res_ivw$study)){ 
            loo_tmp <- calculate_loo_sample_size(study_name)
            loo_sample_sizes <- rbind(loo_sample_sizes, loo_tmp)
            }
        study_loo_res_tmp_formatted <-
            left_join(study_loo_res_tmp_formatted, loo_sample_sizes,
            by = join_by(analysis))
        study_loo_res <- rbind(study_loo_res, study_loo_res_tmp_formatted)

        ### 5 - run maternal MR adjusted for fetal genotype
        duos_outdat_tmp <- generate_outdat_with_proxies(exp_dat, duos_outdat,
            outcome_name, proxies)
        # harmonise
        duos_dat_tmp <- harmonise_data(exp_dat, duos_outdat_tmp)
        # run mr
        print("Running maternal (adjusted)")
        duos_res_tmp <- mr(duos_dat_tmp, method_list = mr_methods) %>% 
            mutate(analysis = "Maternal (adjusted)",
            exposure = exposure_name, study = "duos")
        if(dim(duos_res_tmp)[1] > 0){
        # sample sizes
        duos_sample_sizes <- duos_outdat_tmp %>%
            group_by(outcome) %>%
            summarise(min_n = min(samplesize.outcome),
            mean_n = mean(samplesize.outcome),
            median_n = median(samplesize.outcome),
            max_n = max(samplesize.outcome),
            # case & control data not available for this dataset,
            # but save NAs to allow rbind()
            median_ncases = median(ncase.outcome),
            median_ncontrol = median(ncontrol.outcome)) %>% 
            mutate_if(is.numeric, round) %>% select(-outcome)
        duos_res_tmp <- cbind(duos_res_tmp, duos_sample_sizes)
        mr_res <- rbind(mr_res, duos_res_tmp)
        # cochrane's Q statistic
        q_stats_tmp_duo <- mr_heterogeneity(duos_dat_tmp) %>% 
            mutate(exposure = exposure_name, outcome = outcome_name,
            analysis = "Maternal (adjusted)", study = "duos")
        q_stats <- rbind(q_stats, q_stats_tmp_duo)
        }

        ### 6 - run fetal MR
        fetal_outdat_tmp <- generate_outdat_with_proxies(exp_dat, fetal_outdat,
            outcome_name, proxies)
        # harmonise
        fetal_dat_tmp <- harmonise_data(exp_dat, fetal_outdat_tmp)
        # run mr
        print("Running fetal (unadjusted)")
        fetal_res_tmp <- mr(fetal_dat_tmp, method_list = mr_methods) %>% 
            mutate(analysis = "Fetal (unadjusted)",
            exposure = exposure_name, study = "duos")
        if(dim(fetal_res_tmp)[1] > 0){
            # sample sizes
        fetal_sample_sizes <- fetal_outdat_tmp %>%
            group_by(outcome) %>%
            summarise(min_n = min(samplesize.outcome),
            mean_n = mean(samplesize.outcome),
            median_n = median(samplesize.outcome),
            max_n = max(samplesize.outcome),
            # case & control data not available for this dataset,
            # but save NAs to allow rbind()
            median_ncases = median(ncase.outcome),
            median_ncontrol = median(ncontrol.outcome)) %>% 
            mutate_if(is.numeric, round) %>% select(-outcome)
        fetal_res_tmp <- cbind(fetal_res_tmp, fetal_sample_sizes)
        mr_res <- rbind(mr_res, fetal_res_tmp)
        }

        ### 7 - run fetal MR adjusted for maternal genotype
        fetal_adj_outdat_tmp <- generate_outdat_with_proxies(exp_dat, fetal_adj_outdat,
            outcome_name, proxies)
        # harmonise
        fetal_adj_dat_tmp <- harmonise_data(exp_dat, fetal_adj_outdat_tmp)
        # run mr
        print("Running fetal (adjusted)")
        fetal_adj_res_tmp <- mr(fetal_adj_dat_tmp, method_list = mr_methods) %>% 
            mutate(analysis = "Fetal (adjusted)",
            exposure = exposure_name, study = "duos")
        if(dim(fetal_adj_res_tmp)[1] > 0){
        # sample sizes
        fetal_adj_sample_sizes <- fetal_adj_outdat_tmp %>%
            group_by(outcome) %>%
            summarise(min_n = min(samplesize.outcome),
            mean_n = mean(samplesize.outcome),
            median_n = median(samplesize.outcome),
            max_n = max(samplesize.outcome),
            # case & control data not available for this dataset,
            # but save NAs to allow rbind()
            median_ncases = median(ncase.outcome),
            median_ncontrol = median(ncontrol.outcome)) %>% 
            mutate_if(is.numeric, round) %>% select(-outcome)
        fetal_adj_res_tmp <- cbind(fetal_adj_res_tmp, fetal_adj_sample_sizes)
        mr_res <- rbind(mr_res, fetal_adj_res_tmp)
        }

    }

}

# inspect to check analyses ran successfully
head(mr_res)
table(mr_res$outcome, mr_res$exposure)
dim(study_loo_res)

###############################################################################
#                         Clean and scale MR results                          #
###############################################################################

# Exposure is currently per 1-unit log odds increase in autoimmune disease
# liability
# To make results more easily interpretable, we will rescale this to per
# doubling in autoimmune disease liability
# as per: Burgess, S. & Labrecque, J. A. Mendelian randomization with a binary
# exposure variable: interpretation and presentation of causal estimates.
# Eur J Epidemiol 33, 947–952 (2018).

### Main MR analysis ##########################################################
### Append full outcome name, full exposure name
mr_res <- mr_res %>%
    left_join(outcome_clean) %>%
    left_join(opengwas_list, by = c("exposure" = "phenotype_id"))

### Scale and add ORs
### binary outcomes
mr_res_binary <- mr_res %>%
  # restrict to binary
  filter(!(outcome %in% continuous_outcomes)) %>%
  # calculate CIs
  generate_odds_ratios() %>%
  # scale to per doubling in liability
  mutate(
    b = b * log(2),
    lo_ci = lo_ci * log(2),
    up_ci = up_ci * log(2),
    # add OR CIs
    or = exp(b),
    or_lci95 = exp(lo_ci),
    or_uci95 = exp(up_ci),
  ) %>%
  # add clean labels for plots & tables
  mutate(
    `Cases / controls` = paste0(median_ncases, " / ", median_ncontrol),
    `OR (95% CI)` = sprintf("%.2f (%.2f - %.2f)", or, or_lci95, or_uci95))

### continuous outcomes - re-scaled both exposure and outcome
mr_res_cont <- mr_res %>%
  # restrict to continuous
  filter(outcome %in% continuous_outcomes) %>%
  mutate(
    # 1. standardise continuous outcomes to 1 SD unit
    # zbw birthweight - z-scores so SD = 1 by definition
    # ga gestational age - 1 SD = 1.89 weeks
    b = case_when(
        outcome == "zbw_all" ~ b/1,
        outcome == "ga_all" ~ b/1.89,
        .default = b),
	  se = case_when(
        outcome == "zbw_all" ~ se/1,
        outcome == "ga_all" ~ se/1.89,
        .default = se),
    # add CIs
    lo_ci = b - qnorm(0.975)*se,
    up_ci = b + qnorm(0.975)*se,
    # 2. scale to per doubling in liability of exposure
    b = b * log(2),
    lo_ci = lo_ci * log(2),
    up_ci = up_ci * log(2),
    # add clean labels for plots & tables
    `Sample size` = paste0(max_n),
     `Beta (95% CI)` = sprintf("%.3f (%.3f - %.3f)", b, lo_ci, up_ci)
     )

### Leave one study out sensitivity analysis ##################################
### Append full outcome name, full exposure name
study_loo_res <- study_loo_res %>%
    left_join(outcome_clean) %>%
    left_join(opengwas_list, by = c("exposure" = "phenotype_id"))

### Scale and add ORs
### binary outcomes
study_loo_res_binary <- study_loo_res %>%
  # restrict to binary
  filter(!(outcome %in% continuous_outcomes)) %>%
  # add CIs
  mutate(
    b = estimate,
    lo_ci = (b - qnorm(0.975) * se),
    up_ci = (b + qnorm(0.975) * se)) %>% 
  # scale to per doubling in liability
    mutate(
      b = b * log(2),
      lo_ci = lo_ci * log(2),
      up_ci = up_ci * log(2),
  # add OR CIs
    or = exp(b),
    or_lci95 = exp(lo_ci),
    or_uci95 = exp(up_ci),
  ) %>%
  # add clean labels for plots & tables
  mutate(
    `Cases / controls` = paste0(sum_median_ncases, " / ", sum_median_ncontrol),
    `OR (95% CI)` = sprintf("%.2f (%.2f - %.2f)", or, or_lci95, or_uci95))

### continuous outcomes - re-scaled both exposure and outcome
study_loo_res_cont <- study_loo_res %>%
  # restrict to continuous
  filter(outcome %in% continuous_outcomes) %>%
  mutate(
    # 1. standardise continuous outcomes to 1 SD unit
    # zbw birthweight - z-scores so SD = 1 by definition
    # ga gestational age - 1 SD = 1.89 weeks
    b = estimate,
    b = case_when(
        outcome == "zbw_all" ~ b/1,
        outcome == "ga_all" ~ b/1.89,
        .default = b),
	  se = case_when(
        outcome == "zbw_all" ~ se/1,
        outcome == "ga_all" ~ se/1.89,
        .default = se),
    # add CIs
    lo_ci = b - qnorm(0.975) * se,
    up_ci = b + qnorm(0.975) * se,
    # 2. scale to per doubling in liability of exposure
    b = b * log(2),
    lo_ci = lo_ci * log(2),
    up_ci = up_ci * log(2),
    # add clean labels for plots & tables
    `Sample size` = paste0(sum_max_n),
     `Beta (95% CI)` = sprintf("%.3f (%.3f - %.3f)", b, lo_ci, up_ci)
     )

###############################################################################
#                            Write out MR results                             #
###############################################################################

### Append: full outcome name, full exposure name, OR 95% CI / Beta 95% CI

### Main results
write.csv(mr_res, paste0(out_dir, "/mr_results.csv"), row.names = FALSE)
write.csv(mr_res_binary, paste0(out_dir, "/mr_results_binary.csv"), row.names = FALSE)
write.csv(mr_res_cont, paste0(out_dir, "/mr_results_cont.csv"), row.names = FALSE)

### Leave one study out sensitivity
write.csv(study_loo_res, paste0(out_dir, "/study_loo_results.csv"), row.names = FALSE)
write.csv(study_loo_res_binary, paste0(out_dir, "/study_loo_results_binary.csv"), row.names = FALSE)
write.csv(study_loo_res_cont, paste0(out_dir, "/study_loo_results_cont.csv"), row.names = FALSE)

### MHC SNPs for each exposure record
mhc_record <- mhc_record %>%
    left_join(outcome_clean) %>%
    left_join(opengwas_list, by = c("exposure" = "phenotype_id"))
write.csv(mhc_record, paste0(out_dir, "/mhc_record.csv"), row.names = FALSE)

## Cochran's Q statistic
q_stats <- q_stats %>%
    left_join(outcome_clean) %>%
    left_join(opengwas_list, by = c("exposure" = "phenotype_id"))
write.csv(q_stats, paste0(out_dir, "/q_stats.csv"), row.names = FALSE)

## Egger intercept
egger_intercept <- egger_intercept %>%
    left_join(outcome_clean) %>%
    left_join(opengwas_list, by = c("exposure" = "phenotype_id"))
write.csv(egger_intercept, paste0(out_dir, "/egger_intercept.csv"), row.names = FALSE)

q('no')
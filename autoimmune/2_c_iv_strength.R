###############################################################################
#                                                                             #
#                            Estimate IV strength                             #
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
data_dir <- file.path(home_dir, "data/autoimmune/")
exp_data_dir <- file.path(data_dir, "/exposure_dat/")
out_dir <- paste0(home_dir, "/results/autoimmune/")
setwd(home_dir)

###############################################################################
#                      Load SNP-exposure data & GWAS N                        #
###############################################################################

# Exposures
opengwas_list <- read.csv(paste0(data_dir, "gwas/opengwas_list.csv"))
opengwas_list$opengwas_id[opengwas_list$phenotype_id == "ms"] <- "ms"
all_exposures <- c(opengwas_list$phenotype_id)
#[1] "axsp"      "cd"        "ht"        "ibd"       "ms"        "ps"       
# [7] "ra-eu"     "ra-eu-eas" "sle"       "ss"        "t1d"   

###############################################################################
#                          Estimate IV strength                               #
###############################################################################

# Using clumped instruments - proxies not included since these will vary by outcome

exp_info <- data.frame()

for(exposure_name in all_exposures){

    # Read in instruments
    exp_dat <- read.csv(paste0(exp_data_dir, exposure_name, "_pval_5e_08_clumped.csv")) %>%
        # Append n cases
        left_join(opengwas_list, by = c("id.exposure" = "opengwas_id"))

    # calculate F-statistic for instrument strength
    exp_dat <- mutate(exp_dat,
                    f  = ((beta.exposure/se.exposure)^2),
                    ) %>%
                # AS rs130075 has se = NA, removing:
                filter(mr_keep.exposure == TRUE)

    tmp_exp_info <- group_by(exp_dat, id.exposure) %>%
        summarise(.,
                id = unique(id.exposure),
                Nsnps = n(), 
                mean_F = round(mean(f, na.rm = TRUE), 0),
                max_F = round(max(f, na.rm = TRUE), 0),
                min_F = round(min(f, na.rm = TRUE), 0),
                median_F = round(median(f, na.rm = TRUE), 0),
                eff_N = round(mean(eff_N, na.rm = TRUE), 0)
                )

    exp_info <- rbind(exp_info, tmp_exp_info)
}

exp_info <- exp_info %>% left_join(opengwas_list, by = c("id.exposure" = "opengwas_id")) %>%
    select(phenotype, phenotype_id, Nsnps, mean_F, max_F, min_F, median_F)

write.csv(exp_info, paste0(out_dir, "exp_instrument_info.csv"),
    row.names = FALSE)

q('no')
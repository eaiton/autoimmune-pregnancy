###############################################################################
#                            Select exposure SNPs                             #
###############################################################################

###############################################################################
#                                  Set up                                     #
###############################################################################

# Clear the work environment
rm(list=ls())

# Required libraries
## install.packages("remotes")
## remotes::install_github("MRCIEU/TwoSampleMR")
## remotes::install_github("richardslab/MRutils")
## remotes::install_github("MRCIEU/ieugwasr")
## remotes::install_github("explodecomputer/plinkbinr")
x <- c("dplyr", "tibble", "TwoSampleMR", "MRutils", "ieugwasr", "plinkbinr")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("DRUGTARGET_DIR"), "working/")
gwas_dir <- file.path(home_dir, "data/autoimmune/gwas/")
out_dir <- file.path(home_dir, "/data/autoimmune/exposure_dat/")
proxies_dir <- file.path(home_dir, "/data/autoimmune/proxies/")
setwd(home_dir)

# Read in exposure GWAS information
opengwas_list <- read.csv(paste0(gwas_dir, "/opengwas_list.csv"))

###############################################################################
#                             Select exposure                                 #
###############################################################################

#array_index <- as.numeric((Sys.getenv("SLURM_ARRAY_TASK_ID")))
#print(array_index)

# MRlink doesn't allow concurrent requests, so looping instead:

for(array_index in seq(1:11)){

    phenotype_id <- opengwas_list$phenotype_id[array_index]
    print(phenotype_id)

    opengwas_id <- opengwas_list$opengwas_id[array_index]
    print(opengwas_id)

###############################################################################
#               Find proxies for all exposure instruments                     #
###############################################################################

    instr <- read.csv(paste0(out_dir, "/", phenotype_id, "_pval_5e_08_clumped.csv"))

    proxies <- MRutils::get_proxies(
        rsids = instr$SNP,
        token = Sys.getenv("LDLINK_TOKEN"),
            # request a personal LDlink token here:
            # https://ldlink.nih.gov/?tab=apiaccess
        population = "CEU",
        results_dir = proxies_dir,
        skip_api = FALSE,
        r2_threshold = 0.8)
    head(proxies)
    dim(proxies)

    proxies_path <- paste0(out_dir, phenotype_id, "_proxies.csv")
    write.csv(proxies, proxies_path, row.names = F)
    print(paste("Potential proxy SNPs saved"))

}

q("no")

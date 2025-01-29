###############################################################################
#                             Create rsid_all.txt                             #
###############################################################################
# Last run 8th Jan 2025

###############################################################################
#                                  Set up                                     #
###############################################################################

# Clear the work environment
rm(list=ls())

# Required libraries
x <- c("dplyr", "tibble")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("DRUGTARGET_DIR"), "working/")
gwas_dir <- file.path(home_dir, "data/autoimmune/gwas/")
out_dir <- file.path(home_dir, "/data/autoimmune/exposure_dat/")
setwd(home_dir)

###############################################################################
#                     Read in all instruments and proxies                     #
###############################################################################

all_snps <- c()
all_proxies <- c()

# Read in exposure GWAS information
opengwas_list <- read.csv(file.path(gwas_dir, "opengwas_list.csv"))

for (i in seq_along(opengwas_list$phenotype_id)) {
    phenotype_id <- opengwas_list$phenotype_id[i]
    opengwas_id <- opengwas_list$opengwas_id[i]

    print(paste0("Phenotype: ", phenotype_id, ", ID: ", opengwas_id))

    # Read in clumped genome-wide significant SNPs
    clumped_path <- file.path(out_dir, paste0(phenotype_id, "_pval_5e_08_clumped.csv"))
    tmp_clumped <- read.csv(clumped_path)
    head(tmp_clumped)

    # Read in proxies, r2 > 0.8
    proxies_path <- file.path(out_dir, paste0(phenotype_id, "_proxies.csv"))
    proxies <- read.csv(proxies_path)
    head(proxies)

    all_proxies <- rbind(all_proxies, proxies)

    all_snps <- c(all_snps, tmp_clumped$SNP, proxies$rsid)
}

###############################################################################
#             Save all exposure SNPs & proxies to be extracted                #
###############################################################################

# Save a list of all exposure SNPs and potential proxies which should be extracted
# as rsid_all.txt; the SNP-outcome associations will be extracted next

all_snps_to_extract <- unique(all_snps)

write(all_snps_to_extract, file.path(home_dir, "data/autoimmune/rsid_all.txt"))

write.csv(all_proxies, file.path(out_dir, "proxies.csv"))


###############################################################################
#    Inspecting instruments - which instruments overlap between conditions?   #
###############################################################################

all_snps_info <- data.frame()

for (i in seq_along(opengwas_list$phenotype_id)) {

    phenotype_id <- opengwas_list$phenotype_id[i]
    opengwas_id <- opengwas_list$opengwas_id[i]

    print(paste0("Phenotype: ", phenotype_id, ", ID: ", opengwas_id))

    # Read in clumped genome-wide significant SNPs
    clumped_path <- file.path(out_dir, paste0(phenotype_id, "_pval_5e_08_clumped.csv"))
    tmp_clumped <- read.csv(clumped_path)
    head(tmp_clumped)

    tmp_clumped$id.exposure <- phenotype_id

    tmp_clumped <- tmp_clumped %>% select(chr.exposure, pos.exposure, SNP,
        other_allele.exposure, effect_allele.exposure, beta.exposure,
        se.exposure, eaf.exposure, id.exposure)

    all_snps_info <- rbind(all_snps_info, tmp_clumped)

}

# Precisely overlapping SNPs
# Count of SNPs per exposure
all_snps_info %>% count(id.exposure) %>% arrange(n)

# Counts of conditions for each snp
tab <- table(all_snps_info$SNP, all_snps_info$id.exposure)
head(tab)

# Shared instruments
dup <- all_snps_info %>%
    # remove second RA GWAS for clarity
    filter(id.exposure != "ra-eu-eas") %>%
    group_by(SNP) %>% filter(n()>1) %>% ungroup() %>% data.frame()
dim(dup)
# 34 duplictaed snps, of 605 unique instruments
dup

# Which chromosomes are shared instruments on?
dup %>% count(chr.exposure) %>% arrange(chr.exposure)
### most on chr6

# Which conditions share instruments?
dup %>% count(id.exposure) %>% arrange(n)
# bigger gwas overlap more
# no overlap for ss! this is not the smallest gwas either

tab_dup <- table(dup$SNP, dup$id.exposure)
tab_dup

# instruments shared in
# 3 conditions: rs1611236 (ps, ra, t1d), rs3184504 (ht, ibd, t1d)
# 4 conditions: rs6679677 (ht, ra, sle, t1d)
# ra & t1d seem to have most overlap

# Could consider overlapping regions/LD blocks for more detailed analyses

q('no')
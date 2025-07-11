###############################################################################
#                                                                             #
#                       Retrieve clumped instruments                          #
#                  using MR-Base API and GWAS catalogue                       #
#                                                                             #
###############################################################################

###############################################################################
#                                  Set up                                     #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
## install.packages("remotes")
## remotes::install_github('MRCIEU/TwoSampleMR')
## remotes::install_github("MRCIEU/MRInstruments")
library(TwoSampleMR)
library(MRInstruments)
library(dplyr)
library(tidyr)

# Set directories
home_dir <- paste0(Sys.getenv("AUTOIMMUNE_DIR"), "working/")
gwas_dir <- file.path(home_dir, "data/gwas/")
out_dir <- gwas_dir

###############################################################################
#                              List exposures                                 #
###############################################################################

## Make csv
phenotype <- c("Ankylosing spondylitis", "Coeliac disease",
    "Hashimoto's thyroiditis", "Inflammatory bowel disease",
    "Multiple sclerosis", "Psoriasis", "Rheumatoid arthritis (Eu)",
    "Rheumatoid arthritis (Eu+EAs)", "Systemic lupus erythematosus",
    "Systemic sclerosis", "Type 1 diabetes")
phenotype_id <- c("axsp", "cd", "ht", "ibd", "ms", "ps",
    "ra-eu", "ra-eu-eas", "sle", "ss", "t1d")
opengwas_id <- c("ebi-a-GCST005529", "ieu-a-1058", "ebi-a-GCST90018855",
    "ieu-a-294", "ieu-b-18", "ebi-a-GCST90019017", "ebi-a-GCST002318",
    "ieu-a-833", "ebi-a-GCST003156", "ss", "ebi-a-GCST90014023")
opengwas_list <- data.frame(cbind(phenotype, phenotype_id, opengwas_id))
write.csv(opengwas_list, paste0(out_dir, "/opengwas_list.csv"), row.names = FALSE)

###############################################################################
#                Retrieve clumped instruments from Open GWAS API              #
###############################################################################

# Loading locally so that these instruments can be accessed without relying on
# IEU OpenGWAS API servers, as described here: https://gwas.mrcieu.ac.uk/about/

### Use API to extract instruments
opengwas_id_available <- opengwas_id[opengwas_id != "ss" & opengwas_id != "ms"]
# multiple sclerosis gwas doesn't load,
# systemic sclerosis gwas not on opengwas and no summary statistics,
# so these are clumped separately below

### Loop to extract and save locally
for(i in seq(1:10)){
    id <- opengwas_id_available[i]
    filtered <- opengwas_list %>% filter(opengwas_id == id)
    phenotype_id <- filtered$phenotype_id
    
    instr <- extract_instruments(outcomes = id)
    
    if(exists("instr")){
        write.csv(instr, paste0(Sys.getenv("DRUGTARGET_DIR"),
        "/working/data/autoimmune/exposure_dat/", phenotype_id,
        "_pval_5e_08_clumped.csv"), row.names = FALSE)

        write.csv(instr, paste0(Sys.getenv("DRUGTARGET_DIR"),
        "/working/data/autoimmune/exposure_dat/", phenotype_id,
        "_pval_5e_08_raw.csv"), row.names = FALSE)

        rm(instr)
    }
}

###############################################################################
#              Download OpenGWAS VCF for MS locally and clump                 #
###############################################################################

# Load vcf locally

id <- "ieu-b-18"
file_url = paste0("https://gwas.mrcieu.ac.uk/files/", id, "/", id, ".vcf.gz")
system(paste0("wget -nc ", file_url, " -O ", out_dir, id, ".vcf.gz"))

### Clump

# Convert vcfs to TwoSampleMR format dataframes
phenotype_id <- "ms"
opengwas_id <- "ieu-b-18"

file_path <- paste0(out_dir, phenotype_id)

### Read in and format VCF - nb. computationally intensive so don't re-run
print("Formatting VCF")
tmp <- format_vcf(id = opengwas_id, vcf_file_path = gwas_dir,
        exposure_name = phenotype_id)
head(tmp)

### Filter to genome-wide significant SNPs only
tmp_genomewidesig <- filter(tmp, pval.exposure < 5e-08)
rm(tmp) # to reduce memory burden
if(nrow(tmp_genomewidesig)>0){
    genomewidesig_path <- paste0(Sys.getenv("DRUGTARGET_DIR"),
        "/working/data/autoimmune/exposure_dat/", phenotype_id,
        "_pval_5e_08_raw.csv")
    write.csv(tmp_genomewidesig, genomewidesig_path, row.names = F)
    print(paste("Genome-wide significant SNPs saved:", genomewidesig_path))
    head(tmp_genomewidesig)
}
    
### Clump
if(nrow(tmp_genomewidesig)>0){
print("Clumping SNPs")
# Rename for ld_clump function
    tmp_clumped <- rename(tmp_genomewidesig, pval = pval.exposure, rsid = SNP) # rename for ld_clump
    
    # To run ld_clump through the ieugwas API use the below line however the server can be busy
    # so this command won't always run
    # tmp_clumped <-  ld_clump(tmp_clumped, clump_kb=10000, clump_r2=0.001, pop = "EUR")
    # To run ld_clump locally, use the below line.
    # See: https://mrcieu.github.io/ieugwasr/articles/local_ld.html
    # To download reference panel:
    # wget -nc http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz data/reference_panels/
    tmp_clumped <-  ld_clump(tmp_clumped, clump_kb = 10000, clump_r2 = 0.001,
                      plink_bin = get_plink_exe(), # needs plinkbinr package loaded
                      bfile = european_reference_panel)
    # Rename back to TwoSampleMR format
    tmp_clumped$id <- NULL
    tmp_clumped <- rename(tmp_clumped, pval.exposure = pval, SNP = rsid)
    
    clumped_path <- paste0(Sys.getenv("DRUGTARGET_DIR"),
        "/working/data/autoimmune/exposure_dat/", phenotype_id,
        "_pval_5e_08_clumped.csv")
    write.csv(tmp_clumped, clumped_path, row.names = F)
    print(paste("Clumped SNPs saved:", clumped_path))
    head(tmp_clumped)
}

# Missing clumped instruments for Systemic sclerosis

###############################################################################
#                 Download other GWAS from GWAS catalog locally               #
###############################################################################

# Lopez-Isac 2019 - Systemic sclerosis - full summary stats:
ss_url <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009131/Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt"
system(paste0("wget -nc ", ss_url, " -O ", out_dir, "other/lopez-isac-2019.txt"))
# check download
system(paste0("head ", out_dir, "other/lopez-isac-2019.txt"))
system(paste0("awk '{print $1}' ", out_dir, "other/lopez-isac-2019.txt | uniq -c"))
# missing allele information (column 'A2') so cannot be used

# download GWAS hits manually since some of these have effect allele information
# https://www.ebi.ac.uk/gwas/studies/GCST009131
# gwas-association-downloaded_2024-11-15-accessionId_GCST009131.tsv

ss <- read.table(paste0(out_dir, "gwas-association-downloaded_2024-11-15-accessionId_GCST009131.tsv"),
    sep="\t", header = TRUE)
### clean
ss_clean <- tibble(ss) %>% 
    separate(STRONGEST.SNP.RISK.ALLELE, c("SNP", "effect_allele"), remove = F) %>%
    mutate(
        beta = log(OR.or.BETA),
        ncase = 9095, ncontrol = 17584, samplesize = ncase + ncontrol,
    # estimate standard error from beta and p.value
        se = abs(beta / qnorm(P.VALUE / 2))) %>%
    format_data(
        type = "exposure",
        snps = NULL,
        header = TRUE,
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "RISK.ALLELE.FREQUENCY",
        effect_allele_col = "effect_allele",
        pval_col = "P.VALUE",
        ncase_col = "ncase",
        ncontrol_col = "ncontrol",
        samplesize_col = "samplesize",
        gene_col = "REPORTED.GENE.S.",
        min_pval = 1e-200,
        chr_col = "CHR_ID",
        pos_col = "CHR_POS",
        log_pval = FALSE)
# SNPs without a known effect allele:
ss_clean[ss_clean$mr_keep.exposure == FALSE, ] # 23, so lose 11 SNPs
ss_clean$exposure <- "Systemic sclerosis"
ss_clean$id.exposure <- "ss"

# SNPs given for b38, we use b37
# only analysis this would affect is MHC region exclusion sensitivity
# so changing base pair position from b38 to b37 for SNP on chr6:
# using dbSNP 19th Nov 24: https://www.ncbi.nlm.nih.gov/snp/?term=rs633724
ss_clean$pos.exposure[ss_clean$SNP == "rs633724"] <- "106734040"

### Write out instruments
write.csv(ss_clean, paste0(Sys.getenv("DRUGTARGET_DIR"),
        "/working/data/autoimmune/exposure_dat/ss_pval_5e_08_clumped.csv"), row.names = FALSE)

write.csv(ss_clean, paste0(Sys.getenv("DRUGTARGET_DIR"),
        "/working/data/autoimmune/exposure_dat/ss_pval_5e_08_raw.csv"), row.names = FALSE)

q("no")
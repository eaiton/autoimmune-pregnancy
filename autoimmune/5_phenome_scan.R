############################################################################
#                                                                          #
#                 Phenome-wide scan for influential SNPs                   #
#                    to identify potential pleiotropy               	     #
#                                                                          #
############################################################################
# Last run 25th April 2025

# Influential SNPs identified in leave-one-SNP-out analyses to perform
# phenome-wide association scan (PheWAS) for:
# * rs2857700 was driving effects of multiple sclerosis (MS) on HDP and PTB
# * rs9275183 was driving effects of rheumatoid arthritis (RA) on PTB
# * rs1794269 was driving effects of type 1 diabetes (T1D) on PTB

###############################################################################
#                                 Set up                                      #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
# devtools::install_github('kevinblighe/EnhancedVolcano')
x <- c("TwoSampleMR", "dplyr", "data.table", "purrr", "ggplot2", "openxlsx",
	"ieugwasr", "EnhancedVolcano")
lapply(x, require, character.only = TRUE)

# Set date for plots and results
date <- "20250425"

# Set directories
home_dir <- paste0(Sys.getenv("AUTOIMMUNE_DIR"), "working/")
out_dir <- file.path(home_dir, "/results/pheWAS/")
exp_data_dir <- file.path(home_dir, "data/exposure_dat/")
outcome_data_dir <- file.path(home_dir, "data/outcome/")
setwd(home_dir)

# Input data
# rsids
rsids <-
	c("rs2857700", # MS
	  "rs9275183", # RA
	  "rs1794269") # T1D

# exposure data
exp_dat_ms <- read.csv(paste0(exp_data_dir, "ms_pval_5e_08_clumped.csv")) %>%
  filter(SNP == "rs2857700")
exp_dat_ra <- read.csv(paste0(exp_data_dir, "ra-eu-eas_pval_5e_08_clumped.csv")) %>%
  filter(SNP == "rs9275183") %>%
  mutate(id.exposure = "ra")
exp_dat_t1d <- read.csv(paste0(exp_data_dir, "t1d_pval_5e_08_clumped.csv")) %>%
  filter(SNP == "rs1794269") %>%
  mutate(id.exposure = "t1d")
cols <- intersect(intersect(colnames(exp_dat_ms), colnames(exp_dat_t1d)), colnames(exp_dat_ra))
loci <- rbind(exp_dat_ms[, cols], exp_dat_ra[, cols], exp_dat_t1d[, cols])
# flip one MS SNP so all effects are per liability increasing allele
loci$other_allele.exposure[1] <- "C"
loci$effect_allele.exposure[1] <- "T"
loci$beta.exposure[1] <- loci$beta.exposure[1] * (-1)
loci

# Available outcomes
# Authenticated token for TwoSampleMR required to download catalogue of GWAS
# effects, see:
# https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication
# Lines below commented out since only these need to be run once
#ao <- available_outcomes(opengwas_jwt=Sys.getenv("OPENGWAS_TOKEN"))
#write.table(ao, "data/autoimmune/phenome_scan/avail_outcomes_", date, ".txt", quote = T, col.names = T, row.names = F)
ao <- read.table("data/autoimmune/phenome_scan/avail_outcomes_20250425.txt", header = T)
nrow(ao) # 50043

# Run PheWAS for each selected SNP
phewas <- map(rsids, ~ieugwasr::phewas(., pval = 1e-2, batch = c(),
    opengwas_jwt = token))
# p-value less than this does not run, to minimise server burden

# Format PheWAS output
out_dat <- map(phewas, ~format_data(.,
									type = "outcome", 
									phenotype_col = "trait",
									snp_col = "rsid",
									effect_allele_col = "ea",
									other_allele_col = "nea",
									pval_col = "p",
									samplesize_col = "n",
									id_col = "id"
					 ))

# Harmonise PheWAS results to reflect exposure increasing allele 
dat <- map(out_dat, ~harmonise_data(loci, ., action = 2)) %>% 
		bind_rows %>%
		mutate(z = qnorm(pval.outcome/2, lower.tail = F) * sign(beta.outcome))
head(dat)

# Inspect to check all results are per autoimmune condition liability-increasing allele
loci
dat %>% group_by(SNP) %>% slice_min(pval.outcome, with_ties = FALSE, n = 1) %>% data.frame()
summary(dat$beta.exposure) # should all be > 0
summary(dat$beta.outcome) # any direction

# Make volcano plots (https://github.com/kevinblighe/EnhancedVolcano)
make_volcplot <- function(x) {

			snp <- x
			
			#g <- filter(loci_tophit, SNP == x) %>%
			#		pull(gene.exposure) %>%
			#		as.character 
			
			#target <- filter(drugs, Ensembl.ID == g) %>% pull(Drug.target)
			
			df <- filter(dat, SNP == x)
			
			p_bonf <- 0.05 / nrow(df)
			
			p <- EnhancedVolcano::EnhancedVolcano(df,
						lab = df$outcome,
						x = 'z',
						y = 'pval.outcome',
						xlim = c(-20, 20),
						ylim = c(0, 100),
						xlab = "Z",
						title = paste0(snp),
						subtitle = NULL,
						pCutoff = p_bonf,
						FCcutoff = 1.96,
						labSize	= 3,
						max.overlaps = Inf
						)
				
			ggsave(paste0(out_dir, "/volcano_", snp, "_PheWAS_", date, ".png"),
				plot = p, width = 10, height = 7, dpi = 150, units = "in", device='png')
}

walk(rsids, ~make_volcplot(.))

# Save results
openxlsx::write.xlsx(dat, file = paste0(out_dir, "/pheWAS_IEUGWAS_", date, ".xlsx"), overwrite = T)

q("no")
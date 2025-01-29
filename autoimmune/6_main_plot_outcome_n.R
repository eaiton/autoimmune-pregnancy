###############################################################################
#                                                                             #
#         Main MR results for autoimmune disease - pregnancy outcomes         #
#             Plot forest plots and estimate outcome sample sizes             #
#                                                                             #
###############################################################################

###############################################################################
#                                 Set up                                      #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
## install.packages("remotes")
## remotes::install_github("NightingaleHealth/ggforestplot")
## remotes::install_github("karthik/wesanderson")
## remotes::install_github('MRCIEU/TwoSampleMR')
x <- c("dplyr", "data.table", "patchwork", "purrr", "tibble", "ggplot2",
    "forestploter", "ggforestplot", "wesanderson", "TwoSampleMR", "viridis",
    "plyr", "patchwork")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("DRUGTARGET_DIR"), "working/")
results_dir <- file.path(home_dir, "/results/autoimmune/")
gwas_dir <- file.path(home_dir, "data/autoimmune/gwas/")
outcome_data_dir <- file.path(home_dir, "data/autoimmune/outcome/")
out_dir <- file.path(results_dir, "/forest_plots/main/")
setwd(home_dir)

# outcomes
primary_outcomes <- read.table(paste0(outcome_data_dir, "primary_outcomes.txt"))$V1
secondary_outcomes <- read.table(paste0(outcome_data_dir, "secondary_outcomes.txt"))$V1 
continuous_outcomes <- c("zbw_all", "ga_all")

###############################################################################
#                               Read in results                               #
###############################################################################

### All main analyses and sensitivities results
mr_res_binary <- read.csv(file.path(results_dir, "/mr_results_binary.csv"))
mr_res_binary$`OR (95% CI)` <- mr_res_binary$`OR..95..CI.`
mr_res_binary$`Cases / controls` <- mr_res_binary$`Cases...controls`

mr_res_cont <- read.csv(file.path(results_dir, "/mr_results_cont.csv"))
mr_res_cont$`Beta (95% CI)` <- mr_res_cont $`Beta..95..CI.`
mr_res_cont $`Sample size` <- mr_res_cont $`Sample.size`

exp_info <- read.csv(file.path(results_dir, "/exp_instrument_info.csv"))

###############################################################################
#                                 Clean                                       #
###############################################################################

mr_res_all <- rbind.fill(mr_res_cont, mr_res_binary) %>%
	mutate(
        id.exposure = case_when(
            id.exposure == "ebi-a-GCST90018907" ~ "psv",
            .default = id.exposure),
        outcome_clean =
        case_when(
            outcome == "pe_subsamp" ~ "Preeclampsia",
            outcome == "pretb_all" ~ "Preterm birth",
            .default = outcome_clean))

### Main IVW results to plot
mr_res_all_main <- mr_res_all %>%
    filter(analysis == "Main" & method == "Inverse variance weighted")

### Append median nsnp across outcomes
median_snp <- aggregate(mr_res_all_main$nsnp, by=list(mr_res_all_main$id.exposure), FUN=median)
median_snp$id.exposure <- median_snp$Group.1
median_snp$median_nsnp <- median_snp$x

### Remove RA in Eu onl
### use European +East Asian GWAS (72% Eu ancestry) since
### higher powered and less than 5% of loci were specific to one ancestry (5/101)
mr_res_all_main <- mr_res_all_main %>% filter(exposure != "ra-eu")
mr_res_all_main$phenotype[mr_res_all_main$exposure == "ra-eu-eas"] <- "Rheumatoid arthritis"

table(mr_res_all_main$outcome_clean, mr_res_all_main$phenotype)
# 15 outcomes, 10 conditions

###############################################################################
#                          Main forest plots                                  #
###############################################################################

# outcome as title, exposures as y axis

formatter <- function(...){
  function(x) format(round(x, 2), ...)
}

formatter_continuous <- function(...){
  function(x) format(signif(x, 3), ...)
}

run_fplot <- function(df, prim, outcome_type, nrows, ncols) {
    # df = dataframe to use
    # prim = filter on primary (1) / secondary (2) / or all (3) outcomes
    # outcome_type= "binary" or "continuous" plot
    # ncols = ncols in plot - probably if axsp included 4, if not 3
	# nrows = number of rows in plot
    
    # primary, secondary or all outcomes
    if(prim == 1){
        dt <- filter(df, outcome %in% primary_outcomes)
    } else if(prim == 2) {
        dt <- filter(df, outcome %in% secondary_outcomes)
    } else if(prim == 3){
        dt <- df # no filtering needed
    } else {
		stop("provide 1, 2, or 3 as argument")
	}
    
    dt <- dt %>%
        mutate(outcome_clean = factor(outcome_clean,
            levels = unique(outcome_clean)))

    # create binary or continuous forest plot
    if(outcome_type == "binary"){
    plot <- dt %>%
        filter(!(outcome %in% continuous_outcomes)) %>%
        ggforestplot::forestplot(
		  name = outcome_clean, # swap 
		  estimate = b,
		  se = se,
		  pvalue = pval,
		  psignif = 0.05,
		  xlab = "OR per doubling in disease susceptibility",
		  colour = phenotype, # swapped
		  logodds = TRUE
		) +
        scale_x_continuous(labels = formatter(nsmall = 2), trans = "log")

    } else if(outcome_type == "continuous") {
    plot <- dt %>%
        filter(outcome %in% continuous_outcomes) %>%
        ggforestplot::forestplot(
		  name = outcome_clean,
		  estimate = b,
		  se = se,
		  pvalue = pval,
		  psignif = 0.05,
		  xlab = paste0("Beta per doubling in disease susceptibility"),
		  colour = phenotype,
		  logodds = FALSE,
		) +
		scale_x_continuous(labels = formatter_continuous())
    }

	plot <- plot +
		theme(
		  legend.position = "none",
		  axis.text.x = element_text(size=11),
		  axis.title.x = element_text(size=14),
		  axis.text.y = element_text(size=13),
		  strip.text = element_text(size = 14)
		) +
		facet_wrap(
			facets = ~ phenotype, nrow = nrows, ncol = ncols,
			scales = "free_x" # lets the axis vary for each facet
		) +
		scale_color_manual(values = viridis(14)[1:12])
		# make titles unbold -
		#theme(strip.text.x = element_text(face = "plain"))

    return(plot)

}

### Figure 1
# Primary outcomes - binary
mr_res_all_main %>%
   run_fplot(., 1, "binary", 6, 2) %>%
   ggsave("primary-all.jpeg", plot = .,
    path = "./results/autoimmune/forest_plots/main/",
    width = 35, height = 40, units = "cm")

### Supplementary Figures 1 & 2
# Secondary outcomes - binary
mr_res_all_main %>%
   run_fplot(., 2, "binary", 6, 2) %>%
   ggsave("secondary-binary.jpeg", plot = .,
    path = "./results/autoimmune/forest_plots/main/",
    width = 35, height = 30, units = "cm")
	
# Secondary outcomes - continuous
mr_res_all_main %>%
   run_fplot(., 2, "continuous", 6, 2) %>%
   ggsave("secondary-continuous.jpeg", plot = .,
    path = "./results/autoimmune/forest_plots/main/",
    width = 35, height = 20, units = "cm")

###############################################################################
#                           Inspect results                                   #
###############################################################################

table(mr_res_all_main$exposure)
mr_res_all_main %>% filter(exposure == "axsp") %>% arrange(or)

table(mr_res_all_main$outcome)
mr_res_all_main %>% filter(outcome == "gdm_subsamp")
mr_res_all_main %>% filter(outcome == "pretb_all")
mr_res_all_main %>% filter(outcome == "lowapgar5")
mr_res_all_main %>% filter(outcome == "nicu")
mr_res_all_main %>% filter(outcome == "gh_subsamp")
mr_res_all_main %>% filter(outcome == "pe_subsamp")
mr_res_all_main %>% filter(outcome == "lbw_all")
mr_res_all_main %>% filter(outcome == "zbw_all")
mr_res_all_main %>% filter(outcome == "zbw_all", id.exposure == "ra")
mr_res_all_main %>% filter(outcome == "zbw_all", id.exposure == "sle")
mr_res_all_main %>% filter(outcome == "ga_all")
mr_res_all_main %>% filter(outcome == "ga_all", id.exposure == "ra")
mr_res_all_main %>% filter(outcome == "ga_all", id.exposure == "sle")
mr_res_all_main %>% filter(outcome == "ga_all", id.exposure == "ms")

###############################################################################
#               Main analysis outcome n for Supp Table 2                      #
###############################################################################

outcome_n <- mr_res_all_main %>%
    filter(exposure == "t1d") %>%
    # very little variation between conditions so just picked one at random
    select(outcome_clean, min_n, mean_n, median_n, max_n, median_ncases,
    median_ncontrol, `Sample size`, `Cases / controls`)

write.csv(outcome_n, file.path(results_dir, "outcome_n.csv"), row.names = FALSE)

q('no')

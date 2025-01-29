###############################################################################
#                                                                             #
#   Plot heatmap for MR results for autoimmune disease - pregnancy outcomes   #
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

### Main analyses and sensitivities results
mr_res_binary <- read.csv(file.path(results_dir, "/mr_results_binary.csv"))
mr_res_binary$`OR (95% CI)` <- mr_res_binary$`OR..95..CI.`
mr_res_binary$`Cases / controls` <- mr_res_binary$`Cases...controls`

mr_res_cont <- read.csv(file.path(results_dir, "/mr_results_cont.csv"))
mr_res_cont$`Beta (95% CI)` <- mr_res_cont $`Beta..95..CI.`
mr_res_cont $`Sample size` <- mr_res_cont $`Sample.size`

exp_info <- read.csv(file.path(results_dir, "/exp_instrument_info.csv"))

mr_res_all <- rbind.fill(mr_res_cont, mr_res_binary) %>%
   mutate(stars = case_when(
		pval < 0.05 ~ "*",
		pval < 0.001 ~ "**",
		pval < 0.0001 ~ "***",
		pval < 0.00001 ~ "****")) %>% 
    left_join(exp_info) %>% # by id.exposure
	  mutate(phenotype_full = paste0(phenotype, "\n", Nsnps, " SNPs, mean F statistic ", mean_F))

mr_res_all_main <- mr_res_all %>%
    filter(analysis == "Main" & method == "Inverse variance weighted") %>%
	mutate(
		outcome_clean = factor(outcome_clean, levels = unique(outcome_clean)),
		OR = exp(b))

###############################################################################
#                 Heatmap - binary & continuous combined                      #
###############################################################################

heatmap_df <- mr_res_all_main %>%
    # add new lines into long names so they fit into plots
    mutate(phenotype = case_when(
        phenotype == "Systemic lupus erythematosus" ~ "Systemic lupus\n erythematosus",
        phenotype == "Inflammatory bowel disease" ~ "Inflammatory bowel\n disease",
        .default = phenotype))

### Primary outcomes
### All exposures incl AxSp (dominates)
p1_all <- heatmap_df %>%
  filter(outcome %in% primary_outcomes) %>% 
	ggplot(., aes(phenotype, outcome_clean, fill = OR)) + 
	  geom_tile() +
	  colorspace::scale_fill_continuous_diverging("Blue-Red 3",
	  	mid = 1) + # default mid = 0, for log odds scale
	  theme_minimal() +
	  theme(
			legend.title = element_text(size=20),
			legend.text = element_text(size=20),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size=13, angle = 360), #330
			axis.text.y = element_text(size=17),
			strip.text = element_blank(),
			# to remove grid lines -
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) + 
		scale_y_discrete(limits=rev) +
		geom_text(aes(label=stars), color="black", size=6, angle = 180, na.rm = T) +
		facet_grid(. ~ phenotype, scales = "free_x", space = "free_x")
ggsave("heatmap_all_primary_one_scale.jpeg", plot = p1_all,
    path = "./results/autoimmune/heatmap/", width = 60,  height = 30, units = "cm")

### All exposures except AxSp
p1 <- heatmap_df %>%
  # remove AxSp since betas so high, drowns out other results
  filter(exposure != "axsp") %>%
  filter(outcome %in% primary_outcomes) %>% 
	ggplot(., aes(phenotype, outcome_clean, fill = OR)) + 
	  geom_tile() +
	  colorspace::scale_fill_continuous_diverging("Blue-Red 3",
	  	mid = 1) + # default mid = 0, for log odds scale
	  theme_minimal() +
	  theme(
			legend.title = element_text(size=20),
			legend.text = element_text(size=20),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size=13, angle = 360), #330
			axis.text.y = element_text(size=17),
			strip.text = element_blank(),
			# to remove grid lines -
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) + 
		scale_y_discrete(limits=rev) +
		geom_text(aes(label=stars), color="black", size=6, angle = 180, na.rm = T) +
		facet_grid(. ~ phenotype, scales = "free_x", space = "free_x")
ggsave("heatmap_all_primary_no_axsp.jpeg", plot = p1,
    path = "./results/autoimmune/heatmap/", width = 60,  height = 30, units = "cm")
### AxSp only
p1_axsp <- mr_res_all_main %>%
    filter(exposure == "axsp") %>%
    filter(outcome %in% primary_outcomes) %>%
    mutate(phenotype = "Axial\n spondyloarthritis") %>% 
	ggplot(., aes(phenotype, outcome_clean, fill = OR)) + 
	  geom_tile() +
	  colorspace::scale_fill_continuous_diverging("Blue-Red 3",
	  	mid = 1) +
	  theme_minimal() +
	  theme(
			legend.title = element_text(size=20),
			legend.text = element_text(size=20),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size=13, angle = 360), #330
			axis.text.y = element_blank(), #element_text(size=16),
			strip.text = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) + 
		scale_y_discrete(limits=rev) +
		geom_text(aes(label=stars), color="black", size=6, angle = 180, na.rm = T) +
		facet_grid(. ~ phenotype, scales = "free_x", space = "free_x")
### Combined all exposures for primary outcomes
ggsave("heatmap_all_primary.jpeg",
	plot = (p1 + p1_axsp) + plot_layout(widths = c(45, 8)),
    path = "./results/autoimmune/heatmap/", width = 53,  height = 30, units = "cm")

### All exposures for secondary binary outcomes
p2 <- heatmap_df %>%
    filter(outcome %in% secondary_outcomes) %>% 
	filter(!(outcome %in% continuous_outcomes)) %>%
	ggplot(., aes(phenotype, outcome_clean, fill = OR)) + 
	  geom_tile() +
	  colorspace::scale_fill_continuous_diverging("Blue-Red 3",
	  	mid = 1) +
	  theme_minimal() +
	  theme(
			legend.title = element_text(size=16),
			legend.text = element_text(size=11),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size=14.5, angle = 360), #330
			axis.text.y = element_text(size=16),
			strip.text = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) + 
		scale_y_discrete(limits=rev) +
		geom_text(aes(label=stars), color="black", size=6, angle = 180, na.rm = T) +
		facet_grid(. ~ phenotype, scales = "free_x", space = "free_x")
ggsave("heatmap_all_secondary_binary.jpeg", plot = p2,
    path = "./results/autoimmune/heatmap/", width = 53,  height = 25, units = "cm")

### All exposures for secondary cont outcomes
p3 <- heatmap_df %>%
    filter(outcome %in% secondary_outcomes) %>% 
	filter(outcome %in% continuous_outcomes) %>%
	ggplot(., aes(phenotype, outcome_clean, fill = b)) + 
	  geom_tile() +
	  colorspace::scale_fill_continuous_diverging("Blue-Red 3") +
	  theme_minimal() +
	  theme(
			legend.title = element_text(size=16),
			legend.text = element_text(size=11),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size=14.5, angle = 360), #330
			axis.text.y = element_text(size=16),
			strip.text = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
			) + 
		scale_y_discrete(limits=rev) +
		geom_text(aes(label=stars), color="black", size=6, angle = 180, na.rm = T) +
		facet_grid(. ~ phenotype, scales = "free_x", space = "free_x")
ggsave("heatmap_all_secondary_cont.jpeg", plot = p3,
    path = "./results/autoimmune/heatmap/", width = 53,  height = 10, units = "cm")

q('no')

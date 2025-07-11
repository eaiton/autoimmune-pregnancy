###############################################################################
#                                                                             #
#                     Generate main MR results forest plots                   #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
x <- c("paletteer", "forestplot", "dplyr", "data.table", "tidyr", "grid",
        "svglite", "stringr")
lapply(x, require, character.only = TRUE)

# Set date for plots
date <- "260625"

# Set directories
home_dir <- paste0(Sys.getenv("AUTOIMMUNE_DIR"), "working/")
results_dir <- file.path(home_dir, "results/")
outcome_data_dir <- file.path(home_dir, "data/outcome/")
out_dir <- file.path(home_dir, "results/autoimmune/forest_plots/")
setwd(home_dir)

# outcomes
primary_outcomes <- read.table(paste0(outcome_data_dir, "primary_outcomes.txt"))$V1
secondary_outcomes <- read.table(paste0(outcome_data_dir, "secondary_outcomes.txt"))$V1 
continuous_outcomes <- c("zbw_all", "ga_all")

# phenotype order for plotting, ordered by precision (exposure GWAS N cases)
phenotype_order <-
  c("Multiple sclerosis", "Rheumatoid arthritis", "Type 1 diabetes",
    "Psoriasis", "Inflammatory bowel disease", "Systemic lupus erythematosus",
    "Coeliac disease", "Hashimoto's thyroiditis", "Systemic sclerosis",
    "Ankylosing spondylitis")

###############################################################################
#                               Read in results                               #
###############################################################################

### All main analyses and sensitivities results
mr_res_binary <- read.csv(file.path(results_dir, "/mr_results_binary.csv")) %>%
  mutate(`OR..95..CI.` = str_replace(`OR..95..CI.`, " - ", ", "))

mr_res_cont <- read.csv(file.path(results_dir, "/mr_results_cont.csv")) %>%
  mutate(`Beta..95..CI.` = str_replace(`Beta..95..CI.`, " - ", ", "))

exp_info <- read.csv(file.path(results_dir, "/exp_instrument_info.csv"))

mr_res_all <-
  plyr::rbind.fill(mr_res_cont, mr_res_binary) %>%
  # Results to plot here
  filter(analysis == "Main" & method == "Inverse variance weighted")

### Using rheumatoid arthritis results from the combined European and East Asian
### exposure GWAS (72% Eu ancestry) since this is higher powered than European
### only and over 95% of risk loci were shared across ancestries
mr_res_all <- mr_res_all %>% filter(exposure != "ra-eu")
mr_res_all$phenotype[mr_res_all$exposure == "ra-eu-eas"] <- "Rheumatoid arthritis"

###############################################################################
#                         Multiple testing correction                         #
###############################################################################

### Primary outcomes
primary <- mr_res_all %>% filter((outcome %in% primary_outcomes))
# Apply FDR 5% correction to primary outcomes
primary$pval.fdr <- p.adjust(primary$pval, method = "BH")
# Add asterisks to OR for plots:
primary$OR..95..CI.[primary$pval.fdr < 0.05] <- paste0(primary$OR..95..CI.[primary$pval.fdr < 0.05], "**")
primary$OR..95..CI.[primary$pval.fdr >= 0.05 & primary$pval < 0.05] <- paste0(primary$OR..95..CI.[primary$pval.fdr >= 0.05 & primary$pval < 0.05], "* ")
primary$OR..95..CI.[primary$pval >= 0.05] <- paste0(primary$OR..95..CI.[primary$pval >= 0.05], "  ")
# Select 18 results for follow up:
primary_follow_up <- primary %>%
  filter(pval<0.05 | ((or<0.95 | or>1.05) & (pval>0.05 & pval<0.1)))
write.csv(primary_follow_up, "results/autoimmune/primary_follow_up.csv", row.names = FALSE)

### Secondary outcomes
secondary <- mr_res_all %>% filter((outcome %in% secondary_outcomes))
# Apply FDR 5% correction to secondary outcomes
secondary$pval.fdr <- p.adjust(secondary$pval, method = "BH")
# Add asterisks to binary secondary outcomes (OR) for plots:
secondary_binary <- drop_na(secondary, or)
secondary_binary$OR..95..CI.[secondary_binary$pval.fdr < 0.05] <- paste0(secondary_binary$OR..95..CI.[secondary_binary$pval.fdr < 0.05], "**")
secondary_binary$OR..95..CI.[secondary_binary$pval.fdr >= 0.05 & secondary_binary$pval <= 0.05] <-
  paste0(secondary_binary$OR..95..CI.[secondary_binary$pval.fdr >= 0.05 & secondary_binary$pval < 0.05], "* ")
secondary_binary$OR..95..CI.[secondary_binary$pval >= 0.05] <- paste0(secondary_binary$OR..95..CI.[secondary_binary$pval > 0.05], "  ")
View(secondary_binary)
# Add asterisks to continuous outcomes (Beta) for plots:
secondary_cont <- drop_na(secondary, Sample.size)
secondary_cont$Beta..95..CI.[secondary_cont$pval.fdr < 0.05] <- paste0(secondary_cont$Beta..95..CI.[secondary_cont$pval.fdr < 0.05], "**")
secondary_cont$Beta..95..CI.[secondary_cont$pval.fdr >= 0.05 & secondary_cont$pval < 0.05] <-
  paste0(secondary_cont$Beta..95..CI.[secondary_cont$pval.fdr >= 0.05 & secondary_cont$pval < 0.05], "* ")
secondary_cont$Beta..95..CI.[secondary_cont$pval >= 0.05] <- paste0(secondary_cont$Beta..95..CI.[secondary_cont$pval > 0.05], "  ")
View(secondary_cont)

###############################################################################
#                             Format data for plots                           #
###############################################################################

### function to make alltext
### columns: condition, outcome, method (optional), effect & 95% CI

make_alltext <- function(df, type, labelled_row) {
    # type: 1 binary, 2 continuous
    # labelled_row: outcome where condition label should be

    # same as below -
    # add empty rows and order by condition, outcome, then method
    empty_rows <- data.frame(phenotype = unique(df$phenotype))
    df <- plyr::rbind.fill(empty_rows, df) %>%
        arrange(match(phenotype, rep(phenotype_order, 15)),
            match(outcome_clean, df$outcome_clean),
            match(method, df$method)) %>%
    # replace phenotype with "" except for row where we want the label
        mutate(phenotype =
            case_when(outcome_clean == labelled_row ~ phenotype,
            .default = ""))

    if(type == 1){ # binary
        df <- df %>% select(phenotype, outcome_clean, method, OR..95..CI.)
        # header titles
        df <- rbind(rep(NA, 4), df)
        df[1,] <- c("Condition", "Outcome", "Method", "OR (95% CI)  ")  
    } else if(type == 2){ # continuous
        df <- df %>% select(phenotype, outcome_clean, method, Beta..95..CI.)
        # header titles
        df <- rbind(rep(NA, 4), df)
        df[1,] <- c("Condition", "Outcome", "Method", "Beta (95% CI)  ")
    }
    
    # return alltext
    df
}

#### function to make resforalltext - adapted version of above
### columns: mean, upper, lower

make_resforalltext <- function(df, type, labelled_row) {
    # type: 1 binary, 2 continuous
    # labelled_row: outcome where condition label should be
    
    if(type == 1){ # binary
        df$mean <- df$or
        df$lower <- df$or_lci95
        df$upper <- df$or_uci95
    } else if(type == 2){ # continuous
        df$mean <- df$b
        df$lower <- df$lo_ci
        df$upper <- df$up_ci
    }

    # same as above -
    # add empty rows and order by condition, outcome, then method
    empty_rows <- data.frame(phenotype = unique(df$phenotype))
    df <- plyr::rbind.fill(empty_rows, df) %>%
        arrange(match(phenotype, rep(phenotype_order, 15)),
            match(outcome_clean, df$outcome_clean),
            match(method, df$method)) %>%
    # replace phenotype with "" except for row where we want the label
        mutate(phenotype =
            case_when(outcome_clean == labelled_row ~ phenotype,
            .default = ""))
    
    # select numeric columns, add leading NA
    df <- df %>% select(mean, lower, upper)
    df <- rbind(rep(NA, 3), df)
    df <- sapply(df, as.numeric)
    
    # return alltext
    df
}

###############################################################################
#                              Plot formatting                                #
###############################################################################

# Outcome as title, exposures as y axis

# Define palette
nb.cols <- 16
mycolours <- colorRampPalette(paletteer_d("MoMAColors::Klein"))(nb.cols)

###############################################################################
#                            Plot primary outcomes                            #
###############################################################################

# Format
df <- primary
alltext_1 <- make_alltext(df, 1, "Gestational diabetes mellitus")
resforalltext_1 <- make_resforalltext(df, 1, "Gestational diabetes mellitus")
# Styles
palette <- mycolours[1:9]
styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8]),
      gpar(col = palette[9], fill = palette[9])
    ), 10), # number of conditions
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8]),
    gpar(col=palette[9])
  ), 10)
)
# Plots
### all 10 conditions
ticks <- c(0.50, 0.70, 0.90, 1.0, 1.10, 1.3, 1.5, 1.7, 1.9)
all <- resforalltext_1[1:100, ] %>%
  forestplot(labeltext = alltext_1[1:100, c(1,2,4)], # no method label
             txt_gp = fpTxtGp(
               label = gpar(cex = .8, fontfamily="sans"),
               ticks = gpar(cex = .8, fontfamily="sans"),
               xlab = gpar(cex = .9, fontfamily="sans")),
             is.summary=c(TRUE, rep(FALSE, 99)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             fn.ci_norm = fpDrawCircleCI,
             xlog = TRUE,
             xticks=ticks,
             xlab = "Odds ratio per doubling in log odds of autoimmune condition",
             zero=1,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles)  |>
  fp_decorate_graph(grid = structure(ticks[-c(1, length(ticks))],
                                     gp = gpar(lty = 2, col = "grey50")))
print(all)

# Save
tiff(paste0(out_dir, "1_main_primary_bin_all_", date, ".tiff"),
     units="in", width=11, height=18, res=400)
plot.new()
all
dev.off()

svglite(paste0(out_dir, "1_main_primary_bin_all_", date, ".svg"),
        system_fonts = list(sans = "sans"),
        width = 11, height = 18)
plot.new()
all
dev.off()

###############################################################################
#                            Plot secondary outcomes                          #
###############################################################################
### Secondary binary outcomes

# Format
df <- secondary_binary
alltext_2 <- make_alltext(df, 1, "Preeclampsia")
resforalltext_2 <- make_resforalltext(df, 1, "Preeclampsia")

# Styles
palette <- mycolours[10:14]
styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5])
    ), 10), # number of conditions
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5])
  ), 10)
)
# Plots
ticks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
all <- resforalltext_2[1:60, ] %>%
  forestplot(labeltext = alltext_2[1:60, c(1,2,4)], # no method label
             txt_gp = fpTxtGp(
               label = gpar(cex = .8, fontfamily="sans"),
               ticks = gpar(cex = .8, fontfamily="sans"),
               xlab = gpar(cex = .9, fontfamily="sans")),
             is.summary=c(TRUE, rep(FALSE, 99)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             fn.ci_norm = fpDrawCircleCI,
             xlog = TRUE,
             xticks=ticks,
             xlab = "Odds ratio per doubling in log odds of autoimmune condition",
             zero=1,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles) |>
  fp_decorate_graph(grid = structure(ticks[-c(1, length(ticks))],
                                     gp = gpar(lty = 2, col = "grey50")))
print(all)

# Save
tiff(paste0(out_dir, "1_main_sec_bin_all_", date, ".tiff"),
     units="in", width=11, height=10, res=400)
plot.new()
all
dev.off()

svglite(paste0(out_dir, "1_main_sec_bin_all_", date, ".svg"),
        system_fonts = list(sans = "sans"),
        width = 11, height = 10)
plot.new()
all
dev.off()

### Secondary continuous outcomes

# Format
df <- secondary_cont
alltext_3 <- make_alltext(df, 2, "Birthweight")
resforalltext_3 <- make_resforalltext(df, 2, "Birthweight")

palette <- mycolours[1:2]
styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2])
    ), 10), # number of conditions
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2])
  ), 10)
)
# Plots
### first 5, plot 1
ticks <- c(-0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075)
all <- resforalltext_3[1:30, ] %>%
  forestplot(labeltext = alltext_3[1:30, c(1,2,4)], # no method label
             txt_gp = fpTxtGp(
               label = gpar(cex = .8, fontfamily="sans"),
               ticks = gpar(cex = .8, fontfamily="sans"),
               xlab = gpar(cex = .9, fontfamily="sans")),
             is.summary=c(TRUE, rep(FALSE, 99)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             fn.ci_norm = fpDrawCircleCI,
             xlog = FALSE,
             xticks=ticks,
             xlab = "Beta per doubling in log odds of autoimmune condition",
             zero=0,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles) |>
  fp_decorate_graph(grid = structure(ticks[-c(1, length(ticks))],
                                     gp = gpar(lty = 2, col = "grey50")))
print(all)

tiff(paste0(out_dir, "1_main_sec_cont_all_", date, ".tiff"),
     units="in", width=11, height=5.25, res=400)
plot.new()
all
dev.off()

svglite(paste0(out_dir, "1_main_sec_cont_all_", date, ".svg"),
        system_fonts = list(sans = "sans"),
        width = 11, height = 5.25)
plot.new()
all
dev.off()

q('no')
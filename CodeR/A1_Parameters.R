# A1_Parameters.R / Heat-Related Mortality / Exposure Response"
# Author: Hans-Aloys Wischmann, Institute for Public Health, Charité - Universitätsmedizin Berlin
# Written: 2023-04-12 to 2023-08-31 / Totally revised: 2024-11-07 to 2025-07-04

# ensure consistency across systems, define presets, set default figure size
Sys.setlocale("LC_ALL", 'en_US.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")
knitr::opts_chunk$set(echo = FALSE, dpi = 1200, comment = NA)
knitr::opts_knit$set(root.dir = getwd())

# create a clean slate
rm(list = ls()); start_time <- proc.time(); set.seed(314159)

# configurable parameters: years to use for modeling and prediction, statistical significance cutoff
YEARS_MODEL <- 2000:2023            # include COVID pandemic years, but with exclusions, for 1-year models
COVID_YEARS <- c(2020, 2022)        # exclude two years with COVID mortality impact during summer
LOADS_YEAR  <- min(COVID_YEARS) - 1 # show factor loadings etc. just before COVID
ALPHA       <- 0.05                 # standard 0.05 value

# configurable parameters: seasonal week range to analyze (summer)
FIRST_SUMMER_WEEK <- 15 # as in RKI model: begin of summer period
LAST_SUMMER_WEEK  <- 40 # as in RKI model: end of summer period
DOF_TMPC          <-  8 # as in RKI model: degrees of freedom for dose-response curves
DOF_WEEK          <-  5 # fewer than in RKI model
DOF_YEAR          <- 10 # as in RKI model: about 1 DoF per 2 calendar years
LAG_MAX           <-  3 # maximum number of days to look back for temperature exposure

# define color palettes for categorical data with 12 entries
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
safe_heatmap_palette    <- c("#0066FF", "#00CCFF", "#66FF00", "#FFFF00", "#FF7F00", "#FF0000", "#CF00FF", "#9F00FF")

# function to replace install.packages/library combination
library_wrapper <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) { install.packages(package) }
    library(package, character.only = TRUE, quietly = TRUE)
  }
}
library_wrapper(c("tidyverse", "flextable", "foreach", "doParallel", "RhpcBLASctl", "ISOweek", "sf"))
  
# use all physical cores for threads, up to max_threads, limit number of BLAS per thread
config_parallel <- function(max_threads = NA) {
  registerDoParallel(cores = min(get_num_cores(), max_threads, na.rm = TRUE))
  blas_set_num_threads(get_num_procs() / min(get_num_cores(), max_threads, na.rm = TRUE))
}
config_parallel()

# utility function to plot to *.pdf and *.png file and inline
ggplot_font_size = 10 # font size for all texts except for geom_text, in points
plot_pdf_png <- function(file_name, aspect_ratio, plot_object, plot_width = 6.25, plot_res = 1200) {
  themed_object <- plot_object + theme(text = element_text(size = ggplot_font_size), plot.title = element_text(size = ggplot_font_size))
  png(sprintf("../Plots/%s.png", file_name), width = plot_width, height = plot_width * aspect_ratio, units = "in", res = plot_res)
  print(themed_object)
  ignore <- dev.off()
  pdf(sprintf("../Plots/%s.pdf", file_name), plot_width, plot_width * aspect_ratio, paper = "a4")
  print(themed_object)
  ignore <- dev.off()
  print(themed_object)
}

# create a standard theme
std_theme <- function() {
  theme(
    panel.border     = element_rect(colour = "black", fill = NA),
    panel.background = element_rect(fill   = "white"),
    panel.grid.major = element_line(colour = "gray80", linewidth = 0.2),
    axis.text  = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "right"
  )
}

# free up memory: stop all clusters and perform a full garbage collection
stop_gc <- function() { stopImplicitCluster(); ignore <- gc(full = TRUE) }

# week lookup table
FIRST_WEEK <- sprintf("%d-W01-1", min(YEARS_MODEL))
LAST_WEEK  <- sprintf("%d-W52-1", max(YEARS_MODEL))
week_table <- data.frame(date_from = seq(ISOweek2date(FIRST_WEEK), ISOweek2date(LAST_WEEK), 7)) %>%
  mutate(date_to = date_from + 6, year = isoyear(date_from), week = isoweek(date_from))

# new year and mid year tick marks for multi-year plots, axis labels
MID_YEAR_DATES <- seq(as.Date(sprintf("%d-07-01", min(YEARS_MODEL))), as.Date(sprintf("%d-07-01", max(YEARS_MODEL))), by = "year")
NEW_YEAR_DATES <- seq(as.Date(sprintf("%d-01-01", min(YEARS_MODEL))), as.Date(sprintf("%d-01-01", max(YEARS_MODEL) + 1)), by = "year")
YEAR_BREAKS <- c(MID_YEAR_DATES, NEW_YEAR_DATES)
YEAR_LABELS <- c(YEARS_MODEL, rep("|", length(NEW_YEAR_DATES)))

# age ranges (standard plus coarse to match available data), plus mapping for aggregation
age_5y_standard <- data.frame(age_from = seq(0, 85, 5),    age_to = c(seq(4, 84, 5), 99)) %>% mutate(age = fct_inorder(paste(age_from, age_to, sep = "-")))
age_4cat_coarse <- data.frame(age_from = c(0, 65, 75, 85), age_to = c(64, 74, 84, 99))    %>% mutate(age = fct_inorder(paste(age_from, age_to, sep = "-")))
age_aggregate   <- left_join(age_5y_standard, age_4cat_coarse, join_by("age_from" >= "age_from", "age_to" <= "age_to")) %>% select(age = age.x, age_cat = age.y)

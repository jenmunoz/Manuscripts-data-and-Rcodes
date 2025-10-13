# ###_###_###_###_###_###_###_###_###_###_######_####_###_###_###_###_###_###_###_###_###_######_########_###_###_###_###_###_###_###_###_###_######_####_###_###_###_###_###_###_###_###_###_######_#
## Manuscript 
## Project: Mixed species flcoks and parasites
# Description:
# This script creates plots of the conditional effects for each of the sleected models 
# Author: J. Munoz
# Last updated:2025
#
# Notes:
# - Adapt paths to your working directory as needed.
# ###_###_###_###_###_###_###_###_###_###_######_####_###_###_###_###_###_###_###_###_###_######_########_###_###_###_###_###_###_###_###_###_######_####_###_###_###_###_###_###_###_###_###_######_#

# 0. Libraries --------------------------------------------------------------
# Check R version
R.Version()

# # Libraries for easier manipulation of data
# install.packages("pacman")
# library("pacman")
# install.packages('BiocManager')
# library("BiocManager")
# install.packages("tidyr")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("data.table")
# install.packages("extrafont")
# install.packages("lubridate")  # for dates
# 
# # Data cleaning libraries
# install.packages("janitor")
# install.packages("assertr")
# 
# # Libraries for data analyses and visualizations
# install.packages("vegan")
# install.packages("ggplot2")
# install.packages("devtools")
# install.packages("knitr")
# install.packages("ts")
# install.packages("RColorBrewer")
# install.packages("ggridges")
# install.packages("ggtree")
# install.packages("aplot")
# 
# # Libraries for models
# install.packages("car")  # ANOVA command
# install.packages("lattice")  # Preliminary plots
# install.packages("lme4")  # For glmer (generalized linear mixed models)
# install.packages("visreg")  # Extract confidence intervals and trend lines from GLMMs
# install.packages("lsmeans")  # Least squared means
# install.packages("MuMIn")  # Pseudo R squared for GLMMs
# install.packages("emmeans")
# 
# # Model assumptions
# install.packages("DHARMa")
# 
# # Bayesian models
# install.packages('tidybayes')
# install.packages('bayesplot')
# install.packages("rstan")
# install.packages('brmstools')
# remotes::install_github("Pakillo/DHARMa.helpers")
# devtools::install_github("mvuorre/brmstools")
# install.packages('rstanarm')
# install.packages('loo')
# devtools::install_github("paul-buerkner/brms")  # Use reloo=TRUE to remove pareto
# install.packages('brms')  # Bayesian approach to model phylogenetic data with repeated observations
# 
# # Phylogenetic component
# install.packages("ape")
# install.packages("here")
# install.packages("phytools")
# install.packages("tidyverse")
# install.packages("metafor")
# install.packages("phangorn")  # To reconstruct a maximum clade credibility tree
# install.packages("rr2")
# install.packages("MCMCglmm")
# remotes::install_github("daijiang/phyr")
# install.packages("TreeTools")

# Libraries for data
library(tidyverse)
library(ggplot2)
library(dplyr)
library(data.table)
library(extrafont)
library(lubridate)

# Data cleaning libraries
library(janitor)
library(assertr)

# Data visualization libraries
library(vegan)
library(ggplot2)
library(devtools)
library(knitr)
library(ts)
library(RColorBrewer)
library(ggridges)
library(ggtree)
library(aplot)

# Libraries for models and visualizations
library(lattice)  # Preliminary plots
library(car)  # ANOVA command
library(lsmeans)  # Least squared means
library(lme4)  # For glmer (generalized linear mixed models)
library(visreg)  # Extract confidence intervals and trend lines from GLMMs
library(MuMIn)  # Pseudo R squared for GLMMs
library(emmeans)

# Model assumption checks
library(DHARMa)

# Phylogenetic component libraries
library(ape)
library(here)
library(phyr)
library(phytools)
library(metafor)
library(phangorn)  # To reconstruct a maximum clade credibility tree
library(rr2)
library(MCMCglmm)
library(tidyverse)
library(skimr)
library(TreeTools)

# Libraries for plots
library(gridExtra)
library(ggpubr)
library(grid)

# Libraries for Bayesian models
library(bayesplot)
library(tidybayes)
library(brms)  # Bayesian approach to model phylogenetic data with repeated observations
library(DHARMa.helpers)
library(brmstools)
library(rstan)
library(rstanarm)
library(loo)

# to plot conditional effects
library(brms)  # Make sure you have the brms library loaded
library(ggplot2)  # For plotting
library(ggeffects)  # For plotting conditional effects

# Load necessary libraries
library(brms)
library(gridExtra)
library(cowplot)


?conditional_effects


# Selected models  ------------------------------------------------------------------

# If I want to plot the conditional effects of one variable in particular use this:

# sociality<-conditional_effects(selected_ecto_infection_brms_bayes_no_int, 'sociality_groups')  # 'sociality_groups', 'elevation', 'year_seasonality', 'mass_tidy_species', 'mass_ind_comp'
# 
# # Valid effects are (combinations of): 'sociality_groups', 'elevation', 'year_seasonality', 'mass_tidy_species', 'mass_ind_comp'"sociality_groups"
# plot_infection_sociality<-plot(sociality, plot = FALSE, points = TRUE,line_args = list(color = "black", size = 3,linetype = "solid"), point_args = list(width = .2, shape = 19, col="orange", size=5, alpha=0.7))[[1]] +
#   scale_color_grey() +
#   scale_fill_grey() +
#   xlab("Sociality")+
#   ylab("Presence Absence")+
#   scale_color_grey() +
#   theme_classic2(30)+
#   labs(caption = "Model: infection ~Categories")+
#   theme(plot.caption = element_text(size = 12))

# Models 
# Infection
selected_ecto_infection_brms_bayes_no_int<-readRDS("results/selected_models/1_M1P_model_INFECTION_bernu_brms_phylo_multiple_obs_no_interactions_priors_SELECTED_antfollowers_included.RDS")

# Abundance Lice
selected_zinb_a_lice_brms_bayes_no_int_priors<-readRDS("results/selected_models/3_M1L.model_brms_LICE_ABUNDANCE_zinb_a_lice_brms_bayes_no_int_priors_ind_mass_scaled_SELECTED_antbirds_included.RDS")

# Abundance non-feather mites
selected_zinb_a_nf_mites_brms_bayes_no_int_prior
selected_zinb_a_nf_mites_brms_bayes_no_int_prior<-readRDS("results/selected_models/3_M1MNF.model_prevalence_zinb_brms_ABUNDANCE_nf_MITES_phylo_multiple_obs_no_interactions_prior_SELECTED_antfollowers_included.RDS")

# Abundance all mites
selected_zinb_a_all_mites_brms_bayes_no_int_prior

# Prevalence
ecto_p_brms_bayes_no_int_species_priors_zobi
ecto_p_brms_bayes_no_int_species_priors_zobi<-readRDS("results/selected_models/P2s.model_prevalence_brms_phylo_SPECIES_no_interactions_priors_zobi_antbirds_included.RDS")

#conditional_effects(ecto_p_brms_bayes_no_int_species_priors_zobi, "new") # 'sociality', 'elevation', 'sample_size', 'mass'
# Lice richness
selected_poisson_lice_diversity_sociality_no_int_priors
selected_poisson_lice_diversity_sociality_no_int_priors<-readRDS( "results/selected_models/5_DL.model_lICE_diversity_brms_phylo_multiple_obs_no_interactions_NO_truncated_antfollowers_included.RDS")


####_####_###_###_###
### For networks
####_####_###_###_###


# Prevalence networks
ecto_p_brms_bayes_no_int_species_priors_degree_zobi
# Lice richness
selected_poisson_lice_diversity_sociality_no_int_priors
# Infection networks
selected_ecto_p_brms_bayes_no_int_degree_prior
# Abundance Lice networks
selected_zinb_a_lice_brms_bayes_no_int_degree_prior
# Abundance non-feather mites networks
selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior
# Lice richness networks 
selected_poisson_lice_diversity_degree_no_int_priors


# Plot the conditional effects---------------------------------------------------------
?conditional_effects
            
# Function to plot conditional effects for sociality

plot_conditional_effect_sociality <- function(model, variable, x_label, y_label, point_color, point_colors, y_lim = NULL) {
  cond_effect <- conditional_effects(model, variable)
  plot <- plot(cond_effect, plot = FALSE, points = TRUE, method = "posterior_epred", robust = TRUE,
               line_args = list(color = point_color, size = 3, linetype = "solid"),
               point_args = list(width = .2, shape = 19, col = point_color, size = 5, alpha = 0.7))[[1]] +
    scale_color_manual(point_color) +
    scale_fill_manual(point_color) +
    xlab(x_label) +
    ylab(y_label) +
    theme_classic2(30) +
    theme(axis.title.y = element_text(size = 23),
          plot.title = element_blank())  # Remove title
  if (!is.null(y_lim)) {
    plot <- plot + coord_cartesian(ylim = y_lim)  # Set y-axis limits
  }
  return(plot)
}


# for teh rest of teh variables conditioned to sociality 
plot_conditional_effect <- function(model, variable, x_label, y_label, point_color, y_lim = NULL) {
  cond_effect <- conditional_effects(model, variable, conditions = make_conditions(model, vars = c("sociality_groups")))
  plot <- plot(cond_effect, plot = FALSE, points = TRUE, method = "posterior_epred", robust = TRUE,
               line_args = list(color = point_color, size = 3, linetype = "solid"),
               point_args = list(width = .2, shape = 19, col = point_color, size = 5, alpha = 0.7))[[1]] +
    scale_color_manual(point_color) +
    scale_fill_manual(point_color) +
    xlab(x_label) +
    ylab(y_label) +
    theme_classic2(30) +
    theme(axis.title.y = element_text(size = 23))
  if (!is.null(y_lim)) {
    plot <- plot + coord_cartesian(ylim = y_lim)  # Set y-axis limits
  }
  return(plot)
}

# for prevalence since teh variable is named differently 
plot_conditional_effect_prevalence <- function(model, variable, x_label, y_label, point_color, y_lim = NULL) {
  cond_effect <- conditional_effects(model, variable, conditions = make_conditions(model, vars = c("sociality")))
  plot <- plot(cond_effect, plot = FALSE, points = TRUE, method = "posterior_epred", robust = TRUE,
               line_args = list(color = point_color, size = 3, linetype = "solid"),
               point_args = list(width = .2, shape = 19, col = point_color, size = 5, alpha = 0.7))[[1]] +
    scale_color_manual(point_color) +
    scale_fill_manual(point_color) +
    xlab(x_label) +
    ylab(y_label) +
    theme_classic2(30) +
    theme(axis.title.y = element_text(size = 23),
          plot.title = element_blank())  # Remove title
  if (!is.null(y_lim)) {
    plot <- plot + coord_cartesian(ylim = y_lim)  # Set y-axis limits
  }
  return(plot)
}



# Define a function to remove y-axis labels from all but the first plot
remove_axis_label <- function(plot) {
  plot + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
}

remove_y_axis_label <- function(plot) {
  plot + theme(axis.title.y = element_blank())
}

remove_x_axis_label <- function(plot) {
  plot + theme(axis.title.x = element_blank())
}

# Plotting conditional effects for each variable with respective y-axis limits
plot_infection_sociality <- plot_conditional_effect_sociality(selected_ecto_infection_brms_bayes_no_int, 'sociality_groups', "Sociality", "Infection", "#0098BA", c(0, 1)) %>% remove_x_axis_label()
plot_infection_elevation <- plot_conditional_effect(selected_ecto_infection_brms_bayes_no_int, 'elevation', "Elevation", "Ectoparasite Infection", "darkgray", c(0, 1)) %>% remove_axis_label()
plot_infection_mass_tidy_species <- plot_conditional_effect(selected_ecto_infection_brms_bayes_no_int, 'mass_tidy_species', " Bird Mass (gr)", "Ectoparasite Infection", "darkgray", c(0, 1)) %>% remove_axis_label()
plot_infection_mass_ind_comp <- plot_conditional_effect(selected_ecto_infection_brms_bayes_no_int, 'mass_ind_comp', " Individual Mass", "Ectoparasite Infection", "darkgray", c(0, 1)) %>% remove_axis_label()
plot_infection_year_seasonality <- plot_conditional_effect(selected_ecto_infection_brms_bayes_no_int, 'year_seasonality', "Seasonality", "Ectoparasite Infection", "#800020", c(0, 1)) %>% remove_axis_label()

plot_liceabundance_sociality <- plot_conditional_effect_sociality(selected_zinb_a_lice_brms_bayes_no_int_priors, 'sociality_groups', "Sociality", "Lice abundance", "#0098BA", c(0, 70)) %>% remove_x_axis_label()
plot_liceabundance_elevation <- plot_conditional_effect(selected_zinb_a_lice_brms_bayes_no_int_priors, 'elevation', "Elevation", "Lice abundance", "darkgray", c(0, 70)) %>% remove_axis_label()
plot_liceabundance_mass_tidy_species <- plot_conditional_effect(selected_zinb_a_lice_brms_bayes_no_int_priors, 'mass_tidy_species', "Bird Mass (gr)", "Lice abundance", "#800020", c(0, 70)) %>% remove_axis_label()
plot_liceabundance_mass_ind_comp <- plot_conditional_effect(selected_zinb_a_lice_brms_bayes_no_int_priors, 'mass_ind_comp', " Individual Mass", "Lice abundance", "darkgray", c(0, 70)) %>% remove_axis_label()
plot_liceabundance_year_seasonality <- plot_conditional_effect(selected_zinb_a_lice_brms_bayes_no_int_priors, 'year_seasonality', "Seasonality", "Lice abundance", "darkgray", c(0, 70)) %>% remove_axis_label()

plot_mitesabundance_sociality <- plot_conditional_effect_sociality(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, 'sociality_groups', "Sociality", "Mites abundance", "#0098BA", c(0, 50))
plot_mitesabundance_elevation <- plot_conditional_effect(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, 'elevation', "Elevation", "mites abundance", "#800020", c(0, 50)) %>% remove_y_axis_label()
plot_mitesabundance_mass_tidy_species <- plot_conditional_effect(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, 'mass_tidy_species', "Bird Mass (gr)", "mites abundance", "darkgray", c(0, 50)) %>% remove_y_axis_label()
plot_mitesabundance_mass_ind_comp <- plot_conditional_effect(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, 'mass_ind_comp', " Individual Mass", "Mites abundance", "darkgray", c(0, 50)) %>% remove_y_axis_label()
plot_mitesabundance_year_seasonality <- plot_conditional_effect(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, 'year_seasonality', "Seasonality", "mites abundance", "darkgray", c(0, 50)) %>% remove_y_axis_label()

plot_prevalence_sociality <- plot_conditional_effect_sociality(ecto_p_brms_bayes_no_int_species_priors_zobi, 'sociality', "Sociality", "[Species level] Prevalence", "#0098BA", c(0, 1))%>% remove_x_axis_label()
plot_prevalence_elevation <- plot_conditional_effect_prevalence (ecto_p_brms_bayes_no_int_species_priors_zobi, 'elevation', "Elevation", "[Species level] Prevalence", "darkgray",c(0, 1)) %>% remove_axis_label()
plot_prevalence_mass_tidy_species <- plot_conditional_effect_prevalence (ecto_p_brms_bayes_no_int_species_priors_zobi, 'mass', "Bird Mass (gr)",  "[Species level] Prevalence", "darkgray", c(0, 1)) %>% remove_axis_label()
#plot_prevalence_samplesize <- plot_conditional_effect(ecto_p_brms_bayes_no_int_species_priors_zobi, 'sample_size', " Individual Mass", "mites abundance", "darkgray", c(0, 50)) %>% remove_y_axis_label()
plot_empty<-plot(1, type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')

#conditional_effects(selected_poisson_lice_diversity_sociality_no_int_priors,"new") # 'sociality_groups', 'elevation_midpoint', 'total_sample_size', 'mass_tidy_species'
plot_licerichness_sociality <- plot_conditional_effect_sociality(selected_poisson_lice_diversity_sociality_no_int_priors, 'sociality_groups', "Sociality", "Lice richness", "#0098BA", c(0, 3)) %>% remove_x_axis_label()
plot_licerichness_elevation <- plot_conditional_effect(selected_poisson_lice_diversity_sociality_no_int_priors, 'elevation_midpoint', "Elevation", "Lice richness", "darkgray", c(0, 3)) %>% remove_axis_label()
plot_licerichness_mass_tidy_species <- plot_conditional_effect(selected_poisson_lice_diversity_sociality_no_int_priors, 'mass_tidy_species', "Bird Mass (gr)", "Lice richness", "darkgray", c(0, 3)) %>% remove_axis_label()
plot_empty1<-plot(1, type = "n", xlab = "Individual Bird Mass (gr)", ylab = "", xaxt = 'n', yaxt = 'n')%>% remove_y_axis_label()
plot_empty2<-plot(1, type = "n", xlab = "Seasonality", ylab = "", xaxt = 'n', yaxt = 'n')%>% remove_y_axis_label()

?conditional_effects
#plot_licerichness_sample_size <- plot_conditional_effect(selected_poisson_lice_diversity_sociality_no_int_priors, 'sample_size', "Sample_size", "Lice richness", "darkgray", c(0, 70)) %>% remove_y_axis_label()


# Arrange the plots in a grid: one row, five columns, ensuring they fit within the y-axis limits
combined_plot <- grid.arrange(
  
  plot_infection_sociality,
  plot_infection_elevation,
  plot_infection_mass_tidy_species,
  plot_infection_mass_ind_comp,
  plot_infection_year_seasonality,
  
  plot_prevalence_sociality,
  plot_prevalence_elevation,
  plot_prevalence_mass_tidy_species,
  plot_empty,
  plot_empty,
  
  plot_liceabundance_sociality,
  plot_liceabundance_elevation,
  plot_liceabundance_mass_tidy_species,
  plot_liceabundance_mass_ind_comp,
  plot_liceabundance_year_seasonality,
  
  plot_licerichness_sociality,
  plot_licerichness_elevation, 
  plot_licerichness_mass_tidy_species,
  plot_empty1,
  plot_empty2,
  
  plot_mitesabundance_sociality,
  plot_mitesabundance_elevation, 
  plot_mitesabundance_mass_tidy_species, 
  plot_mitesabundance_mass_ind_comp,
  plot_mitesabundance_year_seasonality,
  ncol= 5,
  nrow = 5)


# Save the final plot

ggsave("Figure3_Conditional_effects_plots_ectoparasites.png", plot = combined_plot, width = 20, height = 15, units = "in", dpi = 300)

# here I conditioned to social species # condition= "social"
ggsave("Figure3a_Conditional_effects_plots_ectoparasites_social_condition.png", plot = combined_plot, width = 20, height = 15, units = "in", dpi = 300)

ggsave("Figure3b_Conditional_effects_plots_ectoparasites_both_conditions.png", plot = combined_plot, width = 30, height = 15, units = "in", dpi = 300)


##3 experimenting 

library(gridExtra)
library(ggplot2)
library(cowplot)

# Create empty plots to use as spacers with dashed lines
empty_plot_with_line <- ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +
  theme_void()

# Combine plots with spacers to add dashed lines between rows and columns
combined_plot <- plot_grid(
  plot_infection_sociality, plot_infection_elevation, plot_infection_mass_tidy_species, plot_infection_mass_ind_comp, plot_infection_year_seasonality,
  plot_prevalence_sociality, plot_prevalence_elevation, plot_prevalence_mass_tidy_species, empty_plot_with_line, empty_plot_with_line,
  plot_liceabundance_sociality, plot_liceabundance_elevation, plot_liceabundance_mass_tidy_species, plot_liceabundance_mass_ind_comp, plot_liceabundance_year_seasonality,
  plot_licerichness_sociality, plot_licerichness_elevation, plot_licerichness_mass_tidy_species, empty_plot_with_line, empty_plot_with_line,
  plot_mitesabundance_sociality, plot_mitesabundance_elevation, plot_mitesabundance_mass_tidy_species, plot_mitesabundance_mass_ind_comp, plot_mitesabundance_year_seasonality,
  ncol = 5,
  nrow = 5,
  align = "h"
)




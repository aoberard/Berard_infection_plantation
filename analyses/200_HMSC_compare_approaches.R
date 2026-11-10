# Global parameters ----


# Remove saved elements (clean workspace)
rm(list = ls())
gc()

# Load libraries
library("dplyr")
library("ggplot2")
if (!requireNamespace("Hmsc", quietly = TRUE)) {
  install.packages("Hmsc", repos = "https://cloud.r-project.org")
}
library("Hmsc")
library("patchwork")

# Make script reproducible
set.seed(1)

#Graphical parameters

#Dual scale for each type of model
pal_prob <- c(low = "#2C7796", mid = "#F7F7F7", high = "#D55763")
pal_ab <- c(low = "#2C7796", mid = "#F7F7F7", high = "#D55763")
# pal_ab   <- c(low = "#3E5641", mid = "#F7F7F7", high = "#BB7E8C")




# Load saved elements ----

#Setting directories
dataDir <- here::here("data", "derived-data", "hmsc_out")
dataDir_ab <- here::here("data", "derived-data", "hmsc_ab_out")

#Loading pairwise comparisons
load(here::here(dataDir, "results", "prediction", "pairwise_comparisons.RData"))
load(here::here(dataDir_ab, "results", "prediction", "pairwise_comparisons.RData"))

#Loading predictive power per species
load(file = here::here(dataDir, "results", "model_fit", "species_order_tjurCV.RData"))
load(file = here::here(dataDir_ab, "results", "model_fit", "species_order_srCV.RData"))

#Loading Varpart plots
load(file = here::here(dataDir, "results", "parameters", "parameters_beta_probit_plot.RData"))
load(file = here::here(dataDir_ab, "results", "parameters", "parameters_beta_ab_plot.RData"))


# (Optional) keeping only taxa explained well enough  ----
#Here based on intersection between the two models

# #Extract not well explained species for each model
# model_name <- "habitat_hostmodel_spat"
# sp_inf <- species_order[[model_name]]$tjur#_cv
# good_sp_inf <- names(sp_inf[sp_inf >= 0.1])
# 
# model_name <- "habitatxseason_spat"
# sp_ab <- species_ab_order[[model_name]]$sr#_cv 
# good_sp_ab  <- names(sp_ab[sp_ab >= 0.1])
# 
# good_sp <- intersect(good_sp_inf, good_sp_ab)
# 
# #Remove them from each dataframe
# df_ab_habitat <- df_ab_habitat %>%
#   filter(species %in% good_sp)
# 
# df_habitat <- df_habitat %>%
#   filter(species %in% good_sp)
# 
# df_ab_season <- df_ab_season %>%
#   filter(species %in% good_sp)
# 
# df_season <- df_season %>%
#   filter(species %in% good_sp)
# 
# df_interaction <- df_interaction %>%
#   filter(species %in% good_sp)
# 
# df_ab_interaction <- df_ab_interaction %>%
#   filter(species %in% good_sp)



# Season comparison heatmap ----

## Generate dataframe ----

#Combine the two modelling approaches into one long dataframe
df_heatmap_season <- bind_rows(
  df_season %>% mutate(model = "probability"),
  df_ab_season %>% mutate(model = "abundance")
)

#Order species (by average difference here)
species_order_season <- df_heatmap_season %>%
  group_by(species) %>%
  summarise(mean_effect = mean(abs(diff_mean))) %>%
  arrange(desc(mean_effect)) %>%
  pull(species)

#Factorize and add labels
df_heatmap_season <- df_heatmap_season %>%
  mutate(species = factor(species, levels = species_order_season),
    model = factor(model, levels = c("probability", "abundance"),
                     labels = c("Probability", "Abundance")),
    habitat = factor(habitat, levels = c("pine_edge", "hedgerows", "broadleaved_forest"),
                     labels = c("Pine edges", "Hedgerows", "Broadleaved forests")))


## Heatmap plotting ----

#Version with dual scale
df_prob_season <- df_heatmap_season %>% filter(model == "Probability")
df_ab_season <- df_heatmap_season %>% filter(model == "Abundance")

lim_prob_season <- max(abs(df_prob_season$diff_mean), na.rm = TRUE)
lim_ab_season <- max(abs(df_ab_season$diff_mean),  na.rm = TRUE)

plot_heatmap_season_dualscale <- ggplot() +
  
  #Infection prob tiles
  geom_tile(data = df_prob_season,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0,
    limits = c(-lim_prob_season, lim_prob_season),
    name = "Δ P(infection)"
  ) +
  geom_text(data = df_prob_season,
            aes(x = model, y = habitat,
                label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  
  #Individual number tiles
  ggnewscale::new_scale_fill() +
  
  geom_tile(data = df_ab_season,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_ab["low"], mid = pal_ab["mid"], high = pal_ab["high"],
    midpoint = 0,
    limits = c(-lim_ab_season, lim_ab_season),
    name = "Δ E(infected individuals)"
  ) +
  geom_text(data = df_ab_season,
            aes(x = model, y = habitat,
                label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  
  facet_wrap(~ species,
             labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL,
       title = "Season effect (Spring - Autumn)") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_season_dualscale




# Habitat comparison heatmap ----

## Generate dataframe ----

#Combine the two modelling approaches into one long dataframe
df_heatmap_habitat <- bind_rows(
  df_habitat %>% 
    mutate(model = "probability"),
  df_ab_habitat %>% 
    mutate(model = "abundance"))

#Create pair label
df_heatmap_habitat <- df_heatmap_habitat %>%
  mutate(pair_label = paste(habitat1, "-", habitat2))

#Order species (by average difference here)
species_order_habitat <- df_heatmap_habitat %>%
  group_by(species) %>%
  summarise(mean_effect = mean(abs(diff_mean))) %>%
  arrange(desc(mean_effect)) %>%
  pull(species)

#Factorize and add labels
df_heatmap_habitat <- df_heatmap_habitat %>%
  mutate(
    species = factor(species, levels = species_order_habitat),
    model = factor(model, levels = c("probability", "abundance"),
                     labels = c("Probability", "Abundance")),
    season = factor(season, levels = c("Spring", "Autumn"))
  )


## Heatmap plotting ----

#Version with dual scale
df_prob_habitat <- df_heatmap_habitat %>% filter(model == "Probability")
df_ab_habitat <- df_heatmap_habitat %>% filter(model == "Abundance")

lim_prob_habitat <- max(abs(df_prob_habitat$diff_mean), na.rm = TRUE)
lim_ab_habitat <- max(abs(df_ab_habitat$diff_mean),  na.rm = TRUE)

plot_heatmap_habitat_dualscale <- ggplot() +
  
  #Infection prob tiles
  geom_tile(data = df_prob_habitat,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0,
    limits = c(-lim_prob_habitat, lim_prob_habitat),
    name = "Δ P(infection)"
  ) +
  geom_text(data = df_prob_habitat,
            aes(x = model, y = pair_label,
                label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  
  #Individual number tiles
  ggnewscale::new_scale_fill() +
  
  geom_tile(data = df_ab_habitat,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_ab["low"], mid = pal_ab["mid"], high = pal_ab["high"],
    midpoint = 0,
    limits = c(-lim_ab_habitat, lim_ab_habitat),
    name = "Δ E(infected individuals)"
  ) +
  geom_text(data = df_ab_habitat,
            aes(x = model, y = pair_label,
                label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  facet_grid(season ~ species, 
             labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Habitat effect across species, season and model") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    strip.text.y = element_text(color = "#F7F7F7"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_habitat_dualscale





# Interaction comparison ----

## Generate dataframe ----

#Combine the two modelling approaches into one long dataframe
df_heatmap_interaction <- bind_rows(
  df_interaction %>% mutate(model = "probability"),
  df_ab_interaction %>% mutate(model = "abundance")
)

#Create pair label
df_heatmap_interaction <- df_heatmap_interaction %>%
  mutate(pair_label = paste(habitat1, "-", habitat2))

#Order species (by average interaction effect here)
species_order_interaction <- df_heatmap_interaction %>%
  group_by(species) %>%
  summarise(mean_effect = mean(abs(interaction_mean))) %>%
  arrange(desc(mean_effect)) %>%
  pull(species)

#Factorize and add labels
df_heatmap_interaction <- df_heatmap_interaction %>%
  mutate(species = factor(species, levels = species_order_interaction),
         model = factor(model, levels = c("probability", "abundance"),
                     labels = c("Probability", "Abundance")))



## Heatmap plotting ----

#Version with dual scale
df_prob_interaction <- df_heatmap_interaction %>% filter(model == "Probability")
df_ab_interaction <- df_heatmap_interaction %>% filter(model == "Abundance")

lim_prob_interaction <- max(abs(df_prob_interaction$interaction_mean), na.rm = TRUE)
lim_ab_interaction <- max(abs(df_ab_interaction$interaction_mean),  na.rm = TRUE)

plot_heatmap_interaction_dualscale <- ggplot() +
  
  #Infection prob tiles
  geom_tile(data = df_prob_interaction,
            aes(x = model, y = pair_label, fill = interaction_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0,
    limits = c(-lim_prob_interaction, lim_prob_interaction),
    name = "ΔΔ P(infection)"
  ) +
  geom_text(data = df_prob_interaction,
            aes(x = model, y = pair_label,
                label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  
  #Individual number tiles
  ggnewscale::new_scale_fill() +
  
  geom_tile(data = df_ab_interaction,
            aes(x = model, y = pair_label, fill = interaction_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_ab["low"], mid = pal_ab["mid"], high = pal_ab["high"],
    midpoint = 0,
    limits = c(-lim_ab_interaction, lim_ab_interaction),
    name = "ΔΔ E(infected individuals)"
  ) +
  geom_text(data = df_ab_interaction,
            aes(x = model, y = pair_label,
                label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season x Habitat interaction across species and model") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_interaction_dualscale










# Combine plots ----

## Heatmaps ----


#Order species 
#This time considering both season mean differences and habitat mean differences (still both approaches)
species_order_combined <- bind_rows(
  df_heatmap_habitat %>% select(species, diff_mean),
  df_heatmap_season %>% select(species, diff_mean)
) %>%
  group_by(species) %>%
  summarise(mean_effect = mean(abs(diff_mean), na.rm = TRUE)) %>%
  arrange(desc(mean_effect)) %>%
  pull(species)

df_heatmap_habitat <- df_heatmap_habitat %>%
  mutate(species = factor(species, levels = species_order_combined))
df_heatmap_season <- df_heatmap_season %>%
  mutate(species = factor(species, levels = species_order_combined))

#Split back dataframes that were previously pooled in previous section of the script
df_prob_habitat <- df_heatmap_habitat %>% filter(model == "Probability")
df_ab_habitat <- df_heatmap_habitat %>% filter(model == "Abundance")
df_prob_season <- df_heatmap_season  %>% filter(model == "Probability")
df_ab_season <- df_heatmap_season  %>% filter(model == "Abundance")


#Choose limits for legend based (in order to make only one shared legend after)
lim_prob_global <- max(abs(df_prob_habitat$diff_mean), abs(df_prob_season$diff_mean), na.rm = TRUE)
lim_ab_global <- max(abs(df_ab_habitat$diff_mean), abs(df_ab_season$diff_mean),na.rm = TRUE)


#Modify species label:
species_labels <- c(
  Bartonella_taylorii = "B. taylorii",
  Bartonella_grahamii = "B. grahamii",
  Bartonella_birtlesii = "B. birtlesii",
  Eperythrozoon_coccoides = "E. coccoides",
  Neoehrlichia_mikurensis= "N. mikurensis",
  Eperythrozoon_muris = "E. muris",
  Eperythrozoon_2 = "Eperythrozoon 2",
  Bartonella_doshiae = "B. doshiae",
  Leptospira = "Leptospira",
  Sarcocystidae =  "Sarcocystidae",
  Eperythrozoon_4 ="Eperythrozoon 4",
  Borrelia_miyamotoi = "B. miyamotoi"
)



#Rebuild heatmaps, but with common limits / legend here :

#First panel : Season effect
plot_heatmap_season_dualscale_sharedlegend <- ggplot() +
  geom_tile(data = df_prob_season,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0, limits = c(-lim_prob_global, lim_prob_global),
    name = "Δ P(infection)"
  ) +
  geom_text(data = df_prob_season,
            aes(x = model, y = habitat, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = df_ab_season,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_ab["low"], mid = pal_ab["mid"], high = pal_ab["high"],
    midpoint = 0, limits = c(-lim_ab_global, lim_ab_global),
    name = "Δ E(infected individuals)"
  ) +
  geom_text(data = df_ab_season,
            aes(x = model, y = habitat, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  facet_grid(. ~ species,
             labeller = labeller(species = species_labels)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    axis.text.x = element_blank(),   # repeated after so masked here
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

#Second panel : Habitat effect across season
plot_heatmap_habitat_dualscale_sharedlegend <- ggplot() +
  geom_tile(data = df_prob_habitat,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0, limits = c(-lim_prob_global, lim_prob_global),
    name = "Δ P(infection)"
  ) +
  geom_text(data = df_prob_habitat,
            aes(x = model, y = pair_label, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = df_ab_habitat,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_ab["low"], mid = pal_ab["mid"], high = pal_ab["high"],
    midpoint = 0, limits = c(-lim_ab_global, lim_ab_global),
    name = "Δ E(infected individuals)"
  ) +
  geom_text(data = df_ab_habitat,
            aes(x = model, y = pair_label, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  facet_grid(season ~ species,
             labeller = labeller(species = species_labels)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    strip.text.y = element_text(color = "#F7F7F7"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


# Assemble both panel
plot_composed <- plot_heatmap_season_dualscale_sharedlegend / plot_heatmap_habitat_dualscale_sharedlegend +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold"))

plot_composed




## Varpartition plot----

# Regenerate no legend version of previously generated plots
p_probit_beta_noleg <- p_probit_beta + theme(legend.position = "none")
p_probit_vp_noleg <- p_probit_vp + theme(legend.position = "none")
p_probit_r2_noleg <- p_probit_r2 + theme(legend.position = "none")
p_ab_beta_noleg <- p_ab_beta + theme(legend.position = "none")
p_ab_vp_noleg <- p_ab_vp + theme(legend.position = "none")
p_ab_r2_noleg <- p_ab_r2 + theme(legend.position = "none")

#Reextract legends (and make them smaller)
leg_probit_beta <- cowplot::get_legend(p_probit_beta + theme(legend.position = "bottom",
                                                                 legend.title = element_text(size = 6),
                                                                 legend.text = element_text(size = 5),
                                                                 legend.key.size = unit(0.3, "cm"),
                                                                 legend.margin = margin(0, 0, 0, 0),
                                                                 legend.box.margin = margin(-2, -2, -2, -2)))
leg_probit_vp <- cowplot::get_legend(p_probit_vp+ theme(legend.position = "bottom",
                                                        legend.title = element_text(size = 6),
                                                        legend.text = element_text(size = 5),
                                                        legend.key.size = unit(0.3, "cm"),
                                                        legend.margin = margin(0, 0, 0, 0),
                                                        legend.box.margin = margin(-2, -2, -2, -2)))
leg_ab_beta <- cowplot::get_legend(p_ab_beta + theme(legend.position = "bottom",
                                                     legend.title = element_text(size = 6),
                                                     legend.text = element_text(size = 5),
                                                     legend.key.size = unit(0.3, "cm"),
                                                     legend.margin = margin(0, 0, 0, 0),
                                                     legend.box.margin = margin(-2, -2, -2, -2)))
# leg_ab_vp <- cowplot::get_legend(p_ab_vp+ theme(legend.position = "bottom",
#                                                 legend.title = element_text(size = 6),
#                                                 legend.text = element_text(size = 5),
#                                                 legend.key.size = unit(0.3, "cm"),
#                                                 legend.margin = margin(0, 0, 0, 0),
#                                                 legend.box.margin = margin(-2, -2, -2, -2)))

#Panel assemblage
plot_full <- p_probit_beta_noleg + coord_fixed(ratio = 1)

plot_full <- plot_full %>%
  aplot::insert_left(p_tree, width = 0.1) %>%
  aplot::insert_right(p_probit_vp_noleg, width = 0.5) %>%
  aplot::insert_right(p_probit_r2_noleg + coord_fixed(ratio = 1), width = 0.05) %>%
  aplot::insert_right(p_ab_beta_noleg +
                        theme(axis.text.y = element_blank()) +
                        coord_fixed(ratio = 1),
                      width = 1) %>%
  aplot::insert_right(p_ab_vp_noleg, width = 0.5) %>%
  aplot::insert_right(p_ab_r2_noleg + coord_fixed(ratio = 1), width = 0.05) %>%
  aplot::set_panel_spacing(spacing = 2) %>%
  ggplotify::as.ggplot()


#Legend section
plot_all_legends <- cowplot::plot_grid(
  leg_probit_beta, leg_ab_beta, leg_probit_vp,
  nrow = 1,
  rel_widths = c(0.5, 0.5, 1)
)

#Assemblage
plot_varpart_final <- cowplot::plot_grid(
  plot_full,
  plot_all_legends,
  ncol = 1,
  rel_heights = c(1, 0.1))

plot_varpart_final




# Supplementary CWM heatmap ----

## Generate dataframes ----

#Season comparison heatmap
##Generate dataframe
#Combine the two modelling approaches into one long dataframe
df_heatmap_season_cwm <- bind_rows(
  df_season_cwm %>% 
    mutate(model = "probability"),
  df_ab_season_cwm %>% 
    mutate(model = "abundance")
)

#Factorize and add labels
df_heatmap_season_cwm <- df_heatmap_season_cwm %>%
  mutate(trait = factor(trait),
         model = factor(model, levels = c("probability", "abundance"),
                        labels = c("Probability", "Abundance")),
         habitat = factor(habitat, levels = c("pine_edge", "hedgerows", "broadleaved_forest"),
                          labels = c("Pine edges", "Hedgerows", "Broadleaved forests")))


# Habitat comparison heatmap 
##Generate dataframe 

#Combine the two modelling approaches into one long dataframe
df_heatmap_habitat_cwm <- bind_rows(
  df_habitat_cwm %>% 
    mutate(model = "probability"),
  df_ab_habitat_cwm %>% 
    mutate(model = "abundance"))

#Create pair label
df_heatmap_habitat_cwm <- df_heatmap_habitat_cwm %>%
  mutate(pair_label = paste(habitat1, "-", habitat2))

#Factorize and add labels
df_heatmap_habitat_cwm <- df_heatmap_habitat_cwm %>%
  mutate(
    trait = factor(trait),
    model = factor(model, levels = c("probability", "abundance"),
                   labels = c("Probability", "Abundance")),
    season = factor(season, levels = c("Spring", "Autumn")))


#Remove intercept from both dataframe (optional)
df_heatmap_season_cwm <- df_heatmap_season_cwm %>% 
  filter(! trait == "(Intercept)")
df_heatmap_habitat_cwm <- df_heatmap_habitat_cwm %>% 
  filter(! trait == "(Intercept)")


#Order trait for COMBINED plot version
#This time considering both season mean differences and habitat mean differences (still both approaches)
trait_order_combined <- bind_rows(
  df_heatmap_habitat_cwm %>% select(trait, diff_mean),
  df_heatmap_season_cwm %>% select(trait, diff_mean)
) %>%
  group_by(trait) %>%
  summarise(mean_effect = mean(abs(diff_mean), na.rm = TRUE)) %>%
  arrange(desc(mean_effect)) %>%
  pull(trait)

df_heatmap_habitat_cwm <- df_heatmap_habitat_cwm %>%
  mutate(trait = factor(trait, levels = trait_order_combined))
df_heatmap_season_cwm <- df_heatmap_season_cwm %>%
  mutate(trait = factor(trait, levels = trait_order_combined))

#Split back dataframes that were previously pooled in previous section of the script
df_prob_habitat_cwm <- df_heatmap_habitat_cwm %>% filter(model == "Probability")
df_ab_habitat_cwm <- df_heatmap_habitat_cwm %>% filter(model == "Abundance")
df_prob_season_cwm <- df_heatmap_season_cwm  %>% filter(model == "Probability")
df_ab_season_cwm <- df_heatmap_season_cwm  %>% filter(model == "Abundance")


#Choose limits for legend based (in order to make only one shared legend after)
lim_global_cwm <- max(
  abs(df_heatmap_season_cwm$diff_mean),
  abs(df_heatmap_habitat_cwm$diff_mean))

#Modify trait label:
trait_labels <- c(
  '(Intercept)' =  "Intercept",
  non_vector_transmission1 = "Without vector",
  vector_transmission1 = "With vector"
)


## Heatmap plotting (combined)----

#Rebuild heatmaps, but with common limits / legend here :

#First panel : Season effect
plot_heatmap_season_dualscale_sharedlegend_cwm <- ggplot() +
  geom_tile(data = df_prob_season_cwm,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "#F7F7F7") +
  geom_text(data = df_prob_season_cwm,
            aes(x = model, y = habitat, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  geom_tile(data = df_ab_season_cwm,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0, limits = c(-lim_global_cwm, lim_global_cwm),
    name = "Δ CWM"
  ) +
  geom_text(data = df_ab_season_cwm,
            aes(x = model, y = habitat, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  facet_grid(. ~ trait,
             labeller = labeller(trait = trait_labels)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    axis.text.x = element_blank(),   # repeated after so masked here
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

#Second panel : Habitat effect across season
plot_heatmap_habitat_dualscale_sharedlegend_cwm <- ggplot() +
  geom_tile(data = df_prob_habitat_cwm,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "#F7F7F7") +
  geom_text(data = df_prob_habitat_cwm,
            aes(x = model, y = pair_label, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  geom_tile(data = df_ab_habitat_cwm,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "#F7F7F7") +
  scale_fill_gradient2(
    low = pal_prob["low"], mid = pal_prob["mid"], high = pal_prob["high"],
    midpoint = 0, limits = c(-lim_global_cwm, lim_global_cwm),
    name = "Δ CWM"
  ) +
  geom_text(data = df_ab_habitat_cwm,
            aes(x = model, y = pair_label, label = if_else(significant, "*", "")),
            size = 5, vjust = 0.75, color = "black") +
  facet_grid(season ~ trait,
             labeller = labeller(trait = trait_labels)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "#F7F7F7"),
    strip.text.y = element_text(color = "#F7F7F7"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


# Assemble both panel
plot_composed_cwm <- plot_heatmap_season_dualscale_sharedlegend_cwm / plot_heatmap_habitat_dualscale_sharedlegend_cwm +
  plot_layout(guides = "collect",
              heights = c(1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold"))

plot_composed_cwm


# Supplementary Omega plot ----

## Generate dataframes ----
omega_combined_df <- bind_rows(
  omega_ab_plot_df %>% mutate(model = "Probability"),
  omega_probit_plot_df %>% mutate(model = "Abundance")
)
  
## Heatmap plotting ----

#Generate heatmaps
p_omega_combined <- omega_combined_df %>%
  filter(species_1 != species_2) %>%  # to remove diagonal
  # filter(as.integer(species_1) > as.integer(species_2)) %>%  # to remove bottom triangle
  ggplot(aes(x = species_2, y = species_1, fill = mean)) +
  geom_tile(aes(fill = mean, alpha = support > 0 | supportNeg > 0), 
            color = "white", linewidth = 0.4) +
  geom_text(
    data = . %>% filter(significant),
    aes(label = "*"),
    color = "black", size = 3, vjust = 0.75
  ) +
  scale_fill_gradient2(
    low = "#2C7796",
    mid = "#F7F7F7",
    high= "#D55763",
    midpoint = 0,
    limits= c(-1, 1),
    name= "Residual association"
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0), guide = "none") + # to make not well supported things transparent
  facet_grid(
    rows = vars(model),
    cols = vars(random_level))+
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x= element_text(angle = 45, hjust = 1, size = 7, color = "grey40"),
    axis.text.y = element_text(size = 7, color = "grey40"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(fill = "grey30"),
    strip.text = element_text(color = "white", size = 7),
    legend.position  = "bottom"
  )

p_omega_combined


# Save plots ----

if(!dir.exists(here::here("figures", "prediction_heatmaps"))) {
  dir.create(here::here("figures", "prediction_heatmaps"), recursive = TRUE)
}

ggsave(filename = here::here("figures", "prediction_heatmaps", "heatmap_season.png"),
       plot = plot_heatmap_season_dualscale + coord_equal(),
       width = 2700, height = 2500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "prediction_heatmaps", "heatmap_habitat.png"),
       plot = plot_heatmap_habitat_dualscale + coord_equal(),
       width = 6300, height = 1750, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "prediction_heatmaps", "heatmap_interaction.png"),
       plot = plot_heatmap_interaction_dualscale + coord_equal(),
       width = 2900, height = 2500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "prediction_heatmaps", "heatmap_composed.pdf"),
       plot = plot_composed & coord_equal(),
       width = 4500, height = 2500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "parameters", "varpart_composed.pdf"),
       plot = plot_varpart_final,
       width = 3000, height = 2000, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "parameters", "heatmap_cwm_composed.pdf"),
       plot = plot_composed_cwm & coord_equal(),
       width = 1500, height = 2500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "parameters", "omega_composed.pdf"),
       plot = p_omega_combined ,
       width = 2500, height = 1500, units = "px", dpi = 300, bg = NULL)


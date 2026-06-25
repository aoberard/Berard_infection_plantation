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

# Make script reproducible
set.seed(1)

#Graphical parameters

#Dual scale for each type of model
pal_prob <- c(low = "#2C7796", mid = "#F7F7F7", high = "#D55763")
pal_ab   <- c(low = "#2C7796", mid = "#F7F7F7", high = "#D55763")
# pal_ab   <- c(low = "#3E5641", mid = "#F7F7F7", high = "#BB7E8C")




# Load saved elements ----

#Setting directories
dataDir <- here::here("data", "derived-data", "hmsc_out")
dataDir_ab <- here::here("data", "derived-data", "hmsc_ab_out")
resultDir <- file.path(dataDir, "results", "prediction")
resultDir_ab <- file.path(dataDir_ab, "results", "prediction")


#Loading pairwise comparisons
load(here::here(resultDir,"pairwise_comparisons.RData"))
load(here::here(resultDir_ab,"pairwise_comparisons.RData"))

#Loading predictive power per species
load(file = here::here(dataDir, "results", "model_fit", "species_order_tjurCV.RData"))
load(file = here::here(dataDir_ab, "results", "model_fit", "species_order_srCV.RData"))


# (Optional) keeping only explained species ----
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

#Add a scaled version of mean differences value (here scaled across species per model)
df_season <- df_season %>%
  mutate(diff_scaled = as.numeric(scale(diff_mean)))
df_ab_season <- df_ab_season %>%
  mutate(diff_scaled = as.numeric(scale(diff_mean)))


#Combine the two modelling approaches into one long dataframe
df_heatmap_season <- bind_rows(
  df_season    %>% mutate(model = "probability"),
  df_ab_season %>% mutate(model = "abundance")
)

#Add support as a signed metric, so that when the effect is mostly negative, supports appears negative
#(because diff_mean per default are not on the same scale for both models)
df_heatmap_season <- df_heatmap_season %>%
  mutate(signed_support = if_else(support_pos > support_neg, support_pos, -support_neg))


#Order species (by average difference here)
species_order_season <- df_heatmap_season %>%
  group_by(species) %>%
  summarise(mean_effect = mean(abs(diff_mean))) %>%
  arrange(desc(mean_effect)) %>%
  pull(species)

#Factorize and add labels
df_heatmap_season <- df_heatmap_season %>%
  mutate(
    species = factor(species, levels = species_order_season),
    model = factor(model, levels = c("probability", "abundance"),
                     labels = c("Probability", "Abundance")),
    habitat = factor(habitat, levels = c("pine_edge", "hedgerows", "broadleaved_forest"),
                     labels = c("Pine edges", "Hedgerows", "Broadleaved forests"))
    )


## Heatmap plotting ----

#Version with diff mean (not on same scale)
plot_heatmap_season_diff <- df_heatmap_season %>%
  ggplot(aes(x = model , y = habitat, fill = diff_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Mean difference"
  ) +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season effect (Spring - Autumn) across species, habitat type and model") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_heatmap_season_diff


#Version with diff mean (initially scaled)
plot_heatmap_season_diffscaled <- df_heatmap_season %>%
  ggplot(aes(x = model , y = habitat, fill = diff_scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Scaled mean difference"
  ) +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season effect (Spring - Autumn) across species, habitat type and model (scaled)") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_heatmap_season_diffscaled


#version with support
plot_heatmap_season_support <- df_heatmap_season %>%
  ggplot(aes(x = model , y = habitat, fill = signed_support)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Signed support"
  ) +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season effect (Spring - Autumn) across species, habitat type and model") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

plot_heatmap_season_support


#Version with dual scale
df_prob_season <- df_heatmap_season %>% filter(model == "Probability")
df_ab_season   <- df_heatmap_season %>% filter(model == "Abundance")

lim_prob_season <- max(abs(df_prob_season$diff_mean), na.rm = TRUE)
lim_ab_season   <- max(abs(df_ab_season$diff_mean),  na.rm = TRUE)

plot_heatmap_season_dualscale <- ggplot() +
  
  #Infection prob tiles
  geom_tile(data = df_prob_season,
            aes(x = model, y = habitat, fill = diff_mean),
            color = "white") +
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
            color = "white") +
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
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_season_dualscale




# Habitat comparison heatmap ----

## Generate dataframe ----

#Add a scaled version of mean differences value (here scaled across species per model)
df_habitat <- df_habitat %>%
  mutate(diff_scaled = as.numeric(scale(diff_mean)))
df_ab_habitat <- df_ab_habitat %>%
  mutate(diff_scaled = as.numeric(scale(diff_mean)))


#Combine the two modelling approaches into one long dataframe
df_heatmap_habitat <- bind_rows(
  df_habitat %>% 
    mutate(model = "probability"),
  df_ab_habitat %>% 
    mutate(model = "abundance")
)

#Add support as a signed metric, so that when the effect is mostly negative, supports appears negative
#(because diff_mean are not on the same scale for both models)
df_heatmap_habitat <- df_heatmap_habitat %>%
  mutate(signed_support = if_else(support_pos > support_neg, support_pos, -support_neg))

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

#Version with diff mean (not on same scale)
plot_heatmap_habitat_diff <- df_heatmap_habitat %>%
  ggplot(aes(x = model, y = pair_label, fill = diff_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Mean difference"
  ) +
  facet_grid(season ~ species, 
             labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Habitat effect across species, season and model") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    strip.text.y = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_habitat_diff


#Version with diff mean (not on same scale)
plot_heatmap_habitat_diffscaled <- df_heatmap_habitat %>%
  ggplot(aes(x = model, y = pair_label, fill = diff_scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Scaled mean difference"
  ) +
  facet_grid(season ~ species, 
             labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Habitat effect across species, season and model (scaled)") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    strip.text.y = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_habitat_diffscaled


#Version with support 
plot_heatmap_habitat_support <- df_heatmap_habitat %>%
  ggplot(aes(x = model, y = pair_label, fill = signed_support)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Signed support"
  ) +
  facet_grid(season ~ species, 
             labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Habitat effect (signed posterior support) across species, season and model") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    strip.text.y = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_habitat_support


#Version with dual scale
df_prob_habitat <- df_heatmap_habitat %>% filter(model == "Probability")
df_ab_habitat   <- df_heatmap_habitat %>% filter(model == "Abundance")

lim_prob_habitat <- max(abs(df_prob_habitat$diff_mean), na.rm = TRUE)
lim_ab_habitat   <- max(abs(df_ab_habitat$diff_mean),  na.rm = TRUE)

plot_heatmap_habitat_dualscale <- ggplot() +
  
  #Infection prob tiles
  geom_tile(data = df_prob_habitat,
            aes(x = model, y = pair_label, fill = diff_mean),
            color = "white") +
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
            color = "white") +
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
    strip.text.x = element_text(color = "white"),
    strip.text.y = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_habitat_dualscale





# Interaction comparison ----

## Generate dataframe ----

#Add a scaled version of mean differences value (here scaled across species per model)
df_interaction <- df_interaction %>%
  mutate(diff_scaled = as.numeric(scale(interaction_mean)))
df_ab_interaction <- df_ab_interaction %>%
  mutate(diff_scaled = as.numeric(scale(interaction_mean)))

#Combine the two modelling approaches into one long dataframe
df_heatmap_interaction <- bind_rows(
  df_interaction    %>% mutate(model = "probability"),
  df_ab_interaction %>% mutate(model = "abundance")
)

#Add support as a signed metric, so that when the effect is mostly negative, supports appears negative
#(because interaction_mean are not on the same scale for both models)
df_heatmap_interaction <- df_heatmap_interaction %>%
  mutate(signed_support = if_else(support_pos > support_neg, support_pos, -support_neg))

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

#Version with interaction mean (not on same scale)
plot_heatmap_interaction_diff <- df_heatmap_interaction %>%
  ggplot(aes(x = model, y = pair_label, fill = interaction_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Difference of difference"
  ) +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season x Habitat interaction across species and model") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_interaction_diff


#Version with interaction mean (scaled)
plot_heatmap_interaction_diffscaled <- df_heatmap_interaction %>%
  ggplot(aes(x = model, y = pair_label, fill = diff_scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Scaled difference of difference"
  ) +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season x Habitat interaction across species and model (scaled)") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_interaction_diffscaled


#Version with support
plot_heatmap_interaction_support <- df_heatmap_interaction %>%
  ggplot(aes(x = model, y = pair_label, fill = signed_support)) +
  geom_tile(color = "white") +
  geom_text(aes(label = if_else(significant, "*", "")), 
            size = 5, vjust = 0.75, color = "black") +
  scale_fill_gradient2(
    low = "#2C7796", mid = "white", high = "#D55763", 
    midpoint = 0,
    name = "Signed support"
  ) +
  facet_wrap(~ species, labeller = labeller(species = function(x) gsub("_", " ", x))) +
  labs(x = NULL, y = NULL, 
       title = "Season x Habitat interaction (signed posterior support) across species and model") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_interaction_support



#Version with dual scale
df_prob_interaction <- df_heatmap_interaction %>% filter(model == "Probability")
df_ab_interaction   <- df_heatmap_interaction %>% filter(model == "Abundance")

lim_prob_interaction <- max(abs(df_prob_interaction$interaction_mean), na.rm = TRUE)
lim_ab_interaction   <- max(abs(df_ab_interaction$interaction_mean),  na.rm = TRUE)

plot_heatmap_interaction_dualscale <- ggplot() +
  
  #Infection prob tiles
  geom_tile(data = df_prob_interaction,
            aes(x = model, y = pair_label, fill = interaction_mean),
            color = "white") +
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
            color = "white") +
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
    strip.text.x = element_text(color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot_heatmap_interaction_dualscale



# Save plots ----
ggsave(filename = here::here("figures", "heatmap_season.png"),
       plot = plot_heatmap_season_dualscale + coord_equal(),
       width = 2700, height = 2500, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures",  "heatmap_habitat.png"),
       plot = plot_heatmap_habitat_dualscale + coord_equal(),
       width = 6300, height = 1750, units = "px", dpi = 300, bg = NULL)

ggsave(filename = here::here("figures", "heatmap_interaction.png"),
       plot = plot_heatmap_interaction_dualscale + coord_equal(),
       width = 2900, height = 2500, units = "px", dpi = 300, bg = NULL)



# Script parameters ----

#Remove saved elements
rm(list = ls())
gc()

# Library
library("dplyr")
library("ggplot2")


#Graphical parameters 
category_palette <- c("pine_edge"          = "#3f826d",
                      "hedgerows"          = "#f2d0a4",
                      "broadleaved_forest" = "#a4a8d1")

season_palette <- c("Spring" = "#7FBB7B",
                    "Autumn" = "#9D695A")


category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
session_order <- c("Juin_2023", "Octobre_2023", "Juin_2024", "Septembre_2024", "Juin_2025")
season_order <- c("Spring", "Autumn")


# Import data ----

#Apodemus sylvaticus infectious data

inf <- data.table::fread(file = here::here("data", "derived-data", "asylvaticus_data", "asylvaticus_infectious.csv"))

names_inf <- readRDS(file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo.rds"))
names_inf10 <- readRDS(file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo10.rds"))


# Generate data ----

#Compute richness for every infectious agent
inf_alpha <- inf %>% 
  mutate(rich = vegan::specnumber(
    inf %>%
      as_tibble() %>%
      select(all_of(names_inf)) %>%
      as.matrix()
  )) %>% 
  mutate(rich10 = vegan::specnumber(
    inf %>%
      as_tibble() %>%
      select(all_of(names_inf10)) %>%
      as.matrix()
  )) %>% 
  select(-all_of(names_inf))


#Dataframes for analysis (no NA for explanatory variables)
df_alpha <- inf_alpha %>% 
  filter(!is.na(sexe)) 

df_alpha2 <- df_alpha %>% 
  mutate(category =  
           case_when(
             category == "pine_edge" ~ "plantation",
             category == "hedgerows" ~ "plantation",
             .default = category
           ))


#Choose reference levels
df_alpha$year <- relevel(as.factor(df_alpha$year), ref = "2024")
df_alpha$season <- relevel(as.factor(df_alpha$season), ref = "Autumn")
df_alpha$category <- relevel(as.factor(df_alpha$category), ref = "hedgerows")

#Dataset with joined category
df_alpha2 <- df_alpha %>% 
  mutate(category =  
           case_when(
             category == "pine_edge" ~ "plantation",
             category == "hedgerows" ~ "plantation",
             .default = category
           ))

# Visualization of distribtion ----

# Calculate percentage
percent_at_least_2 <- mean(df_alpha$rich >= 2) * 100
label_text <- paste0(round(percent_at_least_2, 1), "% with ≥2 infectious agents")

# Create histogram with annotation
richness_plot <- ggplot(df_alpha, aes(x = rich)) +
  geom_histogram(binwidth = 1, fill = "#fddc62", color = "#7a6a55", linewidth = 1.2) +
  labs(
    title = "Distribution of infectious agents richness per individual",
    x = "Number of infectious agents detected",
    y = "Number of individuals"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.1, vjust = 2,
    size = 5,
    fontface = "italic",
    color = "#7a6a55"
  ) 
richness_plot



# GLM : richness variation ----

#Check effectifs
df_alpha %>% 
  group_by(category, season) %>% 
  count()

#Models' fixed effects formula
XFormula.null = rich ~ 1 + (1 | numero_ligne) + (1 | year) + (1|numero_centre)
XFormula.host = rich ~ sexe + poids + maturity + (1 | numero_ligne) + (1 | year) + (1|numero_centre)
XFormula.hab =rich ~ category + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year) + (1|numero_centre)
XFormula.season = rich ~ season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year) + (1|numero_centre)
XFormula.hab_season_add = rich ~  category + season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year) + (1|numero_centre)
XFormula.hab_season_int = rich ~  category * season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year) + (1|numero_centre)


#Specification of GLMMs
m.null <- glmmTMB::glmmTMB(
  XFormula.null,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.null)

m.host <- glmmTMB::glmmTMB(
  XFormula.host,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.host)

m.hab <- glmmTMB::glmmTMB(
  XFormula.hab,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.hab)

m.season <- glmmTMB::glmmTMB(
  XFormula.season,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.season)


m.hab_season_add <- glmmTMB::glmmTMB(
  XFormula.hab_season_add,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.hab_season_add)


m.hab_season_int <- glmmTMB::glmmTMB(
  XFormula.hab_season_int,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.hab_season_int)


# Model selection

#Based on AICc
model_selection <- MuMIn::model.sel(
  m.null,
  m.host,
  m.hab,
  m.season,
  m.hab_season_add,
  m.hab_season_int,
  rank = "AICc"
)
model_selection
model_selection %>% filter(delta <= 2)

# #Choosing between the two best models
# avg <- MuMIn::model.avg(model_selection, subset = delta < 2)
# summary(avg)
# 
# anova(m.hab_season_int, m.season)

DHARMa::simulateResiduals(m.hab_season_int) %>%
  DHARMa::testResiduals()



summary(m.hab_season_int)
confint(m.hab_season_int)

#Looking at marginal means contrast for the complex model
em_inf <- emmeans::emmeans(m.hab_season_int, ~  category | season, type = "response")
emmeans::contrast(em_inf, "pairwise", adjust = "Tukey")


em_summary_inf <- as.data.frame(em_inf)


em_inf_plot <- em_summary_inf %>%
  mutate(
    category = forcats::fct_rev(forcats::fct_relevel(category, category_order))
  ) %>%
  ggplot(aes(x = rate, y = category, color = category)) +
  geom_point(size = 4.5) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = category, yend = category), linewidth = 1.5) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Estimated richness of infectious taxa",
    y = "Habitat category",
  ) +
  theme(
    axis.text = element_text(color = "gray20", size = 11),
    axis.title = element_text(color = "gray20", size = 12),
    strip.text = element_text(color = "gray20", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5)
  ) +
  facet_grid(season ~ .) +
  guides(color = "none")

em_inf_plot 

ggsave(filename = here::here("figures","richness_trend_plot.pdf"),
       plot = em_inf_plot,
       width = 2250, height = 1350, units = "px", dpi = 300, bg = NULL)


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


#Dataframes for GLM (no NA for explanatory variables)
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

# Visualization ----

# Calculate percentage
percent_at_least_2 <- mean(inf_alpha$rich >= 2) * 100
label_text <- paste0(round(percent_at_least_2, 1), "% with ≥2 infectious agents")

# Create histogram with annotation
richness_plot <- ggplot(inf_alpha, aes(x = rich)) +
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

ggstatsplot::gghistostats(data = inf_alpha,
                          x = rich,
                          binwidth = 1)


#Richness across category
inf_alpha %>% 
  ggplot(mapping = aes(x = category, y = rich, colour = category)) +
  geom_boxplot(outliers = FALSE , 
               size = 2.5) +
  scale_color_manual(values = category_palette) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggstatsplot::ggbetweenstats(data = inf_alpha,
                            x = category,
                            y = rich)

#Richness across season
inf_alpha %>% 
  ggplot(mapping = aes(x = season, y = rich, colour = season)) +
  geom_boxplot(outliers = FALSE , 
               size = 2.5) +
  scale_color_manual(values = season_palette) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggstatsplot::ggbetweenstats(data = inf_alpha,
                            x = season,
                            y = rich)


#Richness across category-season
inf_alpha %>% 
  ggplot(mapping = aes(x = category, y = rich, colour = category)) +
  geom_boxplot(outliers = FALSE , 
               size = 2.5) +
  scale_color_manual(values = category_palette) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  facet_wrap(.~season)


inf_alpha %>% 
  ggplot(mapping = aes(x = category, y = rich, colour = category)) +
  geom_violin(trim = TRUE, linewidth = 1.5) +
  geom_boxplot(outliers = FALSE ,
               width = 0.2, size = 1) +
  scale_color_manual(values = category_palette) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  facet_wrap(.~season)



# GLM : richness variation ----

## All category ----

#Check effectifs
df_alpha %>% 
  group_by(category, season) %>% 
  count()

#Global model specification
if (exists("m_richness")) rm(m_richness)
m_richness <- glmmTMB::glmmTMB(
  rich ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne) + (1 | year),
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m_richness)

#Model selection
model_selection <- MuMIn::dredge(m_richness, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Chosen model specification
if (exists("m_richness")) rm(m_richness)
m_richness <- glmmTMB::glmmTMB(
  rich ~ season + poids + sexe + (1 | numero_ligne) + (1 | year),
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)

#Chosen model validation
DHARMa::simulateResiduals(m_richness) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_richness)

performance::check_collinearity(m_richness)


#Chosen model look
summary(m_richness)
confint(m_richness)
gtsummary::tbl_regression(m_richness)
ggstats::ggcoef_model(m_richness)

#Marginal means contrast
# em_rich <- emmeans::emmeans(m_richness, ~  category | season  , type = "response")
# emmeans::contrast(em_rich, "pairwise", adjust = "Tukey")
# emmeans::contrast(em_rich, "pairwise", adjust = "Tukey") %>% confint()
# 
# em_summary_rich <- as.data.frame(em_rich)
# 
# 
# em_summary_rich %>%
#   mutate(
#     category = forcats::fct_rev(forcats::fct_relevel(category, category_order))
#   ) %>%
#   ggplot(aes(x = rate, y = category, color = category)) +
#   geom_point(size = 4.5) +
#   geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = category, yend = category), linewidth = 1.5) +
#   scale_color_manual(values = category_palette) +
#   theme_minimal(base_size = 16) +
#   labs(
#     x = "Estimated number of infectious",
#     y = "Habitat category",
#   ) +
#   theme(
#     axis.text = element_text(color = "gray20", size = 11),
#     axis.title = element_text(color = "gray20", size = 12),
#     strip.text = element_text(color = "gray20", size = 12),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5)
#   ) +
#   facet_grid(season ~ .) +
#   guides(color = "none")




## Regrouped plantations ----

#Check effectifs
df_alpha2 %>% 
  group_by(category, season) %>% 
  count()

#Global model specification
if (exists("m_richness2")) rm(m_richness2)
m_richness2 <- glmmTMB::glmmTMB(
  rich ~ category * season + poids + sexe + maturity + (1 | numero_ligne) + (1 | year),
  data = df_alpha2,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m_richness2)

#Model selection
model_selection <- MuMIn::dredge(m_richness2, rank = "AICc")
model_selection %>% filter(delta <= 2)


#Chosen model specification
if (exists("m_richness2")) rm(m_richness2)
m_richness2 <- glmmTMB::glmmTMB(
  rich ~ category + season + poids + sexe + (1 | numero_ligne) + (1 | year),
  data = df_alpha2,
  family = poisson(link = "log"),
  na.action = "na.fail",
)

#Chosen model validation
DHARMa::simulateResiduals(m_richness2) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_richness2)

performance::check_collinearity(m_richness2)


#Chosen model look
summary(m_richness2)
confint(m_richness2)
gtsummary::tbl_regression(m_richness2)
ggstats::ggcoef_model(m_richness2)




## Inference approach ----

#Models' fixed effects formula
XFormula.null = rich ~ 1 + (1 | numero_ligne) + (1 | year)
XFormula.host = rich ~ sexe + poids + maturity + (1 | numero_ligne) + (1 | year)
XFormula.hab_host = rich ~ category + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year)
XFormula.season_host = rich ~ season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year) 

XFormula.hab_season = rich ~  season * category + nb_apo + (1 | numero_ligne) + (1 | year)
XFormula.hab_season_host = rich ~  season * category + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year)


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

m.hab_host <- glmmTMB::glmmTMB(
  XFormula.hab_host,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.hab_host)

m.season_host <- glmmTMB::glmmTMB(
  XFormula.season_host,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.season_host)


m.hab_season <- glmmTMB::glmmTMB(
  XFormula.hab_season,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.hab_season)


m.hab_season_host <- glmmTMB::glmmTMB(
  XFormula.hab_season_host,
  data = df_alpha,
  family = poisson(link = "log"),
  na.action = "na.fail",
)
performance::check_collinearity(m.hab_season_host)


# Model selection

#Based on AICc
model_selection <- MuMIn::model.sel(
  m.null,
  m.host,
  m.hab_season,
  m.hab_host,
  m.season_host,
  m.hab_season_host,
  rank = "AICc"
)
model_selection
model_selection %>% filter(delta <= 2)

#Choosing between the two best models
avg <- MuMIn::model.avg(model_selection, subset = delta < 2)
summary(avg)

AIC(m.hab_season_host, m.season_host)
BIC(m.hab_season_host, m.season_host)

anova(m.season_host, m.hab_season_host)




#Looking at marginal means contrast for the complex model
em_inf <- emmeans::emmeans(m.hab_season_host, ~  category | season  , type = "response")
emmeans::contrast(em_inf, "pairwise", adjust = "Tukey")

em_summary_inf <- as.data.frame(em_inf)


em_summary_inf %>%
  mutate(
    category = forcats::fct_rev(forcats::fct_relevel(category, category_order))
  ) %>%
  ggplot(aes(x = rate, y = category, color = category)) +
  geom_point(size = 4.5) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = category, yend = category), linewidth = 1.5) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Estimated number of infectious",
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



### version new test ----


XFormula.null <- rich ~ 1 + (1 | numero_ligne) + (1 | year)
XFormula.host <- rich ~ sexe + poids + maturity + (1 | numero_ligne) + (1 | year)
XFormula.hab <- rich ~ category + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year)
XFormula.season <- rich ~ season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year)
XFormula.hab_season_add <- rich ~ category + season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year)
XFormula.hab_season_int <- rich ~ category * season + nb_apo + sexe + poids + maturity + (1 | numero_ligne) + (1 | year)


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

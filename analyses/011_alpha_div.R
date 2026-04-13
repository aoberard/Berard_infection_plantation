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

season_palette <- c("Spring" = "#ED6A5A",
                    "Autumn" = "#6B4B3E")


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


#Dataframe for GLM (no NA for explanatory variables)
df_alpha <- inf_alpha %>% 
  filter(!is.na(sexe)) 



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

#Check effectifs
df_alpha %>% 
  group_by(category, season) %>% 
  count()

#Global model specification
if (exists("m_richness")) rm(m_richness)
m_richness <- glmmTMB::glmmTMB(
  rich ~ category * season + poids + sexe + (1 | numero_ligne) + (1 | year),
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
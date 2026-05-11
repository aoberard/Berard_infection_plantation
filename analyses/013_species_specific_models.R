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

category_shapes <- c("pine_edge" = 17,
                     "hedgerows" = 16,
                     "broadleaved_forest" = 15)

category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
session_order <- c("Juin_2023", "Octobre_2023", "Juin_2024", "Septembre_2024", "Juin_2025")
season_order <- c("Spring", "Autumn")


# Import data ----

#Apodemus sylvaticus infectious data

inf <- data.table::fread(file = here::here("data", "derived-data", "asylvaticus_data", "asylvaticus_infectious.csv"))

names_inf <- readRDS(file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo.rds"))
names_inf10 <- readRDS(file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo10.rds"))


# Generate data ----

#Dataframe for individual infection probability GLM (no NA for explanatory variables)
df_probainf <- inf %>% 
  filter(!is.na(sexe)) 

#Dataframe for  GLM (data for glm regrouped by site-session)
df_nbinf <- df_probainf %>% 
  group_by(numero_ligne, field_session) %>% 
  summarise(across(all_of(names_inf10), ~sum(.x, na.rm = TRUE), 
                   .names = "sum_{.col}"),
            .groups = "drop_last") %>% 
  left_join(df_probainf %>% 
              select(numero_ligne, field_session, category, season, year) %>% 
              distinct())

#Choose reference levels
df_probainf$year <- relevel(as.factor(df_probainf$year), ref = "2024")
df_probainf$season <- relevel(as.factor(df_probainf$season), ref = "Autumn")
df_probainf$category <- relevel(as.factor(df_probainf$category), ref = "hedgerows")

df_nbinf$year <- relevel(as.factor(df_nbinf$year), ref = "2024")
df_nbinf$season <- relevel(as.factor(df_nbinf$season), ref = "Autumn")
df_nbinf$category <- relevel(as.factor(df_nbinf$category), ref = "hedgerows")

#Dataset with joined category
df_probainf2 <- df_probainf %>% 
  mutate(category =  
           case_when(
             category == "pine_edge" ~ "plantation",
             category == "hedgerows" ~ "plantation",
             .default = category
           ))

#Dataset with only broadleaved forests
df_probainf_3 <- df_probainf %>% 
  filter(category == "broadleaved_forest")



# GLM : individual infection probability ----

## Model : Bartonella taylorii  ----

#Check effectifs
df_probainf %>% 
  group_by(category, Bartonella_taylorii) %>% 
  count()

#Global model specification
if (exists("m_barttaylo")) rm(m_barttaylo)
m_barttaylo <- glmmTMB::glmmTMB(
  Bartonella_taylorii ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_barttaylo)
performance::check_convergence(m_barttaylo)


#Model selection
model_selection <- MuMIn::dredge(m_barttaylo, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_barttaylo")) rm(m_barttaylo)
m_barttaylo <- glmmTMB::glmmTMB(
  Bartonella_taylorii ~ season + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_barttaylo) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_barttaylo)
performance::check_convergence(m_barttaylo)
performance::check_collinearity(m_barttaylo)

#Selected model look
summary(m_barttaylo)
confint(m_barttaylo)



## Model : Bartonella birtlesii  ----

#Check effectifs
df_probainf %>% 
  group_by(category, Bartonella_birtlesii) %>% 
  count()

#Global model specification
if (exists("m_bartbirt")) rm(m_bartbirt)
m_bartbirt <- glmmTMB::glmmTMB(
  Bartonella_birtlesii ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_bartbirt)

#Model selection
model_selection <- MuMIn::dredge(m_bartbirt, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_bartbirt")) rm(m_bartbirt)
m_bartbirt <- glmmTMB::glmmTMB(
  Bartonella_birtlesii ~ poids  + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_bartbirt) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_bartbirt)
performance::check_convergence(m_bartbirt)
performance::check_collinearity(m_bartbirt)

#Selected model look
summary(m_bartbirt)
confint(m_bartbirt)



## Model : Eperythrozoon coccoides  ----

#Check effectifs
df_probainf %>% 
  group_by(category, Eperythrozoon_coccoides) %>% 
  count()


#Global model specification
if (exists("m_ecocco")) rm(m_ecocco)
m_ecocco <- glmmTMB::glmmTMB(
  Eperythrozoon_coccoides ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_ecocco)

#Model selection
model_selection <- MuMIn::dredge(m_ecocco, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_ecocco")) rm(m_ecocco)
m_ecocco <- glmmTMB::glmmTMB(
  Eperythrozoon_coccoides ~ season + poids + nb_apo + (1 | numero_ligne), 
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_ecocco) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_ecocco)
performance::check_convergence(m_ecocco)
performance::check_collinearity(m_ecocco)

#Selected model look
summary(m_ecocco)                                      
confint(m_ecocco)



## Model : Bartonella grahamii  ----

#Check effectifs
df_probainf %>% 
  group_by(category, Bartonella_grahamii) %>% 
  count()


#Global model specification
if (exists("m_bartgrah")) rm(m_bartgrah)
m_bartgrah <- glmmTMB::glmmTMB(
  Bartonella_grahamii ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_bartgrah)

#Model selection
model_selection <- MuMIn::dredge(m_bartgrah, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_bartgrah")) rm(m_bartgrah)
m_bartgrah <- glmmTMB::glmmTMB(
  Bartonella_grahamii ~ season + poids + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_bartgrah) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_bartgrah)

performance::check_collinearity(m_bartgrah)

#Selected model look
summary(m_bartgrah)
confint(m_bartgrah)



## Model : Neoehrlichia mikurensis  ----

#Check effectifs
df_probainf %>% 
  group_by(category, Neoehrlichia_mikurensis, season) %>% 
  count()

df_probainf %>% 
  group_by(category, season) %>% 
  summarise(prev = mean(Neoehrlichia_mikurensis), n = n())              ##### TOTAL SEPARATION FOR ONE LEVEL


#Global model specification
if (exists("m_neoehr")) rm(m_neoehr)
m_neoehr <- glmmTMB::glmmTMB(
  Neoehrlichia_mikurensis ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_neoehr)
performance::check_convergence(m_neoehr)

#Model selection
model_selection <- MuMIn::dredge(m_neoehr, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_neoehr")) rm(m_neoehr)
m_neoehr <- glmmTMB::glmmTMB(
  Neoehrlichia_mikurensis ~ category + season + poids + sexe + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_neoehr) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_neoehr)

performance::check_convergence(m_neoehr)

performance::check_collinearity(m_neoehr)

#Selected model look
summary(m_neoehr)
confint(m_neoehr)

#Selected model - Contrasts
em_neoehr <- emmeans::emmeans(m_neoehr, ~ category | season, type = "response")
emmeans::contrast(em_neoehr, "pairwise", adjust = "Tukey")
emmeans::contrast(em_neoehr, "pairwise", adjust = "Tukey") %>% confint()

em_summary_neoehr <- as.data.frame(em_neoehr)


em_summary_neoehr %>%
  mutate(
    category = forcats::fct_rev(forcats::fct_relevel(category, category_order))
  ) %>%
  ggplot(aes(x = prob, y = category, color = category)) +
  geom_point(size = 4.5) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = category, yend = category), linewidth = 1.5) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Estimated infection probability for Neoehrlichia mikurensis",
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



## Model : Haemobartonella muris ----

#Check effectifs
df_probainf %>% 
  group_by(category, Haemobartonella_muris) %>% 
  count()


#Global model specification
if (exists("m_haemomuris")) rm(m_haemomuris)
m_haemomuris <- glmmTMB::glmmTMB(
  Haemobartonella_muris ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_haemomuris)

#Model selection
model_selection <- MuMIn::dredge(m_haemomuris, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_haemomuris")) rm(m_haemomuris)
m_haemomuris <- glmmTMB::glmmTMB(
  Haemobartonella_muris ~ season + poids + (1 | numero_ligne) + (1 | year),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_haemomuris) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_haemomuris)

performance::check_collinearity(m_haemomuris)

#Selected model look
summary(m_haemomuris)
confint(m_haemomuris)


## Model : Haemo 2 ----

#Check effectifs
df_probainf %>% 
  group_by(category, Haemo_2) %>% 
  count()


#Global model specification
if (exists("m_haemo2")) rm(m_haemo2)
m_haemo2 <- glmmTMB::glmmTMB(
  Haemo_2 ~ category * season + poids + sexe + maturity + nb_apo + (1 | numero_ligne),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_haemo2)

#Model selection
model_selection <- MuMIn::dredge(m_haemo2, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_haemo2")) rm(m_haemo2)
m_haemo2 <- glmmTMB::glmmTMB(
  Haemo_2 ~ season + poids + (1 | numero_ligne),
  data = df_probainf,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_haemo2) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_haemo2)

performance::check_collinearity(m_haemo2)


#Selected model look
summary(m_haemo2)
confint(m_haemo2)








# GLM : Number of infected individuals ----


## Model : Bartonella taylorii  ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Bartonella_taylorii,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_barttaylo")) rm(m_nb_barttaylo)
m_nb_barttaylo <- glmmTMB::glmmTMB(
  sum_Bartonella_taylorii ~ category * season + (1 | numero_ligne) + (1 | year) , 
  data = df_nbinf,
  family =  poisson(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_barttaylo)

#Model selection
model_selection <- MuMIn::dredge(m_nb_barttaylo, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_barttaylo")) rm(m_nb_barttaylo)
m_nb_barttaylo <- glmmTMB::glmmTMB(
  sum_Bartonella_taylorii ~ category * season + (1 | numero_ligne) + (1 | year) , 
  data = df_nbinf,
  family =  poisson(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_barttaylo) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_barttaylo)

performance::check_collinearity(m_nb_barttaylo)

#Selected model look
summary(m_nb_barttaylo)
confint(m_nb_barttaylo)


#Selected model - Contrasts
em_nb_barttaylo <- emmeans::emmeans(m_nb_barttaylo, ~ category | season , type = "response")
emmeans::contrast(em_nb_barttaylo, "pairwise", adjust = "Tukey")
emmeans::contrast(em_nb_barttaylo, "pairwise", adjust = "Tukey") %>% confint()

em_summary_nb_barttaylo <- as.data.frame(em_nb_barttaylo)


em_summary_nb_barttaylo %>%
  mutate(
    category = forcats::fct_rev(forcats::fct_relevel(category, category_order))
  ) %>%
  ggplot(aes(x = rate, y = category, color = category)) +
  geom_point(size = 4.5) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = category, yend = category), linewidth = 1.5) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Estimated infection probability for Neoehrlichia mikurensis",
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





## Model : Bartonella_birtlesii  ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Bartonella_birtlesii,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_bartbirt")) rm(m_nb_bartbirt)
m_nb_bartbirt <- glmmTMB::glmmTMB(
  sum_Bartonella_birtlesii ~ category * season + (1 | numero_ligne) + (1 | year),
  data = df_nbinf,
  family = glmmTMB::nbinom2(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_bartbirt)

#Model selection
model_selection <- MuMIn::dredge(m_nb_bartbirt, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_bartbirt")) rm(m_nb_bartbirt)
m_nb_bartbirt <- glmmTMB::glmmTMB(
  sum_Bartonella_birtlesii ~ category + (1 | numero_ligne), 
  data = df_nbinf,
  family = glmmTMB::nbinom2(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_bartbirt) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_bartbirt)

performance::check_collinearity(m_nb_bartbirt)

#Selected model look
summary(m_nb_bartbirt)
confint(m_nb_bartbirt)



## Model : Eperythrozoon_coccoides  ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Eperythrozoon_coccoides,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_ecocco")) rm(m_nb_ecocco)
m_nb_ecocco <- glmmTMB::glmmTMB(
  sum_Eperythrozoon_coccoides ~ category * season + (1 | numero_ligne) + (1|year),
  data = df_nbinf,
  family = glmmTMB::nbinom2(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_ecocco)

#Model selection
model_selection <- MuMIn::dredge(m_nb_ecocco, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_ecocco")) rm(m_nb_ecocco)
m_nb_ecocco <- glmmTMB::glmmTMB(
  sum_Eperythrozoon_coccoides ~ category + (1 | numero_ligne) + (1|year),
  data = df_nbinf,
  family = glmmTMB::nbinom2(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_ecocco) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_ecocco)

performance::check_collinearity(m_nb_ecocco)

#Selected model look
summary(m_nb_ecocco)
confint(m_nb_ecocco)


## Model : Bartonella grahamii  ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Bartonella_grahamii,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_bartgrah")) rm(m_nb_bartgrah)
m_nb_bartgrah <- glmmTMB::glmmTMB(
  sum_Bartonella_grahamii ~ category * season + (1 | numero_ligne) + (1|year),
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_bartgrah)

#Model selection
model_selection <- MuMIn::dredge(m_nb_bartgrah, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_bartgrah")) rm(m_nb_bartgrah)
m_nb_bartgrah <- glmmTMB::glmmTMB(
  sum_Bartonella_grahamii ~ season + (1 | numero_ligne) + (1|year), 
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_bartgrah) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_bartgrah)

performance::check_collinearity(m_nb_bartgrah)

#Selected model look
summary(m_nb_bartgrah)
confint(m_nb_bartgrah)


## Model : Neoehrlichia mikurensis  ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Neoehrlichia_mikurensis,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_neoehr")) rm(m_nb_neoehr)
m_nb_neoehr <- glmmTMB::glmmTMB(
  sum_Neoehrlichia_mikurensis  ~ category * season + (1 | numero_ligne) + (1 | year),        #Interaction generate issues
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_neoehr)
performance::check_convergence(m_nb_neoehr)
performance::check_overdispersion(x = m_nb_neoehr)

#Model selection
model_selection <- MuMIn::dredge(m_nb_neoehr, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_neoehr")) rm(m_nb_neoehr)
m_nb_neoehr <- glmmTMB::glmmTMB(
  sum_Neoehrlichia_mikurensis ~  category * season + (1 | numero_ligne) + (1|year), 
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_neoehr) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_neoehr)
performance::check_collinearity(m_nb_neoehr)

#Selected model look
summary(m_nb_neoehr)
confint(m_nb_neoehr)

#Selected model - Contrasts
em_nb_neoehr <- emmeans::emmeans(m_nb_neoehr, ~ category | season , type = "response")
emmeans::contrast(em_nb_neoehr, "pairwise", adjust = "Tukey")
emmeans::contrast(em_nb_neoehr, "pairwise", adjust = "Tukey") %>% confint()

em_summary_nb_neoehr <- as.data.frame(em_nb_neoehr)


em_summary_nb_neoehr %>%
  mutate(
    category = forcats::fct_rev(forcats::fct_relevel(category, category_order))
  ) %>%
  ggplot(aes(x = rate, y = category, color = category)) +
  geom_point(size = 4.5) +
  geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, y = category, yend = category), linewidth = 1.5) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Estimated infection probability for Neoehrlichia mikurensis",
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


## Model : Haemobartonella muris ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Haemobartonella_muris,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_haemomuris")) rm(m_nb_haemomuris)
m_nb_haemomuris <- glmmTMB::glmmTMB(
  sum_Haemobartonella_muris ~ category * season + (1 | numero_ligne) + (1 | year),
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_haemomuris)

#Model selection
model_selection <- MuMIn::dredge(m_nb_haemomuris, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_haemomuris")) rm(m_nb_haemomuris)
m_nb_haemomuris <- glmmTMB::glmmTMB(
  sum_Haemobartonella_muris ~ season + (1 | numero_ligne),
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_haemomuris) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_haemomuris)

performance::check_collinearity(m_nb_haemomuris)

#Selected model look
summary(m_nb_haemomuris)
confint(m_nb_haemomuris)


## Model : Haemo 2 ----

#Check effectifs
ggstatsplot::gghistostats(data = df_nbinf,
                          x = sum_Haemo_2,
                          binwidth = 1)

#Global model specification
if (exists("m_nb_haemo2")) rm(m_nb_haemo2)
m_nb_haemo2 <- glmmTMB::glmmTMB(
  sum_Haemo_2 ~ category * season + (1 | numero_ligne) + (1| year),
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)
performance::check_collinearity(m_nb_haemo2)

#Model selection
model_selection <- MuMIn::dredge(m_nb_haemo2, rank = "AICc")
model_selection %>% filter(delta <= 2)

#Selected model specification
if (exists("m_nb_haemo2")) rm(m_nb_haemo2)
m_nb_haemo2 <- glmmTMB::glmmTMB(
  sum_Haemo_2 ~ season + (1 | numero_ligne),
  data = df_nbinf,
  family = poisson(link = "log"),
  na.action = "na.fail"
)

#Selected model validation
DHARMa::simulateResiduals(m_nb_haemo2) %>%
  DHARMa::testResiduals()

performance::check_overdispersion(x = m_nb_haemo2)

performance::check_collinearity(m_nb_haemo2)

#Selected model look
summary(m_nb_haemo2)
confint(m_nb_haemo2)










# GLM : Neoehrlichia in broadleaved forests ----

## Model : Neoehrlichia mikurensis  ----

#Check effectifs
df_probainf_3 %>% 
  group_by(category, Neoehrlichia_mikurensis, season) %>% 
  count()

df_probainf_3 %>% 
  group_by(category, season) %>% 
  summarise(prev = mean(Neoehrlichia_mikurensis), n = n())  

df_probainf_3$nb_clet %>% hist()
df_probainf_3$nb_apo %>% hist()


#Global model specification
if (exists("m_neoehr_brd")) rm(m_neoehr_brd)
m_neoehr_brd <- glmmTMB::glmmTMB(
  Neoehrlichia_mikurensis ~ season + poids  + sexe + nb_clet + (1 | numero_ligne) + (1| year),
  data = df_probainf_3,
  family = binomial(link = "logit"),
  na.action = "na.fail"
)
performance::check_collinearity(m_neoehr_brd)

#Model selection
model_selection <- MuMIn::dredge(m_neoehr_brd, rank = "AICc")
model_selection %>% filter(delta <= 2)



# GLM Bayesien : Neoehrlichian mikurensis ----

## Model : Neoehrlichia mikurensis  ----

#Check effectifs
df_probainf %>% 
  group_by(category, Neoehrlichia_mikurensis, season) %>% 
  count()

df_probainf %>% 
  group_by(category, season) %>% 
  summarise(prev = mean(Neoehrlichia_mikurensis), n = n())              ##### TOTAL SEPARATION FOR ONE LEVEL




# Prior faiblement informatif — classique pour régression logistique
# Normal(0, 2.5) sur les β est standard (Gelman et al.)
prior_weakly_informative <- c(
  brms::prior(normal(0, 2.5), class = "b"),
  brms::prior(normal(0, 2.5), class = "Intercept"),
  brms::prior(exponential(1),  class = "sd")   # effets aléatoires
)

# Prior pour le modèle nul — seulement intercept et effets aléatoires
prior_null <- c(
  brms::prior(normal(0, 2.5), class = "Intercept"),
  brms::prior(exponential(1),  class = "sd")
)

m_neoehr_bayes <- brms::brm(
  Neoehrlichia_mikurensis ~ category * season + poids + sexe + maturity + nb_apo +
    (1 | numero_ligne) + (1 | year),
  data    = df_probainf,
  family  = brms::bernoulli(link = "logit"),
  prior   = prior_weakly_informative,
  chains  = 4,
  iter    = 6000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(
    adapt_delta   = 0.99,   
    max_treedepth = 15      
  )
)

summary(m_neoehr_bayes)
plot(m_neoehr_bayes)           # convergence des chaînes
brms::pp_check(m_neoehr_bayes)


#Comparaison de modèles dérivés (pas tous)

# Ajouter le critère LOO à ton modèle global
m_neoehr_bayes <- brms::add_criterion(m_neoehr_bayes, "loo")

# Construire les modèles candidats 
# Modèle sans interaction
m_neoehr_additive <- brms::brm(
  Neoehrlichia_mikurensis ~ category + season + poids + sexe + maturity + nb_apo +
    (1 | numero_ligne) + (1 | year),
  data    = df_probainf,
  family  = brms::bernoulli(link = "logit"),
  prior   = prior_weakly_informative,
  chains  = 4, iter = 6000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
m_neoehr_additive <- brms::add_criterion(m_neoehr_additive, "loo")


# Modèle sans aucun effets hôtes (poids, sexe, maturity, nb_apo)
m_neoehr_env <- brms::brm(
  Neoehrlichia_mikurensis ~ category * season +
    (1 | numero_ligne) + (1 | year),
  data    = df_probainf,
  family  = brms::bernoulli(link = "logit"),
  prior   = prior_weakly_informative,
  chains  = 4, iter = 6000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
m_neoehr_env <- brms::add_criterion(m_neoehr_env, "loo")


# Modèle sans certains effets hôtes (maturity, nb_apo)
m_neoehr_env_host <- brms::brm(
  Neoehrlichia_mikurensis ~ category * season + poids + sexe +
    (1 | numero_ligne) + (1 | year),
  data    = df_probainf,
  family  = brms::bernoulli(link = "logit"),
  prior   = prior_weakly_informative,
  chains  = 4, iter = 6000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
m_neoehr_env_host <- brms::add_criterion(m_neoehr_env_host, "loo")



# Modèle nul (intercept only)
m_neoehr_null <- brms::brm(
  Neoehrlichia_mikurensis ~ 1 + (1 | numero_ligne) + (1 | year),
  data    = df_probainf,
  family  = brms::bernoulli(link = "logit"),
  prior   = prior_null,
  chains  = 4, iter = 6000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
m_neoehr_null <- brms::add_criterion(m_neoehr_null, "loo")


loo::loo_compare(
  m_neoehr_bayes,    
  m_neoehr_additive, 
  m_neoehr_env,      
  m_neoehr_env_host,
  m_neoehr_null      
)

#                    elpd_diff se_diff
# m_neoehr_env_host   0.0       0.0  
# m_neoehr_bayes     -1.2       1.8  
# m_neoehr_additive  -3.8       2.3  
# m_neoehr_env      -22.1       7.3  
# m_neoehr_null     -29.5       8.1 


# Probabilité que l'effet soit dans la direction estimée
result <- bayestestR::describe_posterior(m_neoehr_env_host, ci = 0.95, ci_method = "HDI")
print(result)

# Ou plus direct :
pd <- bayestestR::p_direction(m_neoehr_env_host)
print(pd)

rope_result <- bayestestR::rope(m_neoehr_env_host, range = c(-0.1, 0.1))
print(rope_result)



#Marginal effects

marginaleffects::avg_comparisons(m_neoehr_env_host, variables = "sexe")
marginaleffects::avg_comparisons(m_neoehr_env_host, variables = "category")


# 2. Probabilités prédites par catégorie × saison (le plus lisible pour un article)
pred <- ggeffects::ggpredict(m_neoehr_env_host, terms = c("category", "season"))
plot(pred)

# 3. Posterior predictive check (vérifier que le modèle reproduit tes données)
# Le plus informatif pour du binaire : proportion de 1 simulés vs observés
brms::pp_check(m_neoehr_env_host, type = "stat", stat = "mean", ndraws = 1000)
# → La moyenne observée (prévalence globale) doit tomber au centre de la distribution simulée

# Vérifier les valeurs extrêmes
brms::pp_check(m_neoehr_env_host, type = "stat_2d", stat = c("mean", "sd"), ndraws = 1000)
# → Le point observé (mean, sd) doit être dans le nuage simulé

# Groupé par catégorie — le plus révélateur pour ton design
brms::pp_check(m_neoehr_env_host, type = "stat_grouped", 
               stat = "mean", group = "category", ndraws = 1000)

# Groupé par saison
brms::pp_check(m_neoehr_env_host, type = "stat_grouped", 
               stat = "mean", group = "season", ndraws = 1000)






#Prevalence variable
#selon ?

#nbr d'individus infectés variable (on regarde que la prevalence a chaque fois c'est bon wsh, c'est pas ça le "danger")
# mais attention implique de standardiser data pour meme effort capture : + probleme c'est que c'est pas des vrais
#abondances qu'on regarde mais du relatif - pas forcement grave là aussi on peut dire c'est relatif non ?
#pour faire ça regarder les individus supprimés dans autre script quand retrait 4 eme nuit de piegeage
#peut-être autre truc a check



#modifier pour que plus de mention de infectious mais agents infectieux



#rien a voir mais en plus serait bien d'avoir figure avec tous les diff environnement entre categorie sur une fig
#comme ça permettrait d'une pierre de coup d'extraire valeur axe acp apres serait jsutifié ?



# Script parameters ----

#Remove saved elements
rm(list = ls())
gc()

# Library
library("dplyr")
library("ggplot2")
library("ggraph")

#Graphical parameters 
category_palette <- c("pine_edge"          = "#3f826d",
                      "hedgerows"          = "#f2d0a4",
                      "broadleaved_forest" = "#a4a8d1")


category_order <- c("pine_edge", "hedgerows", "broadleaved_forest")
session_order <- c("Juin_2023", "Octobre_2023", "Juin_2024", "Septembre_2024", "Juin_2025")
season_order <- c("Spring", "Autumn")


# Import data ----

#Samples data
dissected <- read.delim(file = here::here("data", "raw-data", "dissected_indiv.csv"), sep = ",")

#Site environmental data
envi <- data.table::fread(file = here::here("data", "raw-data", "lines_variability.csv"))

#Cleaned infectious agents data
infecti <- data.table::fread(file = here::here("data", "derived-data", "infectious_cleaning", "infectious_data.csv"))

#Apodemus sylvaticus number per site-session (over the three first trapping nights)
effectif <- data.table::fread(file = here::here("data", "raw-data", "apodemus_number.csv"))


# Generate data ----

## Join dataframes

#Infectious agents and individuals
infected <- infecti %>% 
  left_join(dissected)

#Adding environment
infected <- infected %>% 
  left_join(envi %>% 
              #keeping only variables that we will use later on
              select(num, category, herb_cover, shrub_cover, tree_cover, longitude_lambert93, latitude_lambert93,
                     PCA_microclim2025_axis1, PCA_microclim2025_axis2, pct500brd, gauss500brd, code_mission), 
            by = c("numero_ligne" = "num",
                         "code_mission" = "code_mission"))

# Add number of A. sylvaticus  caught for each line-session (in 3 first days of capture)
infected <- infected %>% 
  left_join(effectif, by = join_by("numero_ligne", "field_session")) %>% 
  rename(nb_apo = total)


#Ensure empty information are converted to NA
infected <- infected %>%
  mutate(across(where(is.character), ~na_if(., "")))


# Exploration ----

#Infectious agents name
names_infectious <- setdiff(colnames(infecti), y = "numero_centre")



## Generate essential subdataframes and vector ----

#Keep only Apodemus sylvaticus individuals
infected_apo <- infected %>%
  filter(taxon_mamm == "Apodemus sylvaticus") 

#Filter out empty infectious agents column for Apodemus sylvaticus
infected_apo <- infected_apo %>% 
  select(-names(which(colSums(infected_apo[, ..names_infectious]) == 0)))


#Vector of infectious agents name present in A. sylvaticus
names_infectious_apo <- intersect(names_infectious, colnames(infected_apo))
  
#Vector of infectious agents name present in A. sylvaticus with at least 10% global prevalence 
names_infectious_apo10 <- infected_apo %>%
  summarise(across(all_of(names_infectious_apo), ~ sum(. > 0) / n())) %>%
  unlist()
names_infectious_apo10 <- names(names_infectious_apo10[names_infectious_apo10 >= 0.10])
names_infectious_apo10


## Prevalence visualizations ----

### Prevalence variation across field session

# Heatmap of most prevalent infectiouss across field_sessions-category
df_pathos <- infected_apo %>% 
  tidyr::pivot_longer(
    cols = all_of(names_infectious_apo),
    names_to = "infectious",
    values_to = "infected"
  ) %>% 
  group_by(infectious, field_session) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    prevalence = cases/n,
    .groups = "drop"
  ) %>% 
  filter(infectious %in% names_infectious_apo10) 

matrix_pathos  <- df_pathos %>%
  select(infectious, field_session, prevalence) %>%
  tidyr::pivot_wider(names_from = field_session, values_from = prevalence) %>%
  tibble::column_to_rownames("infectious") %>%
  as.matrix()
matrix_pathos <- matrix_pathos[, session_order, drop = FALSE]

effectif <- df_pathos %>%
  distinct(field_session, n) %>%
  slice(match(session_order, field_session))

col_ha = ComplexHeatmap::HeatmapAnnotation(
  effectif = ComplexHeatmap::anno_barplot(
    effectif$n,
    bar_width = 1,
    baseline = 0,
    add_numbers = TRUE,
    gp = grid::gpar(fill = cividis::cividis(2, begin = 0.2, end = 0.9)),
    border = FALSE
  )
)

ComplexHeatmap::Heatmap(
  matrix_pathos,
  name = "Prévalence",
  top_annotation = col_ha,
  col = cividis_palette <- cividis::cividis(2),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE
)


  
#Line plot of most prevalent infectious across field_sessions-category
infected_apo %>% 
  tidyr::pivot_longer(
    cols = all_of(names_infectious_apo),
    names_to = "infectious",
    values_to = "infected"
  ) %>% 
  group_by(infectious, category, field_session) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    prevalence = cases/n,
    .groups = "drop"
  ) %>% 
  filter(infectious %in% names_infectious_apo10) %>% 
  rowwise() %>%
  mutate(
    ci = list(prevalence::propCI(cases, n, method = "wilson")),
    lower = ci[["lower"]],
    upper = ci[["upper"]]
  ) %>%
  ungroup() %>%
  select(-ci) %>% 
  mutate(category = forcats::fct_relevel(category, category_order),
         field_session = forcats::fct_relevel(field_session, session_order)) %>%
  ggplot(mapping = aes(x = field_session, y = prevalence, colour = category)) +
  geom_line(
    aes(group = category),
    linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~ infectious, scales = "free_y") +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )



#### Prevalence variation across season

# Heatmap of most prevalent infectious across season-category
df_pathos <- infected_apo %>% 
  tidyr::pivot_longer(
    cols = all_of(names_infectious_apo),
    names_to = "infectious",
    values_to = "infected"
  ) %>% 
  group_by(infectious, season) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    prevalence = cases/n,
    .groups = "drop"
  ) %>% 
  filter(infectious %in% names_infectious_apo10) 

matrix_pathos  <- df_pathos %>%
  select(infectious, season, prevalence) %>%
  tidyr::pivot_wider(names_from = season, values_from = prevalence) %>%
  tibble::column_to_rownames("infectious") %>%
  as.matrix()
matrix_pathos <- matrix_pathos[, season_order, drop = FALSE]

effectif <- df_pathos %>%
  distinct(season, n) %>%
  slice(match(season_order, season))

col_ha = ComplexHeatmap::HeatmapAnnotation(
  effectif = ComplexHeatmap::anno_barplot(
    effectif$n,
    bar_width = 1,
    baseline = 0,
    add_numbers = TRUE,
    gp = grid::gpar(fill = cividis::cividis(2, begin = 0.2, end = 0.9)),
    border = FALSE
  )
)

ComplexHeatmap::Heatmap(
  matrix_pathos,
  name = "Prévalence",
  top_annotation = col_ha,
  col = cividis_palette <- cividis::cividis(2),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  cluster_rows = TRUE
)




# Write files ----
dir.create(
  here::here("data", "derived-data", "asylvaticus_data"),
  recursive = TRUE,
  showWarnings = FALSE
)

data.table::fwrite(x = infected_apo, 
                   file = here::here("data", "derived-data", "asylvaticus_data", "asylvaticus_infectious.csv"),
                   na = "NA")

saveRDS(object = names_infectious_apo, file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo.rds"))
saveRDS(object = names_infectious_apo10, file = here::here("data", "derived-data", "asylvaticus_data", "names_infectious_apo10.rds"))




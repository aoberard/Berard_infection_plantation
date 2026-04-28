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

common_rodent_palette <- c(
  "Clethrionomys glareolus" = "#748CAB",  
  "Apodemus sylvaticus"      = "#D0A98F"
)

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


# Generate data ----

## Join dataframes

#Infectious agents and individuals
infected <- infecti %>% 
  left_join(dissected)

#Adding environment
infected <- infected %>% 
  left_join(envi, by = c("numero_ligne" = "num",
                         "code_mission" = "code_mission"))

#Compute number of A. sylvaticus and C. glareolus caught for each line-session (in 3 first days of capture)
infected <- infected %>% 
  left_join(
  infected %>% 
    filter(after_3days == FALSE) %>% 
    group_by(numero_ligne, field_session) %>% 
    summarise(
      nb_apo = sum(taxon_mamm == "Apodemus sylvaticus"),
      nb_clet = sum(taxon_mamm == "Clethrionomys glareolus"),
      .groups = "drop"
    ), by = c("numero_ligne", "field_session")) %>%
  mutate(
    nb_apo = coalesce(nb_apo, 0),
    nb_clet = coalesce(nb_clet, 0)
  )

#Ensure empty information are converted to NA
infected <- infected %>%
  mutate(across(where(is.character), ~na_if(., "")))


# Exploration ----

#Infectious agents name
names_infectious <- setdiff(colnames(infecti), y = "numero_centre")


## Across rodent species ----

#Number of individuals per rodent species
infected %>%
  group_by(taxon_mamm) %>%
  count()

#Number infected per taxa
infected %>%
  mutate(across(all_of(names_infectious), ~ replace(., . > 0, 1))) %>%
  group_by(taxon_mamm) %>%
  summarise(across(all_of(names_infectious), sum))

#List of infectious agents per species
list_infectious_per_species <- list()
for (i in unique(infected$taxon_mamm)) {
  list_infectious_per_species[[i]] <- infected %>%
    filter(taxon_mamm == i) %>%
    select(all_of(names_infectious)) %>%
    select(where(~ any(. > 0))) %>%
    colnames()
}
list_infectious_per_species


### Prevalence per species ----

#observed number of cases
prev_df <- infected %>%
  tidyr::pivot_longer(
    cols = all_of(names_infectious),
    names_to = "infectious",
    values_to = "infected"
  ) %>%
  group_by(taxon_mamm, infectious) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    .groups = "drop"
  )

#prevalence IC
prev_df <- prev_df %>%
  rowwise() %>%
  mutate(
    ci = list(prevalence::propCI(cases, n, method = "wilson")),
    prevalence = ci[, "p"],
    lower = ci[, "lower"],
    upper = ci[, "upper"]
  ) %>%
  ungroup() %>%
  select(-ci)

#remove rodent species without infectious agents
prev_df <- prev_df %>% 
  filter(!taxon_mamm %in% c("Mus musculus", "Micromys minutus"))

#Plot confidence interval of infectious agents across rodent species
prev_df %>%
  mutate(
    infectious = as.factor(reorder(infectious, cases)),
    taxon_mamm = as.factor(taxon_mamm)
  ) %>%
  ggplot(aes(x = prevalence, y = infectious, color = taxon_mamm)) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2
  ) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(width = 0.5),
    width = 0.3
  ) +
  scale_color_manual(values = common_rodent_palette) +
  theme_bw()


### Bipartite network ----

#Choose data to use ▲
bipartite_data <- infected %>% 
  filter(if_any(all_of(names_infectious), ~ . > 0))  

#Generate biadjacency matrix
bimatrix <- bipartite_data %>% 
  tibble::column_to_rownames(var = "numero_centre") %>% 
  select(all_of(names_infectious)) %>% 
  as.matrix

#Order matrix 
rodent_order <- hclust(dist(bimatrix))$order
infectious_order <- hclust(dist(t(bimatrix)))$order
bimatrix <- bimatrix[rodent_order, infectious_order]

#Generate graph from biadjacency matrix using igraph
bipartite_rodent <- bipartite_rodent <- igraph::graph_from_biadjacency_matrix(incidence = bimatrix ,
                                                               directed = FALSE)

#Wrap igraph object in tidygraph object
tidy_bipart <- tidygraph::as_tbl_graph(bipartite_rodent)

id_species <- bipartite_data %>%
  select(numero_centre, taxon_mamm) %>%
  distinct()

tidy_bipart <- tidy_bipart %>% 
  tidygraph::activate(nodes) %>% 
  left_join(id_species, by = c("name" = "numero_centre")) %>% 
  mutate(node_color = if_else(condition = type, true = "infectious", false = taxon_mamm))

#Plot tidygraph object using ggraph
ggraph(tidy_bipart, layout = "kk") +
  geom_edge_link(alpha = 0.7, color = "grey80") +
  geom_node_point(aes(color = node_color), size = 3) +
  geom_node_text(
    aes(label = ifelse(type, name, "")),
    repel = TRUE,
    size = 3
  ) +
  scale_color_manual(values = c(
    infectious = "tomato",
    common_rodent_palette
  )) +
  theme_void()


ggraph(tidy_bipart, layout = "bipartite") +
  geom_edge_link(alpha = 0.7, color = "grey80") +
  geom_node_point(aes(color = node_color), size = 3) +
  geom_node_text(
    aes(label = ifelse(type, name, "")),
    repel = TRUE,
    size = 3
  ) +
  coord_polar() + 
  scale_color_manual(values = c(
    infectious = "tomato",
    common_rodent_palette
  )) +
  theme_void()






tidy_bipart_ord <- tidy_bipart %>%
  tidygraph::activate(nodes) %>%
  dplyr::group_by(type) %>%   
  dplyr::mutate(
    angle = seq(0, 2*pi, length.out = dplyr::n()),
    radius = ifelse(type, 1, 2)  
  ) %>%
  dplyr::ungroup()

ggraph(tidy_bipart_ord,
       layout = "manual",
       x = radius * cos(angle),
       y = radius * sin(angle)) +
  geom_edge_link(alpha = 0.6, color = "grey80") +
  geom_node_point(aes(color = node_color, size = type)) +
  scale_size_manual(values = c(2, 5)) +
  geom_node_text(
    data = tidy_bipart_ord %>%
      tidygraph::activate(nodes) %>%
      as_tibble() %>%
      dplyr::mutate(
        x = radius * cos(angle),
        y = radius * sin(angle)
      ) %>%
      dplyr::filter(type),
    aes(x = x, y = y, label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_color_manual(values = c(
    infectious = "tomato",
    common_rodent_palette
  )) +
  
  theme_void()






#Quick version using bipartite
bimatrix2 <- bipartite_data %>%
  tidyr::pivot_longer(
    cols = all_of(names_infectious),
    names_to = "infectious",
    values_to = "infected"
  ) %>%
  group_by(taxon_mamm, infectious) %>%
  summarise(prevalence = mean(infected > 0), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = infectious,
    values_from = prevalence,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("taxon_mamm") %>%
  as.matrix()

bipartite::plotweb(
  bimatrix2,
  horizontal = TRUE, sorting = "dec"
)


## For A. sylvaticus ----

### Generate essential subdataframes and vector ----

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


### Prevalence visualizations ----

#### Prevalence variation across field session

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


#Bubble plot of most prevalent infectiouss across field_sessions-category
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
  filter(infectious %in% names_infectious_apo10,
         prevalence > 0) %>% 
  mutate(category = forcats::fct_relevel(category, category_order),
         field_session = forcats::fct_relevel(field_session, session_order)) %>%
  ggplot(mapping = aes(x = field_session, y = infectious, color = category)) +
  geom_point(aes(size = prevalence)) +
  scale_size_continuous(range = c(1, 10)) + 
  facet_grid(. ~ category) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
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


#Bubble plot of most prevalent infectiouss across season-category
infected_apo %>% 
  tidyr::pivot_longer(
    cols = all_of(names_infectious_apo),
    names_to = "infectious",
    values_to = "infected"
  ) %>% 
  group_by(infectious, category, season) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    prevalence = cases/n,
    .groups = "drop"
  ) %>% 
  filter(infectious %in% names_infectious_apo10,
         prevalence > 0) %>% 
  mutate(category = forcats::fct_relevel(category, category_order),
         season = forcats::fct_relevel(season, season_order)) %>%
  ggplot(mapping = aes(x = season, y = infectious, color = category)) +
  geom_point(aes(size = prevalence)) +
  scale_size_continuous(range = c(1, 10)) + 
  facet_grid(. ~ category) +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


#Confidence interval of prevalence for most prevalent infectious across season-category
infected_apo %>% 
  tidyr::pivot_longer(
    cols = all_of(names_infectious_apo),
    names_to = "infectious",
    values_to = "infected"
  ) %>% 
  group_by(infectious, category, season) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    prevalence = cases / n,
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
         season = forcats::fct_relevel(season, season_order)) %>%
  ggplot(aes(x = prevalence, y = season, colour = category)) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    position = position_dodge(width = 0.5),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  facet_wrap(~ infectious, scales = "free_y") +
  scale_color_manual(values = category_palette) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )


### Bipartite networks ----

#Choose data to use ▲
bipartite_data_apo <- infected_apo %>% 
  filter(if_any(all_of(names_infectious_apo10), ~ . > 0))  

#Generate biadjacency matrix
bimatrix_apo <- bipartite_data_apo %>% 
  tibble::column_to_rownames(var = "numero_centre") %>% 
  select(all_of(names_infectious_apo10)) %>% 
  as.matrix

#Generate graph from biadjacency matrix using igraph
bipartite_apo <- bipartite_apo <- igraph::graph_from_biadjacency_matrix(incidence = bimatrix_apo ,
                                                                              directed = FALSE)

#Wrap igraph object in tidygraph object
tidy_bipart_apo <- tidygraph::as_tbl_graph(bipartite_apo)

id_category_apo <- bipartite_data_apo %>%
  select(numero_centre, category) %>%
  distinct()

tidy_bipart_apo <- tidy_bipart_apo %>% 
  tidygraph::activate(nodes) %>% 
  left_join(id_category_apo, by = c("name" = "numero_centre")) %>% 
  mutate(node_color = if_else(condition = type, true = "infectious", false = category)) 


#Plot tidygraph object using ggraph
ggraph(tidy_bipart_apo, layout = "auto") +
  geom_edge_link(alpha = 0.7, color = "grey80") +
  geom_node_point(aes(color = node_color), size = 3) +
  geom_node_text(
    aes(label = ifelse(type, name, "")),
    repel = TRUE,
    size = 3
  ) +
  scale_color_manual(values = c(
    infectious = "tomato",
    category_palette
  )) +
  theme_void()


### Projection networks and else ----
projection <- igraph::bipartite_projection(bipartite_apo)

infect_netw <- projection$proj2

indiv_netw <- projection$proj1



#Wrap igraph object in tidygraph object
tidy_bipart_infect <- tidygraph::as_tbl_graph(infect_netw)

prev_apo <- infected_apo %>% 
  tidyr::pivot_longer(
    cols = all_of(names_infectious_apo),
    names_to = "infectious",
    values_to = "infected"
  ) %>% 
  group_by(infectious) %>%
  summarise(
    n = n(),
    cases = sum(infected),
    prevalence = cases/n,
    .groups = "drop"
  )

tidy_bipart_infect <- tidy_bipart_infect %>% 
  tidygraph::activate(nodes) %>% 
  left_join(prev_apo, by = c("name" = "infectious"))

#Plot tidygraph object using ggraph
ggraph(tidy_bipart_infect, layout = "auto") +
  geom_edge_arc(aes(width = weight), alpha = 0.4, color = "grey70", strength = 1) +
  geom_node_point(aes(size = prevalence), color = "tomato")+
  scale_size_continuous(range = c(2, 10)) +
  geom_node_text(aes(label = name),
                 repel = TRUE,
                 size = 4) +
  theme_void()


### SCA screening ----

#First check observed combination effectifs using UpsetR package
UpSetR::upset(infected_apo %>% 
                select(all_of(names_infectious_apo10)))

#Source the SCA function
source(file = here::here("R", "FctTestScreenENV.txt"))

#Generate matrix for use in screening
m_infectious <- infected_apo %>% 
  tibble::column_to_rownames(var = "numero_centre") %>% 
  select(all_of(names_infectious_apo10)) %>% 
  as.matrix()
  
#Indicator showing whether the number of hosts is sufficient for the number of infectious agents combinations considered
ifelse(nrow(m_infectious) < (2^ncol(m_infectious)),
       "Not enough hosts for the number of combinations", 
       "Sufficient number of hosts")

# Run screening
ResSCRENV <- FctTestScreenENV(m_infectious)

### --- Reconstruct combination labels --- ###

agent_names <- colnames(m_infectious)   # infectious agent names
nG <- length(agent_names)               # number of agents
nC <- length(ResSCRENV[[3]])            # number of combinations (2^nG)


combo_codes <- ResSCRENV[[3]]            # numeric codes (decimal representation of binary pattern)

#The original function builds DatXMul: each row is a combination, each column is an agent (1/0)
#We rebuild it here identically to decode each combination correctly
DatXMul <- array(rep(0, nC * nG), c(nC, nG))
for (k in 1:nG) {
  k1 <- 2^(nG - k)
  k2 <- 2^(k - 1)
  DatXMul[, k] <- rep(c(rep(1, k1), rep(0, k1)), k2)
}
#Now decode each row of DatXMul into a readable label
decode_combo_from_matrix <- function(row_bits, agent_names) {
  active <- agent_names[row_bits == 1]
  if (length(active) == 0) return("None")
  paste(active, collapse = " + ")
}
combo_labels <- apply(DatXMul, 1, decode_combo_from_matrix, agent_names = agent_names)

#Build data frame for ggplot
df_plot <- data.frame(
  index       = seq_len(nC),
  combo_code  = combo_codes,
  combo_label = combo_labels,
  yobs        = ResSCRENV[[6]],
  ysup        = ResSCRENV[[7]],
  yinf        = ResSCRENV[[8]],
  outside_env = ResSCRENV[[4]],
  pvalue      = ResSCRENV[[5]]
) %>%
  mutate(
    status = case_when(
      outside_env == 1 & yobs > ysup ~ "Over-represented",
      outside_env == 1 & yobs < yinf ~ "Under-represented",
      TRUE ~ "Within envelope"
    ),
    combo_label = forcats::fct_reorder(combo_label, index)
  )

# Print significant combinations
df_plot %>%
  filter(outside_env == 1) %>%
  select(combo_label, yobs, yinf, ysup, pvalue, status) %>%
  arrange(pvalue) %>%
  print()


colour_status <- c(
  "Over-represented"  = "#E74C3C",   # red
  "Under-represented" = "#3498DB",   # blue
  "Within envelope"   = "#7F8C8D"    # grey
)

p_sca <- ggplot(df_plot, aes(x = as.factor(combo_code))) +
  geom_linerange(
    aes(ymin = yinf, ymax = ysup),
    color = "grey70", linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = yinf, ymax = ysup),
    width = 0.4, color = "grey50", linewidth = 0.5
  ) +
  #Observed points coloured by status
  geom_point(
    aes(y = yobs, color = status, size = status)
  ) +
  
  #Labels only  combinations outside the envelope
  ggrepel::geom_label_repel(
    data = filter(df_plot, outside_env == 1),
    aes(
      y     = yobs,
      label = paste0(combo_label, "\np = ", round(pvalue, 4)),
      color = status
    ),
    size          = 2,
    nudge_y       = max(df_plot$ysup) * 0.12,
    box.padding   = 0.5,
    segment.color = "grey40",
    fontface      = "bold",
    show.legend   = FALSE
  ) +
  
  scale_color_manual(values = colour_status, name = "Status") +
  scale_size_manual(
    values = c("Over-represented" = 3, "Under-represented" = 3, "Within envelope" = 2),
    name   = "Status"
  ) +
  labs(
    title    = "Association Screening",
    x = "Combination code",
    y = "Count"
    ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    plot.title       = element_text(face = "bold"),
    plot.caption     = element_text(size = 8, color = "grey40"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )
p_sca












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

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

## Generate dissimilarity matrices

#Subset community matrix
m_inf <- inf %>%
  tibble::column_to_rownames(var = "numero_centre") %>%
  select(all_of(names_inf10)) %>%                              #Only on 10% infectious agent
  as.matrix()
dim(m_inf)

#Remove rows where all species are 0 (empty rows)
m_inf_clean <- m_inf[rowSums(m_inf) > 0, ]
dim(m_inf_clean)

#Generate dataframe with kept rows
kept_rows <- which(rowSums(m_inf) > 0)
metadata_clean <- inf[kept_rows, ]
rm(kept_rows)

#Generate Jaccard matrix
m_jacc <- vegan::vegdist(m_inf_clean, method = "jaccard", binary = TRUE)



# Statistical analyses ----

## Non-parametric permutation tests ----

### PERMANOVA ----
vegan::adonis2(formula = m_jacc ~ category ,
               data = metadata_clean,
               permutations = 9999,
               by = "margin") 



#Post-hoc tests




### BETADISPER ----
vegan::betadisper(d = m_jacc, 
                  group = metadata_clean$category) %>%
  vegan::permutest(permutations = 9999)












## db-RDA ----

### STANDARDIZE EXPLANATORY VARIABLES IF CONTINUOUS!!!

#Specify null and full models for db-RDA ordistep
mod_dbrda_0 <- vegan::dbrda(
  formula = metadata_clean %>% 
    select(all_of(names_inf10))
  ~ 1,
  data = metadata_clean,
  distance = "jaccard", 
  binary = TRUE
)

mod_dbrda_full <- vegan::dbrda(
  formula = metadata_clean %>% 
    select(all_of(names_inf10))
  ~ category * season + poids + sexe,
  data = metadata_clean,
  distance = "jaccard", 
  binary = TRUE
)

vegan::vif.cca(mod_dbrda_full)

#Execute ordiR2step (best db-RDA selection)
if (exists("dbrda")) rm(dbrda)
dbrda <- vegan::ordiR2step(object = mod_dbrda_0,
                          scope = formula(mod_dbrda_full),
                          direction = "both",
                          R2permutations = 999,
                          trace = FALSE)
#add species scores
vegan::sppscores(dbrda) <- metadata_clean %>% 
  select(all_of(names_inf10))

#Check model
dbrda$call
vegan::vif.cca(dbrda)
dbrda$anova
anova(dbrda, by = "margin")
summary(dbrda)

round(vegan::RsquareAdj(mod_dbrda_full)$adj.r.squared, 2)
round(vegan::RsquareAdj(dbrda)$adj.r.squared, 2)


#Variance partitioning
varpart <- vegan::varpart(Y = m_jacc, 
                      X = ~ season , ~ poids,
                      data = metadata_clean)
plot(varpart, bg=3:4)


# test posthoc capscale : multiconstrained package BiodiversityR









# Representations ----

## nMDS ----

#Jaccard nMDS on cleaned matrix
nmds_jacc <- vegan::metaMDS(
  comm = m_inf_clean,
  distance = "jaccard",
  k = 2,
  autotransform = FALSE,
  binary =
)

cat(paste0(
  "Final NMDS has a stress of ", round(nmds_jacc$stress, 3), " — ",
  dplyr::case_when(
    nmds_jacc$stress < 0.05 ~ "excellent goodness of fit (inferences very reliable)",
    nmds_jacc$stress < 0.1  ~ "good goodness of fit (inferences confident)",
    nmds_jacc$stress < 0.2  ~ "fair goodness of fit (some distances misleading)",
    TRUE                    ~ "poor goodness of fit (risks in interpretation)"
  ), "\n"
))

vegan::stressplot(nmds_jacc)
# Extract scores
nmds_jacc_point <- vegan::scores(nmds_jacc)$sites %>%
  as_tibble(rownames = "numero_centre")
nmds_jacc_variable <- vegan::scores(nmds_jacc)$species %>%
  as_tibble(rownames = "species")

nmds_jacc_plot <- nmds_jacc_point %>% 
  left_join(metadata_clean) %>%
  mutate(category = forcats::fct_relevel(category, category_order)) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = as.factor(category))) +
  geom_jitter(width = 0.05, height = 0.05, size = 1.5, alpha = 0.5) +  
  stat_ellipse(show.legend = FALSE, linewidth = 2) +
  geom_text(data = nmds_jacc_variable, 
            aes(x = NMDS1, y = NMDS2, label = species), 
            inherit.aes = FALSE,
            colour = "grey20") +
  scale_color_manual(values = category_palette)  +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black")
  ) +
  coord_equal()
nmds_jacc_plot


## db-RDA ----

#Extracting scores
scores_dbrda <- vegan::scores(dbrda,
                              scaling = "sites") # BEWARE of scaling chosen

#Community coordinates (individuals infracommunity)
sites <- scores_dbrda$sites %>% 
  as.data.frame()

sites$numero_centre <- metadata_clean$numero_centre
sites$category <- metadata_clean$category
sites$season <- metadata_clean$season

#Species coordinates (infectious agents)
species <- scores_dbrda$species %>% 
  as.data.frame()

species$taxa <- rownames(species)

#Variable coordinates
env <- scores_dbrda$biplot %>% 
  as.data.frame()
env$var <- rownames(env)

#Eigenvalues
eig <- vegan::eigenvals(dbrda)
var_explained <- eig / sum(eig)
xlab <- paste0("dbRDA1 (", round(100 * var_explained[1], 1), "%)")
ylab <- paste0("dbRDA2 (", round(100 * var_explained[2], 1), "%)")


#Plot it
ggplot(sites, aes(x = dbRDA1, y = dbRDA2)) +
  stat_ellipse(
    geom = "polygon",
    aes(color = season, fill = season),
    alpha = 0.1,
    linewidth = 1.5
  ) +  
  geom_point(aes(color = season), size = 3) +
  geom_segment(data = env,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(type = "closed", length = unit(0.3, "cm")),
               color = "black") +
  
  geom_text(data = env,
            aes(x = dbRDA1 * 1.35, y = dbRDA2 * 1.35,
                label = var),
            size = 4, fontface = "bold",
            color = "black") +
  
  # geom_segment(data = species,
  #              aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #              arrow = arrow(type = "closed", length = unit(0.3, "cm")),
  #              color = "red") +
  # geom_text(data = species,
  #           aes(x = dbRDA1 * 1.35, y = dbRDA2 * 1.35,
  #               label = taxa),
  #           size = 4, fontface = "bold",
  #           color = "red") +
  
  labs(x = xlab, y = ylab) +
  
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour="black", fill = NA, linewidth = 2),
    legend.background = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = season_palette)  +
  scale_fill_manual(values = season_palette)  +
  coord_equal()


## PCoA ----

#checker script lines variability à l'ancienne...







#QUICK write matrix distance for MAPi : 
m_jacc %>%
  as.matrix() %>%
  write.csv("m_jacc.csv", row.names = TRUE)









# TEST zone ----



## Test glm ----

#Import dataframe of pairwise distance for loiselle, spatial and time
loiselle_pairwise <- data.table::fread(file = here::here("data", "raw-data", "pairwise_loiselle.csv"))

#Rename loiselle value
loiselle_pairwise <- loiselle_pairwise %>%
  rename(loiselle = value)

#Remove pop id in ind id for loiselle
loiselle_pairwise <- loiselle_pairwise %>%
  mutate(
    ind1 = gsub("^.*_", "", ind1),
    ind2 = gsub("^.*_", "", ind2)
  )


### Co-presence model ----

#Infection matrix
m_inf
#or
m_inf_clean
rownames(m_inf_clean) #Only ones infected by at least one thing

#Generate co_presence matrix
co_presence <- do.call(rbind, lapply(names_inf10, function(p){

  # vecteur pour un pathogène
  v <- m_inf_clean[, p]                                #Matrix with only

  # matrice de co-présence
  copresence_matrix <- outer(v, v, "*")

  # melt en data.frame
  df <- as.data.frame(as.table(copresence_matrix))
  colnames(df) <- c("ind1", "ind2", "shared")

  # ajouter nom du pathogène
  df$infectious_ag <- p

  return(df)
}))

#Merge with Loiselle
copresence_relatedness <- merge(x = loiselle_pairwise, y = co_presence, by = c("ind1", "ind2"))


#Standardize the order of elements individuals in pairs (ind1 is always the "smaller" one from the pair)
copresence_relatedness_strd <- copresence_relatedness %>%
  mutate(
    ind_min = pmin(ind1, ind2),
    ind_max = pmax(ind1, ind2)
  ) %>%
  mutate(
    ind1 = ind_min,
    ind2 = ind_max
  ) %>%
  select(-ind_min, -ind_max)

#Remove self_relatedness
copresence_relatedness_strd <- copresence_relatedness_strd %>%
  filter(ind1 != ind2)

#Remove repeated pairs
copresence_relatedness_strd <- copresence_relatedness_strd %>%
  distinct()

#Small checks
copresence_relatedness_strd %>% count(shared, infectious_ag)


#Doing linear model (which is likely wrong considering non-independance)
#keeping only small distance things
copresence_relatedness_strd_test <- copresence_relatedness_strd %>%
  filter(geo < 300) %>%
  mutate(dyad =  paste(pmin(ind1, ind2), pmax(ind1, ind2), sep = "_"))

#modelling and doing model selection
bad_model <- glmmTMB::glmmTMB(
  shared ~ (loiselle + infectious_ag + geo)^3 +
    (1 | ind1) + (1 | ind2) + (1 | dyad),
  family = binomial,
  data = copresence_relatedness_strd_test,
  na.action = "na.fail"
)

performance::check_collinearity(bad_model)

#Model selection
model_selection <- MuMIn::dredge(bad_model, rank = "AICc")
model_selection %>% filter(delta <= 2)



### with jaccard ----

#Melt jaccard
m_jacc_dt <- m_jacc %>%
  as.matrix() %>%
  data.table::as.data.table(keep.rownames = TRUE)

m_jacc_melted <- data.table::melt(data = m_jacc_dt, id.vars = "rn")

#Prepare jaccard for join
colnames(m_jacc_melted) <- c("ind1", "ind2", "jacc")
m_jacc_melted

#Join to loiselle pairwise data
jacc_relatedness <- merge(x = loiselle_pairwise, y = m_jacc_melted, by = c("ind1", "ind2"))


#Visualization
jacc_relatedness %>% 
  filter(field_session1 == field_session2) %>% 
  filter(geo <100) %>% 
  ggplot(aes(x = jacc, y = loiselle)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm")




## Test mantel ---- 

#Import symmetric loiselle data
loiselle_matrix <- data.table::fread(file = here::here("data", "raw-data", "loiselle_matrix_sym.csv"))

#Make it a matrix object
loiselle_matrix <- loiselle_matrix %>% 
  tibble::column_to_rownames("V1")

loiselle_matrix <- loiselle_matrix %>% 
  as.matrix()

#Modify column and row names for loiselle
colnames(loiselle_matrix) <- sub("^[^_]+_", "", colnames(loiselle_matrix))
rownames(loiselle_matrix) <- sub("^[^_]+_", "", colnames(loiselle_matrix))

#Make jaccard distance object a matrix too
m_jacc_mat <- as.matrix(m_jacc)


#Generate a spatial distance matrix
#Generate data with lambert coordinates for each individuals
samp0g <- sf::st_as_sf(inf,
                       coords = c("longitude_piege", "latitude_piege"),
                       crs = 4326)
samp1g <- sf::st_transform(samp0g, crs = "EPSG:2154")
samp2g <- cbind(samp1g, sf::st_coordinates(samp1g))


#Then generate the spatial distance matrix
dgeo_indiv <- as.matrix(sf::st_distance(samp2g))
dgeo_indiv <- units::drop_units(dgeo_indiv)
colnames(dgeo_indiv) <- samp2g$numero_centre 
rownames(dgeo_indiv) <- samp2g$numero_centre 




#Keep only shared names among matrix
common_ids <- intersect(rownames(loiselle_matrix), rownames(m_jacc_mat)) %>% 
  sort()


loiselle_sub <- loiselle_matrix[common_ids, common_ids]
m_jacc_sub   <- m_jacc_mat[common_ids, common_ids]
dgeo_indiv_sub <- dgeo_indiv[common_ids, common_ids]

all(rownames(loiselle_sub) == rownames(m_jacc_sub))


#Mantel m_jacc - loiselle 
vegan::mantel(m_jacc_sub, loiselle_sub, 
         method = "spearman", 
         permutations = 999)

#Mantel m_jacc - dgeo
vegan::mantel(m_jacc_sub, dgeo_indiv_sub, 
              method = "spearman", 
              permutations = 999)

#Mantel partial test
vegan::mantel.partial(m_jacc_sub, loiselle_sub, dgeo_indiv_sub)


#Mantel correlogram
mc <- vegan::mantel.correlog(D.eco = loiselle_sub, 
                             D.geo = dgeo_indiv_sub,
                             r.type = "spearman")

plot(mc)










## Mantel par classes de distance géographique ----

# Classes de distance en mètres
breaks_m <- c(0, 100, 200, 1000, 2000, Inf)
labels_m <- c("0-100m", "100-200m", "200-1km", "1-2km", ">2km")

# Cor de Spearman par classe — même session uniquement
cor_results <- jacc_relatedness %>%
  #filter(field_session1 == field_session2) %>%                                                          #ONLY WITHIN SESSION HERE
  mutate(dist_class = cut(geo, breaks = breaks_m, labels = labels_m, include.lowest = TRUE)) %>%
  group_by(dist_class) %>%
  filter(n() >= 5) %>%
  summarise(
    r_spearman = cor(jacc, loiselle, method = "spearman", use = "complete.obs"),
    n_pairs    = n(),
    .groups    = "drop"
  ) %>%
  mutate(dist_class = factor(dist_class, levels = labels_m))

print(cor_results)

# Plot
ggplot(cor_results, aes(x = dist_class, y = r_spearman)) +
  geom_col(fill = "#2196F3", width = 0.6) +
  geom_text(aes(label = paste0("n=", n_pairs),
                vjust = ifelse(r_spearman >= 0, -0.4, 1.2)),
            size = 3.5, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x     = "Classe de distance géographique",
    y     = "r de Spearman",
    title = "Corrélation Jaccard ~ Kinship (Loiselle)\npar classe de distance — même session"
  ) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))


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
m_jacc <- vegan::vegdist(m_inf_clean, 
                         method = "jaccard", 
                         binary = TRUE)



# Statistical analyses ----

## Non-parametric permutation tests ----

### PERMANOVA ----
vegan::adonis2(formula = m_jacc ~ category ,
               data = metadata_clean,
               permutations = 999,
               by = "margin") 

vegan::adonis2(formula = m_jacc ~ season ,
               data = metadata_clean,
               permutations = 999,
               by = "margin") 

vegan::adonis2(formula = m_jacc ~ poids  ,
               data = metadata_clean,
               permutations = 999,
               by = "margin") 

vegan::adonis2(formula = m_jacc ~ sexe  ,
               data = metadata_clean,
               permutations = 999,
               by = "margin") 

vegan::adonis2(formula = m_jacc ~ nb_apo ,
               data = metadata_clean,
               permutations = 999,
               by = "margin") 


#Post-hoc tests


### BETADISPER ----
vegan::betadisper(d = m_jacc, 
                  group = metadata_clean$category) %>%
  vegan::permutest(permutations = 999)

vegan::betadisper(d = m_jacc, 
                  group = metadata_clean$season) %>%
  vegan::permutest(permutations = 999)



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
  ~ category * season + poids + sexe + nb_apo,
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


# Posthoc dbrda : multiconstrained package BiodiversityR



# Representations ----

## nMDS ----

#Jaccard nMDS on cleaned matrix
nmds_jacc <- vegan::metaMDS(
  comm = m_inf_clean,
  distance = "jaccard",
  k = 2,
  autotransform = FALSE,
  binary = TRUE
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
nmdspoint <- vegan::scores(nmds_jacc)$sites %>%
  as_tibble(rownames = "numero_centre")
nmdsvariable <- vegan::scores(nmds_jacc)$species %>%
  as_tibble(rownames = "species")


#filter species by contribution to dispersion 
fit <- vegan::envfit(nmds_jacc, m_inf_clean, permutations = 999)
sig_sp <- names(which(fit$vectors$pvals < 0.05))
nmdsvariable <- nmdsvariable %>%
  filter(species %in% sig_sp)

#Generate point data containing metadata
nmdspoint_plot <- nmdspoint %>%
  left_join(metadata_clean) %>%
  mutate(category = factor(category))

#Plot nMDS
nmds_jacc_plot <- ggplot() +
  stat_ellipse(
    data = nmdspoint_plot,
    aes(x = NMDS1, y = NMDS2, color = category),
    linewidth = 1.5, show.legend = FALSE) +
  
  geom_jitter(data = nmdspoint_plot,
              aes(x = NMDS1, y = NMDS2, color = category, shape = category),
              width = 0.05, 
              height = 0.05, 
              size = 1.5, 
              alpha = 0.5) +  
  
  ggrepel::geom_text_repel(
    data = nmdsvariable,
    aes(x = NMDS1, y = NMDS2, label = species),
    color = "grey40", size = 2.5, fontface = "italic",
    box.padding = 0.3, max.overlaps = 20
  ) +
  
  scale_color_manual(values = category_palette) +
  scale_shape_manual(values = category_shapes) +
  labs(
    x = "NMDS1", y = "NMDS2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    panel.background = element_rect(fill = NA),
    plot.background = element_rect(fill = NA, color = NA),
    axis.title = element_text(size = 13, face = "bold", color = "grey30"),
    axis.text = element_text(size = 11, color = "grey30"),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.background = element_blank(),
    aspect.ratio = 1,
    panel.grid = element_blank()
  ) +
  coord_equal()
nmds_jacc_plot






## db-RDA ----

# SCALING 1 : "sites" 
# Interpretation: distances between sites preserved : focus on patterns between communities
# What is represented:
# - Site points: relative position of infracommunities within the constrained space
# - Environmental arrows (biplot): direction and strength of environmental variables
# - Ellipses: dispersion by group (season, category...)

#Extracting scores at sites scale
scores_sites <- vegan::scores(dbrda, scaling = "sites")

#Community coordinates (individuals infracommunity)
sites_s1 <- scores_sites$sites %>%
  as.data.frame() %>%
  mutate(
    numero_centre = metadata_clean$numero_centre,
    category = metadata_clean$category,
    season = metadata_clean$season
  )

#Environnement variable coordinates
env_s1 <- scores_sites$biplot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var")

#Eigenvalues (%explained by axes)
eig <- vegan::eigenvals(dbrda)
var_constr <- eig / sum(abs(eig)) 
xlab_s1    <- paste0("dbRDA1 (", round(100 * var_constr[1], 1), "%)")
ylab_s1    <- paste0("dbRDA2 (", round(100 * var_constr[2], 1), "%)")

#R2ajusted
r2adj <- vegan::RsquareAdj(dbrda)$adj.r.squared

plot_scaling_sites <- ggplot(sites_s1, aes(x = dbRDA1, y = dbRDA2)) +
  stat_ellipse(
    aes(color = season),
    linewidth = 1.5, show.legend = FALSE
  ) +
  geom_point(aes(color = season, shape = category), size = 3) +
  geom_segment(
    data = env_s1,
    aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    color = "black", linewidth = 0.8
  ) +
  geom_label(
    data = env_s1,
    aes(x = dbRDA1 * 1.25, y = dbRDA2 * 1.25, label = var),
    size = 3.5, fontface = "bold", color = "black",
    fill = alpha("white", 0.4), label.size = 0
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.5,
           label = paste0("R²adj = ", round(r2adj, 3)),
           size = 4, fontface = "italic") +
  labs(
    title = "db-RDA (Scaling sites)",
    x = xlab_s1, y = ylab_s1
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour="black", fill = NA, linewidth = 2),
    legend.background = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = season_palette)  +
  scale_fill_manual(values = season_palette)  +
  coord_equal()


# SCALING 2 : "species"
# Interpretation: distances between species preserved : focus on species-environment associations
# What is represented:
# - Species arrows: infectious agents and their association with constrained axes
# - Environmental arrows (biplot): direction of environmental variables
# - Angle species-variable ≈ correlation between the species and the variable
# - Site points: approximate positions only (not directly interpretable)

#Extracting scores at species scale
scores_sp <- vegan::scores(dbrda, scaling = "species")

#Community coordinates (individuals infracommunity) — approximate positions
sites_s2 <- scores_sp$sites %>%
  as.data.frame() %>%
  mutate(
    numero_centre = metadata_clean$numero_centre,
    category = metadata_clean$category,
    season = metadata_clean$season
  )

#Species coordinates (infectious agents)
species_s2 <- scores_sp$species %>%
  as.data.frame() %>%
  tibble::rownames_to_column("taxa")

#Environment variable coordinates
env_s2 <- scores_sp$biplot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var")

#Eigenvalues (%explained by axes) — same as scaling 1
xlab_s2 <- paste0("dbRDA1 (", round(100 * var_constr[1], 1), "%)")
ylab_s2 <- paste0("dbRDA2 (", round(100 * var_constr[2], 1), "%)")

#Rescaling factor to make sites and species readable on the same plot
scale_sp <- max(abs(sites_s2[, c("dbRDA1", "dbRDA2")])) / max(abs(species_s2[, c("dbRDA1", "dbRDA2")])) * 0.7

plot_scaling_species <- ggplot() +
  geom_point(
    data = sites_s2,
    aes(x = dbRDA1, y = dbRDA2, color = season, shape = category),
    size = 3, alpha = 0.3
  ) +
  geom_segment(
    data = species_s2,
    aes(x = 0, y = 0,
        xend = dbRDA1 * scale_sp,
        yend = dbRDA2 * scale_sp),
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    color = "#c0392b", linewidth = 0.8
  ) +
  geom_label(
    data = species_s2,
    aes(x = dbRDA1 * scale_sp * 1.25,
        y = dbRDA2 * scale_sp * 1.25,
        label = taxa),
    size = 3.5, fontface = "bold", color = "#c0392b",
    fill = alpha("white", 0.4), label.size = 0
  ) +
  geom_segment(
    data = env_s2,
    aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
    arrow = arrow(type = "closed", length = unit(0.25, "cm")),
    color = "black", linewidth = 0.8
  ) +
  geom_label(
    data = env_s2,
    aes(x = dbRDA1 * 1.25, y = dbRDA2 * 1.25, label = var),
    size = 3.5, fontface = "bold", color = "black",
    fill = alpha("white", 0.6), label.size = 0
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.5,
           label = paste0("R²adj = ", round(r2adj, 3)),
           size = 4, fontface = "italic") +
  labs(
    title = "db-RDA (Scaling species)",
    x = xlab_s2, y = ylab_s2
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 2),
    legend.background = element_blank(),
    legend.position   = "bottom"
  ) +
  scale_color_manual(values = season_palette) +
  coord_equal()

# Combined plot
library(patchwork)
plot_scaling_sites + plot_scaling_species +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



## PCoA ----

#Generate PCoA
PCOA_jacc_inf<- vegan::wcmdscale(m_jacc, eig = TRUE)

#Check proportion of variance explained by each axis
var_explained <- round(100 * PCOA_jacc_inf$eig / sum(PCOA_jacc_inf$eig), 2)
print(var_explained)

tibble(
  Axis = seq_along(var_explained),
  Variance = var_explained
) %>%
  ggplot(aes(x = Axis, y = Variance)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 3) +
  labs(
    x = "PCoA Axis",
    y = "Variance Explained (%)",
    title = "Scree Plot - PCoA"
  ) +
  theme_minimal(base_size = 12)

#Extract PCoA scores for plotting
pcoapoint <- PCOA_jacc_inf$points %>%
  as_tibble(rownames = "numero_centre") %>% 
  left_join(metadata_clean) %>% 
  rename(PCoA1 = Dim1, PCoA2 = Dim2)  %>%
  mutate(category = factor(category))

#envfit arrows (correlation between axis and their abundances)
fit <- vegan::envfit(
  pcoapoint[, c("PCoA1", "PCoA2")],
  m_inf_clean,
  permutations = 999
)

spec_df <- as.data.frame(fit$vectors$arrows) %>%
  tibble::rownames_to_column(var = "species") %>%
  mutate(
    r    = fit$vectors$r,
    pval = fit$vectors$pvals,
    p_adj = p.adjust(pval, method = "fdr")
  )

#Significativity filter
sig_spec <- spec_df %>%
  filter(p_adj < 0.05) %>%
  mutate(
    PCoA1 = PCoA1 * r,
    PCoA2 = PCoA2 * r
  )

#Plot PCoA
pcoa_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5, linetype = "dashed") +
  stat_ellipse(
    data = pcoapoint,
    aes(x = PCoA1, y = PCoA2, color = category),
    linewidth = 1.5, show.legend = FALSE) +
  geom_point(
    data = pcoapoint,
    aes(x = PCoA1, y = PCoA2, color = category, shape = category),
    size = 3) +
  geom_segment(
    data = sig_spec,
    aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = "grey40", linewidth = 0.6)  +
  ggrepel::geom_text_repel(
    data = sig_spec,
    aes(x = PCoA1, y = PCoA2, label = species),
    color = "grey40", size = 3, fontface = "italic",
    box.padding = 0.3, max.overlaps = 20) +
  scale_color_manual(values = category_palette) +
  scale_shape_manual(values = category_shapes) +
  labs(
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    panel.background = element_rect(fill = NA),
    plot.background = element_rect(fill = NA, color = NA),
    axis.title = element_text(size = 13, face = "bold", color = "grey30"),
    axis.text = element_text(size = 11, color = "grey30"),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.background = element_blank(),
    aspect.ratio = 1,
    panel.grid = element_blank()) +
  coord_equal()

pcoa_plot






# TEST zone ----

#QUICK write matrix distance for MAPi : 
m_jacc %>%
  as.matrix() %>%
  write.csv("m_jacc.csv", row.names = TRUE)




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
# 
# #Infection matrix
# m_inf
# #or
# m_inf_clean
# rownames(m_inf_clean) #Only ones infected by at least one thing
# 
# #Generate co_presence matrix
# co_presence <- do.call(rbind, lapply(names_inf10, function(p){
# 
#   # vecteur pour un pathogène
#   v <- m_inf_clean[, p]                                #Matrix with only
# 
#   # matrice de co-présence
#   copresence_matrix <- outer(v, v, "*")
# 
#   # melt en data.frame
#   df <- as.data.frame(as.table(copresence_matrix))
#   colnames(df) <- c("ind1", "ind2", "shared")
# 
#   # ajouter nom du pathogène
#   df$infectious_ag <- p
# 
#   return(df)
# }))
# 
# #Merge with Loiselle
# copresence_relatedness <- merge(x = loiselle_pairwise, y = co_presence, by = c("ind1", "ind2"))
# 
# 
# #Standardize the order of elements individuals in pairs (ind1 is always the "smaller" one from the pair)
# copresence_relatedness_strd <- copresence_relatedness %>%
#   mutate(
#     ind_min = pmin(ind1, ind2),
#     ind_max = pmax(ind1, ind2)
#   ) %>%
#   mutate(
#     ind1 = ind_min,
#     ind2 = ind_max
#   ) %>%
#   select(-ind_min, -ind_max)
# 
# #Remove self_relatedness
# copresence_relatedness_strd <- copresence_relatedness_strd %>%
#   filter(ind1 != ind2)
# 
# #Remove repeated pairs
# copresence_relatedness_strd <- copresence_relatedness_strd %>%
#   distinct()
# 
# #Small checks
# copresence_relatedness_strd %>% count(shared, infectious_ag)
# 
# 
# #Doing linear model (which is likely wrong considering non-independance)
# #keeping only small distance things
# copresence_relatedness_strd_test <- copresence_relatedness_strd %>%
#   filter(geo < 300) %>%
#   mutate(dyad =  paste(pmin(ind1, ind2), pmax(ind1, ind2), sep = "_"))
# 
# #modelling and doing model selection
# bad_model <- glmmTMB::glmmTMB(
#   shared ~ (loiselle + infectious_ag + geo)^3 +
#     (1 | ind1) + (1 | ind2) + (1 | dyad),
#   family = binomial,
#   data = copresence_relatedness_strd_test,
#   na.action = "na.fail"
# )
# 
# performance::check_collinearity(bad_model)
# 
# #Model selection
# model_selection <- MuMIn::dredge(bad_model, rank = "AICc")
# model_selection %>% filter(delta <= 2)
# 


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


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







#Infection matrix
m_inf
#or
m_inf_clean
rownames(m_inf_clean) #Only ones infected by at least one thing

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
  distinct(ind1, ind2, infectious_ag)






copresence_relatedness %>% count(shared, infectious_ag)










# test modelisation proximité selon proximité uwu 

#Tronche chelou distribution de metrique jaccard :
#soit modeliser variation apparentement en fonction de jacc (en gros faire l'inverse)
#
#soit modeliser part de shared pathogens avec poisson ? soit avec binomial (succes / tentative)
#soit pour l'ensemble des pathogènes

#soit pour l'ensemble des pathogènes x paires avec pathogen en effet aléatoire (OU FIXE?)
#glmer(shared ~ relatedness + distance + (1|pathogen), family = binomial)

#juste faut faire gaffe psk du coup faut bien controler pas trop de non indépendance etc là pfiou ça devient complicado

#est ce que quand meme garder que ceux infecté par au moins un truc ?


























# #Melt jaccard 
# m_jacc_dt <- m_jacc %>%
#   as.matrix() %>% 
#   data.table::as.data.table(keep.rownames = TRUE)
# 
# m_jacc_melted <- data.table::melt(data = m_jacc_dt, id.vars = "rn")
# 
# #Prepare jaccard for join
# colnames(m_jacc_melted) <- c("ind1", "ind2", "jacc")
# m_jacc_melted
# 
# #Join to loiselle pairwise data
# jacc_relatedness <- merge(x = loiselle_pairwise, y = m_jacc_melted, by = c("ind1", "ind2"))
# 
# 
# 
# 
# jacc_relatedness %>%
#   count(ind1) %>% View(  )
# #y'a des ind1 tres peu representé
# jacc_relatedness %>%
#   count(ind2) %>% View(  )
# 
# jacc_relatedness$jacc %>% hist()

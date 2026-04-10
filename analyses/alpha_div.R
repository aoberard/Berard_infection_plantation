









adapter ce script :


# Calculate percentage
percent_at_least_2 <- mean(data_for_m$number_pathos >= 2) * 100
label_text <- paste0(round(percent_at_least_2, 1), "% with ≥2 pathogens")

# Create histogram with annotation
pathos_richness_plot <- ggplot(data_for_m, aes(x = number_pathos)) +
  geom_histogram(binwidth = 1, fill = "#fddc64", color = "#7a6a45", linewidth = 1.2) +
  labs(
    title = "Distribution of Pathogen Richness per Individual",
    x = "Number of Pathogens Detected",
    y = "Number of Individuals"
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
    color = "#7a6a45"
  )
pathos_richness_plot






# Plot richness (hist, boxplots et co)
# plot beta div





#Prevalence variable
#selon ?

#nbr d'individus infectés variable (on regarde que la prevalence a chaque fois c'est bon wsh, c'est pas ça le "danger")
# mais attention implique de standardiser data pour meme effort capture : + probleme c'est que c'est pas des vrais
#abondances qu'on regarde mais du relatif - pas forcement grave là aussi on peut dire c'est relatif non ?
#pour faire ça regarder les individus supprimés dans autre script quand retrait 4 eme nuit de piegeage
#peut-être autre truc a check


#pour les variables environnements : jamais regardé hobo pour l'instant





#modifier pour que plus de mention de infectious mais agents infectieux
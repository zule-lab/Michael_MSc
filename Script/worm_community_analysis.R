# Data and packages
worms <- read.csv("Input/worm_species_data.csv")
worms_env <- read.csv("Input/worm_env_variables.csv")
install.packages(c("tidyverse","vegan","ggrepel"))
my_packages <- c("tidyverse","vegan","ggrepel")
lapply(my_packages, require, character.only = TRUE)

# Creating new column for site ID
worms <- worms %>%
  mutate(site_id = str_extract(complete_id, "[^-]+"), .before = "complete_id")

# Creating new column for plot ID
worms <- worms %>%
  mutate(plot_id = str_extract(complete_id, "\\b\\w+$"), .before = "complete_id")

# Creating community matrix of sites (rows) by species (columns)
worms_comm_matrix <- worms %>%
  group_by(site_id, sp_code) %>% 
  summarise(n = n() )%>%
  pivot_wider(names_from = sp_code, values_from = n, values_fill = 0) %>%
  as.data.frame()

# Calling rownames in dataframe to sites code and removing first column
rownames(worms_comm_matrix) <- worms_comm_matrix$site_id
worms_comm_matrix <- worms_comm_matrix[,-1]

# Set RNG
set.seed(123) # To obtain the same NMDS results at every R session

# Convert absolute abundances into relative abundances 
worms_comm_matrix <- decostand(worms_comm_matrix, method = "total")

# Create distance matrix
worms_comm_dist <- vegdist(worms_comm_matrix, method = "bray") 

# Transform distance matrix into csv
worms_comm_dist <- as.matrix(worms_comm_dist, labels = TRUE)
write.csv(worms_comm_dist, "Output/worms_comm_dist.csv")

# Load distance matrix
worms_comm_dist <- read.csv("Output/worms_comm_dist.csv")
worms_comm_dist <- worms_comm_dist[-1] # Removing first column of distance matrix

# Run NMDS
worms_comm_NMDS <- metaMDS(worms_comm_dist, k = 2, distance = "bray", wascores = TRUE)
worms_comm_NMDS$stress # Indicates a sort of goodness-of-fit
stressplot(worms_comm_NMDS) # Shepard Plot

# Extracting point coordinates
NMDS_coords <- as.data.frame(worms_comm_NMDS$points)

# Calculating species scores
sp <- wascores(x = NMDS_coords, w = worms_comm_matrix, expand = TRUE)
sp_scores <- as.data.frame(sp)
sp_scores$sp <- row.names(sp_scores)

# Adding historical land-use to dataframe
NMDS_coords$site <- row.names(NMDS_coords)
NMDS_coords <- NMDS_coords %>%
  mutate(legacy = case_when(
    startsWith(site, "AGR") ~ "Agricultural",
    startsWith(site, "FOR") ~ "Forested",
    startsWith(site, "IND") ~ "Industrial"),
    .after = "site")

# Graphing
ggplot(data = NMDS_coords, aes(x = MDS1, y = MDS2), colour = legacy) +
  geom_point(aes(colour = legacy), size = 2) +
  #geom_text_repel(aes(colour = legacy, label = rownames(NMDS_coords)), max.overlaps = 25,fontface = "italic") +
  geom_point(data = sp_scores, aes(x = MDS1, y = MDS2), size = 2, colour = "magenta2") +
  #geom_text_repel(data = sp_scores, aes(label = rownames(sp_scores)), hjust = -0.7, vjust = -2.3, fontface = "bold") +
  stat_ellipse(aes(colour = legacy, fill = legacy), level = 0.95, geom = "polygon", alpha = 0.25, lwd = 1) +
  theme_minimal() +
  labs(y = "NMDS2", x = "NMDS1", colour = NULL, fill = NULL) + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom") + 
  scale_color_manual(values = c("goldenrod","green4","blue3"))
ggsave("nmds_plot.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/diversity")

# PERMANOVA
worms_perma <- adonis2(worms_comm_dist ~ legacy + soil_ph + soil_temp + soil_moist, data = worms_env, method = "bray")
worms_perma

# ANOSIM
worms_ano <- anosim(worms_comm_dist, worms_env$legacy, distance = "bray", permutations = 999)
worms_ano


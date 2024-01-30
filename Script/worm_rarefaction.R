# Data and packages
worms <- read.csv("Input/worm_species_data.csv")
install.packages(c("tidyverse","ggrepel","iNEXT"))
my_packages <- c("tidyverse","ggrepel","iNEXT")
lapply(my_packages, require, character.only = TRUE)

# Creating new column for site ID
worms <- worms %>%
  mutate(site_id = str_extract(complete_id, "[^-]+"), .before = "complete_id")

# Creating new column for plot ID
worms <- worms %>%
  mutate(plot_id = str_extract(complete_id, "\\b\\w+$"), .before = "complete_id")

# Creating community matrix of sites by species
worms_comm_matrix <- worms %>%
  group_by(legacy, sp_code) %>% 
  summarize(n = n() )%>%
  pivot_wider(names_from = legacy, values_from = n, values_fill = 0) %>%
  as.data.frame()

# Calling rownames in dataframe to species code and removing first column
rownames(worms_comm_matrix) <- worms_comm_matrix$sp_code
worms_comm_matrix <- worms_comm_matrix[,-1]

# Generate rarefaction/extrapolation with iNEXT
worms_re <- iNEXT(worms_comm_matrix, q = 0, datatype = "abundance")
worms_re$DataInfo # Summarizes data information
worms_re$iNextEst # Shows diversity estimates along with related statistics for a series of rarefied and extrapolated samples
worms_re$AsyEst # Shows asymptotic diversity estimates along with related statistics

# Extract minimum sampling coverage
cov <- min(worms_re$DataInfo$SC)

# Estimate diversity metrics at lowest sampling coverage
div <- estimateD(worms_comm_matrix, datatype = "abundance", base = "coverage", level = cov, conf = 0.95)

# Reformat for easier interpretation
div_wide <- pivot_wider(div, id_cols = Assemblage, names_from = Order.q,
                     names_sep = ".", values_from = c(Method, SC, qD))
div_final <- div_wide %>%
  select(-c(Method.1, Method.2, SC.1, SC.2)) %>%
  rename(species = Assemblage,
         method = Method.0,
         min_samp_cov = SC.0,
         richness = qD.0,
         Shannon = qD.1,
         Simpson = qD.2)

# Saving output
write.csv(div_final,"Output/worms_biodiv_metrics.csv", row.names = FALSE) 

# Plotting estimated curves
ggiNEXT(worms_re, type = 1) # Sample‐size‐based R/E curves
ggiNEXT(worms_re, type = 2) # Sample completeness curve
ggiNEXT(worms_re, type = 3) # Coverage-based R/E curve

# Drawing own sample-size-based R/E curves
df <- fortify(worms_re, type = 1) # Saving curve as dataframe
df.point <- df[which(df$Method == "Observed"),] # Creating object of observed points
df.line <- df[which(df$Method!= "Observed"),] # Creating object of every other point
df.line$Method <- factor(df.line$Method, c("Rarefaction", "Extrapolation")) # Making two linetypes into factor to distinguish rarefied from extrapolated parts

ggplot(df, aes(x = x, y = y, colour = Assemblage)) + 
  geom_point(size = 3, data = df.point) +
  geom_line(aes(linetype = Method), lwd = 1.5, data = df.line, show.legend = FALSE) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), alpha = 0.1, show.legend = FALSE) +
  labs(colour = " ", 
       x = "\nNumber of individuals",
       y = "Species richness\n") +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "top") +
  scale_colour_manual(values = c("goldenrod","forestgreen","blue"),
                      labels = c("Agricultural","Forested","Industrial")) 
ggsave("rarefac_extrap_curves.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Diversity")


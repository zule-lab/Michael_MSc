# Data and packages
soil <- read.csv("Input/bulk_data.csv")
soil$legacy <- as.factor(soil$legacy)
install.packages(c("tidyverse","ggrepel","ggcorrplot","lme4","lmerTest","MuMIn","performance","DHARMa","car"))
my_packages <- c("tidyverse","ggrepel","ggcorrplot","lme4","lmerTest","MuMIn","performance","DHARMa","car")
lapply(my_packages, require, character.only = TRUE)

# Determining which ρb to use in analysis
cor.test(soil$bulk_core, soil$bulk_fine, method = "pearson")
cor.test(soil$bulk_core, soil$bulk_hyb, method = "pearson")
cor.test(soil$bulk_fine, soil$bulk_hyb, method = "pearson")

soil$legacy <- factor(soil$legacy, levels = c("forested", "agriculture", "industrial")) # Reorder levels
p1 <- ggplot(data = soil, aes(x = legacy, y = bulk_core, fill = legacy)) + 
  geom_boxplot() + 
  labs(y = expression(paste(rho)[core])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16)) + 
  scale_fill_manual(values = c("forestgreen","gold","blue"))

p2 <- ggplot(data = soil, aes(x = legacy, y = bulk_hyb, fill = legacy)) + 
  geom_boxplot() +
  labs(y = expression(paste(rho)[hybrid])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16)) + 
  scale_fill_manual(values = c("forestgreen","gold","blue")) 

p3 <- ggplot(data = soil, aes(x = legacy, y = bulk_fine, fill = legacy)) + 
  geom_boxplot() +
  labs(y = expression(paste(rho)[fine])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16)) + 
  scale_fill_manual(values = c("forestgreen","gold","blue"))

require(patchwork)
p <- p1 + p2 + p3
p
ggsave("bulk_comparison.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_soil_bulk")

# Fitting the model
soil$legacy <- relevel(soil$legacy, ref = "forested") # Setting model reference 
contrasts(soil$legacy)
mod_bulk <- lmer(bulk_hyb ~ legacy + age + legacy*age + (1|site_id), data = soil, REML = TRUE)

# Goodness-of-fit
r.squaredGLMM(mod_bulk)

# Test of effects: KR approximation
anova(mod_bulk, type = 3, ddf = "Kenward-Roger")

# Checking for violated model assumptions
check_model(mod_bulk)
simul_output <- (simulateResiduals(fittedModel = mod_bulk))
plot(simul_output) 

# Check for model independence
e <- resid(mod_bulk)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
soil$predict_bulk <- predict(mod_bulk, type = "response", newdata = soil)
ggplot(soil, aes(x = legacy, y = predict_bulk, fill = legacy)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) +
  geom_hline(yintercept = 1.47, linetype = "dashed", color = "red", size = 1) + # Represents the density at which root growth is impeded in clayey soil
  theme_minimal() +
  labs(y = expression(paste(rho)[b]~(gm^-3)),
       x = "\nHistorical land-use") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial")) + 
  scale_fill_manual(values = c("forestgreen","gold","blue"))  
ggsave("glmm_result.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_soil_bulk")

# Visualizing fixed effects
require(visreg)
legs <- c(
  `forested` = "Forested",
  `agriculture` = "Agricultural",
  `industrial` = "Industrial")

visreg(mod_bulk, "age", by = "legacy", gg = TRUE) +
  theme_minimal() + 
  labs(y = expression(paste(rho)[b]~(gm^-3)),
       x = "\nTime since establishment (years)") +
  geom_point(aes(fill = legacy)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "none") +
  facet_wrap(legacy ~ ., labeller = as_labeller(legs)) 
ggsave("visreg_age_contrast_BD.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_soil_bulk")


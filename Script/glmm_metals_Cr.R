# Data and packages
metal <- read.csv("Input/metals_glmm.csv")
metal$legacy <- as.factor(metal$legacy)
install.packages(c("tidyverse","ggrepel","lme4","lmerTest","MuMIn","performance","DHARMa","car"))
my_packages <- c("tidyverse","ggrepel","lme4","lmerTest","MuMIn","performance","DHARMa","car")
lapply(my_packages, require, character.only = TRUE)

# Boxplot of soil chromium
metal$legacy <- factor(metal$legacy, levels = c("forested", "agricultural", "industrial")) # Reorder levels
ggplot(data = metal, aes(x = legacy, y = Cr_conc, fill = legacy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) +
  labs(x = "\nHistorical land-use",
       y = "Concentration (ppm)\n") +
  theme_gray() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +  
  scale_fill_manual(values = c("forestgreen","gold","blue")) +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial"))
ggsave("chromium_boxplot.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_metals/Chromium") 

# Fitting the model
metal <- metal[-116,] # Removing extreme outlier from analysis
mod_Cr <- lmer(Cr_conc ~ legacy + (1|site_id), data = metal, REML = TRUE)

# Goodness-of-fit
r.squaredGLMM(mod_Cr)

# Test of effects: KR approximation
anova(mod_Cr, type = 3, ddf = "Kenward-Roger")

# Checking for violated model assumptions
check_model(mod_Cr)
simul_output <- (simulateResiduals(fittedModel = mod_Cr))
plot(simul_output) 

# Check for model independence
e <- resid(mod_Cr)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
metal$predict_Cr <- predict(mod_Cr, type = "response", newdata = metal)
ggplot(metal, aes(x = legacy, y = predict_Cr, fill = legacy)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) + 
  theme_gray() +
  labs(y = "Predicted concentration (ppm)\n",
       x = "\nHistorical land-use") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") + 
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial")) + 
  scale_fill_manual(values = c("forestgreen","gold","blue"))  
ggsave("glmm_result.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_metals/Chromium")


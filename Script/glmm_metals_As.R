# Data and packages
metal <- read.csv("Input/metals_glmm.csv")
metal$legacy <- as.factor(metal$legacy)
install.packages(c("tidyverse","ggrepel","lme4","lmerTest","MuMIn","performance","DHARMa","car"))
my_packages <- c("tidyverse","ggrepel","lme4","lmerTest","MuMIn","performance","DHARMa","car")
lapply(my_packages, require, character.only = TRUE)

# Boxplot of soil arsenic
metal$legacy <- factor(metal$legacy, levels = c("forested", "agricultural", "industrial")) # Reorder levels
ggplot(data = metal, aes(x = legacy, y = As_conc, fill = legacy)) + 
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
ggsave("arsenic_boxplot.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_metals/Arsenic") 

# Fitting the model
mod_As <- lmer(As_conc ~ legacy + (1|site_id), data = metal, REML = TRUE)

# Goodness-of-fit
r.squaredGLMM(mod_As)

# Test of effects: KR approximation
anova(mod_As, type = 3, ddf = "Kenward-Roger")

# Checking for violated model assumptions
check_model(mod_As)
simul_output <- (simulateResiduals(fittedModel = mod_As))
plot(simul_output) 

# Check for model independence
e <- resid(mod_As)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
metal$predict_As <- predict(mod_As, type = "response", newdata = metal)
ggplot(metal, aes(x = legacy, y = predict_As, fill = legacy)) +
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
ggsave("glmm_result.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_metals/Arsenic")


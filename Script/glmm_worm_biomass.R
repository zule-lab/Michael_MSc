# Data and packages
worm <- read.csv("Input/worm_glmm.csv")
worm$legacy <- as.factor(worm$legacy)
install.packages(c("tidyverse","ggrepel","fitdistrplus","lme4","lmerTest","MuMIn","performance","DHARMa","car"))
my_packages <- c("tidyverse","ggrepel","fitdistrplus","lme4","lmerTest","MuMIn","performance","DHARMa","car")
lapply(my_packages, require, character.only = TRUE)

# Setting model reference 
worm$legacy <- relevel(worm$legacy, ref = "forested") 
contrasts(worm$legacy)

# Fitting model
mod_biom <- lmer(biomass ~ legacy + soil_moist + (1|site_id), data = worm, REML = TRUE)

# Goodness-of-fit
r.squaredGLMM(mod_biom)

# Summary of model statistics
summary(mod_biom)

# Test of effects: KR approximation
anova(mod_biom, type = 3, ddf = "Kenward-Roger")

# Checking for violated model assumptions
check_model(mod_biom)
simul_output <- (simulateResiduals(fittedModel = mod_biom))
plot(simul_output) 

# Check for model independence
e <- resid(mod_biom)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
worm$predict_biom <- predict(mod_biom, type = "response", newdata = worm)
worm$legacy <- factor(worm$legacy, levels = c("forested", "agricultural", "industrial"))

ggplot(worm, aes(x = legacy, y = predict_biom, fill = legacy)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) + 
  theme_minimal() +
  labs(y = expression(Fresh~biomass~(kg~m^-2)),
       x = "\nHistorical land-use") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial")) + 
  scale_fill_manual(values = c("forestgreen","gold","blue")) 
ggsave("glmm_result1.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Biomass")

# Visualizing fixed effects
require(visreg)
legs <- c(
  `forested` = "Forested",
  `agricultural` = "Agriculture",
  `industrial` = "Industrial")

visreg(mod_biom, "soil_moist", by = "legacy", gg = TRUE) +
  theme_minimal() + 
  labs(y = expression(Fresh~biomass~(kg~m^-2)),
       x = "\nSoil moisture") +
  geom_point(aes(fill = legacy)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "none") +
  facet_wrap(legacy ~ ., labeller = as_labeller(legs)) 
ggsave("visreg_moist_contrast.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Biomass")


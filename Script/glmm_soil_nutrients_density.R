# Install and load necessary packages
install.packages(c("tidyverse","ggrepel","lme4","lmerTest","MuMIn","performance","DHARMa","car"))
my_packages <- c("tidyverse","ggrepel","lme4","lmerTest","MuMIn","performance","DHARMa","car")
lapply(my_packages, require, character.only = TRUE)

## ANALYSING CARBON DENSITY ----------------------------------------------------

# Load data
nut <- read.csv("Input/nutrients_data.csv")
nut$legacy <- as.factor(nut$legacy)

# Creating new column for site ID
nut <- nut %>%
  mutate(site_id = str_extract(complete_id, "[^-]+"), .before = "complete_id")

# Creating new column for plot ID
nut <- nut %>%
  mutate(plot_id = str_extract(complete_id, "\\b\\w+$"), .before = "complete_id")

# Removing false zeros
nut <- filter_if(nut, is.numeric, all_vars((.) != 0))

# Setting model reference 
nut$legacy <- relevel(nut$legacy, ref = "forested") 
contrasts(nut$legacy)

# Fitting model
mod_Cden <- lmer(log(C_den) ~ legacy + age + legacy*age + (1|site_id), data = nut, REML = TRUE)

# Goodness-of-fit
r.squaredGLMM(mod_Cden)

# Summary of model statistics
summary(mod_Cden)

# Test of effects: KR approximation
anova(mod_Cden, type = 3, ddf = "Kenward-Roger")

# Post-hoc pairwise comparison
require(emmeans)
Cden_emm <- emmeans(mod_Cden, "legacy")
pairs(Cden_emm, type = "response")
plot(Cden_emm, comparisons = TRUE)
emtrends(mod_Cden, pairwise ~ legacy, var = "age", infer = T)

# Checking for violated model assumptions
check_model(mod_Cden)
simul_output <- (simulateResiduals(fittedModel = mod_Cden))
plot(simul_output) 

# Check for model independence
e <- resid(mod_Cden)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
nut$predict_Cden <- predict(mod_Cden, type = "response", newdata = nut)
nut$predict_Cden <- exp(nut$predict_Cden)
nut$legacy <- factor(nut$legacy, levels = c("forested", "agriculture", "industrial"))

ggplot(nut, aes(x = legacy, y = predict_Cden, fill = legacy)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) + 
  theme_minimal() +
  labs(y = expression(Predicted~carbon~density~(kg~m^-2)),
       x = "\nHistorical land-use") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial")) + 
  scale_fill_manual(values = c("forestgreen","gold","blue")) 
ggsave("C_density_result1.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_carbon")


# Visualizing fixed effects
require(visreg)
legs <- c(
  `forested` = "Forested",
  `agriculture` = "Agricultural",
  `industrial` = "Industrial")

visreg(mod_Cden, "age", by = "legacy", gg = TRUE, trans = exp, rug = FALSE) +
  theme_minimal() + 
  labs(y = "C density\n",
       x = "\nTime since establishment (years)") +
  geom_point(aes(fill = legacy)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "none") +
  facet_wrap(legacy ~ ., labeller = as_labeller(legs)) 
ggsave("visreg_age_contrast_Cden.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_carbon")


## ANALYZING NITROGEN DENSITY --------------------------------------------------

# Load data
nut <- read.csv("Input/nutrients_data.csv")
nut$legacy <- as.factor(nut$legacy)

# Creating new column for site ID
nut <- nut %>%
  mutate(site_id = str_extract(complete_id, "[^-]+"), .before = "complete_id")

# Creating new column for plot ID
nut <- nut %>%
  mutate(plot_id = str_extract(complete_id, "\\b\\w+$"), .before = "complete_id")

# Setting model reference 
nut$legacy <- relevel(nut$legacy, ref = "forested") 
contrasts(nut$legacy)

# Fitting model
mod_Nden <- lmer(log(N_den) ~ legacy + age + legacy*age + (1|site_id), data = nut, REML = TRUE)

# Goodness-of-fit
r.squaredGLMM(mod_Nden)

# Summary of model statistics
summary(mod_Nden)

# Test of effects: KR approximation
anova(mod_Nden, type = 3, ddf = "Kenward-Roger")

# Checking for violated model assumptions
check_model(mod_Nden)
simul_output <- (simulateResiduals(fittedModel = mod_Nden))
plot(simul_output) 

# Check for model independence
e <- resid(mod_Nden)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
nut$predict_Nden <- predict(mod_Nden, type = "response", newdata = nut)
nut$predict_Nden <- exp(nut$predict_Nden)
nut$legacy <- factor(nut$legacy, levels = c("forested", "agriculture", "industrial"))

ggplot(nut, aes(x = legacy, y = predict_Nden, fill = legacy)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) + 
  theme_minimal() +
  labs(y = expression(Predicted~nitrogen~density~(kg~m^-2)),
       x = "\nHistorical land-use") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial")) + 
  scale_fill_manual(values = c("forestgreen","gold","blue")) 
ggsave("N_density_result1.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_nitrogen")

# Visualizing fixed effects
require(visreg)
legs <- c(
  `forested` = "Forested",
  `agriculture` = "Agricultural",
  `industrial` = "Industrial")

visreg(mod_Nden, "age", by = "legacy", gg = TRUE, trans = exp, rug = FALSE) +
  theme_minimal() + 
  labs(y = "N density\n",
       x = "\nTime since establishment (years)") +
  geom_point(aes(fill = legacy)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "none") +
  facet_wrap(legacy ~ ., labeller = as_labeller(legs)) 
ggsave("visreg_age_contrast_Nden.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_nitrogen")


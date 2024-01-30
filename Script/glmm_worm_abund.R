# Data and packages
worm <- read.csv("Input/worm_glmm.csv")
worm$legacy <- as.factor(worm$legacy)
install.packages(c("tidyverse","ggrepel","ggcorrplot","lme4","lmerTest","performance","DHARMa","car"))
my_packages <- c("tidyverse","ggrepel","ggcorrplot","lme4","lmerTest","performance","DHARMa","car")
lapply(my_packages, require, character.only = TRUE)

# Boxplot of earthworm abundance
worm$legacy <- factor(worm$legacy, levels = c("forested", "agricultural", "industrial")) # Reorder levels
ggplot(worm, aes(x = legacy, y = worm_abund, fill = legacy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) +
  labs(x = "\nHistorical land-use",
       y = "Count (# indiv.)\n") +
  theme_gray() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_fill_manual(values = c("forestgreen","gold","blue")) +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial"))
ggsave("abundance_boxplot.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Abundance") 

# Checking collinearity: copper/zinc
met_lm <- lm(Cu_conc ~ Zn_conc, data = worm)
summary(met_lm)
cor(worm$Cu_conc, worm$Zn_conc)
ggplot(worm, aes(x = Zn_conc, y = Cu_conc)) +
  geom_point(shape = 16, size = 2) +
  geom_smooth(method = lm, formula = y ~ x, col = "blue") +
  labs(y = expression(Copper~levels~(mg~kg^-1)),
       x = expression(Zinc~levels~(mg~kg^-1)),
       title = "Heavy metal collinearity",
       subtitle = "Pearson's r = 0.8359") +
  theme(axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 14))

# Checking collinearity: carbon/nitrogen
worm_lm <- lm(carb ~ nitro, data = worm)
summary(worm_lm)
cor(worm$carb, worm$nitro)
ggplot(worm, aes(x = carb, y = nitro)) +
  geom_point(shape = 16, size = 2) +
  geom_smooth(method = lm, formula = y ~ x, col = "blue") +
  labs(y = "Total soil carbon (%)",
       x = "Total soil nitrogen (%)",
       title = "Nutrients collinearity",
       subtitle = "Pearson's r = 0.8947") +
  theme(axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 14))

# Fitting the models
worm$legacy <- relevel(worm$legacy, ref = "forested") # Setting model reference 
contrasts(worm$legacy)

M0 <- glmer.nb(worm_abund ~ legacy + soil_moist + soil_ph + soil_temp + (1|site_id),
               data = worm,
               control = glmerControl(optimizer = "bobyqa"))

M1 <- glmer.nb(worm_abund ~ legacy + soil_moist + soil_ph + soil_temp + mean_palat + (1|site_id),
                      data = worm,
                      control = glmerControl(optimizer = "bobyqa"))

sc <- function(x){(x-min(x))/(max(x)-min(x))} # Scaling heavy metal levels: between 0 and 1
worm <- worm %>%
  mutate(Cu_sc = sc(worm$Cu_conc),
         Zn_sc = sc(worm$Zn_conc))

M2 <- glmer.nb(worm_abund ~ legacy + soil_moist + soil_ph + soil_temp + Cu_sc + (1|site_id),
               data = worm,
               control = glmerControl(optimizer = "bobyqa"))

M3 <- glmer.nb(worm_abund ~ legacy + soil_moist + soil_ph + soil_temp + Zn_sc + (1|site_id),
               data = worm,
               control = glmerControl(optimizer = "bobyqa"))

M4 <- glmer.nb(worm_abund ~ legacy + soil_moist + soil_ph + soil_temp + nitro + (1|site_id),
               data = worm,
               control = glmerControl(optimizer = "bobyqa"))

M5 <- glmer.nb(worm_abund ~ legacy + soil_moist + soil_ph + soil_temp + carb + (1|site_id),
               data = worm,
               control = glmerControl(optimizer = "bobyqa"))

# Model performance using LRT
anova(M0, M1) 
anova(M0, M2) 
anova(M0, M3) 
anova(M0, M4) 
anova(M0, M5) 

# Goodness-of-fit
deviance(M0) 
df.residual(M0) 
1-pchisq(41.69526, 37) 
logLik(M0)

# Test of effects: Type III Wald Chi-square
Anova(M0, type = "III")

# Checking for violated model assumptions
check_model(M0)
simul_output <- (simulateResiduals(fittedModel = M0))
plot(simul_output) 

# Check overdispersion
check_overdispersion(M0)

# Check zero-inflation
check_zeroinflation(M0)

# Check for model independence
e <- resid(M0)
durbinWatsonTest(e)
acf(e)

# Graphing model predictions
worm$predict_abund <- predict(M0, type = "response", newdata = worm)
ggplot(worm, aes(x = legacy, y = predict_abund, fill = legacy)) +
  geom_boxplot() +
  geom_jitter(shape = 16, size = 1.5, position = position_jitter(0.2)) + 
  theme_minimal() +
  labs(y = "Predicted abundance (indvs)\n",
       x = "\nHistorical land-use") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none") +
  scale_x_discrete(labels = c("Forested","Agricultural","Industrial")) + 
  scale_fill_manual(values = c("forestgreen","gold","blue")) 
ggsave("glmm_result1.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Abundance")

# Visualizing fixed effects
require(visreg)
legs <- c(
  `forested` = "Forested",
  `agricultural` = "Agricultural",
  `industrial` = "Industrial")

visreg(M0, "soil_moist", by = "legacy", gg = TRUE) +
  theme_minimal() + 
  labs(y = "log(Abundance)\n",
       x = "\nSoil moisture") +
  geom_point(aes(fill = legacy)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "none") +
  facet_wrap(legacy ~ ., labeller = as_labeller(legs)) 
ggsave("visreg_moist_contrast.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Abundance")

visreg(M0, "soil_ph", by = "legacy", gg = TRUE) +
  theme_minimal() + 
  labs(y = "log(Abundance)\n",
       x = "\nSoil pH") +
  geom_point(aes(fill = legacy)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18, color = "black"),
        legend.position = "none") +
  facet_wrap(legacy ~ ., labeller = as_labeller(legs)) 
ggsave("visreg_ph_contrast.png", plot = last_plot(), units = "px", dpi = 450, path = "Figures/glmm_worms/Abundance")


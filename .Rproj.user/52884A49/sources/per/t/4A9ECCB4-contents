################################################################################
#### Analysis of Simulations
################################################################################
# Description: Analysis of simulated data

# Clear R's brain
rm(list = ls())

# Load required packages
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(raster)        # To handle spatial data
library(amt)           # To fit distributions
library(survival)      # To run conditional logistic regression
library(ggh4x)         # For nested plots

# Load custom functions
# setwd("/home/david/Schreibtisch/SimulatingSSFs")
source("00_Functions.R")

# Reload simulated data and the associated truth
dat   <- read_rds("Simulation.rds")
truth <- read_rds("Truth.rds")

################################################################################
#### Global Model Parameters
################################################################################
n_rsteps  <- 10                                                                                 # Number of random steps to be generated
formula   <- ~ water + elev + dist + water_start:log_sl + water_start:sl + water_start:cos_ta   # Formula for the model (we'll use the same formula we used to simulate data)
alongstep <- F

################################################################################
#### Ensure All Functions Work as Intended using one of the Simulations
################################################################################
# Extract one of the simulated datasets at random
i <- sample(1:nrow(dat), size = 1)
cov <- dat$Covariates[[i]]
mov <- dat$Movement[[i]]

# Identify movement bursts (there is only going to be one as the data is
# perfectly regular, however, the function is needed as the subsequent function
# requires a burst ID)
mov_bursted <- computeBursts(mov, forgiveness = 1)

# Compute step metrics
mov_metrics <- computeMetrics(mov_bursted)

# We should find that movement kernels differ depending on whether a step starts
# in "water" cell or not
mov_metrics$water_start <- as.data.frame(extract(cov, mov[, c("x", "y")]))$water
mov_metrics %>%
  subset(!is.na(sl)) %>%
  ggplot(aes(x = sl, color = as.factor(water_start))) +
    geom_density() +
    theme_minimal() +
    scale_color_manual(values = c("orange", "cornflowerblue"))
mov_metrics %>%
  subset(!is.na(relta)) %>%
  ggplot(aes(x = relta, color = as.factor(water_start))) +
    geom_density() +
    theme_minimal() +
    xlim(c(-pi, pi)) +
    scale_color_manual(values = c("orange", "cornflowerblue"))

# Fit step- and turning-angle distributions to original and rarified data
fit_gamma <- fit_distr(mov_metrics$sl, dist_name = "gamma")
fit_mises <- fit_distr(mov_metrics$relta, dist_name = "vonmises")

# Generate random steps based on those distributions
mov_stepsel <- randomSteps(mov_metrics
  , n_rsteps   = n_rsteps
  , dist_gamma = fit_gamma
  , dist_mises = fit_mises
)

# Extract covariates
mov_covaris <- computeCovars(mov_stepsel
  , covars    = cov
  , alongstep = alongstep
)

# Estimate model parameters
mod <- runModel(mov_covaris, formula)
mod

# Update tentative distribution parameters for "in nonwater"
updateDists(
    shape   = fit_gamma$params$shape
  , scale   = fit_gamma$params$scale
  , kappa   = fit_mises$params$kappa
  , mu      = fit_mises$params$mu
  , beta_sl     = mod$Estimate[mod$Coefficient == "sl"]
  , beta_log_sl = mod$Estimate[mod$Coefficient == "log_sl"]
  , beta_cos_ta = mod$Estimate[mod$Coefficient == "cos_ta"]
)

# Update tentative distribution parameters for "in water"
updateDists(
    shape   = fit_gamma$params$shape
  , scale   = fit_gamma$params$scale
  , kappa   = fit_mises$params$kappa
  , mu      = fit_mises$params$mu
  , beta_sl     = mod$Estimate[mod$Coefficient == "sl"] + mod$Estimate[mod$Coefficient == "sl:water_start"]
  , beta_log_sl = mod$Estimate[mod$Coefficient == "log_sl"] + mod$Estimate[mod$Coefficient == "log_sl:water_start"]
  , beta_cos_ta = mod$Estimate[mod$Coefficient == "cos_ta"] + mod$Estimate[mod$Coefficient == "cos_ta:water_start"]
)

# Now let's clean up a bit
rm(
    cov
  , mov
  , mod
  , mov_bursted
  , mov_metrics
  , mov_stepsel
  , mov_covaris
  , fit_gamma
  , fit_mises
)
gc()

################################################################################
#### Proper Analysis
################################################################################
# Run the analysis for all scenarios
pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
dat$Estimates <- lapply(1:nrow(dat), function(x) {

  # Extract respective simulated data
  mov <- dat$Movement[[x]]
  cov <- dat$Covariates[[x]]

  # Run analysis
  mov_bursted <- computeBursts(mov, forgiveness = 1)
  mov_metrics <- computeMetrics(mov_bursted)
  fit_gamma <- fit_distr(mov_metrics$sl, dist_name = "gamma")
  fit_mises <- fit_distr(mov_metrics$relta, dist_name = "vonmises")
  mov_stepsel <- randomSteps(mov_metrics
    , n_rsteps   = n_rsteps
    , dist_gamma = fit_gamma
    , dist_mises = fit_mises
  )
  mov_covaris <- computeCovars(mov_stepsel, covars = cov, alongstep = alongstep)
  mod <- runModel(mov_covaris, formula = formula)

  # Update tentative distribution parameters for water = 0
  movement_nowater <- updateDists(
      shape   = fit_gamma$params$shape
    , scale   = fit_gamma$params$scale
    , kappa   = fit_mises$params$kappa
    , mu      = fit_mises$params$mu
    , beta_sl     = mod$Estimate[mod$Coefficient == "sl"]
    , beta_log_sl = mod$Estimate[mod$Coefficient == "log_sl"]
    , beta_cos_ta = mod$Estimate[mod$Coefficient == "cos_ta"]
  ) %>% setNames(c("Parameter", "Estimate")) %>% mutate(Scenario = "NoWater")

  # Update tentative distribution parameters for water = 1
  movement_water <- updateDists(
      shape   = fit_gamma$params$shape
    , scale   = fit_gamma$params$scale
    , kappa   = fit_mises$params$kappa
    , mu      = fit_mises$params$mu
    , beta_sl     = mod$Estimate[mod$Coefficient == "sl"] + mod$Estimate[mod$Coefficient == "sl:water_start"]
    , beta_log_sl = mod$Estimate[mod$Coefficient == "log_sl"] + mod$Estimate[mod$Coefficient == "log_sl:water_start"]
    , beta_cos_ta = mod$Estimate[mod$Coefficient == "cos_ta"] + mod$Estimate[mod$Coefficient == "cos_ta:water_start"]
  ) %>% setNames(c("Parameter", "Estimate")) %>% mutate(Scenario = "Water")

  # Let's also beautify results on the habitat-selection
  result_habitat <- subset(mod, Coefficient %in% c("water", "elev", "dist")) %>%
    dplyr::select(Parameter = Coefficient, Estimate) %>%
    expand_grid(., Scenario = c("NoWater", "Water"))

  # Put everything together
  results <- rbind(movement_nowater, movement_water, result_habitat) %>%
    arrange(Parameter)

  # Return the model results
  setTxtProgressBar(pb, x)
  return(results)
})

# Join the results with the truth and compute the bias
final <- dat %>%
  select(AutocorrRange, Replicate, Estimates) %>%
  unnest(Estimates) %>%
  left_join(., truth, by = c("Parameter", "Scenario")) %>%
  mutate(Bias = Truth - Estimate) %>%
  subset(Parameter != "mu") %>%
  mutate(AutocorrRange = as.factor(AutocorrRange)) %>%
  mutate(Kernel = ifelse(Parameter %in% c("shape", "scale", "kappa"), "Movement-Kernel", "Habitat-Selection"))

# Remove overly biased entries, those are due to convergence issues
final <- subset(final, Bias < 50)

# Plot
p <- ggplot(final, aes(x = AutocorrRange, y = Bias, fill = AutocorrRange, color = AutocorrRange)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_fill_viridis_d(alpha = 0.5, begin = 0.5) +
  scale_color_viridis_d(begin = 0.5) +
  facet_nested_wrap(~ Scenario + Kernel + Parameter, scales = "free", nrow = 2) +
  theme_minimal() +
  theme(
      legend.position  = "bottom"
    , strip.background = element_rect(fill = "gray95", color = "white")
  )
p

# Store plot to file
ggsave("FinalPlot.png"
  , plot   = p
  , device = png
  , width  = 6
  , height = 6
  , bg     = "white"
)

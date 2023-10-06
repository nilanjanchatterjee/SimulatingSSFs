################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Simulate movement with known preferences on virtual landscape

# Clear R's brain
rm(list = ls())

# Load required packages
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(raster)        # To handle spatial data
library(sf)            # For plotting spatial features
library(lubridate)     # To handle dates and times

# Load custom functions
# setwd("/home/david/Schreibtisch/SimulatingSSFs")
source("00_Functions.R")

################################################################################
#### Simulation Parameters
################################################################################
n              <- 300                    # Resolution of the covariate layers
n_rsteps       <- 10                     # Number of random steps to be generated
n_steps        <- 1000                   # Number of consecutive steps to be simulated
formula        <- ~ water + elev + dist  # Formula used to predict step-selection scores (movement metrics will automatically be added later)
prefs          <- c(-0.5, 0.5, -15)        # Preferences of the simulated individuals (the order needs to match the formula above)
stop           <- F                      # Should the simulation terminate at boundaries?
autocorr_range <- c(5, 10, 20)           # Autocorrelation range for the simulated layers
n_replicates   <- 100                    # Number of simulated landscapes/individuals per configuration
alongstep      <- F                      # Should covariate values be extracted along steps? Otherwise at their ends

# Define step distributions that depend on water cover. The first column of each
# matrix is going to represent the "no-water" scenario, the second column the
# "water" scenario
sl_dist <- cbind(c(3, 1.5), c(1.5, 1)) # Step-length parameters. First value = shape, second value = scale
ta_dist <- cbind(c(0.5, 0), c(0.2, 0)) # Turning-angle parameters. First value = kappa, second value = center

# Let's plot the step-length and turning-angle distributions for the two
# scenarios
toshow <- tibble(Scenario = c("nowater", "water")) %>%
  mutate(Values = map(Scenario, function(x) {
    index <- as.numeric(x == "water") + 1
    ta  <- seq(-pi, pi, length.out = 100)
    sl  <- seq(0, 20, length.out = 100)
    prob_ta <- dvonmises(ta, kappa = ta_dist[1, index], mu = ta_dist[2, index])
    prob_sl <- dgamma(sl, shape = sl_dist[1, index], scale = sl_dist[2, index])
    all <- tibble(
        Metric = c(rep("sl", 100), rep("ta", 100))
      , x = c(ta, sl)
      , y = c(prob_ta, prob_sl)
    )
    return(all)
  })) %>% unnest(Values)

# Plot. You'll see that individuals in water move less directed and with shorter
# steps
ggplot(toshow, aes(x = x, y = y, col = Scenario)) +
  geom_line() +
  facet_wrap(~ Metric, scales = "free") +
  scale_color_manual(values = c("orange", "cornflowerblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Let's put the true simulation parameters together and store them
truth_habitat <- data.frame(
    Scenario  = rep(c("NoWater", "Water"), each = length(c(attr(terms(formula), "term.labels"))))
  , Parameter = rep(c(attr(terms(formula), "term.labels")), 2)
  , Truth     = rep(c(prefs), 2)
)
truth_movement <- data.frame(
    Scenario  = rep(c("NoWater", "Water"), each = 4)
  , Parameter = c("shape", "scale", "kappa", "mu")
  , Truth     = c(sl_dist[, 1], ta_dist[, 1], sl_dist[, 2], ta_dist[, 2])
)
truth <- rbind(truth_habitat, truth_movement) %>% arrange(Scenario)
write_rds(truth, "Truth.rds")

# Let's specify the extent on which animals are allowed to move
ext <- extent(0, n, 0, n)
ext <- as(ext, "SpatialPolygons")

# Let's also specify an extent within which individuals will be released
ext2 <- extent(ext) - 100
ext2 <- as(ext2, "SpatialPolygons")

# And an extent to which we will expand covariates before extracting covariate
# values in the SSF (this is only to avoid numerical issues in case random steps
# leave the simulated area, which probably never happens, but anyhow)
ext3 <- extent(c(-50, 350, -50, 350))
ext3 <- as(ext3, "SpatialPolygons")

################################################################################
#### Main Simulation Functions
################################################################################
# Function to simulate center of attraction (using a distance)
simDistance <- function(n, x, y) {
  r    <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
  cent <- SpatialPoints(t(c(x, y)))
  dist <- distanceFromPoints(r, cent)
  dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))
  return(dist)
}

# Function to simulate a CONTINUOUS (e.g. elevation) layer
simElevation <- function(n, autocorr_range, seed = NULL) {
  r <- raster(res = 1, xmn = 0, xmx = n, ymn = 0, ymx = n)
  w <- focalWeight(r, d = autocorr_range, type = "circle")
  add_rows <- (nrow(w) - 1) / 2
  add_cols <- (ncol(w) - 1) / 2
  r <- raster(res = 1
    , xmn = 0 - add_cols
    , xmx = n + add_cols
    , ymn = 0 - add_rows
    , ymx = n + add_rows
  )
  if (is.null(seed)) {
      set.seed(round(runif(1, 0, 1e8)))
    } else {
      set.seed(seed)
  }
  r[] <- rnorm(n = ncell(r))
  r   <- focal(r, w = w, fum = "sum")
  r   <- trim(r)
  r <- (r - cellStats(r, min)) / (cellStats(r, max) - cellStats(r, min))
  return(r)
}

# Function to simulate a BINARY (e.g. water) layer
simWater <- function(n, autocorr_range, prop = 0.5, seed = NULL) {
  if (is.null(seed)) {
    seed <- round(runif(1, 0, 1e8))
  }
  r      <- simElevation(n, autocorr_range, seed)
  cutoff <- quantile(r, 1 - prop)
  r      <- r > cutoff
  return(r)
}

# Wrap those functions into a single "landscape simulation" function. Note that
# we're ensuring the seeds for forest and elevation to be different!!!
simLandscape <- function(n, autocorr_range, proportion_water, seed = NULL) {
  landscape <- stack(list(
      simDistance(n = n, x = n / 2, y = n / 2)
    , simElevation(n = n, autocorr_range = autocorr_range, seed = seed)
    , simWater(n = n, autocorr_range = autocorr_range, prop = proportion_water, seed = seed + 1)
  )) %>% setNames(c("dist", "elev", "water"))
  return(landscape)
}

# Try it
covars <- simLandscape(n = 300, autocorr_range = 20, proportion_water = 0.5, seed = 123)
plot(covars, col = hcl.colors(100))

# Wrapper function for the "move" function to simulate movement on a simulated
# landscape. Returns "observed" GPS data
simMove <- function(covars, messages = F, dropmetrics = T) {

  # Simulate movement
  sim <- move(
      xy        = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
    , covars    = covars
    , formula   = formula
    , prefs     = prefs
    , sl_dist   = sl_dist
    , ta_dist   = ta_dist
    , ext       = ext
    , n_steps   = n_steps
    , n_rsteps  = n_rsteps
    , stop      = stop
    , alongstep = alongstep
    , messages  = messages
  )

  # Add an artifial timestamp, the ID, as well as a unique step_id to the
  # data
  sim$timestamp <- lubridate::ymd_hms("2000-01-01 00:00:00") + hours(1:nrow(sim))

  # In reality, we don't observe step lengths etc.
  if (dropmetrics) {
    sim$absta <- NULL
    sim$ta    <- NULL
    sim$sl    <- NULL
  }
  return(sim)
}

################################################################################
#### Ensure All Functions Work as Intended
################################################################################
# Simulate covariates and movement across the layers
cat("Simulating example landscape\n")
cov <- simLandscape(n = n, autocorr_range = 10, proportion_water = 0.5, seed = 1)
cat("Simulating example trajectory\n")
mov <- simMove(cov, messages = T, dropmetrics = F)

# Visualize
as.data.frame(cov, xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext2), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    geom_path(data = mov, aes(x = x, y = y), inherit.aes = F) +
    geom_point(data = mov[1, ], aes(x = x, y = y), inherit.aes = F, col = "red", shape = 0) +
    geom_point(data = mov[nrow(mov), ], aes(x = x, y = y), inherit.aes = F, col = "green", shape = 2) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(
        axis.title.y    = element_text(angle = 0, vjust = 0.5)
      , legend.position = "none"
    )

################################################################################
#### Simulation
################################################################################
# Specify the different design combinations for which we want to simulate
dat <- expand_grid(
    AutocorrRange = autocorr_range
  , Replicate     = 1:n_replicates
)

# Simulate the covariate layers
cat("Simulating covariate layers...\n")
dat$Covariates <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    cov <- simLandscape(
        n                 = n
      , proportion_water  = 0.5
      , autocorr_range    = dat$AutocorrRange[x]
      , seed              = dat$Replicate[x]
    )
    return(cov)
})

# Now simulate movement across the covariates
cat("Simulating movement...\n")
dat$Movement <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    sim         <- simMove(covars = dat$Covariates[[x]], messages = F)
    sim$ID      <- x
    sim$step_id <- (x * n_steps) + (1:n_steps)
    return(sim)
})

# Let's expand the covariates slightly so that later we won't have random
# steps leaving the study area
cat("Extending covariate layers...\n")
dat$Covariates <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    cov <- extendRaster(dat$Covariates[[x]], ext3)
    return(cov)
})

# Check out structure of this -> contains simulation parameters and then the
# associated simulated covariates and movement.
print(dat)

# Visualize some (red = start, green = end). Individuals tend to move towards
# the center. Rerun the function below multiple times to show different
# simulations. The bottom right indicates the landscape simulation parameters.
plotSimulation(dat, index = sample(nrow(dat), size = 1))

# Store the entire simulation to file
write_rds(dat, "Simulation.rds")

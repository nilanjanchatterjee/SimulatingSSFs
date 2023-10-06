################################################################################
#### Functions to Compute Step Metrics
################################################################################
# Function to compute the absolute turning angle
absTA <- function(dx, dy) {
  absta <- atan2(dy, dx)
  absta <- (absta - pi / 2) * (-1)
  absta <- ifelse(absta < 0, 2 * pi + absta, absta)
  return(absta)
}

# Function to compute the relative turning angle
relTA <- function(absta) {
  relta <- (absta[-1] - absta[-length(absta)])
  relta <- c(NA, relta)
  relta <- ifelse(relta > +pi, relta - 2 * pi, relta)
  relta <- ifelse(relta < -pi, 2 * pi + relta, relta)
  return(relta)
}

# Function to compute step metrics
stepMet <- function(x, y) {

    # Compute distances moved in x and y direction
    dx <- c(x[-1], NA) - x
    dy <- c(y[-1], NA) - y

    # Calculate step length
    sl <- sqrt(dx ** 2 + dy ** 2)

    # Compute absolute turn angle
    absta <- absTA(dx, dy)

    # Compute relative turn angle
    relta <- relTA(absta)

    # Put metrics into data.frame
    metrics <- data.frame(sl = sl, absta = absta, relta = relta)

    # Return the metrics
    return(metrics)
}

################################################################################
#### Function to Interpolate Between Two Coordinates
################################################################################
# Extracting covariate values along spatial lines can be really really slow. To
# improve efficiency, you can therefore extract covariates along points that are
# distributed on the line instead. Note that you may want to modify this
# function if you work with non-projected data!

# Function to interpolate coordinates between two points
interpolatePoints <- function(x1, x2, y1, y2, by = 1){

  # Calculate length of line between points
  length <- sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

  # Calculate how many segments we need
  nsegs <- max(ceiling(length / by), 1)

  # Interpolate between points
  x <- seq(x1, x2, length.out = nsegs + 1)
  y <- seq(y1, y2, length.out = nsegs + 1)
  return(cbind(x, y))
}

# # Try it
# interpolatePoints(x1 = 0, x2 = 10, y1 = 0, y2 = 10, by = 2)
#
# # Let's generate a random line
# start <- runif(n = 2, min = 0, max = 100)
# end <- runif(n = 2, min = 0, max = 100)
# line <- spLines(SpatialPoints(rbind(start, end)))
#
# # Now distribute points on that line using our custom function
# # The smaller the value for "by", the more points are put on the line
# line_inter <- interpolatePoints(
#     x1 = start[1]
#   , x2 = end[1]
#   , y1 = start[2]
#   , y2 = end[2]
#   , by = 0.1
# )
#
# # Plot the line and interpolated coordinates
# plot(line)
# points(line_inter, col = "red", pch = 20, cex = 0.5)
#
# # Are results the same?
# a <- raster::extract(cov, line, fun = "mean")
# b <- colMeans(raster::extract(cov, line_inter))
# rbind(along_line = as.vector(a), at_points = b)
#
# # Let's see how extraction speeds compare
# benchmark <- microbenchmark(
#     AlongLine = raster::extract(cov, line, fun = "mean")
#   , AtPoints  = colMeans(raster::extract(cov, line_inter, fun = "mean"))
#   , times     = 20
# )
# summary(benchmark)

################################################################################
#### Von Mises Distribution Functions
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonmises <- function(x, kappa, mu){
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonmises <- function(n, kappa, mu, by = 0.01){
  x <- seq(-pi, +pi, by = by)
  probs <- dvonmises(x, kappa = kappa, mu = mu)
  random <- sample(x, size = n, prob = probs, replace = T)
  return(random)
}

# # Let's ensure that the function works as expected
# x <- seq(-pi, +pi, by = 0.01)
# y1 <- dvonmises(x, mu = 0, kappa = 0)
# y2 <- dvonmises(x, mu = 0, kappa = 1)
# y3 <- dvonmises(x, mu = 0, kappa = 2)
#
# # Visualize the density for different kappas
# plot(NA, xlim = c(-pi, +pi), ylim = c(0, 0.6), xlab = "Turning Angle", ylab = "Prob. Density")
# abline(v = 0, lty = 2, col = "gray80")
# lines(y1 ~ x, col = "blue")
# lines(y2 ~ x, col = "purple")
# lines(y3 ~ x, col = "red")
# legend("topright"
#   , lty    = 1
#   , col    = c("blue", "purple", "red")
#   , legend = c("kappa = 0", "kappa = 1", "kappa = 2")
# )

################################################################################
#### Function to Simulate Movement from a fitted SSF
################################################################################
# Function to simulate movement
move <- function(
      xy        = NULL    # Source point (in matrix form -> n * 2)
    , covars    = NULL    # Stack of covariate layers (names need to match formula!)
    , formula   = NULL    # Model formula used to predict selection score
    , prefs     = NULL    # Preferences used to predict selection score
    , sl_dist   = NULL    # Parameters describing the step length distribution
    , ta_dist   = NULL    # Parameters describing the turning angle distribution
    , ext       = NULL    # Extent on which animals are allowed to move
    , n_steps   = 10      # Number of steps simulated
    , n_rsteps  = 25      # Number of random steps proposed at each step
    , alongstep = T       # Should covariates be extracted along steps? Otherwise at their end
    , stop      = TRUE    # Should the simulation stop at boundaries?
    , messages  = F       # Would you like to print update messages?
  ) {

  # For testing only
  # xy        <- coordinates(spsample(ext2, type = "random", n = 1))
  # covars    <- covars
  # formula   <- ~ water + elev + dist  # Formula used to predict step-selection scores
  # prefs     <- c(-1, 0.5, -10)        # Preferences of the simulated individuals
  # sl_dist   <- cbind(c(2, 0.5), c(3, 1))   # Step-length parameters (gamma)
  # ta_dist   <- cbind(c(2, 0), c(2, 0)) # Turning-angle parameters (von mises)
  # n_rsteps  <- 25
  # n_steps   <- 200
  # stop      <- F
  # messages  <- T
  # alongstep <- F

  # Create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = NA
    , sl    = NA
  )

  # Simulate random steps
  if (messages) {
    pb <- txtProgressBar(0, n_steps, style = 3)
  }
  for (i in 2:n_steps) {

    # # For testing only
    # i <- 2

    # Figure out if individual is currently in water or not. This will determine
    # the distributions from which we sample step-lengths and turning-angles. 1
    # = nowater, 2 = water.
    inwater <- extract(covars, track[i, c("x", "y")])[, c("water")] + 1

    # Draw random turning angles
    ta_new <- rvonmises(n_rsteps
      , kappa = ta_dist[1, inwater]
      , mu    = ta_dist[2, inwater]
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist[1, inwater]
      , scale = sl_dist[2, inwater]
    )

    # Make sure that the steps cover at least a minimal distance (this is
    # relevant if we need to compute the log of it)
    sl_new[sl_new < 0.0001] <- 0.0001

    # Put the step lengths and turning angles into a new dataframe. These are
    # our proposed random steps.
    rand <- data.frame(
        absta  = track$absta[i - 1] + ta_new
      , ta     = ta_new
      , sl     = sl_new
    )

    # We need to make sure that the absolute turning angle ranges from 0 to 2 *
    # pi
    rand$absta[rand$absta > 2 * pi] <-
      rand$absta[rand$absta > 2 * pi] - 2 * pi
    rand$absta[rand$absta < 0] <-
      rand$absta[rand$absta < 0] + 2 * pi

    # Calculate new endpoints
    rand$x <- track$x[i] + sin(rand$absta) * rand$sl
    rand$y <- track$y[i] + cos(rand$absta) * rand$sl

    # Create spatial points from endpoints
    coordinates(rand) <- c("x", "y")

    # Depending on the answer in the beginning, the loop breaks if one of the
    # new coordinates is outside the map boundaries
    if (stop){
      if (nrow(rand[ext, ]) != n_rsteps){
        break
      }
    } else {
      rand <- rand[ext, ]
    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand, xy = T)
    rand$xy <- NULL

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")]
    endcoords   <- rand[, c("x", "y")]
    if (alongstep) {
      extracted <- sapply(1:nrow(endcoords), function(x){
        line <- interpolatePoints(
            x1 = begincoords[1, 1]
          , x2 = endcoords[x, 1]
          , y1 = begincoords[1, 2]
          , y2 = endcoords[x, 2]
          , by = 0.1
        )
        extr <- raster::extract(covars, line)
        extr <- colMeans(extr)
        return(extr)
      })
      extracted <- t(extracted)
    } else {
      extracted <- raster::extract(covars, endcoords)
    }

    # Bind with data on random steps
    rand <- cbind(rand, extracted)

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta)
    rand$log_sl <- log(rand$sl)

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(formula, rand)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% prefs)

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Sample one of the steps based on predicted probabilities
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track$absta[i] <- rand$absta
    track$ta[i] <- rand$ta
    track$sl[i] <- rand$sl
    track[i + 1, "x"] <- rand$x
    track[i + 1, "y"] <- rand$y

    # Print update
    if (messages) {
      setTxtProgressBar(pb, i)
    }
  }

  # Assign step numbers
  track$step_number <- 0:(nrow(track) - 1)

  # Return track, yet remove initial pseudo-fix
  return(track[-1, ])
}


################################################################################
#### Function to Extend a Covariate Layer
################################################################################
# Function that extends a covariate layer and fills the added border with values
# sampled from the layer
extendRaster <- function(x, y) {
  x_mask  <- is.na(x)
  x_mask  <- extend(x_mask, y, value = 0)
  x_large <- extend(x, y, value = NA)
  for (i in 1:nlayers(x_large)) {
    indices <- which(is.na(x_large[[i]][]))
    n       <- length(indices)
    values  <- na.omit(x[[i]][])
    x_large[[i]][indices] <- sample(values, size = n, replace = T)
    x_large[[i]] <- mask(x_large[[i]], x_mask[[i]], maskvalues = 1, updatevalue = NA)
  }
  names(x_large) <- names(x)
  return(x_large)
}

################################################################################
#### Function to Visualize Simulations
################################################################################
# Helper function to quickly visualize some simulations
plotSimulation <- function(data, index = 1) {

  # Extract covariates and movement at the respective index
  cov <- data$Covariates[[index]]
  mov <- data$Movement[[index]]

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
      ) +
      labs(caption = paste0(
          "AutocorrRange = ", data$AutocorrRange[[index]]
        , ", Replicate: ", data$Replicate[[index]])
      )
}

# Function to compute bursts (per ID, depending on the fogriveness). A new burst
# always starts if the step-duration exceeds the forgiveness.
computeBursts <- function(data, forgiveness) {

  # Nest data by id (a burst cannot expand across multiple ids)
  data_bursted <- data %>%
    group_by(ID) %>%
    nest() %>%

    # Compute bursts by id. A new burst is defined if the duration is >
    # forgiveness
    mutate(data = map(data, function(x) {
      x$duration <- lead(x$step_number) - x$step_number
      x$irregular <- x$duration > forgiveness
      x$burst <- NA
      x$burst[1] <- 1
      x$burst[2:nrow(x)] <- lag(cumsum(x$irregular) + 1)[-1]
    return(x)

  # Unnest the data again
  })) %>% unnest(data)

  # Return the "bursted" data
  return(data_bursted)
}

# Function to compute step metrics (step length, relative turning angle,
# absolute turning angle). Step metrics are calculated on bursts.
computeMetrics <- function(data) {

  # Nest by ID and burst
  data_metrics <- data %>%
    group_by(ID, burst) %>%
    nest() %>%

    # Compute step metrics
    mutate(data = map(data, function(z) {
      metrics <- stepMet(x = z$x, y = z$y)
      metrics <- cbind(z, metrics)
      return(metrics)
    })) %>%

    # Unnest again and tidy up
    unnest(data) %>%
    dplyr::select(burst, step_number, step_id, everything()) %>%
    ungroup()

  # Return the data containing step metrics
  return(data_metrics)
}

# Function to generate random steps. The approach parameter is used to specify
# the approach that should be used to generate random steps
randomSteps <- function(data, n_rsteps, dist_gamma, dist_mises) {

  # Generate a new column that indicates that the steps are "observed" steps
  data$case <- 1

  # Cannot work with steps that have no turning angle, so remove them
  data <- subset(data, !is.na(relta))

  # Create a new dataframe into which we can put alternative/random steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that these steps are random steps (case = 0)
  rand$case <- 0

  # Create random steps
  # Step lengths sampled from "minimal" step-duration distributions
  rand$sl <- rgamma(n = nrow(rand)
    , scale = dist_gamma$params$scale
    , shape = dist_gamma$params$shape
  )
  rand$relta_new <- rvonmises(n = nrow(rand)
    , kappa = dist_mises$params$kappa
    , mu    = dist_mises$params$mu
    , by    = 0.01
  )

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta_new - rand$relta
  rand$absta <- rand$absta + rand$relta_diff
  rand$absta <- ifelse(rand$absta < 0, 2 * pi + rand$absta, rand$absta)
  rand$absta <- ifelse(rand$absta > 2 * pi, rand$absta - 2 * pi, rand$absta)
  rand$relta <- rand$relta_new

  # Remove undesired stuff
  rand$relta_new <- NULL
  rand$relta_diff <- NULL

  # Put steps together
  all <- rbind(data, rand)
  all <- arrange(all, ID, step_id, desc(case))

  # Calculate new endpoints
  all$x_to <- all$x + sin(all$absta) * all$sl
  all$y_to <- all$y + cos(all$absta) * all$sl

  # Sort and ungroup
  all <- dplyr::select(
      all
    , ID
    , step_number
    , step_id
    , x
    , y
    , x_to
    , y_to
    , everything()
  )
  all <- ungroup(all)
  return(all)
}

# Function to extract covariates along steps and compute step covariates
computeCovars <- function(data, covars, multicore = F, alongstep = T) {

  # Extract
  if (alongstep) {
    extracted <- sapply(1:nrow(data), function(x){
      line <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 0.1
      )
      extr <- raster::extract(covars, line)
      extr <- colMeans(extr)
      return(extr)
    })
    extracted <- t(extracted)
  } else {
    extracted <- raster::extract(covars, cbind(data$x_to, data$y_to))
  }

  # We also want to get the values at the current location
  extracted_start <- as.data.frame(raster::extract(covars, cbind(data$x, data$y)))
  names(extracted_start) <- paste0(names(extracted_start), "_start")

  # Put all together
  extracted <- cbind(as.data.frame(extracted), extracted_start)
  data <- cbind(data, extracted)

  # Ensure that step lengths cover a minimal distance
  data$sl[data$sl == 0] <- min(data$sl[data$sl != 0])

  # Calculate derived movement metrics
  data$log_sl <- log(data$sl)
  data$cos_ta <- cos(data$relta)

  # Return the results
  return(data)
}

# Function to run the step selection model using two different approaches
runModel <- function(data, formula) {

  # Adjust model formula
  formula <- update(formula, case ~ . + sl + log_sl + cos_ta + strata(step_id))

  # Run cond. logistic regression model
  mod <- clogit(formula, data = data)

  # Extract model coefficients and put them into a dataframe. Also, compute
  # confidence intervals.
  ci <- confint(mod, level = 0.95)
  coefs <- summary(mod)$coefficients
  coefs <- data.frame(
      Coefficient = rownames(coefs)
    , Estimate    = coefs[, "coef"]
    , SE          = coefs[, "se(coef)"]
    , Z           = coefs[, "z"]
    , LCI         = ci[, 1]
    , UCI         = ci[, 2]
  )
  rownames(coefs) <- NULL

  # Return them
  return(coefs)
}

################################################################################
#### Function to update Step Distirbutions given Parameters
################################################################################
# Function to udpate distributions
updateDists <- function(shape, scale, kappa, mu, beta_log_sl, beta_sl, beta_cos_ta) {

  # Update distribution parameters
  shape_upd <- shape + beta_log_sl
  scale_upd <- 1 / (1 / scale - beta_sl)
  kappa_upd <- kappa + beta_cos_ta
  mu_upd    <- mu

  # Return the updated distribution parameters
  updated <- data.frame(
      Parameter = c("shape", "scale", "kappa", "mu")
    , Value     = c(shape_upd, scale_upd, kappa_upd, mu_upd)
  )
  return(updated)
}

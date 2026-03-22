
# Function for simulating seed data along an elevation gradient.
simulate_seed_data_nb <- function(
    elevation_seq = seq(1100, 2000, by = 100),
    n_per_type = 50,
    species = c("sp1", "sp2"),

    # Fixed effects (on log scale)
    intercept_log = log(100),
    elevation_slope = -0.001,   # effect per unit elevation
    treatment_effect = log(0.95), # hand vs natural (multiplicative)
    reference_elevation = 1100,

    # Random effects (SDs on log scale)
    sd_species_intercept = 0.2,
    sd_species_slope = 0.0005,
    sd_site_intercept = 0.15,

    # Dispersion (negative binomial: smaller = more variance)
    size = 10,

    # Output control
    return_list = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  datalist <- list()
  counter <- 1

  # Species random effects
  species_intercepts <- stats::rnorm(length(species), 0, sd_species_intercept)
  species_slopes <- stats::rnorm(length(species), 0, sd_species_slope)
  names(species_intercepts) <- species
  names(species_slopes) <- species

  for (elev in elevation_seq) {

    # Site random intercept
    site_intercept <- stats::rnorm(1, 0, sd_site_intercept)

    for (sp in species) {

      elev_centered <- elev - reference_elevation

      # Linear predictors (log scale)
      eta_natural <- intercept_log +
        (elevation_slope + species_slopes[sp]) * elev_centered +
        species_intercepts[sp] +
        site_intercept

      eta_hand <- eta_natural + treatment_effect

      # Convert to means
      mu_natural <- exp(eta_natural)
      mu_hand <- exp(eta_hand)

      # Simulate counts
      seed_natural <- stats::rnbinom(n_per_type, mu = mu_natural, size = size)
      seed_hand <- stats::rnbinom(n_per_type, mu = mu_hand, size = size)

      datalist[[counter]] <- data.frame(
        Seed_count = c(seed_natural, seed_hand),
        Elevation = elev,
        Pollination = rep(c("natural", "hand"), each = n_per_type),
        Species = sp,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }

  if (return_list) {
    return(datalist)
  } else {
    return(do.call(rbind, datalist))
  }
}


# Function for simulating binomial survival data
simulate_binomial_data <- function(
    n_groups = 10,                 # e.g. sites
    n_per_group = 50,              # individuals OR trials per group
    species = c("sp1", "sp2"),

    # Fixed effects (logit scale)
    intercept_logit = qlogis(0.6),   # baseline survival probability
    treatment_effect = qlogis(0.7) - qlogis(0.6),  # effect of treatment
    env_slope = 0,                  # continuous predictor effect

    # Optional environmental gradient
    env_values = seq_len(n_groups),

    # Random effects (SDs on logit scale)
    sd_species_intercept = 0.5,
    sd_species_slope = 0.2,
    sd_group_intercept = 0.5,

    # Structure
    treatment_labels = c("control", "treatment"),
    aggregated = FALSE,   # TRUE = binomial counts, FALSE = individual rows

    # Output
    return_list = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  datalist <- list()
  counter <- 1

  # Species random effects
  species_intercepts <- stats::rnorm(length(species), 0, sd_species_intercept)
  species_slopes <- stats::rnorm(length(species), 0, sd_species_slope)
  names(species_intercepts) <- species
  names(species_slopes) <- species

  for (g in seq_len(n_groups)) {

    group_intercept <- stats::rnorm(1, 0, sd_group_intercept)
    env_val <- env_values[g]

    for (sp in species) {

      for (trt in treatment_labels) {

        # Linear predictor
        eta <- intercept_logit +
          env_slope * env_val +
          species_intercepts[sp] +
          species_slopes[sp] * env_val +
          group_intercept

        if (trt == treatment_labels[2]) {
          eta <- eta + treatment_effect
        }

        p <- stats::plogis(eta)

        if (aggregated) {

          # Binomial counts
          successes <- stats::rbinom(1, size = n_per_group, prob = p)

          datalist[[counter]] <- data.frame(
            Successes = successes,
            Failures = n_per_group - successes,
            Trials = n_per_group,
            Proportion = successes / n_per_group,
            Group = g,
            Env = env_val,
            Treatment = trt,
            Species = sp,
            stringsAsFactors = FALSE
          )

        } else {

          # Individual-level data
          outcome <- stats::rbinom(n_per_group, size = 1, prob = p)

          datalist[[counter]] <- data.frame(
            Outcome = outcome,
            Group = g,
            Env = env_val,
            Treatment = trt,
            Species = sp,
            stringsAsFactors = FALSE
          )
        }

        counter <- counter + 1
      }
    }
  }

  if (return_list) {
    return(datalist)
  } else {
    return(do.call(rbind, datalist))
  }
}


# Function for simulating longitude associated datasets
simulate_longitudinal_data <- function(
    n_individuals = 100,
    n_timepoints = 5,
    time_values = seq_len(n_timepoints),

    # Optional grouping (e.g. site)
    n_groups = 1,

    # Fixed effects
    intercept = 10,
    time_slope = 0.5,
    treatment_effect = 2,
    time_treatment_interaction = 0,

    # Treatment structure
    treatment_labels = c("control", "treatment"),

    # Random effects (SDs)
    sd_individual_intercept = 2,
    sd_individual_slope = 0.3,
    cor_intercept_slope = 0,   # correlation between intercept & slope
    sd_group_intercept = 1,

    # Residual error
    sd_residual = 1,

    # Output control
    return_list = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  datalist <- list()
  counter <- 1

  # Assign individuals to groups
  group_ids <- rep(seq_len(n_groups), length.out = n_individuals)

  # Assign treatments (balanced)
  treatment_assign <- rep(treatment_labels, length.out = n_individuals)

  # Covariance matrix for random effects
  cov_matrix <- matrix(
    c(
      sd_individual_intercept^2,
      cor_intercept_slope * sd_individual_intercept * sd_individual_slope,
      cor_intercept_slope * sd_individual_intercept * sd_individual_slope,
      sd_individual_slope^2
    ),
    nrow = 2
  )

  # Draw individual random effects
  re_individual <- MASS::mvrnorm(
    n = n_individuals,
    mu = c(0, 0),
    Sigma = cov_matrix
  )

  colnames(re_individual) <- c("b_intercept", "b_slope")

  # Group-level random intercepts
  group_intercepts <- stats::rnorm(n_groups, 0, sd_group_intercept)

  for (i in seq_len(n_individuals)) {

    for (t in seq_along(time_values)) {

      time_val <- time_values[t]
      group_id <- group_ids[i]
      trt <- treatment_assign[i]

      # Linear predictor
      eta <- intercept +
        time_slope * time_val +
        re_individual[i, "b_intercept"] +
        re_individual[i, "b_slope"] * time_val +
        group_intercepts[group_id]

      if (trt == treatment_labels[2]) {
        eta <- eta +
          treatment_effect +
          time_treatment_interaction * time_val
      }

      # Add residual error
      y <- stats::rnorm(1, mean = eta, sd = sd_residual)

      datalist[[counter]] <- data.frame(
        ID = i,
        Group = group_id,
        Time = time_val,
        Treatment = trt,
        Response = y,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }

  if (return_list) {
    return(datalist)
  } else {
    return(do.call(rbind, datalist))
  }
}


# Function for generating phylogenetic trait data.
simulate_phylo_traits <- function(
    n_species = 50,
    tree = NULL,

    # Model choice
    model = c("BM", "OU"),

    # Brownian motion parameter
    sigma2 = 1,

    # OU parameters
    alpha = 1,     # strength of attraction
    theta = 0,     # optimum

    # Trait options
    n_traits = 1,
    trait_names = NULL,

    # Measurement error
    sd_error = 0,

    # Tree simulation
    tree_method = "pbtree",  # or "rtree"

    # Output
    return_tree = TRUE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  model <- match.arg(model)

  # Load required packages
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is required.")
  }

  # Generate tree if not supplied
  if (is.null(tree)) {
    if (tree_method == "pbtree") {
      tree <- phytools::pbtree(n = n_species)
    } else {
      tree <- ape::rtree(n = n_species)
    }
  }

  n_species <- length(tree$tip.label)

  # Trait names
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  trait_data <- matrix(NA, nrow = n_species, ncol = n_traits)
  colnames(trait_data) <- trait_names
  rownames(trait_data) <- tree$tip.label

  for (j in seq_len(n_traits)) {

    if (model == "BM") {

      trait <- phytools::fastBM(tree, sig2 = sigma2)

    } else if (model == "OU") {

      trait <- phytools::fastBM(
        tree,
        sig2 = sigma2,
        alpha = alpha,
        theta = theta
      )
    }

    # Add measurement error if requested
    if (sd_error > 0) {
      trait <- trait + stats::rnorm(n_species, 0, sd_error)
    }

    trait_data[, j] <- trait
  }

  # Convert to data frame
  df <- data.frame(
    Species = tree$tip.label,
    trait_data,
    stringsAsFactors = FALSE
  )

  # Return structure
  if (return_tree) {
    return(list(
      data = df,
      tree = tree
    ))
  } else {
    return(df)
  }
}

# A function for simulating community matrixes
simulate_community_matrix <- function(
    n_sites = 20,
    n_species = 30,

    # Environmental gradient
    env_values = seq(0, 1, length.out = n_sites),

    # Species niche parameters
    niche_optima = NULL,      # species optima along gradient
    niche_breadth = 0.2,      # tolerance (sd of response curve)

    # Abundance parameters
    max_abundance = 100,
    sd_species_intercept = 1,

    # Distribution choice
    family = c("poisson", "nbinom"),
    size = 10,   # for negative binomial

    # Zero inflation
    zero_inflation = 0,

    # Output
    return_long = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  family <- match.arg(family)

  # Environmental gradient
  env <- env_values

  # Species-specific parameters
  if (is.null(niche_optima)) {
    niche_optima <- stats::runif(n_species, min(env), max(env))
  }

  species_intercepts <- stats::rnorm(n_species, 0, sd_species_intercept)

  # Storage matrix
  comm_matrix <- matrix(0, nrow = n_sites, ncol = n_species)
  rownames(comm_matrix) <- paste0("Site_", seq_len(n_sites))
  colnames(comm_matrix) <- paste0("Sp_", seq_len(n_species))

  for (i in seq_len(n_sites)) {
    for (j in seq_len(n_species)) {

      # Gaussian niche response
      mu <- max_abundance * exp(
        - (env[i] - niche_optima[j])^2 / (2 * niche_breadth^2)
      )

      # Add species-specific abundance variation
      mu <- mu * exp(species_intercepts[j])

      # Draw counts
      if (family == "poisson") {
        abundance <- stats::rpois(1, lambda = mu)
      } else {
        abundance <- stats::rnbinom(1, mu = mu, size = size)
      }

      # Zero inflation
      if (zero_inflation > 0) {
        if (stats::runif(1) < zero_inflation) {
          abundance <- 0
        }
      }

      comm_matrix[i, j] <- abundance
    }
  }

  if (return_long) {

    df <- as.data.frame(as.table(comm_matrix))
    colnames(df) <- c("Site", "Species", "Abundance")
    df$Env <- rep(env, times = n_species)

    return(df)

  } else {
    return(comm_matrix)
  }
}


# Function for simulating multivariate traits
simulate_multivariate_traits <- function(
    n_obs = 100,
    n_traits = 3,

    # Trait names
    trait_names = NULL,

    # Mean structure
    intercepts = NULL,        # vector length n_traits
    env_effects = NULL,       # vector length n_traits
    treatment_effects = NULL, # vector length n_traits

    # Predictors
    env = rnorm(n_obs),
    treatment = rep(c("control", "treatment"), length.out = n_obs),

    # Covariance structure
    cov_matrix = NULL,
    sd_traits = NULL,
    cor_matrix = NULL,

    # Random effects (optional grouping)
    n_groups = 1,
    sd_group = 0,

    # Output
    return_cov = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # Trait names
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  # Defaults for fixed effects
  if (is.null(intercepts)) {
    intercepts <- rep(0, n_traits)
  }
  if (is.null(env_effects)) {
    env_effects <- rep(0, n_traits)
  }
  if (is.null(treatment_effects)) {
    treatment_effects <- rep(0, n_traits)
  }

  # Build covariance matrix
  if (is.null(cov_matrix)) {

    if (is.null(sd_traits)) {
      sd_traits <- rep(1, n_traits)
    }

    if (is.null(cor_matrix)) {
      cor_matrix <- diag(n_traits)
    }

    cov_matrix <- diag(sd_traits) %*% cor_matrix %*% diag(sd_traits)
  }

  # Assign groups
  group <- rep(seq_len(n_groups), length.out = n_obs)
  group_effects <- matrix(
    stats::rnorm(n_groups * n_traits, 0, sd_group),
    nrow = n_groups,
    ncol = n_traits
  )

  # Storage
  Y <- matrix(NA, nrow = n_obs, ncol = n_traits)

  for (i in seq_len(n_obs)) {

    trt_indicator <- ifelse(treatment[i] == "treatment", 1, 0)

    # Mean vector for observation i
    mu <- intercepts +
      env_effects * env[i] +
      treatment_effects * trt_indicator +
      group_effects[group[i], ]

    # Draw multivariate observation
    Y[i, ] <- MASS::mvrnorm(1, mu = mu, Sigma = cov_matrix)
  }

  colnames(Y) <- trait_names

  df <- data.frame(
    ID = seq_len(n_obs),
    Group = group,
    Env = env,
    Treatment = treatment,
    Y,
    stringsAsFactors = FALSE
  )

  if (return_cov) {
    return(list(
      data = df,
      covariance = cov_matrix
    ))
  } else {
    return(df)
  }
}

# Function that simulates spatial data
simulate_spatial_data <- function(
    n_points = 100,

    # Spatial domain
    x_range = c(0, 1),
    y_range = c(0, 1),

    # Generate coordinates or supply them
    coords = NULL,

    # Fixed effects
    intercept = 0,
    env_effect = 1,
    treatment_effect = 0,

    # Environmental covariate
    env = NULL,

    # Treatment
    treatment = rep(c("control", "treatment"), length.out = n_points),

    # Spatial autocorrelation parameters
    spatial_sd = 1,        # variance of spatial effect
    spatial_range = 0.2,   # how quickly correlation decays

    # Residual error
    sd_residual = 1,

    # Optional grouping
    n_groups = 1,
    sd_group = 0,

    # Output
    return_components = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # Coordinates
  if (is.null(coords)) {
    x <- stats::runif(n_points, x_range[1], x_range[2])
    y <- stats::runif(n_points, y_range[1], y_range[2])
    coords <- cbind(x, y)
  } else {
    n_points <- nrow(coords)
    x <- coords[, 1]
    y <- coords[, 2]
  }

  # Environmental variable
  if (is.null(env)) {
    env <- scale(x + y)  # simple spatial gradient
  }

  # Distance matrix
  dist_mat <- as.matrix(stats::dist(coords))

  # Exponential covariance function
  cov_mat <- spatial_sd^2 * exp(-dist_mat / spatial_range)

  # Spatial random effect
  spatial_effect <- as.numeric(
    MASS::mvrnorm(1, mu = rep(0, n_points), Sigma = cov_mat)
  )

  # Group effects
  group <- rep(seq_len(n_groups), length.out = n_points)
  group_effects <- stats::rnorm(n_groups, 0, sd_group)

  # Treatment indicator
  trt_indicator <- ifelse(treatment == "treatment", 1, 0)

  # Linear predictor
  eta <- intercept +
    env_effect * env +
    treatment_effect * trt_indicator +
    spatial_effect +
    group_effects[group]

  # Response
  y_resp <- stats::rnorm(n_points, mean = eta, sd = sd_residual)

  df <- data.frame(
    X = x,
    Y = y,
    Env = as.numeric(env),
    Treatment = treatment,
    Group = group,
    Response = y_resp,
    stringsAsFactors = FALSE
  )

  if (return_components) {
    return(list(
      data = df,
      spatial_effect = spatial_effect,
      covariance_matrix = cov_mat
    ))
  } else {
    return(df)
  }
}

# Function for simulating datasets with zero-inflated data
simulate_zero_inflated_data <- function(
    n = 100,

    # Predictors
    env = NULL,
    treatment = rep(c("control", "treatment"), length.out = n),

    # Linear predictor (count model, log scale)
    count_intercept = 1,
    env_effect_count = 0.5,
    treatment_effect_count = 0,

    # Linear predictor (zero-inflation model, logit scale)
    zi_intercept = -1,
    env_effect_zi = 0,
    treatment_effect_zi = 0,

    # Distribution for counts
    family = c("poisson", "nbinom"),
    size = 10,

    # Optional grouping
    n_groups = 1,
    sd_group_count = 0,
    sd_group_zi = 0,

    # Residual/env
    sd_env = 1,

    # Output
    return_components = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  family <- match.arg(family)

  # Environmental covariate
  if (is.null(env)) {
    env <- rnorm(n)
  }

  # Group assignment
  group <- rep(seq_len(n_groups), length.out = n)

  # Group-level random effects
  group_count_effect <- stats::rnorm(n_groups, 0, sd_group_count)
  group_zi_effect <- stats::rnorm(n_groups, 0, sd_group_zi)

  # Treatment indicator
  trt_indicator <- ifelse(treatment == "treatment", 1, 0)

  # Linear predictors
  eta_count <- count_intercept +
    env_effect_count * env +
    treatment_effect_count * trt_indicator +
    group_count_effect[group]

  eta_zi <- zi_intercept +
    env_effect_zi * env +
    treatment_effect_zi * trt_indicator +
    group_zi_effect[group]

  # Convert to probabilities
  mu <- exp(eta_count)
  zi_prob <- plogis(eta_zi)

  # Storage
  y <- numeric(n)
  is_zero_inflated <- logical(n)

  for (i in seq_len(n)) {

    # Structural zero?
    if (runif(1) < zi_prob[i]) {
      y[i] <- 0
      is_zero_inflated[i] <- TRUE
    } else {

      # Count component
      if (family == "poisson") {
        y[i] <- rpois(1, lambda = mu[i])
      } else {
        y[i] <- rnbinom(1, mu = mu[i], size = size)
      }

      is_zero_inflated[i] <- FALSE
    }
  }

  df <- data.frame(
    Y = y,
    Env = env,
    Treatment = treatment,
    Group = group,
    ZeroInflated = is_zero_inflated,
    stringsAsFactors = FALSE
  )

  if (return_components) {
    return(list(
      data = df,
      mu = mu,
      zi_prob = zi_prob
    ))
  } else {
    return(df)
  }
}

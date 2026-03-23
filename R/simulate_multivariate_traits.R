#' Simulate multivariate trait data
#'
#' Generates multivariate normally distributed traits with optional
#' environmental, treatment, and group-level effects, and a user-defined
#' covariance structure.
#'
#' @param n_obs Number of observations.
#' @param n_traits Number of traits.
#' @param trait_names Optional trait names.
#' @param intercepts Trait-specific intercepts.
#' @param env_effects Environmental effects per trait.
#' @param treatment_effects Treatment effects per trait.
#' @param env Environmental predictor values.
#' @param treatment Treatment labels.
#' @param cov_matrix Covariance matrix.
#' @param sd_traits Trait standard deviations.
#' @param cor_matrix Correlation matrix.
#' @param n_groups Number of groups.
#' @param sd_group SD of group-level effects.
#' @param return_cov Logical; return covariance matrix.
#' @param seed Optional random seed.
#'
#' @return A data frame or list containing simulated traits.
#' @export
#'
#' @examples
#' dat <- simulate_multivariate_traits(seed = 1)


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

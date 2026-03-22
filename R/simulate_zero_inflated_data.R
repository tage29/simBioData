
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

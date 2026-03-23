#' Simulate zero-inflated count data
#'
#' Generates count data with excess zeros using a two-part model:
#' a count process (Poisson or negative binomial) and a zero-inflation process.
#'
#' @param n Number of observations.
#' @param env Environmental predictor.
#' @param treatment Treatment labels.
#' @param count_intercept Intercept for count model (log scale).
#' @param env_effect_count Environmental effect on count model.
#' @param treatment_effect_count Treatment effect on count model.
#' @param zi_intercept Intercept for zero-inflation model (logit scale).
#' @param env_effect_zi Environmental effect on zero inflation.
#' @param treatment_effect_zi Treatment effect on zero inflation.
#' @param family Distribution ("poisson" or "nbinom").
#' @param size Dispersion parameter for negative binomial.
#' @param n_groups Number of groups.
#' @param sd_group_count SD of group effects (count model).
#' @param sd_group_zi SD of group effects (zero-inflation model).
#' @param sd_env SD of environmental variation.
#' @param return_components Logical; return model components.
#' @param seed Optional random seed.
#'
#' @return A data frame or list containing simulated zero-inflated data.
#' @export
#'
#' @examples
#' dat <- simulate_zero_inflated_data(seed = 1)


simulate_zero_inflated_data <- function(
    n = 100,

    env = NULL,
    treatment = rep(c("control", "treatment"), length.out = n),

    count_intercept = 1,
    env_effect_count = 0.5,
    treatment_effect_count = 0,

    zi_intercept = -1,
    env_effect_zi = 0,
    treatment_effect_zi = 0,

    family = c("poisson", "nbinom"),
    size = 10,

    n_groups = 1,
    sd_group_count = 0,
    sd_group_zi = 0,

    sd_env = 1,

    return_components = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  family <- match.arg(family)

  if (is.null(env)) {
    env <- stats::rnorm(n)
  }

  group <- rep(seq_len(n_groups), length.out = n)

  group_count_effect <- stats::rnorm(n_groups, 0, sd_group_count)
  group_zi_effect <- stats::rnorm(n_groups, 0, sd_group_zi)

  trt_indicator <- ifelse(treatment == "treatment", 1, 0)

  eta_count <- count_intercept +
    env_effect_count * env +
    treatment_effect_count * trt_indicator +
    group_count_effect[group]

  eta_zi <- zi_intercept +
    env_effect_zi * env +
    treatment_effect_zi * trt_indicator +
    group_zi_effect[group]

  mu <- exp(eta_count)
  zi_prob <- stats::plogis(eta_zi)

  y <- numeric(n)
  is_zero_inflated <- logical(n)

  for (i in seq_len(n)) {

    if (stats::runif(1) < zi_prob[i]) {
      y[i] <- 0
      is_zero_inflated[i] <- TRUE
    } else {

      if (family == "poisson") {
        y[i] <- stats::rpois(1, lambda = mu[i])
      } else {
        y[i] <- stats::rnbinom(1, mu = mu[i], size = size)
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

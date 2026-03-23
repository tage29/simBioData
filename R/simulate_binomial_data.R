#' Simulate binomial (success/failure) data
#'
#' Generates binomial data representing outcomes such as survival or success,
#' including treatment effects, environmental gradients, and random effects
#' at species and group levels.
#'
#' @param n_groups Number of groups (e.g. sites).
#' @param n_per_group Number of observations or trials per group.
#' @param species Character vector of species names.
#' @param intercept_logit Baseline intercept on the logit scale.
#' @param treatment_effect Treatment effect on the logit scale.
#' @param env_slope Effect of environmental predictor.
#' @param env_values Numeric vector of environmental values.
#' @param sd_species_intercept Standard deviation of species intercepts.
#' @param sd_species_slope Standard deviation of species slopes.
#' @param sd_group_intercept Standard deviation of group intercepts.
#' @param treatment_labels Labels for treatment levels.
#' @param aggregated Logical; return aggregated binomial counts if TRUE.
#' @param return_list Logical; return a list instead of a data frame.
#' @param seed Optional random seed.
#'
#' @return A data frame (or list) of simulated binomial data.
#' @export
#'
#' @examples
#' dat <- simulate_binomial_data(seed = 1)


simulate_binomial_data <- function(
    n_groups = 10,
    n_per_group = 50,
    species = c("sp1", "sp2"),

    intercept_logit = stats::qlogis(0.6),
    treatment_effect = stats::qlogis(0.7) - stats::qlogis(0.6),
    env_slope = 0,

    env_values = seq_len(n_groups),

    sd_species_intercept = 0.5,
    sd_species_slope = 0.2,
    sd_group_intercept = 0.5,

    treatment_labels = c("control", "treatment"),
    aggregated = FALSE,

    return_list = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  datalist <- list()
  counter <- 1

  species_intercepts <- stats::rnorm(length(species), 0, sd_species_intercept)
  species_slopes <- stats::rnorm(length(species), 0, sd_species_slope)
  names(species_intercepts) <- species
  names(species_slopes) <- species

  for (g in seq_len(n_groups)) {

    group_intercept <- stats::rnorm(1, 0, sd_group_intercept)
    env_val <- env_values[g]

    for (sp in species) {
      for (trt in treatment_labels) {

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

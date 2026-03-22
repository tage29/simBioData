

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

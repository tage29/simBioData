
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

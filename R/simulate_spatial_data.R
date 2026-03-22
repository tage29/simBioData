
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

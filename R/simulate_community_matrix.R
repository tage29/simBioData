
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

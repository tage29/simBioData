
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

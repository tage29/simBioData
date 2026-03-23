#' Simulate seed production along an elevation gradient (negative binomial)
#'
#' Generates overdispersed count data representing seed production across an
#' elevation gradient, with species- and site-level random effects and an
#' optional treatment effect (e.g. hand vs natural pollination).
#'
#' @param elevation_seq Numeric vector of elevation values.
#' @param n_per_type Number of observations per treatment per species and elevation.
#' @param species Character vector of species names.
#' @param intercept_log Baseline intercept on the log scale.
#' @param elevation_slope Effect of elevation on the log scale.
#' @param treatment_effect Treatment effect (log scale).
#' @param reference_elevation Reference value for centering elevation.
#' @param sd_species_intercept Standard deviation of species random intercepts.
#' @param sd_species_slope Standard deviation of species-specific slopes.
#' @param sd_site_intercept Standard deviation of site-level random intercepts.
#' @param size Dispersion parameter for the negative binomial distribution.
#' @param return_list Logical; return a list instead of a single data frame.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A data frame or list of data frames containing simulated seed counts.
#' @export
#'
#' @examples
#' dat <- simulate_seed_data_nb(seed = 1)


# Function for simulating seed data along an elevation gradient.
simulate_seed_data_nb <- function(
    elevation_seq = seq(1100, 2000, by = 100),
    n_per_type = 50,
    species = c("sp1", "sp2"),

    # Fixed effects (on log scale)
    intercept_log = log(100),
    elevation_slope = -0.001,
    treatment_effect = log(0.95),
    reference_elevation = 1100,

    # Random effects (SDs on log scale)
    sd_species_intercept = 0.2,
    sd_species_slope = 0.0005,
    sd_site_intercept = 0.15,

    # Dispersion
    size = 10,

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

  for (elev in elevation_seq) {

    site_intercept <- stats::rnorm(1, 0, sd_site_intercept)

    for (sp in species) {

      elev_centered <- elev - reference_elevation

      eta_natural <- intercept_log +
        (elevation_slope + species_slopes[sp]) * elev_centered +
        species_intercepts[sp] +
        site_intercept

      eta_hand <- eta_natural + treatment_effect

      mu_natural <- exp(eta_natural)
      mu_hand <- exp(eta_hand)

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

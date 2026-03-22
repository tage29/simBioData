

# Function for generating phylogenetic trait data.
simulate_phylo_traits <- function(
    n_species = 50,
    tree = NULL,

    # Model choice
    model = c("BM", "OU"),

    # Brownian motion parameter
    sigma2 = 1,

    # OU parameters
    alpha = 1,     # strength of attraction
    theta = 0,     # optimum

    # Trait options
    n_traits = 1,
    trait_names = NULL,

    # Measurement error
    sd_error = 0,

    # Tree simulation
    tree_method = "pbtree",  # or "rtree"

    # Output
    return_tree = TRUE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  model <- match.arg(model)

  # Load required packages
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is required.")
  }

  # Generate tree if not supplied
  if (is.null(tree)) {
    if (tree_method == "pbtree") {
      tree <- phytools::pbtree(n = n_species)
    } else {
      tree <- ape::rtree(n = n_species)
    }
  }

  n_species <- length(tree$tip.label)

  # Trait names
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  trait_data <- matrix(NA, nrow = n_species, ncol = n_traits)
  colnames(trait_data) <- trait_names
  rownames(trait_data) <- tree$tip.label

  for (j in seq_len(n_traits)) {

    if (model == "BM") {

      trait <- phytools::fastBM(tree, sig2 = sigma2)

    } else if (model == "OU") {

      trait <- phytools::fastBM(
        tree,
        sig2 = sigma2,
        alpha = alpha,
        theta = theta
      )
    }

    # Add measurement error if requested
    if (sd_error > 0) {
      trait <- trait + stats::rnorm(n_species, 0, sd_error)
    }

    trait_data[, j] <- trait
  }

  # Convert to data frame
  df <- data.frame(
    Species = tree$tip.label,
    trait_data,
    stringsAsFactors = FALSE
  )

  # Return structure
  if (return_tree) {
    return(list(
      data = df,
      tree = tree
    ))
  } else {
    return(df)
  }
}

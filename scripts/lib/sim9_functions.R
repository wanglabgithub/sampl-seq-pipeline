# Calculate the binary distance between every unique pair of ASVs
pairwise_binary_distance <- function(binary_matrix, threads = 1L) {
  # Transpose the Particle-by-ASV matrix and calculate distances between ASVs
  distance_matrix <- as.matrix(parallelDist::parDist(
    t(binary_matrix), method = "binary", threads = threads
  ))

  # Select the lower triangle to keep each unordered ASV pair only once
  positions <- which(lower.tri(distance_matrix), arr.ind = TRUE)

  # Return the ASV pair names and their binary distance as a tidy table
  tibble::tibble(
    ASVset1 = rownames(distance_matrix)[positions[, 1]],
    ASVset2 = colnames(distance_matrix)[positions[, 2]],
    Distance = distance_matrix[positions]
  )
}

# Run the SIM9 null model and summarize the result for every ASV pair
run_sim9 <- function(binary_matrix, random_communities = 50L,
                     swaps = 25000L, cores = 16L, seed = 123L) {
  # Confirm that the matrix and SIM9 settings have valid values
  stopifnot(is.matrix(binary_matrix), random_communities > 1L, swaps > 0L, cores > 0L)

  # Calculate the observed binary distance for each ASV pair
  observed <- pairwise_binary_distance(binary_matrix, threads = cores) %>%
    dplyr::rename(ObservedDistance = Distance)

  # Create and register a parallel worker cluster for random communities
  cluster <- parallel::makeCluster(min(cores, random_communities))
  doParallel::registerDoParallel(cluster)

  # Stop the worker cluster when this function finishes or encounters an error
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  # Generate SIM9 random communities reproducibly and combine their distances
  simulations <- foreach::foreach(
    iteration = seq_len(random_communities),
    .combine = dplyr::bind_rows,
    .multicombine = TRUE,
    .packages = c("dplyr", "EcoSimR", "parallelDist", "tibble"),
    .export = "pairwise_binary_distance",
    .options.RNG = seed
  ) %dorng% {
    # Start each randomization from the observed binary matrix
    randomized <- binary_matrix

    # Randomize the matrix while preserving its row and column totals
    for (swap_index in seq_len(swaps)) {
      randomized <- EcoSimR::sim9_single(randomized)
    }

    # Calculate ASV-pair distances and record the random-community iteration
    pairwise_binary_distance(randomized, threads = 1L) %>%
      dplyr::mutate(Iteration = iteration)
  }

  # Summarize simulated distances, join observed distances, and calculate statistics
  simulations %>%
    dplyr::group_by(ASVset1, ASVset2) %>%
    dplyr::summarise(
      MeanSimulatedDistance = mean(Distance),
      SDSimulatedDistance = stats::sd(Distance),
      MedianSimulatedDistance = stats::median(Distance),
      .groups = "drop"
    ) %>%
    # Add the distance calculated from the observed community
    dplyr::left_join(observed, by = c("ASVset1", "ASVset2")) %>%
    # Calculate the Z-score, two-sided P-value, BH FDR, and significance flag
    dplyr::mutate(
      Zscore = dplyr::if_else(
        SDSimulatedDistance > 0,
        (MeanSimulatedDistance - ObservedDistance) / SDSimulatedDistance,
        NA_real_
      ),
      PValue = 2 * stats::pnorm(abs(Zscore), lower.tail = FALSE),
      FDR = stats::p.adjust(PValue, method = "BH"),
      Significant = !is.na(FDR) & FDR < 0.05
    ) %>%
    # Return ASV pairs in a deterministic order
    dplyr::arrange(ASVset1, ASVset2)
}

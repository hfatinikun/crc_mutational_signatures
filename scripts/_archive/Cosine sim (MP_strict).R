source("utils.R")

## === CHOOSE WHICH COSMIC FILE TO USE (local file, same folder as script) ===
## Just change the filename on the next line, or uncomment one of the examples.
cosmic_file <- "COSMIC_v3.2_SBS_GRCh38.txt"
#cosmic_file <- "COSMIC_v3_SBS_GRCh38.txt"
#cosmic_file <- "COSMIC_v3.2_SBS_GRCh37.txt"
#cosmic_file <- "COSMIC_v3.4_SBS_GRCh38.txt"
#cosmic_file <- "COSMIC_v2_SBS_GRCh38.txt"
#cosmic_file <- "COSMIC_v2_SBS_GRCh37.txt"
#cosmic_file <- "COSMIC_v1_SBS_GRCh38.txt"

## Load the selected COSMIC signatures file
cosmic_data <- load_cosmic_signatures(cosmic_file, cosmic_dir, base_dir)


## Your three input data frames
# Load WES data
CRC_WES <- load_crc_data("WES", input_dir, base_dir)
# Load WGS data
CRC_WGS <- load_crc_data("WGS", input_dir, base_dir)
# Load TSO-500 data
CRC_TSO_500 <- load_crc_data("TSO-500", input_dir, base_dir)

# Strict refit using MutationalPatterns (backwards by default)
fit_crc_with_MutationalPatterns_strict <- function(crc_data, cosmic_data,
                                                   dataset_type = c("WES","WGS","TSO-500"),
                                                   max_delta = 0.004,
                                                   method = c("backwards","best_subset")) {
  dataset_type <- match.arg(dataset_type)
  method <- match.arg(method)
  message(sprintf("[MutationalPatterns strict] Fitting %s (method=%s, max_delta=%.3f)",
                  dataset_type, method, max_delta))
  # mut_mat: contexts x samples, signatures: contexts x signatures
  strict_refit <- MutationalPatterns::fit_to_signatures_strict(
    mut_mat   = as.matrix(crc_data),
    signatures = as.matrix(cosmic_data),
    max_delta = max_delta,
    method    = method
  )
  # Return the same shape you expect elsewhere
  list(
    contribution  = strict_refit$fit_res$contribution,   # signatures × samples
    reconstructed = strict_refit$fit_res$reconstructed,  # contexts × samples
    method        = paste0("MutationalPatterns_strict_", method),
    dataset_type  = dataset_type,
    sim_decay_fig = strict_refit$sim_decay_fig            # list of ggplots
  )
}

print("---> Mutational Signatures Analysis with MutationalPatterns <---")
# Fit the CRC_WES data
fitting_result_WES <- fit_crc_with_MutationalPatterns_strict(CRC_WES, cosmic_data, "WES")
# Fit the CRC_WGS data
fitting_result_WGS <- fit_crc_with_MutationalPatterns_strict(CRC_WGS, cosmic_data, "WGS")
# Fit the CRC_TSO_500 data
fitting_result_TSO_500 <- fit_crc_with_MutationalPatterns_strict(CRC_TSO_500, cosmic_data, "TSO-500")


# Normalize the contributions
CRC_WES_cosmic_MP_norm <- normalize_signature_contributions(fitting_result_WES)
CRC_WGS_cosmic_MP_norm <- normalize_signature_contributions(fitting_result_WGS)
CRC_TSO_500_cosmic_MP_norm <- normalize_signature_contributions(fitting_result_TSO_500)


print("--> Cosine similarity Analysis")
# Calculate cosine similarities
cosine_data_type <- calculate_cosine_similarity(fitting_result_WES, 
                                                fitting_result_TSO_500, 
                                                fitting_result_WGS, 
                                                CRC_WES, CRC_TSO_500, CRC_WGS)


print("--> Figure 2A (display)")
# Plot the Cosine Similarity data
p <- plot_cosine_similarity(cosine_data_type)

# Display the plot
print(p)


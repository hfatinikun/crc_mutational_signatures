## ---------- utils.R ----------
## Helper library for mutational signature workflows
## Dependencies used by main Rscripts
#cat("\014")

# Loading fittting tools libraries
library(MutationalPatterns)
library(deconstructSigs)
library(signature.tools.lib)
library(siglasso)
library(sigfit)
library(sigminer)
library(SigsPack)
library(YAPSA)

# Loading other dependecies
library(ggplot2)             # plotting
library(ggh4x)               # panel size helper
library(grid)                # unit()
library(dplyr)
library(ggpubr)

#Name of files and directory for inputs
base_dir   <- getwd()
cosmic_dir <- "cosmic_files"
input_dir  <- "input_files"

#Function to load COSMIC files
load_cosmic_signatures <- function(cosmic_file, cosmic_dir, base_dir) {
  # Create full path using file.path() for safety (handles spaces and platform differences)
  full_path <- file.path(base_dir, cosmic_dir, cosmic_file)
  # Read the file
  cosmic_data <- read.table(full_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # Convert to numeric matrix
  cosmic_data <- as.matrix(cosmic_data)
  return(cosmic_data)
}


#Function to load CRC files
load_crc_data <- function(dataset_type = c("WES", "WGS", "TSO-500"),
                          input_dir, base_dir) {
  # Match the dataset type to the correct filename
  dataset_type <- match.arg(dataset_type)  # Ensures valid input
  file_map <- list(
    WES = "Preclinical_Dataset_WES.txt",
    WGS = "Preclinical_Dataset_WGS.txt",
    `TSO-500`  = "Preclinical_Dataset_TSO-500.txt"
  )
  # Build the full path to the file
  full_path <- file.path(base_dir, input_dir, file_map[[dataset_type]])
  # Read the file
  crc_data <- read.table(full_path, header = TRUE, row.names = 1)
  return(crc_data)
}


#Function to fit CRC data with MutationalPatterns
fit_crc_with_MutationalPatterns <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[MutationalPatterns] Fitting", dataset_type))
  
  fitting_result <- fit_to_signatures(crc_data, cosmic_data)
  return(fitting_result)
}


#Function to normalize the contributions
normalize_signature_contributions <- function(fitting_result) {
  contribution <- fitting_result$contribution
  normalized <- contribution / colSums(contribution)[col(contribution)]
  return(as.data.frame(normalized))
}


#Function for Cosine similarity Analysis
calculate_cosine_similarity <- function(fitting_result_WES, fitting_result_TSO_500, fitting_result_WGS,
                                        CRC_WES, CRC_TSO_500, CRC_WGS) {
  
  # Calculate cosine similarities for each dataset
  cosine_sim_CRC_WES_cosmic_MP <- diag(cos_sim_matrix(fitting_result_WES$reconstructed, CRC_WES))
  cosine_sim_CRC_TSO500_cosmic_MP <- diag(cos_sim_matrix(fitting_result_TSO_500$reconstructed, CRC_TSO_500))
  cosine_sim_CRC_WGS_cosmic_MP    <- diag(cos_sim_matrix(fitting_result_WGS$reconstructed, CRC_WGS))
  
  # Create data frames for each dataset type
  cosine_data_type_WES <- data.frame(cosine_sim_CRC_WES_cosmic_MP, "WES")
  colnames(cosine_data_type_WES) <- c("Cosine_Similarity", "Data_Type")
  
  cosine_data_type_TSO500 <- data.frame(cosine_sim_CRC_TSO500_cosmic_MP, "TSO-500")
  colnames(cosine_data_type_TSO500) <- c("Cosine_Similarity", "Data_Type")
  
  cosine_data_type_WGS <- data.frame(cosine_sim_CRC_WGS_cosmic_MP, "WGS")
  colnames(cosine_data_type_WGS) <- c("Cosine_Similarity", "Data_Type")
  
  # Combine all the data frames into one
  cosine_data_type <- rbind(cosine_data_type_WES, cosine_data_type_TSO500, cosine_data_type_WGS)
  colnames(cosine_data_type) <- c("Cosine_Similarity", "Data_Type")
  
  return(cosine_data_type)
}

# Extended cosine similarity function for merging results
# Improved version with Sample column
calculate_cosine_similarity_extended <- function(fitting_result_WES, fitting_result_TSO_500, fitting_result_WGS,
                                                 CRC_WES, CRC_TSO_500, CRC_WGS, tool_name, cosmic_version) {
  
  make_cos_df <- function(reconstructed, observed, data_type) {
    
    reconstructed <- as.matrix(reconstructed)
    observed      <- as.matrix(observed)
    
    # 1) Align contexts (rows)
    common_ctx <- intersect(rownames(reconstructed), rownames(observed))
    if (length(common_ctx) == 0L) stop(paste("No overlapping contexts for", data_type))
    reconstructed <- reconstructed[common_ctx, , drop = FALSE]
    observed      <- observed[common_ctx, , drop = FALSE]
    
    # 2) Align samples (columns)
    common_samples <- intersect(colnames(reconstructed), colnames(observed))
    if (length(common_samples) == 0L) stop(paste("No overlapping samples for", data_type))
    reconstructed <- reconstructed[, common_samples, drop = FALSE]
    observed      <- observed[, common_samples, drop = FALSE]
    
    # 3) Cosine similarities (per-sample diagonal)
    cm <- cos_sim_matrix(reconstructed, observed)
    cos_vec <- diag(cm)
    
    # 4) IMPORTANT: attach sample names so data.frame() works
    names(cos_vec) <- colnames(observed)
    
    data.frame(
      Sample = names(cos_vec),
      Cosine_Similarity = as.numeric(cos_vec),
      Data_Type = data_type,
      Tool = tool_name,
      COSMIC_Version = cosmic_version,
      stringsAsFactors = FALSE
    )
  }
  
  df_WES <- make_cos_df(fitting_result_WES$reconstructed,     CRC_WES,     "WES")
  df_TSO <- make_cos_df(fitting_result_TSO_500$reconstructed, CRC_TSO_500, "TSO-500")
  df_WGS <- make_cos_df(fitting_result_WGS$reconstructed,     CRC_WGS,     "WGS")
  
  bind_rows(df_WES, df_TSO, df_WGS)
}



#Function to Create the Boxplot for Cosine Similarity
plot_cosine_similarity <- function(cosine_data, title = "Cosine Similarity v3.2 (sigfit)") {
  
  # Create the ggplot for the boxplot
  p <- ggplot(cosine_data,
              aes(x = factor(Data_Type, levels = c("WGS", "WES", "TSO-500")),
                  y = Cosine_Similarity)) +
    geom_boxplot(color = "black", fill = "grey") +
    coord_flip() +
    labs(x = "NGS Workflow", y = "Cosine Similarity") +
    ggtitle(title) +
    scale_y_continuous(breaks = seq(0.9, 1, 0.025),
                       limits = c(0.89, 1),
                       expand = c(0, 0)) +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    force_panelsizes(rows = unit(6, "cm"), cols = unit(8, "cm")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text     = element_text(size = 10),
          panel.border  = element_rect(color = "black", linewidth = 0.4))
  
  return(p)
}


#--------------DeconstructSigs-------------------------

#Function to fit CRC data with deconstructSigs
# Fit CRC data with deconstructSigs and return contributions + reconstructed profiles
fit_crc_with_deconstructSigs <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[deconstructSigs] Fitting", dataset_type))
  
  crc_t    <- as.data.frame(t(crc_data))      # samples × contexts
  cosmic_t <- as.data.frame(t(cosmic_data))   # signatures × contexts
  common   <- intersect(colnames(crc_t), colnames(cosmic_t))
  crc_t    <- crc_t[, common]
  cosmic_t <- cosmic_t[, common]
  
  res_list <- lapply(rownames(crc_t), function(s) {
    cat("Processing sample", s, "\n")
    out <- deconstructSigs::whichSignatures(
      tumor.ref        = crc_t,
      sample.id        = s,
      signatures.ref   = cosmic_t,
      contexts.needed  = TRUE,
      signature.cutoff = 0
    )
    w <- out$weights[1, , drop = FALSE]
    rownames(w) <- s
    w
  })
  
  weights <- do.call(rbind, res_list)         # samples × signatures
  contrib <- t(weights)                       # signatures × samples
  recon   <- t(cosmic_t) %*% contrib          # contexts × samples
  
  list(contribution = contrib,
       reconstructed = recon,
       method = "deconstructSigs",
       dataset_type = dataset_type)
}

#----------------------------signature.lib.tools-----------------------------------
#Function to fit CRC data with signature.lib.tools
fit_crc_with_signature_tools_lib <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[signature.tools.lib] Fitting", dataset_type))
  
  # crc_data: contexts × samples
  # cosmic_data: contexts × signatures
  res <- signature.tools.lib::Fit(
    catalogues          = as.matrix(crc_data),
    signatures          = as.matrix(cosmic_data),
    exposureFilterType  = "giniScaledThreshold",
    useBootstrap        = TRUE,
    nboot               = 100,
    nparallel           = 1
  )
  
  # Exposures -> contributions (signatures × samples), drop "unassigned"
  w <- t(res$exposures)  # signatures × samples
  if ("unassigned" %in% rownames(w)) {
    w <- w[setdiff(rownames(w), "unassigned"), , drop = FALSE]
  }
  
  # Ensure the signature order matches between cosmic_data (columns) and w (rows)
  common_sigs <- intersect(colnames(cosmic_data), rownames(w))
  if (length(common_sigs) == 0) {
    stop("No overlapping signatures between cosmic_data and exposures.")
  }
  cosmic_use <- as.matrix(cosmic_data[, common_sigs, drop = FALSE])      # contexts × signatures
  contribution <- as.matrix(w[common_sigs, , drop = FALSE])              # signatures × samples
  
  # Reconstruct catalogues (contexts × samples)
  reconstructed <- cosmic_use %*% contribution
  
  list(
    contribution  = contribution,      # signatures × samples
    reconstructed = reconstructed,     # contexts × samples
    method        = "signature.tools.lib",
    dataset_type  = dataset_type
  )
}


#----------------sigLASSO----------------
fit_crc_with_sigLASSO <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  message(paste("[sigLASSO] Fitting", dataset_type))
  
  # Ensure matrices
  spectrum  <- as.matrix(crc_data)    # contexts × samples
  signature <- as.matrix(cosmic_data) # contexts × signatures
  
  # 1) Align mutation contexts (rows)
  common_ctx <- intersect(rownames(spectrum), rownames(signature))
  if (length(common_ctx) == 0L) stop("No overlapping mutation contexts between CRC and COSMIC.")
  spectrum  <- spectrum [common_ctx, , drop = FALSE]
  signature <- signature[common_ctx, , drop = FALSE]
  
  # 2) Run sigLASSO
  # (set.seed if you want reproducibility)
  res <- siglasso::siglasso(spectrum, signature = signature)  # returns signatures × samples
  
  # 3) Coerce to numeric matrix (no transpose!)
  contribution <- as.matrix(res)                # signatures × samples
  
  # 4) Align signatures by name (defensive)
  common_sigs <- intersect(colnames(signature), rownames(contribution))
  if (length(common_sigs) == 0L) stop("No overlapping signature names between result and COSMIC.")
  signature    <- signature[, common_sigs, drop = FALSE]
  contribution <- contribution[common_sigs, , drop = FALSE]
  
  # 5) Reconstruct catalogues (contexts × samples)
  reconstructed <- signature %*% contribution
  
  list(
    contribution  = contribution,   # signatures × samples
    reconstructed = reconstructed,  # contexts × samples
    method        = "sigLASSO",
    dataset_type  = dataset_type
  )
}


# ============================
# Extra fitters: sigfit / sigminer / SigsPack / YAPSA
# ============================

# Helper: align contexts (rows) between CRC and COSMIC matrices
align_contexts <- function(crc_data, cosmic_data) {
  crc <- as.matrix(crc_data)          # contexts × samples
  cos <- as.matrix(cosmic_data)       # contexts × signatures
  common_ctx <- intersect(rownames(crc), rownames(cos))
  if (length(common_ctx) == 0L) stop("No overlapping mutation contexts between CRC and COSMIC.")
  crc <- crc[common_ctx, , drop = FALSE]
  cos <- cos[common_ctx, , drop = FALSE]
  list(crc = crc, cos = cos)
}

# ---------------- sigminer (NNLS) ----------------
fit_crc_with_sigminer <- function(crc_data, cosmic_data,
                                  dataset_type = c("WES", "WGS", "TSO-500"),
                                  method = "NNLS") {
  dataset_type <- match.arg(dataset_type)
  message(paste("[sigminer]", method, "Fitting", dataset_type))
  
  # sigminer expects mutation matrix with features in rows, samples in columns (96 × N)
  aligned <- align_contexts(crc_data, cosmic_data)
  mut_mat <- aligned$crc
  cosmic  <- aligned$cos
  
  # sigminer::sig_fit returns a matrix/data.frame of exposures (signatures × samples)
  res <- sigminer::sig_fit(mut_mat, cosmic, method = method)
  contribution <- as.matrix(res)  # signatures × samples
  
  # Defensive alignment by signature names (if present)
  common_sigs <- intersect(colnames(cosmic), rownames(contribution))
  if (length(common_sigs) > 0) {
    cosmic       <- cosmic[, common_sigs, drop = FALSE]
    contribution <- contribution[common_sigs, , drop = FALSE]
  }
  
  reconstructed <- cosmic %*% contribution  # contexts × samples
  
  list(
    contribution  = contribution,
    reconstructed = reconstructed,
    method        = paste0("sigminer_", method),
    dataset_type  = dataset_type
  )
}

# ---------------- SigsPack ----------------
fit_crc_with_SigsPack <- function(crc_data, cosmic_data,
                                  dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  message(paste("[SigsPack] Fitting", dataset_type))
  
  aligned <- align_contexts(crc_data, cosmic_data)
  mut_mat <- aligned$crc          # contexts × samples
  cosmic  <- aligned$cos          # contexts × signatures
  
  res <- SigsPack::signature_exposure(mut_mat, P = cosmic)
  exp_raw <- as.matrix(res$exposures)
  
  # Coerce to signatures × samples
  if (nrow(exp_raw) == ncol(mut_mat) && ncol(exp_raw) == ncol(cosmic)) {
    contribution <- t(exp_raw)
  } else if (ncol(exp_raw) == ncol(mut_mat) && nrow(exp_raw) == ncol(cosmic)) {
    contribution <- exp_raw
  } else {
    stop("SigsPack exposures dims don't match expected (signatures×samples).")
  }
  
  # FORCE names (this is the key)
  rownames(contribution) <- colnames(cosmic)      # signatures
  colnames(contribution) <- colnames(mut_mat)     # samples
  
  reconstructed <- cosmic %*% contribution
  
  # FORCE reconstructed sample names too
  colnames(reconstructed) <- colnames(mut_mat)
  rownames(reconstructed) <- rownames(mut_mat)
  
  list(
    contribution  = contribution,
    reconstructed = reconstructed,
    method        = "SigsPack",
    dataset_type  = dataset_type
  )
}


# ---------------- YAPSA (LCD) ----------------
fit_crc_with_YAPSA <- function(crc_data, cosmic_data,
                               dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  message(paste("[YAPSA] Fitting", dataset_type))
  
  aligned <- align_contexts(crc_data, cosmic_data)
  mut_mat <- aligned$crc
  cosmic  <- aligned$cos
  
  # YAPSA::LCD expects (mutation spectra) and (signature profiles)
  # Often expects: mut_mat = 96×N, cosmic = 96×K
  res <- YAPSA::LCD(mut_mat, cosmic)
  
  # LCD returns exposures-like output (format can vary); try common patterns:
  contribution <- as.matrix(res)
  if (!all(dim(contribution) == c(ncol(cosmic), ncol(mut_mat))) &&
      all(dim(contribution) == c(ncol(mut_mat), ncol(cosmic)))) {
    contribution <- t(contribution)
  }
  
  # Try to add names if missing
  if (is.null(colnames(contribution))) colnames(contribution) <- colnames(mut_mat)
  if (is.null(rownames(contribution))) rownames(contribution) <- colnames(cosmic)
  
  common_sigs <- intersect(colnames(cosmic), rownames(contribution))
  if (length(common_sigs) > 0) {
    cosmic       <- cosmic[, common_sigs, drop = FALSE]
    contribution <- contribution[common_sigs, , drop = FALSE]
  }
  
  reconstructed <- cosmic %*% contribution
  
  list(
    contribution  = contribution,
    reconstructed = reconstructed,
    method        = "YAPSA_LCD",
    dataset_type  = dataset_type
  )
}

# ---------------- sigfit ----------------
fit_crc_with_sigfit <- function(context_file, cosmic_data,
                                dataset_type = c("WES", "WGS", "TSO-500"),
                                chains = 1, seed = 0, hpd_prob = 0.90) {
  dataset_type <- match.arg(dataset_type)
  message(paste("[sigfit] Fitting", dataset_type))
  
  # context_file should be the "long" table your template uses (sigfit::build_catalogues needs it)
  counts <- sigfit::build_catalogues(context_file)  # contexts × samples
  
  # Align contexts between counts and COSMIC
  aligned <- align_contexts(counts, cosmic_data)
  counts <- aligned$crc
  cosmic <- aligned$cos
  
  fit <- sigfit::fit_signatures(counts = counts, signatures = cosmic, chains = chains, seed = seed)
  
  exposures <- sigfit::retrieve_pars(fit, par = "exposures", hpd_prob = hpd_prob)
  
  # exposures$mean is typically samples × signatures OR signatures × samples depending on sigfit version.
  exp_mean <- as.matrix(exposures$mean)
  
  # We want signatures × samples
  if (nrow(exp_mean) == ncol(counts) && ncol(exp_mean) == ncol(cosmic)) {
    contribution <- t(exp_mean)
    colnames(contribution) <- colnames(counts)
    rownames(contribution) <- colnames(cosmic)
  } else {
    contribution <- exp_mean
  }
  
  # Align by signature names if possible
  common_sigs <- intersect(colnames(cosmic), rownames(contribution))
  if (length(common_sigs) > 0) {
    cosmic       <- cosmic[, common_sigs, drop = FALSE]
    contribution <- contribution[common_sigs, , drop = FALSE]
  }
  
  reconstructed <- cosmic %*% contribution
  
  list(
    contribution  = contribution,
    reconstructed = reconstructed,
    method        = "sigfit",
    dataset_type  = dataset_type,
    sigfit_fit    = fit,
    sigfit_hpd    = exposures
  )
}

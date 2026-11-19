#Function to fit CRC data with deconstructSigs
# deconstructSigs fitter (simple + compatible with your pipeline)
fit_crc_with_deconstructSigs <- function(crc_data, cosmic_data,
                                         dataset_type = c("WES","WGS","TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  message(sprintf("[deconstructSigs] Fitting %s", dataset_type))
  
  # tumor.ref: samples x contexts
  tumor_df <- as.data.frame(t(crc_data), check.names = FALSE)
  
  # signatures.ref: signatures x contexts (transpose + relabel)
  signatures_ref <- as.matrix(t(cosmic_data))
  colnames(signatures_ref) <- rownames(cosmic_data)
  
  # align by contexts
  common_ctx <- intersect(colnames(tumor_df), colnames(signatures_ref))
  if (length(common_ctx) == 0) {
    stop("No overlapping mutation contexts between CRC and COSMIC. Ensure both are SBS96.")
  }
  common_ctx <- sort(common_ctx)
  tumor_df       <- tumor_df[, common_ctx, drop = FALSE]
  signatures_ref <- signatures_ref[, common_ctx, drop = FALSE]
  
  # outputs
  exposure <- matrix(0, nrow = nrow(signatures_ref), ncol = nrow(tumor_df),
                     dimnames = list(rownames(signatures_ref), rownames(tumor_df)))
  reconstructed <- matrix(0, nrow = nrow(crc_data), ncol = ncol(crc_data),
                          dimnames = list(rownames(crc_data), colnames(crc_data)))
  
  for (s in rownames(tumor_df)) {
    res <- deconstructSigs::whichSignatures(
      tumor.ref       = tumor_df,
      sample.id       = s,
      signatures.ref  = signatures_ref,
      contexts.needed = FALSE,
      signature.cutoff = 0
    )
    w <- res$weights
    if (is.data.frame(w)) w <- w[1, , drop = TRUE]
    w <- w[rownames(signatures_ref)]
    exposure[, s] <- as.numeric(w)
    
    recon_common <- t(signatures_ref[, common_ctx, drop = FALSE]) %*% exposure[, s]
    reconstructed[common_ctx, s] <- as.numeric(recon_common)
  }
  
  list(contribution = exposure, reconstructed = reconstructed)
}






#Function to fit CRC data with signature.tools.lib
fit_crc_with_signaturetools <- function(crc_data, cosmic_data,
                                        dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[signature.tools.lib] Fitting", dataset_type))
  
  fit <- signature.tools.lib::SignatureFit(
    catalogues = as.matrix(crc_data),
    signatures = as.matrix(cosmic_data),
    ncpu = 1, useBootstrap = FALSE
  )
  exposures <- fit$Exposure      # signatures x samples
  reconstructed <- as.matrix(cosmic_data) %*% exposures
  
  return(list(contribution = exposures, reconstructed = reconstructed))
}


#Function to fit CRC data with SigProfilerAssignment
fit_crc_with_SigProfilerAssignment <- function(crc_data, cosmic_data,
                                               dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[SigProfilerAssignment] Fitting", dataset_type))
  
  exposures <- SigProfilerAssignment::fit_to_signatures(
    catalogue_matrix = as.matrix(crc_data),
    signature_matrix = as.matrix(cosmic_data)
  )
  reconstructed <- as.matrix(cosmic_data) %*% exposures
  
  return(list(contribution = exposures, reconstructed = reconstructed))
}






fit_crc_with_SignatureAnalyzer <- function(crc_data, cosmic_data,
                                           dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[SignatureAnalyzer] Fitting (NNLS assignment)", dataset_type))
  
  exposures <- nnls::nnls(as.matrix(cosmic_data), as.matrix(crc_data))$x
  exposures <- matrix(exposures, nrow = ncol(cosmic_data), ncol = ncol(crc_data))
  colnames(exposures) <- colnames(crc_data)
  rownames(exposures) <- colnames(cosmic_data)
  
  reconstructed <- as.matrix(cosmic_data) %*% exposures
  
  return(list(contribution = exposures, reconstructed = reconstructed))
}
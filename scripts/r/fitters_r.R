## Uniform R fitters: each returns list(contribution, reconstructed)
## contribution: signatures × samples
## reconstructed: contexts × samples

align_contexts_local <- function(crc_data, cosmic_data) {
  crc <- as.matrix(crc_data)
  cos <- as.matrix(cosmic_data)
  common_ctx <- rownames(cos)[rownames(cos) %in% rownames(crc)]
  if (length(common_ctx) == 0L) stop("No overlapping mutation contexts.")
  list(crc = crc[common_ctx, , drop = FALSE], cos = cos[common_ctx, , drop = FALSE])
}

fit_MutationalPatterns <- function(catalogue, cosmic) {
  if (!requireNamespace("MutationalPatterns", quietly = TRUE)) {
    stop("Package MutationalPatterns not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  res <- MutationalPatterns::fit_to_signatures(al$crc, al$cos)
  list(contribution = as.matrix(res$contribution),
       reconstructed = as.matrix(res$reconstructed))
}

fit_deconstructSigs <- function(catalogue, cosmic) {
  if (!requireNamespace("deconstructSigs", quietly = TRUE)) {
    stop("Package deconstructSigs not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  crc_t <- as.data.frame(t(al$crc))
  cosmic_t <- as.data.frame(t(al$cos))
  weights <- lapply(rownames(crc_t), function(s) {
    out <- deconstructSigs::whichSignatures(
      tumor.ref = crc_t,
      sample.id = s,
      signatures.ref = cosmic_t,
      contexts.needed = TRUE,
      signature.cutoff = 0
    )
    w <- out$weights[1, , drop = FALSE]
    rownames(w) <- s
    w
  })
  wmat <- do.call(rbind, weights)   # samples × signatures
  contrib <- t(wmat)                # signatures × samples
  recon <- al$cos %*% contrib
  list(contribution = contrib, reconstructed = recon)
}

fit_signature_tools_lib <- function(catalogue, cosmic) {
  if (!requireNamespace("signature.tools.lib", quietly = TRUE)) {
    stop("Package signature.tools.lib not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  # Bootstrap can be slow; use lighter defaults for full grid reproducibility
  res <- signature.tools.lib::Fit(
    catalogues = al$crc,
    signatures = al$cos,
    exposureFilterType = "giniScaledThreshold",
    useBootstrap = FALSE,
    nparallel = 1
  )
  w <- as.matrix(t(res$exposures))
  if ("unassigned" %in% rownames(w)) {
    w <- w[setdiff(rownames(w), "unassigned"), , drop = FALSE]
  }
  common_sigs <- intersect(colnames(al$cos), rownames(w))
  w <- w[common_sigs, , drop = FALSE]
  list(contribution = w, reconstructed = al$cos[, common_sigs, drop = FALSE] %*% w)
}

fit_sigLASSO <- function(catalogue, cosmic) {
  if (!requireNamespace("siglasso", quietly = TRUE)) {
    stop("Package siglasso not installed")
  }
  # siglasso calls nnls()/glmnet()/brewer.pal() unqualified
  for (pkg in c("nnls", "glmnet", "RColorBrewer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package ", pkg, " required by sigLASSO is not installed")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  al <- align_contexts_local(catalogue, cosmic)
  res <- siglasso::siglasso(al$crc, signature = al$cos)
  contrib <- as.matrix(res)
  common_sigs <- intersect(colnames(al$cos), rownames(contrib))
  if (!length(common_sigs)) {
    if (nrow(contrib) == ncol(al$cos)) {
      rownames(contrib) <- colnames(al$cos)
      common_sigs <- colnames(al$cos)
    } else {
      stop("sigLASSO contribution signature names do not match COSMIC")
    }
  }
  contrib <- contrib[common_sigs, , drop = FALSE]
  list(
    contribution = contrib,
    reconstructed = al$cos[, common_sigs, drop = FALSE] %*% contrib
  )
}

fit_sigminer <- function(catalogue, cosmic) {
  if (!requireNamespace("sigminer", quietly = TRUE)) {
    stop("Package sigminer not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  res <- sigminer::sig_fit(al$crc, al$cos, method = "NNLS")
  contrib <- as.matrix(res)
  list(contribution = contrib, reconstructed = al$cos %*% contrib)
}

fit_SigsPack <- function(catalogue, cosmic) {
  if (!requireNamespace("SigsPack", quietly = TRUE)) {
    stop("Package SigsPack not installed")
  }
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("Package nnls required for SigsPack fallback is not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  # Default decomposeQP fails when signatures are linearly dependent
  # (e.g. COSMIC v3.6 has K>96). Fall back to NNLS in that case.
  decompose_nnls <- function(m, P, ...) {
    as.numeric(nnls::nnls(as.matrix(P), as.numeric(m))$x)
  }
  res <- tryCatch(
    SigsPack::signature_exposure(al$crc, P = al$cos),
    error = function(e) {
      message("  [SigsPack] QP failed (", conditionMessage(e), "); using NNLS fallback")
      SigsPack::signature_exposure(al$crc, P = al$cos, FUN = decompose_nnls)
    }
  )
  exp_raw <- as.matrix(res$exposures)
  if (nrow(exp_raw) == ncol(al$crc) && ncol(exp_raw) == ncol(al$cos)) {
    contrib <- t(exp_raw)
  } else if (ncol(exp_raw) == ncol(al$crc) && nrow(exp_raw) == ncol(al$cos)) {
    contrib <- exp_raw
  } else {
    stop("SigsPack exposures dims unexpected")
  }
  rownames(contrib) <- colnames(al$cos)
  colnames(contrib) <- colnames(al$crc)
  list(contribution = contrib, reconstructed = al$cos %*% contrib)
}

fit_YAPSA <- function(catalogue, cosmic) {
  if (!requireNamespace("YAPSA", quietly = TRUE)) {
    stop("Package YAPSA not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  res <- YAPSA::LCD(al$crc, al$cos)
  contrib <- as.matrix(res)
  if (nrow(contrib) == ncol(al$crc) && ncol(contrib) == ncol(al$cos)) {
    contrib <- t(contrib)
  }
  if (is.null(rownames(contrib))) rownames(contrib) <- colnames(al$cos)
  if (is.null(colnames(contrib))) colnames(contrib) <- colnames(al$crc)
  list(contribution = contrib, reconstructed = al$cos %*% contrib)
}

fit_sigfit <- function(catalogue, cosmic, chains = 1L, iter = 1000L, seed = 1L) {
  if (!requireNamespace("sigfit", quietly = TRUE)) {
    stop("Package sigfit not installed")
  }
  al <- align_contexts_local(catalogue, cosmic)
  # sigfit expects counts: samples × contexts, signatures: signatures × contexts
  counts_sx <- t(al$crc)
  sigs_sx <- t(al$cos)
  fit <- sigfit::fit_signatures(
    counts = counts_sx,
    signatures = sigs_sx,
    iter = iter,
    chains = chains,
    seed = seed
  )
  exposures <- sigfit::retrieve_pars(fit, par = "exposures")
  exp_mean <- as.matrix(exposures$mean)
  # typically samples × signatures
  if (nrow(exp_mean) == nrow(counts_sx) && ncol(exp_mean) == nrow(sigs_sx)) {
    contrib <- t(exp_mean)
    rownames(contrib) <- rownames(sigs_sx)
    colnames(contrib) <- rownames(counts_sx)
  } else if (nrow(exp_mean) == nrow(sigs_sx) && ncol(exp_mean) == nrow(counts_sx)) {
    contrib <- exp_mean
    rownames(contrib) <- rownames(sigs_sx)
    colnames(contrib) <- rownames(counts_sx)
  } else {
    stop("sigfit exposures dims unexpected: ", paste(dim(exp_mean), collapse = "x"))
  }
  list(contribution = contrib, reconstructed = al$cos %*% contrib)
}

FITTERS_R <- list(
  MutationalPatterns = fit_MutationalPatterns,
  deconstructSigs = fit_deconstructSigs,
  `signature.tools.lib` = fit_signature_tools_lib,
  sigLASSO = fit_sigLASSO,
  sigfit = fit_sigfit,
  sigminer = fit_sigminer,
  SigsPack = fit_SigsPack,
  YAPSA = fit_YAPSA
)

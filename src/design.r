# ==============================
# Archivo: design.R
# ==============================

# Generador de primers qPCR

generar_qPCR <- function(seqs, fs, fe, rs, re, kr, gmin, gmax, mode, tmn, tmx) {
  resb <- list(); reso <- list(); N <- length(seqs)
  pb <- txtProgressBar(min = 0, max = N, style = 3)

  for (i in seq_along(seqs)) {
    nm <- names(seqs)[i]; s <- as.character(seqs[[nm]])
    if (nchar(s) < re) { setTxtProgressBar(pb, i); next }

    fwR <- substr(s, fs, fe); rvR <- substr(s, rs, re)
    fw_c <- obtener_candidatos(fwR, kr, gmin, gmax)
    rv_c <- obtener_candidatos(rvR, kr, gmin, gmax)

    best <- fallback <- NULL; candidates <- list()
    for (fw in fw_c) {
      fchk <- verificar_primer(fw, seqs)
      for (rv in rv_c) {
        size <- re - fs + nchar(rv)
        if (size < 80 || size > 150) next
        rchk <- verificar_primer(rv, seqs)
        valid <- (fchk$total == 1 && grepl(extraer_iso(nm), fchk$isoformas) &&
                  rchk$total == 1 && grepl(extraer_iso(nm), rchk$isoformas))
        primer <- list(forward = fw, reverse = rv, product_size = size)
        if (mode == "specificity" && valid) { best <- primer; break }
        if (is.null(fallback)) fallback <- primer
        if (length(candidates) < 10) candidates[[length(candidates) + 1]] <- primer
      }
      if (!is.null(best)) break
    }

    resb[[nm]] <- best %||% fallback
    reso[[nm]] <- candidates
    setTxtProgressBar(pb, i)
  }
  close(pb)
  list(best = resb, other = reso)
}

# Generador de primers RACE

generar_RACE <- function(seqs, w, kr, gmin, gmax, mode, tmn, tmx) {
  resb <- list(); reso <- list(); N <- length(seqs)
  pb <- txtProgressBar(min = 0, max = N, style = 3)

  for (i in seq_along(seqs)) {
    nm <- names(seqs)[i]; s <- as.character(seqs[[nm]])
    if (nchar(s) < 2 * w) { setTxtProgressBar(pb, i); next }

    fwR <- substr(s, 1, w); rvR <- substr(s, nchar(s) - w + 1, nchar(s))
    fw_c <- obtener_candidatos(fwR, kr, gmin, gmax)
    rv_c <- obtener_candidatos(rvR, kr, gmin, gmax)

    best <- fallback <- NULL; candidates <- list()
    for (fw in fw_c) {
      fchk <- verificar_primer(fw, seqs)
      for (rv in rv_c) {
        rchk <- verificar_primer(rv, seqs)
        valid <- (fchk$total == 1 && grepl(extraer_iso(nm), fchk$isoformas) &&
                  rchk$total == 1 && grepl(extraer_iso(nm), rchk$isoformas))
        primer <- list(forward = fw, reverse = rv, product_size = NA)
        if (mode == "specificity" && valid) { best <- primer; break }
        if (is.null(fallback)) fallback <- primer
        if (length(candidates) < 10) candidates[[length(candidates) + 1]] <- primer
      }
      if (!is.null(best)) break
    }

    resb[[nm]] <- best %||% fallback
    reso[[nm]] <- candidates
    setTxtProgressBar(pb, i)
  }
  close(pb)
  list(best = resb, other = reso)
}

# Diseño por cobertura usando set-covering greedy

generar_coverage <- function(seqs, params_qPCR, params_RACE, primer_type) {
  iso_names <- names(seqs)
  if (primer_type == "qPCR") {
    fw_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_qPCR$k_range, function(k) extract_kmers(substr(s, params_qPCR$fw_start, params_qPCR$fw_end), k)))))
    rv_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_qPCR$k_range, function(k) extract_kmers(substr(s, params_qPCR$rev_start, params_qPCR$rev_end), k)))))
    gc_min <- params_qPCR$gc_min; gc_max <- params_qPCR$gc_max
  } else {
    fw_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_RACE$k_range, function(k) extract_kmers(substr(s, 1, params_RACE$window_size), k)))))
    rv_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_RACE$k_range, function(k) extract_kmers(substr(s, nchar(s) - params_RACE$window_size + 1, nchar(s)), k)))))
    gc_min <- params_RACE$gc_min; gc_max <- params_RACE$gc_max
  }

  fw_kmers <- unique(fw_pool); rv_kmers <- unique(rv_pool)
  fw_kmers <- fw_kmers[sapply(fw_kmers, function(k) {
    g <- calcular_gc(k); g >= gc_min && g <= gc_max
  })]
  rv_kmers <- rv_kmers[sapply(rv_kmers, function(k) {
    g <- calcular_gc(k); g >= gc_min && g <= gc_max
  })]

  fw_hits <- setNames(lapply(fw_kmers, function(k) names(seqs)[vcountPattern(k, seqs) > 0]), fw_kmers)
  rv_hits <- setNames(lapply(rv_kmers, function(k) names(seqs)[vcountPattern(k, seqs) > 0]), rv_kmers)

  candidates <- expand.grid(forward = fw_kmers, reverse = rv_kmers, stringsAsFactors = FALSE)
  pending <- iso_names; selected <- list()

  for (i in seq_len(nrow(candidates))) {
    pair <- candidates[i, ]
    common <- intersect(fw_hits[[pair$forward]], rv_hits[[pair$reverse]])
    cover <- intersect(common, pending)
    if (length(cover) > 0) {
      selected[[length(selected)+1]] <- list(
        forward  = pair$forward,
        reverse  = pair$reverse,
        isoforms = cover
      )
      pending <- setdiff(pending, cover)
      if (length(pending) == 0) break
    }
  }

  list(selected = selected, pending = pending)
}


# Fallback operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
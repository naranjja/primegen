# ======================================================
# PRIMEGEN
# ======================================================

# -------------------------------
# Cargar librerías
# -------------------------------
library(Biostrings)
library(openxlsx)

# -------------------------------
# CONFIGURACIÓN
# -------------------------------
cat("Seleccione el archivo FASTA con las isoformas:\n")
archivo_fasta <- "./avengers.fa"
seq_set <- readDNAStringSet(archivo_fasta)
cat("Se han leído", length(seq_set), "secuencias.\n")

# Tipo de primers: "qPCR", "RACE" o "coverage"
primer_type <- "RACE"   # Opciones: "qPCR", "RACE"

# Modo de generación: "specificity" o "coverage"
mode_text <- "coverage"   # Opciones: "specificity", "coverage"

# Parámetros para qPCR
autos_params_qPCR <- list(
  fw_start = 50, fw_end = 100,
  rev_start = 150, rev_end = 200,
  k_range = 18:25,
  gc_min = 50, gc_max = 70,
  prod_min = 80, prod_max = 150,
  tm_min = 57, tm_max = 62
)

# Parámetros para RACE
params_RACE <- list(
  window_size = 200,
  k_range     = 23:28,
  gc_min      = 50, gc_max = 70,
  tm_min      = 57, tm_max = 62
)

# -------------------------------
# Funciones auxiliares y de diseño
# -------------------------------

# Cálculo de GC (%)
calcular_gc <- function(s) {
  v <- strsplit(s, NULL)[[1]]
  round(sum(v %in% c("G", "C")) / length(v) * 100, 2)
}

# Extraer kmers de longitud k
extract_kmers <- function(s, k) {
  if (nchar(s) < k) return(character(0))
  substring(s, 1:(nchar(s) - k + 1), k:nchar(s))
}

# Funciones de especificidad (reverso-complemento, dimerización)
is_complement <- function(a, b) {
  (a == "A" && b == "T") || (a == "T" && b == "A") ||
    (a == "C" && b == "G") || (a == "G" && b == "C")
}

longest_complementary <- function(s, t) {
  sc <- strsplit(s, NULL)[[1]]; tc <- strsplit(t, NULL)[[1]]
  m <- length(sc); n <- length(tc)
  L <- matrix(0, m+1, n+1); ml <- 0
  for (i in seq_len(m)) for (j in seq_len(n)) {
    if (is_complement(sc[i], tc[j])) {
      L[i+1, j+1] <- L[i, j] + 1
      ml <- max(ml, L[i+1, j+1])
    }
  }
  ml
}

detect_hairpin <- function(p) {
  longest_complementary(p, as.character(reverseComplement(DNAString(p))))
}

detect_self_dimer <- function(p) longest_complementary(p, p)

detect_cross_dimer <- function(p1, p2) {
  longest_complementary(p1, as.character(reverseComplement(DNAString(p2))))
}

# Cálculo de Tm usando parámetros de vecindario dinucleótido
calcular_Tm <- function(primer, conc = 250e-9, Na = 50e-3) {
  nn <- list(
    "AA" = c(-7.9, -22.2), "TT" = c(-7.9, -22.2),
    "AT" = c(-7.2, -20.4), "TA" = c(-7.2, -21.3),
    "CA" = c(-8.5, -22.7), "TG" = c(-8.5, -22.7),
    "GT" = c(-8.4, -22.4), "AC" = c(-8.4, -22.4),
    "CT" = c(-7.8, -21.0), "AG" = c(-7.8, -21.0),
    "GA" = c(-8.2, -22.2), "TC" = c(-8.2, -22.2),
    "CG" = c(-10.6, -27.2), "GC" = c(-9.8, -24.4),
    "GG" = c(-8.0, -19.9)
  )
  p <- toupper(primer); n <- nchar(p)
  if (n < 2) stop("El primer debe tener al menos 2 nt")
  dp <- substring(p, 1:(n - 1), 2:n)
  vals <- sapply(dp, function(x) if (is.null(nn[[x]])) c(0, 0) else nn[[x]])
  dH <- sum(vals[1, ]) * 1000
  dS <- sum(vals[2, ])
  R  <- 1.987
  TmK <- dH / (dS + R * log(conc / 4))
  TmC <- TmK - 273.15 + 16.6 * log10(Na)
  round(TmC, 1)
}

# Funciones de especificidad adicionales
extraer_iso <- function(nm) sub(".*-([0-9]+).*", "\\1", nm)
verificar_primer <- function(primer, seqs) {
  tots <- vcountPattern(primer, DNAStringSet(seqs), fixed = TRUE)
  hits <- names(seqs)[tots > 0]
  list(total = sum(tots), isoformas = paste(hits, collapse = ","))
}
obtener_candidatos <- function(region, kr, gmin, gmax) {
  km <- unique(unlist(lapply(kr, function(k) extract_kmers(region, k))))
  gcv <- sapply(km, calcular_gc)
  km[gcv >= gmin & gcv <= gmax]
}

# Diseño por especificidad: generar_qPCR() y generar_RACE()
generar_qPCR <- function(seqs, fs, fe, rs, re, kr, gmin, gmax,
                         pmin, pmax, mode, tmn, tmx) {
  resb <- list(); reso <- list(); N <- length(seqs)
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in seq_along(seqs)) {
    nm <- names(seqs)[i]; s <- as.character(seqs[[nm]])
    if (nchar(s) < re) { resb[[nm]] <- NULL; setTxtProgressBar(pb, i); next }
    fwR <- substr(s, fs, fe); rvR <- substr(s, rs, re)
    fw0 <- obtener_candidatos(fwR, kr, gmin, gmax)
    rv0 <- obtener_candidatos(rvR, kr, gmin, gmax)
    fw_c <- fw0; rv_c <- rv0
    best <- NULL; fallback <- NULL
    for (fw in fw_c) {
      fchk <- verificar_primer(fw, seqs)
      for (rv in rv_c) {
        size <- re - fs + nchar(rv)
        if (size < pmin || size > pmax) next
        rchk <- verificar_primer(rv, seqs)
        if (mode == "specificity" && fchk$total == 1 &&
            grepl(extraer_iso(nm), fchk$isoformas) &&
            rchk$total == 1 && grepl(extraer_iso(nm), rchk$isoformas)) {
          best <- list(forward = fw, reverse = rv, product_size = size); break
        } else if (is.null(fallback)) {
          fallback <- list(forward = fw, reverse = rv, product_size = size)
        }
      }
      if (!is.null(best)) break
    }
    if (is.null(best)) best <- fallback
    resb[[nm]] <- best
    o <- list(); cnt <- 0
    for (fw in fw_c) {
      if (cnt >= 10) break
      for (rv in rv_c) {
        if (cnt >= 10) break
        size <- re - fs + nchar(rv)
        if (size < pmin || size > pmax) next
        o[[length(o) + 1]] <- list(forward = fw, reverse = rv, product_size = size)
        cnt <- cnt + 1
      }
    }
    reso[[nm]] <- o; setTxtProgressBar(pb, i)
  }
  close(pb)
  list(best = resb, other = reso)
}

generar_RACE <- function(seqs, w, kr, gmin, gmax, mode, tmn, tmx) {
  resb <- list(); reso <- list(); N <- length(seqs)
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in seq_along(seqs)) {
    nm <- names(seqs)[i]; s <- as.character(seqs[[nm]])
    if (nchar(s) < 2 * w) { resb[[nm]] <- NULL; setTxtProgressBar(pb, i); next }
    fwR <- substr(s, 1, w); rvR <- substr(s, nchar(s) - w + 1, nchar(s))
    fw0 <- obtener_candidatos(fwR, kr, gmin, gmax)
    rv0 <- obtener_candidatos(rvR, kr, gmin, gmax)
    fw_c <- fw0; rv_c <- rv0
    best <- NULL; fallback <- NULL
    for (fw in fw_c) {
      for (rv in rv_c) {
        fchk <- verificar_primer(fw, seqs)
        rchk <- verificar_primer(rv, seqs)
        if (mode == "specificity" && fchk$total == 1 &&
            grepl(extraer_iso(nm), fchk$isoformas) &&
            rchk$total == 1 && grepl(extraer_iso(nm), rchk$isoformas)) {
          best <- list(forward = fw, reverse = rv, product_size = NA); break
        } else if (is.null(fallback)) {
          fallback <- list(forward = fw, reverse = rv, product_size = NA)
        }
      }
      if (!is.null(best)) break
    }
    if (is.null(best)) best <- fallback
    resb[[nm]] <- best
    o <- list(); cnt <- 0
    for (fw in fw_c) {
      if (cnt >= 10) break
      for (rv in rv_c) {
        if (cnt >= 10) break
        o[[length(o) + 1]] <- list(forward = fw, reverse = rv, product_size = NA)
        cnt <- cnt + 1
      }
    }
    reso[[nm]] <- o; setTxtProgressBar(pb, i)
  }
  close(pb)
  list(best = resb, other = reso)
}

# Generar coverage (greedy set-cover)
generar_coverage <- function(seqs, params_qPCR, params_RACE, primer_type) {
  iso_names <- names(seqs)
  if (primer_type == "qPCR") {
    fw_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_qPCR$k_range, function(k) extract_kmers(
        substr(s, params_qPCR$fw_start, params_qPCR$fw_end), k)))))
    rv_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_qPCR$k_range, function(k) extract_kmers(
        substr(s, params_qPCR$rev_start, params_qPCR$rev_end), k)))))
    gc_min <- params_qPCR$gc_min; gc_max <- params_qPCR$gc_max
  } else {
    fw_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_RACE$k_range, function(k) extract_kmers(
        substr(s, 1, params_RACE$window_size), k)))))
    rv_pool <- unlist(lapply(as.character(seqs), function(s) unlist(
      lapply(params_RACE$k_range, function(k) extract_kmers(
        substr(s, nchar(s) - params_RACE$window_size + 1, nchar(s)), k)))))
    gc_min <- params_RACE$gc_min; gc_max <- params_RACE$gc_max
  }
  fw_kmers <- unique(fw_pool); rv_kmers <- unique(rv_pool)
  fw_kmers <- fw_kmers[sapply(fw_kmers, calcular_gc) >= gc_min & sapply(fw_kmers, calcular_gc) <= gc_max]
  rv_kmers <- rv_kmers[sapply(rv_kmers, calcular_gc) >= gc_min & sapply(rv_kmers, calcular_gc) <= gc_max]
  fw_hits <- setNames(lapply(fw_kmers, function(k) names(seqs)[vcountPattern(k, seqs) > 0]), fw_kmers)
  rv_hits <- setNames(lapply(rv_kmers, function(k) names(seqs)[vcountPattern(k, seqs) > 0]), rv_kmers)
  candidates <- expand.grid(forward=fw_kmers, reverse=rv_kmers, stringsAsFactors=FALSE)
  pending <- iso_names; selected <- list()
  for (i in seq_len(nrow(candidates))) {
    pair <- candidates[i,]
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

# -------------------------------
# Ejecución principal según tipo de diseño y modo
# -------------------------------
# Inicializar variables
archivo_salida <- NULL
resultados <- NULL
pendientes <- NULL

if (primer_type == "qPCR") {
  if (mode_text == "specificity") {
    cat("
--- Diseño qPCR (especificidad) ---
")
    resultados <- generar_qPCR(
      seq_set,
      autos_params_qPCR$fw_start, autos_params_qPCR$fw_end,
      autos_params_qPCR$rev_start, autos_params_qPCR$rev_end,
      autos_params_qPCR$k_range,
      autos_params_qPCR$gc_min, autos_params_qPCR$gc_max,
      mode_text,
      autos_params_qPCR$tm_min, autos_params_qPCR$tm_max
    )
    archivo_salida <- "resultados_qPCR_specificity.xlsx"
  } else if (mode_text == "coverage") {
    cat("
--- Diseño qPCR (cobertura) ---
")
    cov_res <- generar_coverage(
      seq_set,
      autos_params_qPCR,
      params_RACE,
      primer_type
    )
    resultados <- cov_res$selected
    pendientes <- cov_res$pending
    archivo_salida <- "resultados_qPCR_coverage.xlsx"
  }
} else if (primer_type == "RACE") {
  if (mode_text == "specificity") {
    cat("
--- Diseño RACE (especificidad) ---
")
    resultados <- generar_RACE(
      seq_set,
      params_RACE$window_size,
      params_RACE$k_range,
      params_RACE$gc_min, params_RACE$gc_max,
      mode_text,
      params_RACE$tm_min, params_RACE$tm_max
    )
    archivo_salida <- "resultados_RACE_specificity.xlsx"
  } else if (mode_text == "coverage") {
    cat("
--- Diseño RACE (cobertura) ---
")
    cov_res <- generar_coverage(
      seq_set,
      autos_params_qPCR,
      params_RACE,
      primer_type
    )
    resultados <- cov_res$selected
    pendientes <- cov_res$pending
    archivo_salida <- "resultados_RACE_coverage.xlsx"
  }
}

# -------------------------------
# Guardar resultados en Excel
# -------------------------------
wb <- createWorkbook()
if (mode_text == "specificity") {
  # resultado especificidad es lista con best y other
  tab_best <- do.call(rbind, lapply(names(resultados$best), function(nm) {
    p <- resultados$best[[nm]]
    if (is.null(p)) return(NULL)
    data.frame(
      Isoforma         = nm,
      Primer_Forward   = p$forward,
      Primer_Reverse   = p$reverse,
      GC_Forward       = calcular_gc(p$forward),
      GC_Reverse       = calcular_gc(p$reverse),
      Tm_Forward       = calcular_Tm(p$forward),
      Tm_Reverse       = calcular_Tm(p$reverse),
      Hairpin_Forward  = detect_hairpin(p$forward),
      Hairpin_Reverse  = detect_hairpin(p$reverse),
      SelfDimer_Forward= detect_self_dimer(p$forward),
      SelfDimer_Reverse= detect_self_dimer(p$reverse),
      CrossDimer       = detect_cross_dimer(p$forward, p$reverse),
      stringsAsFactors = FALSE
    )
  }))
  addWorksheet(wb, "Best_Primers")
  writeData(wb, "Best_Primers", tab_best)
  
  tab_other <- do.call(rbind, lapply(names(resultados$other), function(nm) {
    lst <- resultados$other[[nm]]
    if (length(lst) == 0) return(NULL)
    do.call(rbind, lapply(lst, function(p) {
      data.frame(
        Isoforma         = nm,
        Primer_Forward   = p$forward,
        Primer_Reverse   = p$reverse,
        GC_Forward       = calcular_gc(p$forward),
        GC_Reverse       = calcular_gc(p$reverse),
        Tm_Forward       = calcular_Tm(p$forward),
        Tm_Reverse       = calcular_Tm(p$reverse),
        Hairpin_Forward  = detect_hairpin(p$forward),
        Hairpin_Reverse  = detect_hairpin(p$reverse),
        SelfDimer_Forward= detect_self_dimer(p$forward),
        SelfDimer_Reverse= detect_self_dimer(p$reverse),
        CrossDimer       = detect_cross_dimer(p$forward, p$reverse),
        stringsAsFactors = FALSE
      )
    }))
  }))
  addWorksheet(wb, "Other_Candidates")
  writeData(wb, "Other_Candidates", tab_other)
} else if (mode_text == "coverage") {
  addWorksheet(wb, "Coverage_Primers")
  out <- do.call(rbind, lapply(resultados, function(x) {
    data.frame(
      Primer_Forward = x$forward,
      Primer_Reverse = x$reverse,
      Isoformas      = paste(x$isoforms, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  writeData(wb, "Coverage_Primers", out)
}

saveWorkbook(wb, archivo_salida, overwrite = TRUE)
cat("✅ Resultados guardados en:", archivo_salida, "")

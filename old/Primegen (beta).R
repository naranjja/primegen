# ======================================================
# PRIMEGEN - Generador de Primers para qPCR y RACE
# ======================================================

# -------------------------------
# Cargar librerías
# -------------------------------
library(Biostrings)
library(openxlsx)

# -------------------------------
# Configuración inicial
# -------------------------------
cat("Seleccione el archivo FASTA con las isoformas:\n")
archivo_fasta <- "./avengers-trimmed.fa"
seq_set <- readDNAStringSet(archivo_fasta)
cat("Se han leído", length(seq_set), "secuencias.\n")

# Tipo de primers: "qPCR", "RACE"
primer_type <- "RACE"

# Modo de generación: "specificity" o "coverage"
mode_text <- "specificity"

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
# Funciones auxiliares
# -------------------------------
calcular_gc <- function(s) {
  v <- strsplit(s, NULL)[[1]]
  round(sum(v %in% c("G", "C")) / length(v) * 100, 2)
}

extract_kmers <- function(s, k) {
  if (nchar(s) < k) return(character(0))
  substring(s, 1:(nchar(s) - k + 1), k:nchar(s))
}

is_complement <- function(a, b) {
  (a == "A" && b == "T") || (a == "T" && b == "A") ||
  (a == "C" && b == "G") || (a == "G" && b == "C")
}

longest_complementary <- function(s, t) {
  sc <- strsplit(s, NULL)[[1]]
  tc <- strsplit(t, NULL)[[1]]
  L <- matrix(0, length(sc)+1, length(tc)+1)
  ml <- 0
  for (i in seq_along(sc)) for (j in seq_along(tc)) {
    if (is_complement(sc[i], tc[j])) {
      L[i+1, j+1] <- L[i, j] + 1
      ml <- max(ml, L[i+1, j+1])
    }
  }
  ml
}

detect_hairpin <- function(p) longest_complementary(p, as.character(reverseComplement(DNAString(p))))
detect_self_dimer <- function(p) longest_complementary(p, p)
detect_cross_dimer <- function(p1, p2) longest_complementary(p1, as.character(reverseComplement(DNAString(p2))))

calcular_Tm <- function(primer, conc = 250e-9, Na = 50e-3) {
  nn <- list(
    "AA"=c(-7.9, -22.2), "TT"=c(-7.9, -22.2), "AT"=c(-7.2, -20.4), "TA"=c(-7.2, -21.3),
    "CA"=c(-8.5, -22.7), "TG"=c(-8.5, -22.7), "GT"=c(-8.4, -22.4), "AC"=c(-8.4, -22.4),
    "CT"=c(-7.8, -21.0), "AG"=c(-7.8, -21.0), "GA"=c(-8.2, -22.2), "TC"=c(-8.2, -22.2),
    "CG"=c(-10.6, -27.2), "GC"=c(-9.8, -24.4), "GG"=c(-8.0, -19.9)
  )
  p <- toupper(primer); n <- nchar(p)
  if (n < 2) stop("El primer debe tener al menos 2 nt")
  dp <- substring(p, 1:(n - 1), 2:n)
  vals <- sapply(dp, function(x) nn[[x]] %||% c(0, 0))
  dH <- sum(vals[1, ]) * 1000
  dS <- sum(vals[2, ])
  R <- 1.987
  TmK <- dH / (dS + R * log(conc / 4))
  round(TmK - 273.15 + 16.6 * log10(Na), 1)
}

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

# -------------------------------
# Funciones de diseño de primers
# -------------------------------
# (No changes made to generar_qPCR, generar_RACE, generar_coverage)
# Due to length, these functions are assumed unchanged and already structured.

# -------------------------------
# Ejecución según tipo y modo
# -------------------------------
archivo_salida <- NULL
resultados <- NULL
pendientes <- NULL

if (primer_type == "qPCR") {
  if (mode_text == "specificity") {
    cat("\n--- Diseño qPCR (especificidad) ---\n")
    resultados <- generar_qPCR(seq_set, autos_params_qPCR$fw_start, autos_params_qPCR$fw_end,
                               autos_params_qPCR$rev_start, autos_params_qPCR$rev_end,
                               autos_params_qPCR$k_range, autos_params_qPCR$gc_min, autos_params_qPCR$gc_max,
                               mode_text, autos_params_qPCR$tm_min, autos_params_qPCR$tm_max)
    archivo_salida <- "resultados_qPCR_specificity.xlsx"
  } else if (mode_text == "coverage") {
    cat("\n--- Diseño qPCR (cobertura) ---\n")
    cov_res <- generar_coverage(seq_set, autos_params_qPCR, params_RACE, primer_type)
    resultados <- cov_res$selected
    pendientes <- cov_res$pending
    archivo_salida <- "resultados_qPCR_coverage.xlsx"
  }
} else if (primer_type == "RACE") {
  if (mode_text == "specificity") {
    cat("\n--- Diseño RACE (especificidad) ---\n")
    resultados <- generar_RACE(seq_set, params_RACE$window_size, params_RACE$k_range,
                               params_RACE$gc_min, params_RACE$gc_max, mode_text,
                               params_RACE$tm_min, params_RACE$tm_max)
    archivo_salida <- "resultados_RACE_specificity.xlsx"
  } else if (mode_text == "coverage") {
    cat("\n--- Diseño RACE (cobertura) ---\n")
    cov_res <- generar_coverage(seq_set, autos_params_qPCR, params_RACE, primer_type)
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
  # Pestaña: Mejores primers
  tab_best <- do.call(rbind, lapply(names(resultados$best), function(nm) {
    p <- resultados$best[[nm]]
    if (is.null(p)) return(NULL)
    data.frame(
      Isoforma          = nm,
      Primer_Forward    = p$forward,
      Primer_Reverse    = p$reverse,
      GC_Forward        = calcular_gc(p$forward),
      GC_Reverse        = calcular_gc(p$reverse),
      Tm_Forward        = calcular_Tm(p$forward),
      Tm_Reverse        = calcular_Tm(p$reverse),
      Hairpin_Forward   = detect_hairpin(p$forward),
      Hairpin_Reverse   = detect_hairpin(p$reverse),
      SelfDimer_Forward = detect_self_dimer(p$forward),
      SelfDimer_Reverse = detect_self_dimer(p$reverse),
      CrossDimer        = detect_cross_dimer(p$forward, p$reverse),
      stringsAsFactors  = FALSE
    )
  }))
  addWorksheet(wb, "Best_Primers")
  writeData(wb, "Best_Primers", tab_best)

  # Pestaña: Otros candidatos
  tab_other <- do.call(rbind, lapply(names(resultados$other), function(nm) {
    lst <- resultados$other[[nm]]
    if (length(lst) == 0) return(NULL)
    do.call(rbind, lapply(lst, function(p) {
      data.frame(
        Isoforma          = nm,
        Primer_Forward    = p$forward,
        Primer_Reverse    = p$reverse,
        GC_Forward        = calcular_gc(p$forward),
        GC_Reverse        = calcular_gc(p$reverse),
        Tm_Forward        = calcular_Tm(p$forward),
        Tm_Reverse        = calcular_Tm(p$reverse),
        Hairpin_Forward   = detect_hairpin(p$forward),
        Hairpin_Reverse   = detect_hairpin(p$reverse),
        SelfDimer_Forward = detect_self_dimer(p$forward),
        SelfDimer_Reverse = detect_self_dimer(p$reverse),
        CrossDimer        = detect_cross_dimer(p$forward, p$reverse),
        stringsAsFactors  = FALSE
      )
    }))
  }))
  addWorksheet(wb, "Other_Candidates")
  writeData(wb, "Other_Candidates", tab_other)

} else if (mode_text == "coverage") {
  addWorksheet(wb, "Coverage_Primers")
  tab_coverage <- do.call(rbind, lapply(resultados, function(x) {
    data.frame(
      Primer_Forward = x$forward,
      Primer_Reverse = x$reverse,
      Isoformas      = paste(x$isoforms, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  writeData(wb, "Coverage_Primers", tab_coverage)
}

saveWorkbook(wb, archivo_salida, overwrite = TRUE)
cat("✅ Resultados guardados en:", archivo_salida, "\n")
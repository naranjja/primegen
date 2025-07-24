calcular_gc <- function(s) sum(strsplit(s, NULL)[[1]] %in% c("G", "C")) / nchar(s) * 100

extract_kmers <- function(s, k) {
  if (nchar(s) < k) return(character(0))
  substring(s, 1:(nchar(s) - k + 1), k:nchar(s))
}

is_complement <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")

longest_complementary <- function(s, t) {
  sc <- strsplit(s, NULL)[[1]]; tc <- strsplit(t, NULL)[[1]]
  L <- matrix(0, length(sc)+1, length(tc)+1)
  for (i in seq_along(sc)) for (j in seq_along(tc)) {
    if (is_complement(sc[i], tc[j])) L[i+1, j+1] <- L[i, j] + 1
  }
  max(L)
}

calcular_Tm <- function(primer, conc = 250e-9, Na = 50e-3) {
  nn <- list(
    AA = c(-7.9, -22.2), TT = c(-7.9, -22.2), AT = c(-7.2, -20.4), TA = c(-7.2, -21.3),
    CA = c(-8.5, -22.7), TG = c(-8.5, -22.7), GT = c(-8.4, -22.4), AC = c(-8.4, -22.4),
    CT = c(-7.8, -21.0), AG = c(-7.8, -21.0), GA = c(-8.2, -22.2), TC = c(-8.2, -22.2),
    CG = c(-10.6, -27.2), GC = c(-9.8, -24.4), GG = c(-8.0, -19.9)
  )
  p <- toupper(primer); n <- nchar(p)
  dp <- substring(p, 1:(n - 1), 2:n)
  vals <- sapply(dp, function(x) if (!is.null(nn[[x]])) nn[[x]] else c(0, 0))
  dH <- sum(vals[1, ]) * 1000; dS <- sum(vals[2, ])
  TmK <- dH / (dS + 1.987 * log(conc / 4))
  round(TmK - 273.15 + 16.6 * log10(Na), 1)
}

extraer_iso <- function(nm) sub(".*-([0-9]+).*", "\\1", nm)

verificar_primer <- function(primer, seqs) {
  tots <- vcountPattern(primer, seqs, fixed = TRUE)
  hits <- names(seqs)[tots > 0]
  list(total = sum(tots), isoformas = paste(hits, collapse = ","))
}

obtener_candidatos <- function(region, kr, gmin, gmax) {
  km <- unique(unlist(lapply(kr, function(k) extract_kmers(region, k))))
  gcv <- sapply(km, calcular_gc)
  km[gcv >= gmin & gcv <= gmax]
}
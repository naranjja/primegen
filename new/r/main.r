# Cargar librerías
library(Biostrings)
library(openxlsx)
source("config.R")
source("utils.R")
source("design.R")
source("output.R")

# Leer FASTA
seq_set <- readDNAStringSet(config$archivo_fasta)
cat("Se han leído", length(seq_set), "secuencias.\n")

# Ejecutar diseño
params_used <- config$params[[config$primer_type]]
cat("\n--- Diseño", config$primer_type, "(", config$mode_text, ") ---\n")

if (config$mode_text == "specificity") {
  resultados <- switch(config$primer_type,
    qPCR = generar_qPCR(seq_set, params_used$fw_start, params_used$fw_end,
                        params_used$rev_start, params_used$rev_end,
                        params_used$k_range, params_used$gc_min, params_used$gc_max,
                        config$mode_text, params_used$tm_min, params_used$tm_max),
    RACE = generar_RACE(seq_set, params_used$window_size, params_used$k_range,
                        params_used$gc_min, params_used$gc_max,
                        config$mode_text, params_used$tm_min, params_used$tm_max)
  )
  archivo_salida <- paste0("./output/resultados_", config$primer_type, "_specificity.xlsx")
} else {
  cov_res <- generar_coverage(seq_set, config$params$qPCR, config$params$RACE, config$primer_type)
  resultados <- cov_res$selected
  pendientes <- cov_res$pending
  archivo_salida <- paste0("./output/resultados_", config$primer_type, "_coverage.xlsx")
}

# Guardar resultados
guardar_resultados_excel(resultados, config$mode_text, archivo_salida)
guardar_resultados_excel <- function(resultados, mode_text, archivo_salida) {
  wb <- createWorkbook()
  if (mode_text == "specificity") {
    tabular <- function(res) do.call(rbind, lapply(names(res), function(nm) {
      p <- res[[nm]]; if (is.null(p)) return(NULL)
      data.frame(
        Isoforma = nm,
        Primer_Forward = p$forward, Primer_Reverse = p$reverse,
        GC_Forward = calcular_gc(p$forward), GC_Reverse = calcular_gc(p$reverse),
        Tm_Forward = calcular_Tm(p$forward), Tm_Reverse = calcular_Tm(p$reverse),
        Hairpin_Forward = longest_complementary(p$forward, as.character(reverseComplement(DNAString(p$forward)))),
        Hairpin_Reverse = longest_complementary(p$reverse, as.character(reverseComplement(DNAString(p$reverse)))),
        SelfDimer_Forward = longest_complementary(p$forward, p$forward),
        SelfDimer_Reverse = longest_complementary(p$reverse, p$reverse),
        CrossDimer = longest_complementary(p$forward, as.character(reverseComplement(DNAString(p$reverse)))),
        stringsAsFactors = FALSE
      )
    }))
    addWorksheet(wb, "Best_Primers")
    writeData(wb, "Best_Primers", tabular(resultados$best))
    addWorksheet(wb, "Other_Candidates")
    writeData(wb, "Other_Candidates", tabular(resultados$other))
  } else {
    tab_coverage <- do.call(rbind, lapply(resultados, function(x) {
      data.frame(
        Primer_Forward = x$forward,
        Primer_Reverse = x$reverse,
        Isoformas      = paste(x$isoforms, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))
    addWorksheet(wb, "Coverage_Primers")
    writeData(wb, "Coverage_Primers", tab_coverage)
  }
  saveWorkbook(wb, archivo_salida, overwrite = TRUE)
  cat("✅ Resultados guardados en:", archivo_salida, "\n")
}
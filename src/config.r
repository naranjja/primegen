config <- list(
  archivo_fasta = "./data/avengers-3.fa",
  primer_type = "RACE",       # "qPCR" o "RACE"
  mode_text = "coverage",  # "specificity" o "coverage"
  params = list(
    qPCR = list(
      fw_start = 50, fw_end = 100, rev_start = 150, rev_end = 200,
      k_range = 18:25, gc_min = 50, gc_max = 70, prod_min = 80, prod_max = 150,
      tm_min = 57, tm_max = 62
    ),
    RACE = list(
      window_size = 200, k_range = 23:28, gc_min = 50, gc_max = 70,
      tm_min = 57, tm_max = 62
    )
  )
)
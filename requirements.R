# requirements.R
install.packages("openxlsx", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("Biostrings")
# Install R dependencies for TCR_Evaluation
# Usage: Rscript R/install_deps.R

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c("openPrimeR", "Biostrings")
cran_packages <- c("ggplot2", "tidyr")

for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
}

for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

cat("\nAll R dependencies installed successfully.\n")
cat("Packages:", paste(bioc_packages, collapse = ", "),
    paste(cran_packages, collapse = ", "), "\n")

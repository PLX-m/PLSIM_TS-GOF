# Install the R packages required by the project scripts.
# Only missing packages will be installed; existing versions are left unchanged.

packages <- c(
  "splines2", "MASS", "openxlsx", "caret", "glmnet",
  "randomForest", "reticulate", "doParallel", "foreach",
  "doRNG", "ggplot2", "patchwork", "VGAMextra"
)

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}

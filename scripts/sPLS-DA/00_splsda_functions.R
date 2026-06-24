```r
############################################################
# 00_splsda_functions.R
#
# Helper functions for sPLS-DA analyses.
#
# This script does NOT run analyses by itself.
# It should be loaded with:
#
#   source("scripts/splsda/00_splsda_functions.R")
#
# from analysis-specific scripts, such as:
#
#   03_run_splsda_set.R
#
# Main functions:
#   - load_splsda_packages()
#   - normalize_text()
#   - prepare_X()
#   - add_multiclass_group()
#   - run_simple_splsda()
############################################################


############################################################
# 1. Load required packages
############################################################

load_splsda_packages <- function() {
  
  cran_packages <- c("ggplot2")
  
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install("mixOmics")
  }
  
  suppressPackageStartupMessages({
    library(mixOmics)
    library(ggplot2)
  })
}


############################################################
# 2. Normalize text
############################################################
# Convert text to lowercase and remove accents.
# Useful for converting labels such as "Cachaça" into "cachaca".

normalize_text <- function(x) {
  x <- tolower(as.character(x))
  x <- chartr(
    "áàãâéèêíìîóòõôúùûç",
    "aaaaeeeiiioooouuuc",
    x
  )
  return(x)
}


############################################################
# 3. Prepare genotype matrix X
############################################################

prepare_X <- function(dat) {
  
  metadata_cols <- c(
    "IID",
    "group_original",
    "class_binary",
    "group_multiclass"
  )
  
  snp_cols <- setdiff(names(dat), metadata_cols)
  
  X <- dat[, snp_cols, drop = FALSE]
  
  for (j in seq_along(X)) {
    X[[j]] <- as.numeric(X[[j]])
  }
  
  X <- as.matrix(X)
  
  has_variation <- apply(X, 2, function(x) {
    length(unique(na.omit(x))) > 1
  })
  
  if (any(!has_variation)) {
    message("Removing SNPs without variation: ", sum(!has_variation))
    X <- X[, has_variation, drop = FALSE]
  }
  
  if (any(is.na(X))) {
    message("Imputing missing values using the mean genotype of each SNP.")
    
    for (j in seq_len(ncol(X))) {
      missing <- is.na(X[, j])
      if (any(missing)) {
        X[missing, j] <- mean(X[, j], na.rm = TRUE)
      }
    }
  }
  
  return(X)
}


############################################################
# 4. Create multiclass group labels
############################################################

add_multiclass_group <- function(dat) {
  
  if (!"group_original" %in% names(dat)) {
    stop("Column 'group_original' was not found in the input table.")
  }
  
  group_clean <- normalize_text(dat$group_original)
  
  # Standardize separators.
  group_clean <- gsub("[^a-z0-9]+", "_", group_clean)
  
  dat$group_multiclass <- NA_character_
  
  dat$group_multiclass[grepl("cachaca", group_clean)] <- "Cachaca"
  dat$group_multiclass[grepl("beer|cerveja", group_clean)] <- "Beer"
  dat$group_multiclass[grepl("bioethanol|bioetanol|ethanol|etanol", group_clean)] <- "Bioethanol"
  dat$group_multiclass[grepl("wine|vinho", group_clean)] <- "Wine"
  
  message("\nOriginal values in group_original:")
  print(table(dat$group_original, useNA = "ifany"))
  
  message("\nRecognized values in group_multiclass:")
  print(table(dat$group_multiclass, useNA = "ifany"))
  
  return(dat)
}


############################################################
# 5. Run a simple sPLS-DA analysis
############################################################

run_simple_splsda <- function(
    dat,
    Y,
    set_name,
    analysis_name,
    analysis_outdir,
    ncomp_use = 2,
    keepX_value = 5,
    folds_max = 5,
    nrepeat = 10,
    dist_method = "max.dist"
) {
  
  message("\n============================================================")
  message("Running analysis: ", analysis_name)
  message("SNP set: ", set_name)
  message("============================================================")
  
  dir.create(analysis_outdir, showWarnings = FALSE, recursive = TRUE)
  
  Y <- factor(Y)
  
  message("\nClass distribution:")
  print(table(Y))
  
  X <- prepare_X(dat)
  
  message("\nNumber of samples: ", nrow(X))
  message("Number of SNPs: ", ncol(X))
  
  if (ncol(X) < 2) {
    stop("Fewer than 2 informative SNPs. sPLS-DA cannot be run.")
  }
  
  ncomp_use <- min(ncomp_use, ncol(X))
  keepX_use <- rep(min(keepX_value, ncol(X)), ncomp_use)
  
  message("\nModel parameters:")
  message("ncomp = ", ncomp_use)
  message("keepX = ", paste(keepX_use, collapse = ", "))
  
  ############################################################
  # Run model
  ############################################################
  
  model <- mixOmics::splsda(
    X = X,
    Y = Y,
    ncomp = ncomp_use,
    keepX = keepX_use
  )
  
  
  ############################################################
  # Extract selected SNPs
  ############################################################
  
  selected_snps <- data.frame()
  
  for (comp in seq_len(ncomp_use)) {
    
    selected <- mixOmics::selectVar(model, comp = comp)$name
    
    temp <- data.frame(
      set = set_name,
      analysis = analysis_name,
      component = comp,
      SNP = selected,
      stringsAsFactors = FALSE
    )
    
    selected_snps <- rbind(selected_snps, temp)
  }
  
  write.table(
    selected_snps,
    file.path(analysis_outdir, paste0("selected_snps_", analysis_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("\nSelected SNPs:")
  print(selected_snps)
  
  
  ############################################################
  # Predict classes in the training dataset
  ############################################################
  
  pred <- predict(
    model,
    X,
    dist = dist_method
  )
  
  predicted_class <- pred$class[[dist_method]][, ncomp_use]
  
  predictions <- data.frame(
    IID = dat$IID,
    group_original = dat$group_original,
    true_class = as.character(Y),
    predicted_class = as.character(predicted_class),
    correct = as.character(Y) == as.character(predicted_class),
    stringsAsFactors = FALSE
  )
  
  write.table(
    predictions,
    file.path(analysis_outdir, paste0("predictions_", analysis_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  confusion_matrix <- as.data.frame(
    table(
      true_class = predictions$true_class,
      predicted_class = predictions$predicted_class
    )
  )
  
  write.table(
    confusion_matrix,
    file.path(analysis_outdir, paste0("confusion_matrix_", analysis_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("\nConfusion matrix:")
  print(confusion_matrix)
  
  
  ############################################################
  # Cross-validation
  ############################################################
  
  folds_use <- min(folds_max, min(table(Y)))
  
  if (folds_use < 2) {
    stop("At least one class has fewer than 2 samples. Cross-validation cannot be run.")
  }
  
  message("\nRunning cross-validation...")
  message("folds = ", folds_use)
  
  perf_model <- mixOmics::perf(
    model,
    validation = "Mfold",
    folds = folds_use,
    nrepeat = nrepeat,
    dist = dist_method,
    progressBar = FALSE
  )
  
  BER <- perf_model$error.rate$BER[ncomp_use, dist_method]
  overall_error <- perf_model$error.rate$overall[ncomp_use, dist_method]
  
  performance <- data.frame(
    set = set_name,
    analysis = analysis_name,
    n_samples = nrow(X),
    n_classes = length(levels(Y)),
    n_snps = ncol(X),
    ncomp = ncomp_use,
    keepX = paste(keepX_use, collapse = ","),
    folds = folds_use,
    nrepeat = nrepeat,
    BER = BER,
    overall_error = overall_error,
    stringsAsFactors = FALSE
  )
  
  write.table(
    performance,
    file.path(analysis_outdir, paste0("performance_", analysis_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("\nPerformance:")
  print(performance)
  
  
  ############################################################
  # Score plot
  ############################################################
  
  scores <- as.data.frame(model$variates$X)
  scores$IID <- dat$IID
  scores$group_original <- dat$group_original
  scores$class <- Y
  
  if (ncomp_use >= 2) {
    
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(x = comp1, y = comp2, color = class)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(
        title = paste0("sPLS-DA — ", set_name, " — ", analysis_name),
        x = "Component 1",
        y = "Component 2",
        color = "Class"
      )
    
  } else {
    
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(x = comp1, y = class, color = class)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(
        title = paste0("sPLS-DA — ", set_name, " — ", analysis_name),
        x = "Component 1",
        y = "Class",
        color = "Class"
      )
  }
  
  ggplot2::ggsave(
    filename = file.path(analysis_outdir, paste0("scores_", analysis_name, ".png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  ############################################################
  # Return results
  ############################################################
  
  return(
    list(
      model = model,
      performance = performance,
      confusion_matrix = confusion_matrix,
      selected_snps = selected_snps,
      predictions = predictions
    )
  )
}
```

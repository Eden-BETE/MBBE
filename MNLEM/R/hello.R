# ============================================================================
# UTILITY FUNCTIONS AND VALIDATION
# ============================================================================

#' Validate required packages
#' @keywords internal
validate_packages <- function() {
  required_packages <- c("mvtnorm", "saemix", "car")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.")
    }
  }

  library(mvtnorm)
  library(saemix)
  library(car)
}

#' Validate and preprocess input data
#' @param data_path Path to the dataset file
#' @return Preprocessed data frame
#' @keywords internal
validate_and_load_data <- function(data_path) {
  if (!file.exists(data_path)) {
    stop("Data file not found at path: ", data_path)
  }

  # Fonction pour tester différents séparateurs
  try_read <- function(sep) {
    tryCatch({
      read.table(data_path,
                 header = TRUE,
                 sep = sep,
                 dec = ".",
                 na.strings = c("", ".", "NA"),
                 stringsAsFactors = FALSE,
                 row.names = NULL,
                 fill = TRUE)
    }, error = function(e) NULL)
  }

  # Tester plusieurs séparateurs
  possible_separators <- c("\t", " ", ",", ";")
  for (sep in possible_separators) {
    tab <- try_read(sep)
    if (!is.null(tab) && ncol(tab) >= 5) break
  }

  if (is.null(tab) || ncol(tab) < 5) {
    stop("Unable to read the dataset. Please check the file format or structure.")
  }

  # Standardiser les noms de colonnes
  names(tab) <- tolower(names(tab))

  # Colonnes attendues
  required_cols <- c("id", "dose", "time", "treatment", "concentration")
  missing_cols <- setdiff(required_cols, names(tab))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Réorganiser les colonnes dans l’ordre attendu
  tab <- tab[, required_cols]

  # Renommer pour uniformiser les noms utilisés plus tard
  names(tab) <- c("Id", "Dose", "Time", "Tr", "Concentration")

  return(tab)
}

#' Validate data content and structure
#' @param tab Data frame with standardized column names
#' @return List with validated data and metadata
#' @keywords internal
validate_data_content <- function(tab) {
  # Check treatment values
  unique_treatments <- unique(tab$Tr)
  if (!all(c("R", "T") %in% unique_treatments)) {
    stop("Treatment column 'treatment' must contain exactly 'R' (reference) and 'T' (test) values")
  }
  if (length(unique_treatments) != 2) {
    warning("Treatment column contains values other than 'R' and 'T': ",
            paste(setdiff(unique_treatments, c("R", "T")), collapse = ", "))
  }

  # Handle dose validation and imputation
  unique_doses <- unique(tab$Dose[!is.na(tab$Dose)])
  if (length(unique_doses) > 1) {
    stop("Multiple dose values found: ", paste(unique_doses, collapse = ", "),
         ". All subjects must have the same dose for bioequivalence analysis.")
  }

  if (any(is.na(tab$Dose))) {
    actual_dose <- unique_doses[1]
    cat("Some dose values are missing. Imputing with dose found in data: ", actual_dose, "\n")
    tab$Dose[is.na(tab$Dose)] <- actual_dose
  }

  # Check sampling times consistency
  subjects <- unique(tab$Id)
  n_subjects <- length(subjects)

  sampling_counts <- sapply(subjects, function(id) sum(tab$Id == id))
  if (length(unique(sampling_counts)) > 1) {
    stop("All subjects must have the same number of sampling times for Gallant correction. ",
         "Current counts: ", paste(unique(sampling_counts), collapse = ", "))
  }

  nb_t <- sampling_counts[1]

  # Final data preparation
  tab$Treat <- ifelse(tab$Tr == "R", 0, 1)
  tab <- tab[, c("Id", "Dose", "Time", "Treat", "Concentration")]

  cat("Data validation completed successfully:\n")
  cat("- Number of subjects:", n_subjects, "\n")
  cat("- Number of sampling times per subject:", nb_t, "\n")
  cat("- Total observations:", nrow(tab), "\n\n")

  return(list(
    data = tab,
    n_subjects = n_subjects,
    nb_timepoints = nb_t
  ))
}

# ============================================================================
# PHARMACOKINETIC MODELS
# ============================================================================

#' Stable difference function to prevent numerical issues
#' @keywords internal
stable_diff <- function(a, b, min_diff = 1e-6) {
  diff <- a - b
  ifelse(abs(diff) < min_diff, sign(diff) * min_diff, diff)
}

#' One-compartment model with first-order absorption
#' @keywords internal
model_1cpt <- function(psi, id, xidep) {
  dose <- xidep[, 1]
  time <- xidep[, 2]
  ka <- psi[id, 1]
  V <- psi[id, 2]
  CL <- psi[id, 3]
  k <- CL / V
  denom <- stable_diff(ka, k)
  ypred <- dose * ka / (V * denom) * (exp(-k * time) - exp(-ka * time))
  return(ypred)
}

#' One-compartment model with first-order absorption and lag time
#' @keywords internal
model_1cpt_tlag <- function(psi, id, xidep) {
  dose <- xidep[,1]
  tim <- xidep[,2]
  ka <- psi[id,1]
  V <- psi[id,2]
  CL <- psi[id,3]
  tlag <- psi[id,4]
  k <- CL / V
  denom <- stable_diff(ka, k)
  tim_adj <- pmax(0, tim - tlag)
  ypred <- ifelse(tim <= tlag, 0,
                  dose * ka / (V * denom) * (exp(-k * tim_adj) - exp(-ka * tim_adj)))
  return(ypred)
}

#' Two-compartment model with first-order absorption
#' @keywords internal
model_2cpt <- function(psi, id, xidep) {
  dose <- xidep[,1]
  tim <- xidep[,2]
  ka <- psi[id,1]
  V1 <- psi[id,2]
  CL <- psi[id,3]
  Q <- psi[id,4]
  V2 <- psi[id,5]

  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2

  alpha <- 0.5 * (k12 + k21 + k10 + sqrt((k12 + k21 + k10)^2 - 4 * k21 * k10))
  beta <- 0.5 * (k12 + k21 + k10 - sqrt((k12 + k21 + k10)^2 - 4 * k21 * k10))

  denom_ka_alpha <- stable_diff(ka, alpha)
  denom_ka_beta <- stable_diff(ka, beta)
  denom_beta_alpha <- stable_diff(beta, alpha)
  denom_alpha_beta <- stable_diff(alpha, beta)

  A <- dose * ka * (k21 - alpha) / (V1 * denom_ka_alpha * denom_beta_alpha)
  B <- dose * ka * (k21 - beta) / (V1 * denom_ka_beta * denom_alpha_beta)
  C <- dose * ka / (V1 * denom_ka_alpha * denom_ka_beta)

  ypred <- A * exp(-alpha * tim) + B * exp(-beta * tim) - C * exp(-ka * tim)
  return(ypred)
}

#' Two-compartment model with first-order absorption and lag time
#' @keywords internal
model_2cpt_tlag <- function(psi, id, xidep) {
  dose <- xidep[,1]
  tim <- xidep[,2]
  ka <- psi[id,1]
  V1 <- psi[id,2]
  CL <- psi[id,3]
  Q <- psi[id,4]
  V2 <- psi[id,5]
  tlag <- psi[id,6]

  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2

  alpha <- 0.5 * (k12 + k21 + k10 + sqrt((k12 + k21 + k10)^2 - 4 * k21 * k10))
  beta <- 0.5 * (k12 + k21 + k10 - sqrt((k12 + k21 + k10)^2 - 4 * k21 * k10))

  tim_adj <- pmax(0, tim - tlag)

  denom_ka_alpha <- stable_diff(ka, alpha)
  denom_ka_beta <- stable_diff(ka, beta)
  denom_beta_alpha <- stable_diff(beta, alpha)
  denom_alpha_beta <- stable_diff(alpha, beta)

  A <- dose * ka * (k21 - alpha) / (V1 * denom_ka_alpha * denom_beta_alpha)
  B <- dose * ka * (k21 - beta) / (V1 * denom_ka_beta * denom_alpha_beta)
  C <- dose * ka / (V1 * denom_ka_alpha * denom_ka_beta)

  ypred <- ifelse(tim <= tlag, 0,
                  A * exp(-alpha * tim_adj) + B * exp(-beta * tim_adj) - C * exp(-ka * tim_adj))
  return(ypred)
}

#' Get model configuration
#' @param model_type Character string specifying the PK model
#' @return List with model configuration
#' @keywords internal
get_model_config <- function(model_type) {
  valid_models <- c("1cpt", "1cpt_tlag", "2cpt", "2cpt_tlag")
  if (!model_type %in% valid_models) {
    stop("Invalid model_type. Choose from: ", paste(valid_models, collapse = ", "))
  }

  model_configs <- list(
    "1cpt" = list(
      model = model_1cpt,
      param_names = c("ka", "V", "CL"),
      nb_param = 3,
      psi0 = matrix(c(1, 20, 1), ncol = 3, byrow = TRUE),
      omega_init = diag(c(0.1, 0.1, 0.1)),
      cmax_formula = "(exp(beta_CL) / (exp(beta_ka + beta_V) - exp(beta_CL))) *
                      log(exp(beta_ka + beta_V) / exp(beta_CL))"
    ),

    "1cpt_tlag" = list(
      model = model_1cpt_tlag,
      param_names = c("ka", "V", "CL", "tlag"),
      nb_param = 4,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0.5, 0, 0, 0, 0), ncol = 4, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05, 0.1)),
      cmax_formula = "(exp(beta_CL) / (exp(beta_ka + beta_V) - exp(beta_CL))) *
                      exp(-exp(beta_CL - beta_V) * exp(beta_tlag)) *
                      ((exp(beta_ka + beta_V) / exp(beta_CL))^
                      (exp(beta_CL - beta_V) / (exp(beta_ka) - exp(beta_CL))))"
    ),

    "2cpt" = list(
      model = model_2cpt,
      param_names = c("ka", "V1", "CL", "Q", "V2"),
      nb_param = 5,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0.1, 0.2, 0, 0, 0, 0, 0), ncol = 5, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05, 0.1, 0.1)),
      cmax_formula = "(exp(beta_ka) * exp(beta_CL)) /
                      (exp(beta_V1) * (exp(beta_Q) + exp(beta_CL)))"
    ),

    "2cpt_tlag" = list(
      model = model_2cpt_tlag,
      param_names = c("ka", "V1", "CL", "Q", "V2", "tlag"),
      nb_param = 6,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0.1, 0.2, 0.5, 0, 0, 0, 0, 0, 0), ncol = 6, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05, 0.1, 0.1, 0.1)),
      cmax_formula = "((exp(beta_ka) * exp(beta_CL)) /
                       (exp(beta_V1) * (exp(beta_Q) + exp(beta_CL)))) *
                       exp(-exp(beta_tlag))"
    )
  )

  return(model_configs[[model_type]])
}

# ============================================================================
# SAEMIX MODEL FITTING
# ============================================================================

#' Fit SAEMIX model
#' @param data Preprocessed data frame
#' @param config Model configuration
#' @param nb_chains Number of MCMC chains
#' @param nb_iter Vector of burn-in and sampling iterations
#' @return Fitted SAEMIX object
#' @keywords internal
fit_saemix_model <- function(data, config, nb_chains, nb_iter) {
  # Remove NA concentrations
  data <- data[!is.na(data$Concentration), ]

  saemix.data <- saemixData(
    name.data = data,
    header = TRUE,
    sep = " ",
    na = NA,
    name.group = c("Id"),
    name.predictors = c("Dose", "Time"),
    name.response = c("Concentration"),
    name.X = "Time",
    units = list(x = "hr", y = "mg/L")
  )

  saemix.modelb <- saemixModel(
    model = config$model,
    description = paste(config$nb_param, "-compartment model"),
    psi0 = config$psi0,
    transform.par = rep(1, config$nb_param),
    covariance.model = diag(config$nb_param),
    omega.init = config$omega_init,
    error.model = "combined",
    error.init = c(0.1, 0.1),
    name.modpar = config$param_names
  )

  cat("Fitting SAEM model...\n")
  mod1 <- saemix(saemix.modelb, saemix.data, list(
    nb.chains = nb_chains,
    nbiter.saemix = nb_iter,
    seed = 1234,
    map = TRUE,
    fim = TRUE,
    ll.is = TRUE,
    save = FALSE,
    print = FALSE
  ))
  cat("Model fitting completed.\n\n")

  return(mod1)
}

# ============================================================================
# TOST ANALYSIS CALCULATIONS
# ============================================================================

#' Calculate delta method estimates for Cmax
#' @param mod1 Fitted SAEMIX object
#' @param config Model configuration
#' @return Data frame with delta method results
#' @keywords internal
calculate_delta_method <- function(mod1, config) {
  estimates <- mod1@results@fixed.effects
  fim <- mod1@results@fim

  if (is.null(estimates) || is.null(fim)) {
    warning("Missing model results: estimates or FIM.")
    return(data.frame(Estimate = NA, SE = NA))
  }

  n_params <- config$nb_param

  # Prendre les effets fixes de population uniquement
  mu_estimates <- estimates[1:n_params]
  names(mu_estimates) <- paste0("beta_", config$param_names)

  fim_mu_block <- fim[1:n_params, 1:n_params]

  if (any(is.na(fim_mu_block)) || any(!is.finite(fim_mu_block))) {
    warning("FIM block for mu contains NA or Inf.")
    return(data.frame(Estimate = NA, SE = NA))
  }

  varcov <- tryCatch(solve(fim_mu_block), error = function(e) {
    warning("FIM is not invertible for mu block: ", e$message)
    return(matrix(NA, nrow = n_params, ncol = n_params))
  })

  if (any(is.na(varcov))) {
    return(data.frame(Estimate = NA, SE = NA))
  }

  formula <- config$cmax_formula

  # Debug output
  cat("=== DEBUG DELTA METHOD ===\n")
  print(mu_estimates)
  print(varcov)
  cat("Formula:\n", formula, "\n")

  result <- tryCatch({
    deltaMethod(mu_estimates, formula, vcov. = varcov)
  }, error = function(e) {
    warning("Delta method failed: ", e$message)
    return(data.frame(Estimate = NA, SE = NA))
  })

  return(result)
}


#' Calculate quantiles for different TOST approaches
#' @param n Number of subjects
#' @param nb_t Number of timepoints
#' @param nb_param Number of parameters
#' @return List with quantiles
#' @keywords internal
calculate_quantiles <- function(n, nb_t, nb_param) {
  list(
    normal = qnorm(0.95),
    student = qt(0.95, df = (n * nb_t - 2 * nb_param - 1)),
    gallant = qt(0.95, df = n - nb_param)
  )
}

#' Perform TOST tests for AUC
#' @param mod1 Fitted SAEMIX object
#' @param config Model configuration
#' @param quantiles List of quantiles
#' @param n Number of subjects
#' @param nb_param Number of parameters
#' @return List with AUC TOST results
#' @keywords internal
perform_tost_auc <- function(mod1, config, quantiles, n, nb_param) {
  delta <- log(1.25)

  # Indice du paramètre CL
  cl_idx <- which(config$param_names == "CL")
  if (length(cl_idx) == 0) {
    warning("CL parameter not found in model")
    return(list(
      normal = list(be = NA, W1 = NA, W2 = NA),
      student = list(be = NA, W1 = NA, W2 = NA),
      gallant = list(be = NA, W1 = NA, W2 = NA)
    ))
  }

  # EXTRACTION CORRIGÉE depuis les mu
  beta_cl <- mod1@results@fixed.effects[cl_idx]
  se_cl <- sqrt(solve(mod1@results@fim)[cl_idx, cl_idx])

  # DEBUG
  cat("AUC TOST calculations:\n")
  cat("mu_CL:", beta_cl, "\n")
  cat("SE_CL:", se_cl, "\n")

  W1 <- (-beta_cl + delta) / se_cl
  W2 <- (-beta_cl - delta) / se_cl
  normal_be <- (W1 >= quantiles$normal) & (W2 <= -quantiles$normal)

  student_be <- (W1 >= quantiles$student) & (W2 <= -quantiles$student)

  corrected_se <- se_cl * sqrt(n / (n - nb_param))
  W1_gallant <- (-beta_cl + delta) / corrected_se
  W2_gallant <- (-beta_cl - delta) / corrected_se
  gallant_be <- (W1_gallant >= quantiles$gallant) & (W2_gallant <= -quantiles$gallant)

  return(list(
    normal = list(be = normal_be, W1 = W1, W2 = W2),
    student = list(be = student_be, W1 = W1, W2 = W2),
    gallant = list(be = gallant_be, W1 = W1_gallant, W2 = W2_gallant)
  ))
}


#' Perform TOST tests for Cmax
#' @param SE_deltam Delta method results
#' @param quantiles List of quantiles
#' @param n Number of subjects
#' @param nb_param Number of parameters
#' @return List with Cmax TOST results
#' @keywords internal
perform_tost_cmax <- function(SE_deltam, quantiles, n, nb_param) {
  if (is.na(SE_deltam[1,1]) || is.na(SE_deltam[1,2])) {
    return(list(
      normal = list(be = NA, W1 = NA, W2 = NA),
      student = list(be = NA, W1 = NA, W2 = NA),
      gallant = list(be = NA, W1 = NA, W2 = NA)
    ))
  }

  delta <- log(1.25)
  beta_cmax_estim <- SE_deltam[1,1]
  secmax <- SE_deltam[1,2]

  # Debug information
  cat("Cmax TOST calculations:\n")
  cat("beta_Cmax:", beta_cmax_estim, "\n")
  cat("SE_Cmax:", secmax, "\n")

  # Normal approach
  W1cmax <- (beta_cmax_estim + delta) / secmax
  W2cmax <- (beta_cmax_estim - delta) / secmax
  normal_be <- (W1cmax >= quantiles$normal) & (W2cmax <= -quantiles$normal)

  # Student approach
  student_be <- (W1cmax >= quantiles$student) & (W2cmax <= -quantiles$student)

  # Gallant approach
  corrected_se <- secmax * sqrt(n / (n - nb_param))
  W1_gallant <- (beta_cmax_estim + delta) / corrected_se
  W2_gallant <- (beta_cmax_estim - delta) / corrected_se
  gallant_be <- (W1_gallant >= quantiles$gallant) & (W2_gallant <= -quantiles$gallant)

  return(list(
    normal = list(be = normal_be, W1 = W1cmax, W2 = W2cmax),
    student = list(be = student_be, W1 = W1cmax, W2 = W2cmax),
    gallant = list(be = gallant_be, W1 = W1_gallant, W2 = W2_gallant)
  ))
}

#' Print analysis summary
#' @param results TOST results
#' @param model_type Model type
#' @param data_path Data path
#' @param n Number of subjects
#' @param nb_t Number of timepoints
#' @keywords internal
print_analysis_summary <- function(results, model_type, data_path, n, nb_t) {
  cat("=== BIOEQUIVALENCE ANALYSIS RESULTS ===\n")
  cat("Model:", model_type, "\n")
  cat("Dataset:", basename(data_path), "\n")
  cat("Subjects:", n, "| Sampling times per subject:", nb_t, "\n\n")

  cat("AUC Bioequivalence Results:\n")
  cat("- Normal approach:  ", ifelse(is.na(results$auc$normal$be), "NA", ifelse(results$auc$normal$be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n")
  cat("- Student approach: ", ifelse(is.na(results$auc$student$be), "NA", ifelse(results$auc$student$be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n")
  cat("- Gallant approach: ", ifelse(is.na(results$auc$gallant$be), "NA", ifelse(results$auc$gallant$be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n\n")

  cat("Cmax Bioequivalence Results:\n")
  cat("- Normal approach:  ", ifelse(is.na(results$cmax$normal$be), "NA", ifelse(results$cmax$normal$be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n")
  cat("- Student approach: ", ifelse(is.na(results$cmax$student$be), "NA", ifelse(results$cmax$student$be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n")
  cat("- Gallant approach: ", ifelse(is.na(results$cmax$gallant$be), "NA", ifelse(results$cmax$gallant$be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n\n")

  overall_be <- if (is.na(results$auc$gallant$be) || is.na(results$cmax$gallant$be)) {
    NA
  } else {
    all(c(results$auc$gallant$be, results$cmax$gallant$be))
  }

  cat("Overall Bioequivalence (Gallant method):", ifelse(is.na(overall_be), "NA", ifelse(overall_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")), "\n")
  cat("========================================\n")
}

# ============================================================================
# MAIN REFACTORED FUNCTION
# ============================================================================

#' Perform TOST Analysis on Pharmacokinetic Data
#'
#' This function performs TOST (Two One-Sided Tests) for bioequivalence analysis on pharmacokinetic parameters (AUC and Cmax),
#' using asymptotic standard errors (Normal and Student) and Gallant correction.
#'
#' @param data_path Path to the dataset file (must be a `.txt` file with tabular structure).
#'   The file should contain the following **columns** with headers **in this exact order**:
#'   \itemize{
#'     \item \code{id}: Subject identifier (numeric or character).
#'     \item \code{dose}: Dose administered (numeric).
#'     \item \code{time}: Time point of concentration measurement (numeric).
#'     \item \code{treatment}: Treatment identifier ("R" for reference treatment and "T" for the test treatment).
#'     \item \code{concentration}: Measured drug concentration (numeric).
#'   }
#' @param model_type Character string specifying the PK model to use. Options:
#'   \itemize{
#'     \item \code{"1cpt"}: One-compartment model with first-order absorption (default)
#'     \item \code{"1cpt_tlag"}: One-compartment model with first-order absorption and lag time
#'     \item \code{"2cpt"}: Two-compartment model with first-order absorption
#'     \item \code{"2cpt_tlag"}: Two-compartment model with first-order absorption and lag time
#'   }
#' @param nb_chains Number of MCMC chains used by SAEM algorithm (default: 10).
#' @param nb_iter A numeric vector of two integers: number of burn-in iterations and sampling iterations for SAEM (default: c(300, 100)).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{saemix_object}: Fitted SAEM model object.
#'     \item \code{res_tost}: Results from the TOST procedure (normal, student, Gallant).
#'     \item \code{res_delta}: Delta method estimates.
#'     \item \code{summary}: Summary of bioequivalence results.
#'   }
#'
#' @examples
#' \dontrun{
#' perform_tost_analysis("path/to/data.txt", model_type = "1cpt")
#' }
#'
#' @export
perform_tost_analysis <- function(data_path, model_type = "1cpt", nb_chains = 10, nb_iter = c(300, 100)) {

  # Step 1: Validate packages and load data
  validate_packages()

  # Step 2: Load and validate data
  raw_data <- validate_and_load_data(data_path)
  validated_data <- validate_data_content(raw_data)

  # Step 3: Get model configuration
  config <- get_model_config(model_type)

  cat("Model selected:", model_type, "with", config$nb_param, "parameters\n\n")

  # Step 4: Fit SAEMIX model
  mod1 <- fit_saemix_model(validated_data$data, config, nb_chains, nb_iter)

  # Step 5: Calculate delta method for Cmax
  SE_deltam <- calculate_delta_method(mod1, config)

  # Step 6: Calculate quantiles
  quantiles <- calculate_quantiles(validated_data$n_subjects, validated_data$nb_timepoints, config$nb_param)

  # Step 7: Perform TOST tests
  auc_results <- perform_tost_auc(mod1, config, quantiles, validated_data$n_subjects, config$nb_param)
  cmax_results <- perform_tost_cmax(SE_deltam, quantiles, validated_data$n_subjects, config$nb_param)

  # Step 8: Compile results

  # Step 8: Compile results
  results <- list(auc = auc_results, cmax = cmax_results)

  # Step 9: Print summary
  print_analysis_summary(results, model_type, data_path, validated_data$n_subjects, validated_data$nb_timepoints)

  # Step 10: Return comprehensive results
  overall_be <- all(c(auc_results$gallant$be, cmax_results$gallant$be))

  return(list(
    saemix_object = mod1,
    res_tost = results,
    res_delta = SE_deltam,
    summary = list(
      model_type = model_type,
      n_subjects = validated_data$n_subjects,
      n_timepoints = validated_data$nb_timepoints,
      n_parameters = config$nb_param,
      bioequivalence = list(
        auc = list(normal = auc_results$normal$be, student = auc_results$student$be, gallant = auc_results$gallant$be),
        cmax = list(normal = cmax_results$normal$be, student = cmax_results$student$be, gallant = cmax_results$gallant$be),
        overall_gallant = overall_be
      )
    )
  ))
}

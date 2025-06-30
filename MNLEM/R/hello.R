#' Perform TOST Analysis on Pharmacokinetic Data
#'
#' This function performs TOST (Two One-Sided Tests) for bioequivalence analysis on pharmacokinetic parameters (AUC and Cmax),
#' using asymptotic standard errors (Normal and Student) and Gallant correction.
#'
#' @param data_path Path to the dataset file (must be a `.txt` file with tabular structure).
#'   The file should contain the following **columns** with headers **in this exact order**:
#'   \itemize{
#'     \item \code{id}: Subject identifier (numeric or character).
#'     \item \code{dose}: Dose administered (numeric). If missing, it will be imputed with \code{dose_value}.
#'     \item \code{time}: Time point of concentration measurement (numeric).
#'     \item \code{treatment}: Treatment identifier ("R" for reference treatment and "T" for the test treatment).
#'     \item \code{concentration}: Measured drug concentration (numeric).
#'   }
#'
#' @param model_type Character string specifying the PK model to use. Options:
#'   \itemize{
#'     \item \code{"1cpt"}: One-compartment model with first-order absorption (default)
#'     \item \code{"1cpt_tlag"}: One-compartment model with first-order absorption and lag time
#'     \item \code{"2cpt"}: Two-compartment model with first-order absorption
#'     \item \code{"2cpt_tlag"}: Two-compartment model with first-order absorption and lag time
#'   }
#' @param dose_value Numeric value of the dose administered, used for imputation only if all dose values are missing (default: 4).
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

perform_tost_analysis <- function(data_path, model_type = "1cpt", dose_value = 4, nb_chains = 10, nb_iter = c(300, 100)) {

  # Load required packages
  if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("Package 'mvtnorm' is required.")
  if (!requireNamespace("saemix", quietly = TRUE)) stop("Package 'saemix' is required.")
  if (!requireNamespace("car", quietly = TRUE)) stop("Package 'car' is required.")

  library(mvtnorm)
  library(saemix)
  library(car)

  # Constants
  delta <- log(1.25)
  quant <- qnorm(0.95)

  # Validate model_type
  valid_models <- c("1cpt", "1cpt_tlag", "2cpt", "2cpt_tlag")
  if (!model_type %in% valid_models) {
    stop("Invalid model_type. Choose from: ", paste(valid_models, collapse = ", "))
  }

  # Define PK models
  model_1cpt <- function(psi, id, xidep) {
    dose <- xidep[,1]
    tim <- xidep[,2]
    ka <- psi[id,1]
    V <- psi[id,2]
    CL <- psi[id,3]
    k <- CL / V
    ypred <- dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    return(ypred)
  }

  model_1cpt_tlag <- function(psi, id, xidep) {
    dose <- xidep[,1]
    tim <- xidep[,2]
    ka <- psi[id,1]
    V <- psi[id,2]
    CL <- psi[id,3]
    tlag <- psi[id,4]
    k <- CL / V
    tim_adj <- pmax(0, tim - tlag)
    ypred <- ifelse(tim <= tlag, 0,
                    dose * ka / (V * (ka - k)) * (exp(-k * tim_adj) - exp(-ka * tim_adj)))
    return(ypred)
  }

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

    A <- dose * ka * (k21 - alpha) / (V1 * (ka - alpha) * (beta - alpha))
    B <- dose * ka * (k21 - beta) / (V1 * (ka - beta) * (alpha - beta))
    C <- dose * ka / (V1 * (ka - alpha) * (ka - beta))

    ypred <- A * exp(-alpha * tim) + B * exp(-beta * tim) - C * exp(-ka * tim)
    return(ypred)
  }

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

    A <- dose * ka * (k21 - alpha) / (V1 * (ka - alpha) * (beta - alpha))
    B <- dose * ka * (k21 - beta) / (V1 * (ka - beta) * (alpha - beta))
    C <- dose * ka / (V1 * (ka - alpha) * (ka - beta))

    ypred <- ifelse(tim <= tlag, 0,
                    A * exp(-alpha * tim_adj) + B * exp(-beta * tim_adj) - C * exp(-ka * tim_adj))
    return(ypred)
  }

  # Select model and parameters
  model_configs <- list(
    "1cpt" = list(
      model = model_1cpt,
      param_names = c("ka", "V", "CL"),
      nb_param = 3,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0, 0, 0), ncol = 3, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05)),
      cmax_formula = "-beta_V-(log((Ka*V)/Cl)+beta_Ka+beta_V-beta_Cl)*Cl*exp(beta_Cl)/(Ka*V*exp(beta_Ka+beta_V)-Cl*exp(beta_Cl))+(Cl/(Ka*V-Cl))*log(Ka*V/Cl)"
    ),
    "1cpt_tlag" = list(
      model = model_1cpt_tlag,
      param_names = c("ka", "V", "CL", "tlag"),
      nb_param = 4,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0.5, 0, 0, 0, 0), ncol = 4, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05, 0.1)),
      cmax_formula = "-beta_V-(log((Ka*V)/Cl)+beta_Ka+beta_V-beta_Cl)*Cl*exp(beta_Cl)/(Ka*V*exp(beta_Ka+beta_V)-Cl*exp(beta_Cl))+(Cl/(Ka*V-Cl))*log(Ka*V/Cl)"
    ),
    "2cpt" = list(
      model = model_2cpt,
      param_names = c("ka", "V1", "CL", "Q", "V2"),
      nb_param = 5,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0.1, 0.2, 0, 0, 0, 0, 0), ncol = 5, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05, 0.1, 0.1)),
      cmax_formula = "-beta_V1"  # Simplified for 2-compartment
    ),
    "2cpt_tlag" = list(
      model = model_2cpt_tlag,
      param_names = c("ka", "V1", "CL", "Q", "V2", "tlag"),
      nb_param = 6,
      psi0 = matrix(c(1.5, 0.5, 0.04, 0.1, 0.2, 0.5, 0, 0, 0, 0, 0, 0), ncol = 6, nrow = 2, byrow = TRUE),
      omega_init = diag(c(0.05, 0.0125, 0.05, 0.1, 0.1, 0.1)),
      cmax_formula = "-beta_V1"  # Simplified for 2-compartment
    )
  )

  config <- model_configs[[model_type]]
  nb_param <- config$nb_param

  # Load and validate data
  if (!file.exists(data_path)) {
    stop("Data file not found at path: ", data_path)
  }

  tab <- read.table(data_path, header = TRUE, sep = " ", dec = ".", na.strings = ".")

  # Data validation and warnings
  required_cols <- c("id", "dose", "time", "treatment", "concentration")
  expected_order <- c("id", "dose", "time", "treatment", "concentration")

  # Check if required columns exist
  missing_cols <- setdiff(required_cols, names(tab))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check column order
  if (!identical(names(tab)[1:5], expected_order)) {
    warning("Column order should be: ", paste(expected_order, collapse = ", "),
            ". Current order: ", paste(names(tab)[1:5], collapse = ", "),
            ". Reordering columns automatically.")
    tab <- tab[, expected_order]
  }

  # Standardize column names for internal processing
  names(tab)[1:5] <- c("Id", "Dose", "Time", "Tr", "Concentration")

  # Check treatment values
  unique_treatments <- unique(tab$Tr)
  if (!all(c("R", "T") %in% unique_treatments)) {
    stop("Treatment column 'treatment' must contain exactly 'R' (reference) and 'T' (test) values")
  }
  if (length(unique_treatments) != 2) {
    warning("Treatment column contains values other than 'R' and 'T': ",
            paste(setdiff(unique_treatments, c("R", "T")), collapse = ", "))
  }

  # Handle dose column
  if (all(is.na(tab$Dose))) {
    # If all doses are missing, use dose_value for imputation
    warning("All dose values are missing. Imputing with dose_value = ", dose_value)
    tab$Dose <- dose_value
  } else {
    # Check if all non-missing doses are identical
    unique_doses <- unique(tab$Dose[!is.na(tab$Dose)])
    if (length(unique_doses) > 1) {
      stop("Multiple dose values found: ", paste(unique_doses, collapse = ", "),
           ". All subjects must have the same dose for bioequivalence analysis.")
    }

    # If some doses are missing, impute with the unique dose found in data
    if (any(is.na(tab$Dose))) {
      actual_dose <- unique_doses[1]
      cat("Some dose values are missing. Imputing with dose found in data: ", actual_dose, "\n")
      tab$Dose[is.na(tab$Dose)] <- actual_dose

      # Warn if dose_value parameter differs from actual dose
      if (dose_value != actual_dose) {
        warning("dose_value parameter (", dose_value, ") differs from actual dose in data (",
                actual_dose, "). Using actual dose from data.")
      }
    }
  }

  # Check sampling times consistency (Gallant correction requirement)
  subjects <- unique(tab$Id)
  n_subjects <- length(subjects)

  sampling_counts <- sapply(subjects, function(id) {
    sum(tab$Id == id)
  })

  if (length(unique(sampling_counts)) > 1) {
    stop("All subjects must have the same number of sampling times for Gallant correction. ",
         "Current counts: ", paste(unique(sampling_counts), collapse = ", "))
  }

  nb_t <- sampling_counts[1]
  n <- n_subjects

  cat("Data validation completed successfully:\n")
  cat("- Number of subjects:", n, "\n")
  cat("- Number of sampling times per subject:", nb_t, "\n")
  cat("- Total observations:", nrow(tab), "\n")
  cat("- Model selected:", model_type, "with", nb_param, "parameters\n\n")

  # Prepare data
  tab <- tab[, c("Id", "Dose", "Time", "Concentration", "Tr")]
  tab$Treat <- ifelse(tab$Tr == "R", 0, 1)
  tab$Tr <- NULL
  tab <- tab[, c("Id", "Dose", "Time", "Treat", "Concentration")]

  # Calculate quantiles
  quant_stud <- qt(1 - 0.05, df = (n * nb_t - 2 * nb_param - 1))
  quant_gallant <- qt(1 - 0.05, df = n - nb_param)

  # SAEM input
  saemix.data <- saemixData(
    name.data = tab,
    header = TRUE,
    sep = " ",
    na = NA,
    name.group = c("Id"),
    name.predictors = c("Dose", "Time"),
    name.covariates = c("Treat"),
    name.response = c("Concentration"),
    name.X = "Time",
    units = list(x = "hr", y = "mg/L")
  )

  # Create covariate matrix
  cov_matrix <- matrix(1, ncol = nb_param)

  saemix.modelb <- saemixModel(
    model = config$model,
    description = paste(model_type, "model"),
    psi0 = config$psi0,
    transform.par = rep(1, nb_param),
    covariate.model = cov_matrix,
    covariance.model = diag(nb_param),
    omega.init = config$omega_init,
    error.model = "combined",
    error.init = c(0.1, 0.1)
  )

  # Fit model
  cat("Fitting SAEM model...\n")
  mod1 <- saemix(saemix.modelb, saemix.data, list(nb.chains = nb_chains, nbiter.saemix = nb_iter))
  cat("Model fitting completed.\n\n")

  # Extract estimates
  n_fixed <- 2 * nb_param  # Parameters + covariate effects
  estimates <- mod1@results@fixed.effects[1:n_fixed]
  fim_fixed <- mod1@results@fim[1:n_fixed, 1:n_fixed]
  varcov <- solve(fim_fixed)

  # Delta method for beta_Cmax (simplified for complex models)
  if (model_type %in% c("1cpt", "1cpt_tlag")) {
    param_names <- paste0(c("Ka", "V", "CL"), rep(c("", "_beta"), each = 3))[1:n_fixed]
    names(estimates) <- param_names

    SE_deltam <- tryCatch({
      deltaMethod(estimates, config$cmax_formula, vcov. = varcov)
    }, error = function(e) {
      warning("Delta method for Cmax failed. Using simplified approach.")
      data.frame(Estimate = -estimates[4], SE = sqrt(varcov[4,4]))  # beta_V
    })
  } else {
    # For 2-compartment models, use simplified approach
    SE_deltam <- data.frame(
      Estimate = -estimates[4],  # beta_V1
      SE = sqrt(varcov[4,4])
    )
  }

  beta_cmax_estim <- SE_deltam[1,1]
  secmax <- SE_deltam[1,2]

  # TOST Results Storage
  results <- list()

  # TOST AUC - Normal
  cl_param_idx <- ifelse(model_type %in% c("2cpt", "2cpt_tlag"), 6, 6)  # beta_CL position
  W1 <- (-mod1@results@fixed.effects[cl_param_idx] + delta) / mod1@results@se.fixed[cl_param_idx]
  W2 <- (-mod1@results@fixed.effects[cl_param_idx] - delta) / mod1@results@se.fixed[cl_param_idx]
  auc_normal_be <- (W1 >= quant) & (W2 <= -quant)

  # TOST AUC - Student
  auc_student_be <- (W1 >= quant_stud) & (W2 <= -quant_stud)

  # TOST AUC - Gallant
  corrected_se_beta_cl <- mod1@results@se.fixed[cl_param_idx] * sqrt(n / (n - nb_param))
  W1_gallant <- (-mod1@results@fixed.effects[cl_param_idx] + delta) / corrected_se_beta_cl
  W2_gallant <- (-mod1@results@fixed.effects[cl_param_idx] - delta) / corrected_se_beta_cl
  auc_gallant_be <- (W1_gallant >= quant_gallant) & (W2_gallant <= -quant_gallant)

  # TOST Cmax - Normal
  W1cmax <- (beta_cmax_estim + delta) / secmax
  W2cmax <- (beta_cmax_estim - delta) / secmax
  cmax_normal_be <- (W1cmax >= quant) & (W2cmax <= -quant)

  # TOST Cmax - Student
  cmax_student_be <- (W1cmax >= quant_stud) & (W2cmax <= -quant_stud)

  # TOST Cmax - Gallant
  corrected_se_beta_Cmax <- secmax * sqrt(n / (n - nb_param))
  W1_gallant_Cmax <- (beta_cmax_estim + delta) / corrected_se_beta_Cmax
  W2_gallant_Cmax <- (beta_cmax_estim - delta) / corrected_se_beta_Cmax
  cmax_gallant_be <- (W1_gallant_Cmax >= quant_gallant) & (W2_gallant_Cmax <= -quant_gallant)

  # Store results
  results$auc <- list(
    normal = list(be = auc_normal_be, W1 = W1, W2 = W2),
    student = list(be = auc_student_be, W1 = W1, W2 = W2),
    gallant = list(be = auc_gallant_be, W1 = W1_gallant, W2 = W2_gallant)
  )

  results$cmax <- list(
    normal = list(be = cmax_normal_be, W1 = W1cmax, W2 = W2cmax),
    student = list(be = cmax_student_be, W1 = W1cmax, W2 = W2cmax),
    gallant = list(be = cmax_gallant_be, W1 = W1_gallant_Cmax, W2 = W2_gallant_Cmax)
  )

  # Print informative summary
  cat("=== BIOEQUIVALENCE ANALYSIS RESULTS ===\n")
  cat("Model:", model_type, "\n")
  cat("Dataset:", basename(data_path), "\n")
  cat("Subjects:", n, "| Sampling times per subject:", nb_t, "\n\n")

  cat("AUC Bioequivalence Results:\n")
  cat("- Normal approach:  ", ifelse(auc_normal_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n")
  cat("- Student approach: ", ifelse(auc_student_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n")
  cat("- Gallant approach: ", ifelse(auc_gallant_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n\n")

  cat("Cmax Bioequivalence Results:\n")
  cat("- Normal approach:  ", ifelse(cmax_normal_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n")
  cat("- Student approach: ", ifelse(cmax_student_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n")
  cat("- Gallant approach: ", ifelse(cmax_gallant_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n\n")

  overall_be <- all(c(auc_gallant_be, cmax_gallant_be))
  cat("Overall Bioequivalence (Gallant method):", ifelse(overall_be, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"), "\n")
  cat("========================================\n")

  # Return comprehensive results
  return(list(
    saemix_object = mod1,
    res_tost = results,
    res_delta = SE_deltam,
    summary = list(
      model_type = model_type,
      n_subjects = n,
      n_timepoints = nb_t,
      n_parameters = nb_param,
      bioequivalence = list(
        auc = list(normal = auc_normal_be, student = auc_student_be, gallant = auc_gallant_be),
        cmax = list(normal = cmax_normal_be, student = cmax_student_be, gallant = cmax_gallant_be),
        overall_gallant = overall_be
      )
    )
  ))
}

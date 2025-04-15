#' Perform TOST Analysis on Pharmacokinetic Data
#'
#' This function performs TOST (Two One-Sided Tests) for bioequivalence analysis on pharmacokinetic data (AUC and Cmax),
#' using asymptotic SE (normal and student) and Gallant correction.
#'
#' @param data_path Path to the dataset (must be a `.txt` file).
#' @param dose_value Numeric value of the dose administered (default: 4).
#' @param nb_chains Number of MCMC chains for SAEM (default: 10).
#' @param nb_iter Vector of number of iterations for SAEM (default: c(300, 100)).
#' @return Prints TOST results for AUC and Cmax (Normal, Student, Gallant).
#' @export
perform_tost_analysis <- function(data_path, dose_value = 4, nb_chains = 10, nb_iter = c(300, 100)) {
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

  # One-compartment model with first-order absorption
  model1cpt <- function(psi, id, xidep) {
    dose <- xidep[,1]
    tim <- xidep[,2]
    ka <- psi[id,1]
    V <- psi[id,2]
    CL <- psi[id,3]
    k <- CL / V
    ypred <- dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    return(ypred)
  }

  # Load and prepare data
  tab <- read.table(data_path, header = TRUE, sep = " ", dec = ".", na.strings = ".")

  tab <- subset(tab, is.na(tab$Dose))
  tab$Dose <- dose_value
  tab <- tab[, c(1, 5, 2, 3, 4)]  # Reorder columns

  # Treatment coding
  tab$Treat <- ifelse(tab$Tr == "R", 0, 1)
  tab$Tr <- NULL
  tab <- tab[, c(1, 2, 3, 5, 4)]  # Reorder again

  nb_t <- sum(tab$Id[tab$Id == 1])
  n <- nrow(tab) / nb_t
  quant_stud <- qt(1 - 0.05, df = (n * nb_t - 5))
  quant_gallant <- qt(1 - 0.05, df = n - 3)

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

  saemix.modelb <- saemixModel(
    model = model1cpt,
    description = "One-compartment model with first-order absorption",
    psi0 = matrix(c(1.5, 0.5, 0.04, 0, 0, 0), ncol = 3, nrow = 2, byrow = TRUE,
                  dimnames = list(NULL, c("ka", "V", "CL"))),
    transform.par = c(1, 1, 1),
    covariate.model = matrix(c(1, 1, 1), ncol = 3),
    covariance.model = diag(3),
    omega.init = diag(c(0.05, 0.0125, 0.05)),
    error.model = "combined",
    error.init = c(0.1, 0.1)
  )

  mod1 <- saemix(saemix.modelb, saemix.data, list(nb.chains = nb_chains, nbiter.saemix = nb_iter))

  # Delta method for beta_Cmax
  estimates <- c(mod1@results@fixed.effects[1:6])
  names(estimates) <- c("Ka", "beta_Ka", "V", "beta_V", "Cl", "beta_Cl")
  fim_fixed <- mod1@results@fim[1:6, 1:6]
  varcov <- solve(fim_fixed)

  SE_deltam <- deltaMethod(estimates,
                           "-beta_V-(log((Ka*V)/Cl)+beta_Ka+beta_V-beta_Cl)*Cl*exp(beta_Cl)/(Ka*V*exp(beta_Ka+beta_V)-Cl*exp(beta_Cl))+(Cl/(Ka*V-Cl))*log(Ka*V/Cl)",
                           vcov. = varcov)
  beta_cmax_estim <- SE_deltam[,1]
  secmax <- SE_deltam[,2]

  # TOST AUC - Normal
  W1 <- (-mod1@results@fixed.effects[6] + delta) / mod1@results@se.fixed[6]
  W2 <- (-mod1@results@fixed.effects[6] - delta) / mod1@results@se.fixed[6]
  cat("TOST Normal sur AUC :", ifelse((W1 >= quant) & (W2 <= -quant), "Rejet de H0 : BE", "Non rejet de H0"), "\n")

  # TOST AUC - Student
  cat("TOST Student sur AUC :", ifelse((W1 >= quant_stud) & (W2 <= -quant_stud), "Rejet de H0 : BE", "Non rejet de H0"), "\n")

  # TOST AUC - Gallant
  corrected_se_beta_cl <- mod1@results@se.fixed[6] * sqrt(n / quant_gallant)
  W1_gallant <- (-mod1@results@fixed.effects[6] + delta) / corrected_se_beta_cl
  W2_gallant <- (-mod1@results@fixed.effects[6] - delta) / corrected_se_beta_cl
  cat("TOST Gallant sur AUC :", ifelse((W1_gallant >= quant_gallant) & (W2_gallant <= -quant_gallant), "Rejet de H0 : BE", "Non rejet de H0"), "\n")

  # TOST Cmax - Normal
  W1cmax <- (beta_cmax_estim + delta) / secmax
  W2cmax <- (beta_cmax_estim - delta) / secmax
  cat("TOST Normal sur Cmax :", ifelse((W1cmax >= quant) & (W2cmax <= -quant), "Rejet de H0 : BE", "Non rejet de H0"), "\n")

  # TOST Cmax - Student
  cat("TOST Student sur Cmax :", ifelse((W1cmax >= quant_stud) & (W2cmax <= -quant_stud), "Rejet de H0 : BE", "Non rejet de H0"), "\n")

  # TOST Cmax - Gallant
  corrected_se_beta_Cmax <- secmax * sqrt(n / quant_gallant)
  W1_gallant_Cmax <- (beta_cmax_estim + delta) / corrected_se_beta_Cmax
  W2_gallant_Cmax <- (beta_cmax_estim - delta) / corrected_se_beta_Cmax
  cat("TOST Gallant sur Cmax :", ifelse((W1_gallant_Cmax >= quant_gallant) & (W2_gallant_Cmax <= -quant_gallant), "Rejet de H0 : BE", "Non rejet de H0"), "\n")
}

#' Fonction principale pour l'analyse de bioéquivalence (TOST) sur AUC et Cmax
#'
#' Cette fonction effectue les tests de bioéquivalence (TOST) sur les paramètres pharmacocinétiques log(AUC) et log(Cmax)
#' pour une étude de bioéquivalence à deux bras parallèles. Le modèle utilisé est un modèle pharmacocinétique à un
#' compartiment ajusté avec l'algorithme SAEM via le package `saemix`.
#'
#' @param filepath Chemin vers le fichier de données contenant les observations (par défaut `NULL`).
#' @param delta Seuil d'équivalence (logarithmique), généralement \code{log(1.25)}. Ce seuil est utilisé pour définir
#'              les marges de bioéquivalence.
#' @param nb_chains Nombre de chaînes pour l'algorithme SAEM (par défaut \code{10}).
#' @param nb_iter Nombre d'itérations pour l'algorithme SAEM, spécifié sous forme de vecteur (par défaut \code{c(300, 100)}).
#'
#' @return Invisible. Une liste contenant le modèle et les résultats des tests TOST.
#' @import mvtnorm saemix car
#' @export
#'
#' @examples
#' \dontrun{
#' # Exécuter l'analyse avec le dataset par défaut
#' tost_bioequivalence_analysis()
#'
#' # Ou avec un chemin spécifique vers le fichier de données
#' # tost_bioequivalence_analysis("chemin/vers/donnees.txt")
#' }
tost_bioequivalence_analysis <- function(filepath = NULL, delta = log(1.25),
                                         nb_chains = 10, nb_iter = c(300, 100)) {
  # Chargement des bibliothèques nécessaires
  requireNamespace("mvtnorm", quietly = TRUE)
  requireNamespace("saemix", quietly = TRUE)
  requireNamespace("car", quietly = TRUE)

  # Si aucun fichier n'est fourni, utiliser le dataset par défaut dans le package
  if (is.null(filepath)) {
    # Chercher dans plusieurs emplacements standards
    potential_paths <- c(
      system.file("extdata", "dataset_1.txt", package = "MNLEM"),
      system.file("dataset_1.txt", package = "MNLEM"),
      system.file("data", "dataset_1.txt", package = "MNLEM"),
      system.file("inst/extdata", "dataset_1.txt", package = "MNLEM")
    )

    for (path in potential_paths) {
      if (path != "" && file.exists(path)) {
        filepath <- path
        break
      }
    }

    if (filepath == "" || is.null(filepath)) {
      stop("Le dataset par défaut 'dataset_1.txt' n'a pas été trouvé dans le package. ",
           "Veuillez spécifier le chemin complet vers votre fichier de données.")
    }

    message("Utilisation du fichier: ", filepath)
  } else {
    # Vérifier que le fichier spécifié existe
    if (!file.exists(filepath)) {
      stop("Le fichier spécifié '", filepath, "' n'existe pas.")
    }
  }

  # Définition du modèle PK à un compartiment avec absorption de premier ordre
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

  # Chargement et préparation des données
  message("Chargement des données...")
  tab <- read.table(filepath, header = TRUE, sep = " ", dec = ".", na.strings = ".")

  # Filtrage et préparation des données
  tab <- subset(tab, is.na(tab$Dose))
  tab$Dose <- 4
  tab <- tab[, c(1,5,2,3,4)]
  tab$Treat <- as.numeric(ifelse(tab$Tr == "R", 0, 1))
  tab <- tab[, c(1,2,3,5,4)]

  # Calcul des quantiles pour les tests statistiques
  nb_t <- sum(tab$Id == 1)
  n <- nrow(tab) / nb_t
  quant <- qnorm(0.95)
  quant_stud <- qt(1 - 0.05, df = (n * nb_t - 5))
  quant_gallant <- qt(1 - 0.05, df = n - 3)

  # Construction des objets saemix
  message("Préparation du modèle SAEM...")
  saemix.data <- saemix::saemixData(
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

  saemix.model <- saemix::saemixModel(
    model = model1cpt,
    description = "Modèle pharmacocinétique à 1 compartiment",
    psi0 = matrix(c(1.5, 0.5, 0.04, 0, 0, 0), ncol = 3, nrow = 2, byrow = TRUE,
                  dimnames = list(NULL, c("ka", "V", "CL"))),
    transform.par = c(1, 1, 1),
    covariate.model = matrix(c(1, 1, 1), ncol = 3),
    covariance.model = diag(3),
    omega.init = diag(c(0.05, 0.0125, 0.05)),
    error.model = "combined",
    error.init = c(0.1, 0.1)
  )

  # Ajustement du modèle
  message("\nAjustement du modèle en cours (peut prendre un moment)...")
  mod1 <- saemix::saemix(saemix.model, saemix.data,
                         list(nb.chains = nb_chains, nbiter.saemix = nb_iter))

  # Extraction des paramètres estimés et de la matrice d'information de Fisher
  message("Extraction des résultats...")
  estimates <- mod1@results@fixed.effects[1:6]
  names(estimates) <- c("Ka", "beta_Ka", "V", "beta_V", "Cl", "beta_Cl")
  fim_fixed <- mod1@results@fim[1:6, 1:6]
  varcov <- solve(fim_fixed)

  # Calcul de beta_Cmax avec la méthode delta
  message("Calcul des statistiques pour Cmax...")
  SE_deltam <- car::deltaMethod(estimates,
                                "-beta_V-(log((Ka*V)/Cl)+beta_Ka+beta_V-beta_Cl)*Cl*exp(beta_Cl)/(Ka*V*exp(beta_Ka+beta_V)-Cl*exp(beta_Cl))+(Cl/(Ka*V-Cl))*log(Ka*V/Cl)",
                                vcov. = varcov)
  beta_cmax_estim <- SE_deltam[,1]
  secmax <- SE_deltam[,2]

  # Exécution des tests TOST
  message("\nRésultats des tests TOST:")
  message("=======================================")

  # Sur AUC (beta_Cl)
  beta_cl <- estimates["beta_Cl"]
  se_beta_cl <- mod1@results@se.fixed[6]

  results_auc <- run_tost_tests(beta_cl, se_beta_cl, delta, quant, quant_stud, quant_gallant, n, "AUC")

  # Sur Cmax
  results_cmax <- run_tost_tests(beta_cmax_estim, secmax, delta, quant, quant_stud, quant_gallant, n, "Cmax")

  # Renvoyer les résultats de manière invisible
  return(invisible(list(
    model = mod1,
    auc_results = results_auc,
    cmax_results = results_cmax,
    estimates = estimates,
    beta_auc = beta_cl,
    se_auc = se_beta_cl,
    beta_cmax = beta_cmax_estim,
    se_cmax = secmax
  )))
}

#' Exécute les tests TOST pour un paramètre PK
#'
#' Cette fonction interne calcule les tests Two One-Sided Tests (TOST) pour évaluer
#' la bioéquivalence d'un paramètre pharmacocinétique en utilisant trois approches statistiques:
#' loi normale, t de Student et approximation de Gallant.
#'
#' @param beta Estimation du paramètre d'effet de traitement (sur échelle logarithmique).
#' @param se Erreur standard de l'estimation.
#' @param delta Marge d'équivalence (log(1.25) typiquement).
#' @param quant_norm Quantile de la loi normale.
#' @param quant_stud Quantile de la loi de Student.
#' @param quant_gallant Quantile de l'approximation de Gallant.
#' @param n Nombre de sujets.
#' @param param_name Nom du paramètre pour l'affichage des résultats.
#'
#' @return Liste des résultats des tests.
#' @keywords internal
run_tost_tests <- function(beta, se, delta, quant_norm, quant_stud, quant_gallant, n, param_name) {
  # Test TOST avec loi normale
  tost_norm_low <- (beta + delta) / se
  tost_norm_up <- (beta - delta) / se
  pval_norm_low <- pnorm(tost_norm_low)
  pval_norm_up <- 1 - pnorm(tost_norm_up)
  pval_norm <- max(pval_norm_low, pval_norm_up)

  # Test TOST avec distribution t de Student
  tost_stud_low <- (beta + delta) / se
  tost_stud_up <- (beta - delta) / se
  pval_stud_low <- pt(tost_stud_low, df = (n - 2))
  pval_stud_up <- 1 - pt(tost_stud_up, df = (n - 2))
  pval_stud <- max(pval_stud_low, pval_stud_up)

  # Test TOST avec distribution de Gallant
  tost_gallant_low <- (beta + delta) / se
  tost_gallant_up <- (beta - delta) / se
  pval_gallant_low <- pt(tost_gallant_low, df = (n - 3))
  pval_gallant_up <- 1 - pt(tost_gallant_up, df = (n - 3))
  pval_gallant <- max(pval_gallant_low, pval_gallant_up)

  # Affichage des résultats
  message("\n----- Test TOST pour ", param_name, " -----")
  message("Estimation de beta: ", round(beta, 6))
  message("Erreur standard: ", round(se, 6))
  message("Delta (seuil): ", round(delta, 6))

  message("\nTest avec loi normale:")
  message("Statistique de test basse: ", round(tost_norm_low, 4))
  message("Statistique de test haute: ", round(tost_norm_up, 4))
  message("P-valeur: ", round(pval_norm, 6))
  message("Conclusion: ", ifelse(pval_norm < 0.05, "Bioéquivalence démontrée", "Bioéquivalence non démontrée"))

  message("\nTest avec loi de Student:")
  message("Statistique de test basse: ", round(tost_stud_low, 4))
  message("Statistique de test haute: ", round(tost_stud_up, 4))
  message("P-valeur: ", round(pval_stud, 6))
  message("Conclusion: ", ifelse(pval_stud < 0.05, "Bioéquivalence démontrée", "Bioéquivalence non démontrée"))

  message("\nTest avec loi de Gallant:")
  message("Statistique de test basse: ", round(tost_gallant_low, 4))
  message("Statistique de test haute: ", round(tost_gallant_up, 4))
  message("P-valeur: ", round(pval_gallant, 6))
  message("Conclusion: ", ifelse(pval_gallant < 0.05, "Bioéquivalence démontrée", "Bioéquivalence non démontrée"))

  # Retourner les résultats sous forme de liste
  return(list(
    param = param_name,
    beta = beta,
    se = se,
    delta = delta,
    norm = list(
      tost_low = tost_norm_low,
      tost_up = tost_norm_up,
      pval = pval_norm,
      conclusion = ifelse(pval_norm < 0.05, "Bioéquivalence démontrée", "Bioéquivalence non démontrée")
    ),
    student = list(
      tost_low = tost_stud_low,
      tost_up = tost_stud_up,
      pval = pval_stud,
      conclusion = ifelse(pval_stud < 0.05, "Bioéquivalence démontrée", "Bioéquivalence non démontrée")
    ),
    gallant = list(
      tost_low = tost_gallant_low,
      tost_up = tost_gallant_up,
      pval = pval_gallant,
      conclusion = ifelse(pval_gallant < 0.05, "Bioéquivalence démontrée", "Bioéquivalence non démontrée")
    )
  ))
}

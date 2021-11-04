#' SMFOLL: Solution to Modeling the Odd Log-Logistic Family.
#'
#' The SMFOLL package provides important
#' functions and datasets
#' for survival analysis.
#' Some distributions of the Odd Log-Logistic Family.
#' (mara.caroline.torres@@uel.br)
#'
#' @section SMFOLL functions:
#' \code{\link[SMFOLL]{BetaBurrXII}}
#'
#' \code{\link[SMFOLL]{BXII}}
#'
#' \code{\link[SMFOLL]{descritiva}}
#'
#' \code{\link[SMFOLL]{Frechet}}
#'
#' \code{\link[SMFOLL]{HalfNormal}}
#'
#' \code{\link[SMFOLL]{KwBXII}}
#'
#' \code{\link[SMFOLL]{LogLogistic}}
#'
#' \code{\link[SMFOLL]{OLLBXII}}
#'
#' \code{\link[SMFOLL]{OLLexpo}}
#'
#' \code{\link[SMFOLL]{OllFit}}
#'
#' \code{\link[SMFOLL]{OLLFrechet}}
#'
#' \code{\link[SMFOLL]{OLLgamma}}
#'
#' \code{\link[SMFOLL]{OLLHN}}
#'
#' \code{\link[SMFOLL]{OLLLL}}
#'
#' \code{\link[SMFOLL]{OLLLN}}
#'
#' \code{\link[SMFOLL]{OLLnorm}}
#'
#' \code{\link[SMFOLL]{OllReg}}
#'
#' \code{\link[SMFOLL]{OLLWeibull}}
#'
#' \code{\link[SMFOLL]{test}}
#'
#' \code{\link[SMFOLL]{Weibull}}
#'
#'
#' @docType package
#' @name SMFOLL
NULL

#' Atuaria data
#'
#' Time data, the data set refers to the time of retired women with temporary disabilities incorporated into the Mexican public insurance system until they died in 2004, corresponding to 280 lives.
#'
#' The variable is as follow:
#'
#' \itemize{
#' \item x Variable of time
#' }
#'
#' @docType data
#' @keywords datasets
#' @name atuaria
#' @usage data(atuaria)
#' @format A data frame with 280 rows
NULL

#' Concrete data
#'
#' Time data
#'
#' The variable is as follow:
#'
#' \itemize{
#' \item x Variable of time
#' }
#'
#' @docType data
#' @keywords datasets
#' @name concrete
#' @usage data(concrete)
#' @format A data frame with 14 rows
NULL

#' Melanoma data
#'
#' Time data, the data are on cutaneous melanoma (a type of malignant cancer) that assesses the performance of a post-operative treatment with a high dose of the drug interferon alfa-2b to prevent a recurrence. Patients were included in the study from 1991 to 1995 and follow-up was carried out through 1998. Survival time X is referred to as the time to patient death.
#'
#' The variable is as follow:
#'
#' \itemize{
#' \item x Variable of time
#' \item cens Censure
#' \item x0 Covariate
#' \item x1 Covariate
#' \item x2 Covariate
#' \item x3 Covariate
#' \item x4 Covariate
#' \item x5 Covariate
#' \item x6 Covariate
#' \item obs Covariate
#' }
#'
#' @docType data
#' @keywords datasets
#' @name melanoma
#' @usage data(melanoma)
#' @format A data frame with 417 rows
NULL

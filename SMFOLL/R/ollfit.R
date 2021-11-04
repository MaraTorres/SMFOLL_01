#' Fit the data with Odd Log-Logistic distributions
#'
#' The function oll.fit shows some fitted results
#'
#' @param dist The name of the studied distributions may be 'fr', 'bxii', 'wei', 'll', 'expo', 'norm', 'gamma', 'ln', 'hn', 'ollfr', 'ollbxii', 'ollwei', 'ollll', 'ollexpo', 'ollnorm', 'ollgamma', 'ollln', 'ollhn', 'bbxii' or 'kwbxii'.
#' @param starts A vector with initial gess.
#' @param data,cens The data set.
#' @param method The method to be used, by default 'B'.
#' @param ...
#'
#' @name OllFit
#' @aliases oll.fit
#'
#' @return Fit data of a distribution.
#' @seealso \url{http://www.r-project.org}
#'
#' @examples
#' data(atuaria)
#' x=atuaria$x
#' oll.fit(dist = 'fr',starts = c(1,1),data = x)
#'
#' data(melanoma)
#' cens=melanoma$cens
#' x=melanoma$x
#' oll.fit(dist = 'ollwei',starts = c(1,1,1),x, cens = cens)

oll.fit = function(dist, starts, data, method = "B", domain = c(0, Inf), mle = NULL, cens = NULL) {
    if (is.null(cens) == TRUE) {
        if (dist == "fr")
            likelihood = vfr
        if (dist == "bxii")
            likelihood = vbxii
        if (dist == "wei")
            likelihood = vwei
        if (dist == "ll")
            likelihood = vll
        if (dist == "expo")
            likelihood = vexpo
        if (dist == "norm")
            likelihood = vnorm
        if (dist == "gamma")
            likelihood = vgamma
        if (dist == "ln")
            likelihood = vln
        if (dist == "hn")
            likelihood = vhn
        if (dist == "ollfr")
            likelihood = vollfr
        if (dist == "ollbxii")
            likelihood = vollbxii
        if (dist == "ollwei")
            likelihood = vollwei
        if (dist == "ollll")
            likelihood = vollll
        if (dist == "ollexpo")
            likelihood = vollexpo
        if (dist == "ollnorm")
            likelihood = vollnorm
        if (dist == "ollgamma")
            likelihood = vollgamma
        if (dist == "ollln")
            likelihood = vollln
        if (dist == "ollhn")
            likelihood = vollhn
        if (dist == "bbxii")
            likelihood = vbbxii
        if (dist == "kwbxii")
            likelihood = vkwbxii

        if (missingArg(data) == TRUE)
            stop("Database missing!")
        if (TRUE %in% is.nan(data) == TRUE)
            warning("The data have missing information!")
        if (length(domain) != 2)
            stop("The domain must have two arguments!")
        if (is.null(mle) == TRUE) {
            if (missingArg(starts) == TRUE)
                stop("The initial shots were not informed!")
        } else {
            starts = mle
        }
        if (is.null(mle) == TRUE) {
#            if (method == "PSO" || method == "P") {
#                result = pso(func = likelihood, ...)
#            }
            if (method == "Nelder-Mead" || method == "N") {
                result = optim(par = starts, fn = likelihood, x = data, method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if (method == "CG" || method == "C") {
                result = optim(par = starts, fn = likelihood, x = data, method = "CG", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if (method == "SANN" || method == "S") {
                result = optim(par = starts, fn = likelihood, x = data, method = "SANN", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if (method == "BFGS" || method == "B") {
                result = optim(par = starts, fn = likelihood, x = data, method = "BFGS", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if ((FALSE %in% (method != c("L", "BFGS", "B", "Nelder-Mead", "N", "SANN", "S", "CG", "C"))) == FALSE) {
                stop("Valid options are: BFGS or B, Nelder-Mead or N, SANN or S, CG or C.")
            }

            return(result)

        }
    } else {
        if (dist == "fr")
            likelihood = sfr
        if (dist == "bxii")
            likelihood = sbxii
        if (dist == "wei")
            likelihood = swei
        if (dist == "ll")
            likelihood = sll
        if (dist == "expo")
            likelihood = sexpo
        if (dist == "norm")
            likelihood = snorm
        if (dist == "gamma")
            likelihood = sgamma
        if (dist == "ln")
            likelihood = sln
        if (dist == "hn")
            likelihood = shn
        if (dist == "ollfr")
            likelihood = sollfr
        if (dist == "ollbxii")
            likelihood = sollbxii
        if (dist == "ollwei")
            likelihood = sollwei
        if (dist == "ollll")
            likelihood = sollll
        if (dist == "ollexpo")
            likelihood = sollexpo
        if (dist == "ollnorm")
            likelihood = sollnorm
        if (dist == "ollgamma")
            likelihood = sollgamma
        if (dist == "ollln")
            likelihood = sollln
        if (dist == "ollhn")
            likelihood = sollhn
        if (dist == "bbxii")
            likelihood = sbbxii
        if (dist == "kwbxii")
            likelihood = skwbxii

        if (missingArg(data) == TRUE)
            stop("Database missing!")
        if (TRUE %in% is.nan(data) == TRUE)
            warning("The data have missing information!")
        if (TRUE %in% is.nan(cens) == TRUE)
            warning("The cens have missing information!")
        if (any(cens > 1))
            stop("Status must be a vector of 0 and 1!")
        if (length(domain) != 2)
            stop("The domain must have two arguments!")
        if (is.null(mle) == TRUE) {
            if (missingArg(starts) == TRUE)
                stop("The initial shots were not informed!")
        } else {
            starts = mle
        }
        if (is.null(mle) == TRUE) {
#            if (method == "PSO" || method == "P") {
#                result = pso(func = likelihood, ...)
#            }
            if (method == "Nelder-Mead" || method == "N") {
                result = optim(par = starts, fn = likelihood, x = data, cens = cens, method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if (method == "CG" || method == "C") {
                result = optim(par = starts, fn = likelihood, x = data, cens = cens, method = "CG", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if (method == "SANN" || method == "S") {
                result = optim(par = starts, fn = likelihood, x = data, cens = cens, method = "SANN", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if (method == "BFGS" || method == "B") {
                result = optim(par = starts, fn = likelihood, x = data, cens = cens, method = "BFGS", hessian = TRUE, control = list(maxit = 1e+06))
            }
            if ((FALSE %in% (method != c("L", "BFGS", "B", "Nelder-Mead", "N", "SANN", "S", "CG", "C"))) == FALSE) {
                stop("Valid options are: BFGS or B, Nelder-Mead or N, SANN or S, CG or C.")
            }
            return(result)
        }
    }
}

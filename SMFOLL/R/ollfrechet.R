#' The Odd Log-Logistic Frechet Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Odd log-logistic Frechet distribution.
#'
#' @param param Vector with three parameters, first shape parameter, second scale parameter and third odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLFrechet
#' @aliases dollfr
#' @aliases pollfr
#' @aliases qollfr
#' @aliases rollfr
#'
#' @return dollfr gives the density, pollfr gives the distribution function, qollfr quantile function, rollfr random generation function.
#'
#' @examples
#' pollfr(c(2,2,3),c(5.1,5.6))
#'
#' dollfr(c(1,2,3),c(4,0.6))
#'
#' x <- rollfr(c(.91,6,2),1e+4)
#' plot(ecdf(x), lwd = 3)
#' curve(pollfr(c(.91,6,2),x), add = TRUE, col = 2, from = 0, lwd = 2)
#' legend('right', legend = c('Random', 'Theoretical'), lty = 1, col = 1:2, bty = 'n', lwd = 2)


dollfr <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {

        lambda <- param[1]
        sigma <- param[2]
        alpha <- param[3]

        u = (alpha * lambda * (sigma^lambda) * ((exp(-(sigma/x)^lambda))^alpha) * (1 - exp(-(sigma/x)^lambda))^(alpha - 1))/(x^(lambda +
            1) * (((exp(-(sigma/x)^lambda))^alpha) + (1 - exp(-(sigma/x)^lambda))^alpha)^2)

        if (log) {
            u = log(u)
        }

        return(u)

    }
}

#' @rdname OLLFrechet
#' @export


pollfr <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        sigma <- param[2]
        alpha <- param[3]

        x = (exp(-(sigma/q)^lambda)^alpha)/((exp(-(sigma/q)^lambda)^alpha) + (1 - exp(-(sigma/q)^lambda))^alpha)

        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)

    }
}

#' @rdname OLLFrechet
#' @export

qollfr <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {

        lambda <- param[1]
        sigma <- param[2]
        alpha <- param[3]

        if (lower.tail == FALSE) {
            p = 1 - p
        }

        x = sigma * (log(((1 - p)^(1/alpha) + p^(1/alpha))/(p^(1/alpha))))^(-(1/lambda))
        if (log) {
            x = log(x)
        }
        return(x)
    }
}

#' @rdname OLLFrechet
#' @export

rollfr <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        sigma <- param[2]
        alpha <- param[3]

        u = runif(n)
        x = sigma * (log(((1 - u)^(1/alpha) + u^(1/alpha))/(u^(1/alpha))))^(-(1/lambda))

        return(x)

    }
}

vollfr <- function(param, x) {
    lambda <- param[1]
    sigma <- param[2]
    alpha <- param[3]

    if (any(lambda < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- (alpha * lambda * (sigma^lambda) * ((exp(-(sigma/x)^lambda))^alpha) * (1 - exp(-(sigma/x)^lambda))^(alpha - 1))/(x^(lambda +
        1) * (((exp(-(sigma/x)^lambda))^alpha) + (1 - exp(-(sigma/x)^lambda))^alpha)^2)

    lv <- log(f)

    sum(-lv)
}

sollfr <- function(param, x, cens) {
    lambda <- param[1]
    sigma <- param[2]
    alpha <- param[3]

    if (any(lambda < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- (alpha * lambda * (sigma^lambda) * ((exp(-(sigma/x)^lambda))^alpha) * (1 - exp(-(sigma/x)^lambda))^(alpha - 1))/(x^(lambda +
        1) * (((exp(-(sigma/x)^lambda))^alpha) + (1 - exp(-(sigma/x)^lambda))^alpha)^2)

    S <- ((1 - exp(-(sigma/x)^lambda))^alpha)/((exp(-(sigma/x)^lambda)^alpha) + (1 - exp(-(sigma/x)^lambda))^alpha)

    lv <- cens * log(f) + (1 - cens) * log(S)

    sum(-lv)
}

#' The Weibull Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Weibull distribution.
#'
#' @param param Vector with scale parameter first and shape parameter second.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name Weibull
#' @aliases pwei
#' @aliases dwei
#' @aliases qwei
#' @aliases rwei
#'
#' @return  dwei gives the density, pwei gives the distribution function, qwei quantile function, rwei random generation function.
#'
#' @examples
#'
#' pwei(c(1,2),1)
#'
#' data(atuaria)
#' pwei(c(0.4,0.2),atuaria)
#'
#' dwei(c(0.1,0.99),0.6)
#'
#' x <- rwei(c(1,1),1e+3)
#' plot(ecdf(x), lwd = 3)
#' curve(pweibull(x,1,1), add = TRUE, col = 2, from = 0, lwd = 2)
#' legend('right', legend = c('Random', 'Theoretical'), lty = 1, col = 1:2, bty = 'n', lwd = 2)


dwei <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]

        u = (gamma/(lambda^gamma)) * (x^(gamma - 1)) * exp(-((x/lambda)^gamma))

        if (log) {
            u = log(u)
        }
        return(u)
    }

}

#' @rdname Weibull
#' @export

pwei <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]

        x = 1 - exp(-(q/lambda)^gamma)

        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
}

#' @rdname Weibull
#' @export

qwei <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]

        if (lower.tail == FALSE) {
            p = 1 - p
        }

        x = lambda * (-log(1 - p))^(1/gamma)
        if (log) {
            x = log(x)
        }
        return(x)
    }

}

#' @rdname Weibull
#' @export

rwei <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]

        u = runif(n)
        x = lambda * (-log(1 - u))^(1/gamma)
        return(x)
    }

}

vwei <- function(param, x) {
    lambda <- param[1]
    gamma <- param[2]

    if (any(lambda < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- (gamma/(lambda^gamma)) * (x^(gamma - 1)) * exp(-((x/lambda)^gamma))

    lv <- log(f)

    sum(-lv)

}

swei <- function(param, x, cens) {
    lambda <- param[1]
    gamma <- param[2]

    if (any(lambda < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20))
        return(.Machine$double.xmax^0.5)

    lv <- cens * log(gamma) - cens * gamma * log(lambda) + (gamma - 1) * cens * log(x) - (1/(lambda^gamma)) * (x^gamma)

    sum(-lv)

}

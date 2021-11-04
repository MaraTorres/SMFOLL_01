#' The Frechet Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Frechet distribution.
#'
#' @param param Vector with shape parameter first and scale parameter second.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name Frechet
#' @aliases dfr
#' @aliases pfr
#' @aliases qfr
#' @aliases rfr
#'
#' @return dfr gives the density, pfr gives the distribution function, qfr quantile function, rfr random generation function.
#'
#' @examples
#' pfr(c(4,2.1,3),5)
#'
#' dfr(c(2,3,3),c(2.1,1.2,0.6))

dfr <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        sigma <- param[2]

        u = (lambda/sigma) * (sigma/x)^(lambda + 1) * exp(-(sigma/x)^lambda)

        if (log) {
            u = log(u)
        }
        return(u)
    }

}

#' @rdname Frechet
#' @export

pfr <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        sigma <- param[2]

        x = exp(-(sigma/q)^lambda)

        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }

}

#' @rdname Frechet
#' @export

qfr <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        sigma <- param[2]

        if (lower.tail == FALSE) {
            p = 1 - p
        }

        x = sigma * (-log(p))^(-1/lambda)
        if (log) {
            x = log(x)
        }
        return(x)
    }

}

#' @rdname Frechet
#' @export

rfr <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        sigma <- param[2]

        u = runif(n)
        x = sigma * (-log(u))^(-1/lambda)
        return(x)
    }
}

vfr <- function(param, x) {
    lambda <- param[1]
    sigma <- param[2]

    if (any(lambda < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- (lambda/sigma) * (sigma/x)^(lambda + 1) * exp(-(sigma/x)^lambda)

    lv <- log(f)

    sum(-lv)
}

sfr <- function(param, x, cens) {
    lambda <- param[1]
    sigma <- param[2]

    if (any(lambda < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)

    lv <- cens * log((lambda/sigma) * (sigma/x)^(lambda + 1) * exp(-(sigma/x)^lambda)) + (1 - cens) * log(1 - exp(-(sigma/x)^lambda))

    sum(-lv)
}

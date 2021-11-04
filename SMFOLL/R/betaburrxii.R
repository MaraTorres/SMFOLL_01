#' The Beta Burr XII Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Beta Burr XII distribution.
#'
#' @param param Vector with five parameters, first scale parameter, the other ones are shape parameters.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name BetaBurrXII
#' @aliases dbbxii
#' @aliases pbbxii
#' @aliases qbbxii
#' @aliases rbbxii
#'
#' @return dbbxii gives the density, pbbxii gives the distribution function, qbbxii quantile function, rbbxii random generation function.
#'
#' @examples
#' pbbxii(c(1,1,1,1,1),5)
#'
#' dbbxii(c(5,5,5,5,5),c(2.1,3,2.5))

dbbxii <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {

        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]

        u = ((c * k * x^(c - 1))/((s^c) * beta(a, b))) * (1 + (x/s)^c)^(-(k * b + 1)) * (1 - (1 + (x/s)^c)^(-k))^(a - 1)


        if (log) {
            u = log(u)
        }
        return(u)
    }
}

#' @rdname BetaBurrXII
#' @export

pbbxii <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {

        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]

        x = pbeta((1 - (1 + (q/s)^c)^(-k)), a, b)
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }

}

#' @rdname BetaBurrXII
#' @export

qbbxii <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]

        if (a == 1 & b == 1) {

            x = qbxii(param = c(s, k, c), p = p, lower.tail = lower.tail, log = log)

        } else {

            if (lower.tail == FALSE) {
                p = 1 - p
            }

            V = 1 - dbeta(p, shape1 = a, shape2 = b)
            x = s * (((1/(V^(1/k))) - 1)^(1/c))

            if (log) {
                x = log(x)
            }
        }
        return(x)
    }
}

#' @rdname BetaBurrXII
#' @export

rbbxii <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]

        V = rbeta(n = n, shape1 = a, shape2 = b)
        x = s * (((1 - V)^(-1/k) - 1)^(1/c))
        return(x)
    }

}

vbbxii <- function(param, x) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    a <- param[4]
    b <- param[5]

    if (any(s < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(a < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(b < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- ((c * k * x^(c - 1))/((s^c) * beta(a, b))) * (1 + (x/s)^c)^(-(k * b + 1)) * (1 - (1 + (x/s)^c)^(-k))^(a - 1)

    lv <- log(f)
    sum(-lv)
}

sbbxii <- function(param, x, cens) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    a <- param[4]
    b <- param[5]

    if (any(s < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(a < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(b < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- (((c * k * x^(c - 1))/((s^c) * beta(a, b))) * (1 + (x/s)^c)^(-(k * b + 1)) * (1 - (1 + (x/s)^c)^(-k))^(a - 1))

    S <- (1 - pbeta((1 - (1 + (x/s)^c)^(-k)), a, b))

    lv <- cens * log(f) + (1 - cens) * log(S)

    sum(-lv)
}

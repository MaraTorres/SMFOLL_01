#' The Half-Normal Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the half-normal distribution.
#'
#' @param param Vector with shape parameter first and scale parameter second.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name HalfNormal
#' @aliases dhn
#' @aliases phn
#' @aliases qhn
#' @aliases rhn
#'
#' @return dhn gives the density, phn gives the distribution function, qhn quantile function, rhn random generation function.
#'
#' @examples
#' phn(c(.4,2.1,3),5)
#'
#' dhn(c(1,1),1)

dhn <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        
        u = exp(-0.5 * (x/theta)^(2 * lambda)) * sqrt(2/pi) * (lambda/x) * (x/theta)^lambda
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname HalfNormal
#' @export

phn <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        
        x = 2 * pnorm(q = (q/theta)^lambda, mean = 0, sd = 1) - 1
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname HalfNormal
#' @export

qhn <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = theta * (qnorm(p = ((p + 1)/2), mean = 0, sd = 1))^(1/lambda)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname HalfNormal
#' @export

rhn <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        
        u = runif(n)
        x = theta * (qnorm(p = ((u + 1)/2), mean = 0, sd = 1))^(1/lambda)
        return(x)
    }
    
}

vhn <- function(param, x) {
    lambda <- param[1]
    theta <- param[2]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(theta < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- exp(-0.5 * (x/theta)^(2 * lambda)) * sqrt(2/pi) * (lambda/x) * (x/theta)^lambda
    
    lv <- log(f)
    
    sum(-lv)
}

shn <- function(param, x, cens) {
    lambda <- param[1]
    theta <- param[2]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(theta < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (exp(-0.5 * (x/theta)^(2 * lambda)) * sqrt(2/pi) * (lambda/x) * (x/theta)^lambda)
    
    S <- 2 * (1 - pnorm(q = (x/theta)^lambda, mean = 0, sd = 1))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

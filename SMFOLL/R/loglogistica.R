#' The Log-Logistic Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Log-Logistic distribution.
#'
#' @param param Vector with scale parameter first and shape parameter second.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name LogLogistic
#' @aliases dll
#' @aliases pll
#' @aliases qll
#' @aliases rll
#'
#' @return dll gives the density, pll gives the distribution function, qll quantile function, rll random generation function.
#'
#' @examples
#' pll(c(4,2.1),5)
#'
#' dll(c(2,3),c(2.1,1.2,0.6))

dll <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        
        u = (gamma/lambda) * (x/lambda)^(gamma - 1) * (1 + (x/lambda)^gamma)^(-2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname LogLogistic
#' @export


pll <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        
        x = (1 + (q/lambda)^(-gamma))^(-1)
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname LogLogistic
#' @export

qll <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = lambda * ((1/p) - 1)^(-(1/gamma))
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname LogLogistic
#' @export

rll <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        
        u = runif(n)
        x = lambda * ((1/u) - 1)^(-(1/gamma))
        return(x)
    }
    
}

vll <- function(param, x) {
    lambda <- param[1]
    gamma <- param[2]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (gamma/lambda) * (x/lambda)^(gamma - 1) * (1 + (x/lambda)^gamma)^(-2)
    
    lv <- log(f)
    
    sum(-lv)
}

sll <- function(param, x, cens) {
    lambda <- param[1]
    gamma <- param[2]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    lv <- cens * log((gamma/lambda) * (x/lambda)^(gamma - 1) * (1 + (x/lambda)^gamma)^(-2)) + (1 - cens) * log(1 - ((1 + (x/lambda)^(-gamma))^(-1)))
    
    sum(-lv)
}

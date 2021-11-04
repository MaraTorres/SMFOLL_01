#' The Odd Log-Logistic exponential Distribution.
#'
#'Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic exponential distribution.
#'
#' @param param Vector with two parameters, first scale parameter and second the odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLexpo
#' @aliases dollexpo
#' @aliases pollexpo
#' @aliases qollexpo
#' @aliases rollexpo
#'
#' @return dollexpo gives the density, pollexpo gives the distribution function, qollexpo quantile function, rollexpo random generation function.
#'
#' @examples
#' pollexpo(c(1,2),c(4,0.6))
#'
#' dollexpo(c(1,2),c(4,0.6))

dollexpo <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        alpha <- param[2]
        
        u = (alpha * lambda * exp(lambda * x) * (exp(lambda * x) - 1)^(alpha - 1))/(((exp(lambda * x) - 1)^(alpha) + 1)^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname OLLexpo
#' @export

pollexpo <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        alpha <- param[2]
        
        x = ((exp(lambda * q) - 1)^alpha)/(((exp(lambda * q) - 1)^alpha) + 1)
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLexpo
#' @export

qollexpo <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        alpha <- param[2]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        x = 1 - exp(-lambda * ((p^(1/alpha))/(((1 - p)^(1/alpha) + (p^(1/alpha))))))
        
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLexpo
#' @export

rollexpo <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        alpha <- param[2]
        
        u = runif(n)
        x = 1 - exp(-lambda * ((u^(1/alpha))/(((1 - u)^(1/alpha) + (u^(1/alpha))))))
        return(x)
    }
    
}

vollexpo <- function(param, x) {
    lambda <- param[1]
    alpha <- param[2]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * lambda * exp(lambda * x) * (exp(lambda * x) - 1)^(alpha - 1))/(((exp(lambda * x) - 1)^(alpha) + 1)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollexpo <- function(param, x, cens) {
    lambda <- param[1]
    alpha <- param[2]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- ((alpha * lambda * exp(lambda * x) * (exp(lambda * x) - 1)^(alpha - 1))/(((exp(lambda * x) - 1)^(alpha) + 1)^2))
    
    S <- (1 - (((exp(lambda * x) - 1)^alpha)/(((exp(lambda * x) - 1)^alpha) + 1)))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

#' The Odd Log-Logistic Log-Logistic Distribution.
#'
#'Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic Log-Logistic distribution.
#'
#' @param param Vector with three parameters, first scale parameter, second shape parameter and third odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLLL
#' @aliases dollll
#' @aliases pollll
#' @aliases qollll
#' @aliases rollll
#'
#' @return dollll gives the density, pollll gives the distribution function, qollll quantile function, rollll random generation function.
#'
#' @examples
#' pollll(c(1,2,3),c(4,0.6))
#'
#' dollll(c(1,2,3),c(4,0.6))

dollll <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        u = (alpha * ((gamma/lambda) * (x/lambda)^(gamma - 1) * (1 + (x/lambda)^gamma)^(-2)) * ((1 + (x/lambda)^(-gamma))^(-1) * 
            (1 - (1 + (x/lambda)^(-gamma))^(-1)))^(alpha - 1))/(((1 + (x/lambda)^(-gamma))^(-alpha) + (1 - (1 + (x/lambda)^(-gamma))^(-1))^alpha)^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
}

#' @rdname OLLLL
#' @export

pollll <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        x = ((1 + (q/lambda)^(-gamma))^(-alpha))/(((1 + (q/lambda)^(-gamma))^(-alpha)) + (1 - (1 + (q/lambda)^(-gamma))^(-1))^alpha)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
}

#' @rdname OLLLL
#' @export

qollll <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = lambda * ((p^(1/alpha))/((1 - p)^(1/alpha)))^(1/gamma)
        
        if (log) {
            x = log(x)
        }
        return(x)
    }
}

#' @rdname OLLLL
#' @export

rollll <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        u = runif(n)
        x = lambda * ((u^(1/alpha))/((1 - u)^(1/alpha)))^(1/gamma)
        return(x)
    }
    
}

vollll <- function(param, x) {
    lambda <- param[1]
    gamma <- param[2]
    alpha <- param[3]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * ((gamma/lambda) * (x/lambda)^(gamma - 1) * (1 + (x/lambda)^gamma)^(-2)) * ((1 + (x/lambda)^(-gamma))^(-1) * (1 - 
        (1 + (x/lambda)^(-gamma))^(-1)))^(alpha - 1))/(((1 + (x/lambda)^(-gamma))^(-alpha) + (1 - (1 + (x/lambda)^(-gamma))^(-1))^alpha)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollll <- function(param, x, cens) {
    lambda <- param[1]
    gamma <- param[2]
    alpha <- param[3]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * ((gamma/lambda) * (x/lambda)^(gamma - 1) * (1 + (x/lambda)^gamma)^(-2)) * ((1 + (x/lambda)^(-gamma))^(-1) * (1 - 
        (1 + (x/lambda)^(-gamma))^(-1)))^(alpha - 1))/(((1 + (x/lambda)^(-gamma))^(-alpha) + (1 - (1 + (x/lambda)^(-gamma))^(-1))^alpha)^2)
    
    S <- ((1 - ((1 + (x/lambda)^(-gamma))^(-1)))^alpha)/(((1 + (x/lambda)^(-gamma))^(-alpha)) + (1 - ((1 + (x/lambda)^(-gamma))^(-1)))^alpha)
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

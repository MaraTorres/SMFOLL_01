#' The Odd Log-Logistic normal Distribution.
#'
#'Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic normal distribution.
#'
#' @param param Vector with three parameters, first scale parameter, second shape parameter and third odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLnorm
#' @aliases dollnorm
#' @aliases pollnorm
#' @aliases qollnorm
#' @aliases rollnorm
#'
#' @return dollnorm gives the density, pollnorm gives the distribution function, qollnorm quantile function, rollnorm random generation function.
#'
#' @examples
#' pollnorm(c(1,2,3),c(4,0.6))
#'
#' dollnorm(c(1,2,3),c(4,0.6))

dollnorm <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        mi <- param[1]
        sigma <- param[2]
        alpha <- param[3]
        
        u = (alpha * dnorm(x = x, mean = mi, sd = sigma) * (pnorm(q = x, mean = mi, sd = sigma) * (1 - pnorm(q = x, mean = mi, sd = sigma)))^(alpha - 
            1))/(((pnorm(q = x, mean = mi, sd = sigma)^alpha) + (1 - pnorm(q = x, mean = mi, sd = sigma))^alpha)^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname OLLnorm
#' @export

pollnorm <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        mi <- param[1]
        sigma <- param[2]
        alpha <- param[3]
        
        x = (pnorm(q = q, mean = mi, sd = sigma)^alpha)/((pnorm(q = q, mean = mi, sd = sigma)^alpha) + (1 - pnorm(q = q, mean = mi, 
            sd = sigma))^alpha)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLnorm
#' @export

qollnorm <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        mi <- param[1]
        sigma <- param[2]
        alpha <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = qnorm(((p^(1/alpha))/(((1 - p)^(1/alpha)) + (p^(1/alpha)))), mean = mi, sd = sigma)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLnorm
#' @export

rollnorm <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        mi <- param[1]
        sigma <- param[2]
        alpha <- param[3]
        
        u = runif(n)
        x = qnorm(((u^(1/alpha))/(((1 - u)^(1/alpha)) + (u^(1/alpha)))), mean = mi, sd = sigma)
        return(x)
    }
    
}

vollnorm <- function(param, x) {
    mi <- param[1]
    sigma <- param[2]
    alpha <- param[3]
    
    if (any(mi < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * dnorm(x = x, mean = mi, sd = sigma) * (pnorm(q = x, mean = mi, sd = sigma) * (1 - pnorm(q = x, mean = mi, sd = sigma)))^(alpha - 
        1))/(((pnorm(q = x, mean = mi, sd = sigma)^alpha) + (1 - pnorm(q = x, mean = mi, sd = sigma))^alpha)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollnorm <- function(param, x, cens) {
    mi <- param[1]
    sigma <- param[2]
    alpha <- param[3]
    
    if (any(mi < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * dnorm(x = x, mean = mi, sd = sigma) * (pnorm(q = x, mean = mi, sd = sigma) * (1 - pnorm(q = x, mean = mi, sd = sigma)))^(alpha - 
        1))/(((pnorm(q = x, mean = mi, sd = sigma)^alpha) + (1 - pnorm(q = x, mean = mi, sd = sigma))^alpha)^2)
    
    S <- (1 - ((pnorm(q = x, mean = mi, sd = sigma)^alpha)/((pnorm(q = x, mean = mi, sd = sigma)^alpha) + (1 - pnorm(q = x, mean = mi, 
        sd = sigma))^alpha)))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}


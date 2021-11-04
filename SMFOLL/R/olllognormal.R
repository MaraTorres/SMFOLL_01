#' The Odd Log-Logistic Log-Normal Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic Log-Normal distribution.
#'
#' @param param Vector with three parameters, first scale parameter, second shape parameter and third odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLLN
#' @aliases dollln
#' @aliases pollln
#' @aliases qollln
#' @aliases rollln
#'
#' @return dollln gives the density, pollln gives the distribution function, qollln quantile function, rollln random generation function.
#'
#' @examples
#' pollln(c(1,2,3),c(4,0.6))
#'
#' dollln(c(1,2,3),c(4,0.6))

dollln <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        meanlog <- param[1]
        sdlog <- param[2]
        alpha <- param[3]
        
        u = (alpha * dlnorm(x = x, meanlog = meanlog, sdlog = sdlog) * (plnorm(q = x, meanlog = meanlog, sdlog = sdlog) * (1 - plnorm(q = x, 
            meanlog = meanlog, sdlog = sdlog)))^(alpha - 1))/(((plnorm(q = x, meanlog = meanlog, sdlog = sdlog)^alpha) + (1 - plnorm(q = x, 
            meanlog = meanlog, sdlog = sdlog))^alpha)^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname OLLLN
#' @export

pollln <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        meanlog <- param[1]
        sdlog <- param[2]
        alpha <- param[3]
        
        x = (plnorm(q = q, meanlog = meanlog, sdlog = sdlog)^alpha)/((plnorm(q = q, meanlog = meanlog, sdlog = sdlog)^alpha) + (1 - 
            plnorm(q = q, meanlog = meanlog, sdlog = sdlog))^alpha)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLLN
#' @export

qollln <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        meanlog <- param[1]
        sdlog <- param[2]
        alpha <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = qlnorm(((p^(1/alpha))/(((1 - p)^(1/alpha)) + (p^(1/alpha)))), meanlog = meanlog, sdlog = sdlog)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLLN
#' @export

rollln <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        meanlog <- param[1]
        sdlog <- param[2]
        alpha <- param[3]
        
        
        u = runif(n)
        x = qlnorm(((u^(1/alpha))/(((1 - u)^(1/alpha)) + (u^(1/alpha)))), meanlog = meanlog, sdlog = sdlog)
        return(x)
    }
    
}

vollln <- function(param, x) {
    meanlog <- param[1]
    sdlog <- param[2]
    alpha <- param[3]
    
    if (any(meanlog < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(sdlog < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * dlnorm(x = x, meanlog = meanlog, sdlog = sdlog) * (plnorm(q = x, meanlog = meanlog, sdlog = sdlog) * (1 - plnorm(q = x, 
        meanlog = meanlog, sdlog = sdlog)))^(alpha - 1))/(((plnorm(q = x, meanlog = meanlog, sdlog = sdlog)^alpha) + (1 - plnorm(q = x, 
        meanlog = meanlog, sdlog = sdlog))^alpha)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollln <- function(param, x, cens) {
    meanlog <- param[1]
    sdlog <- param[2]
    alpha <- param[3]
    
    if (any(meanlog < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(sdlog < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- ((alpha * dlnorm(x = x, meanlog = meanlog, sdlog = sdlog) * (plnorm(q = x, meanlog = meanlog, sdlog = sdlog) * (1 - plnorm(q = x, 
        meanlog = meanlog, sdlog = sdlog)))^(alpha - 1))/(((plnorm(q = x, meanlog = meanlog, sdlog = sdlog)^alpha) + (1 - plnorm(q = x, 
        meanlog = meanlog, sdlog = sdlog))^alpha)^2))
    
    S <- (1 - ((plnorm(q = x, meanlog = meanlog, sdlog = sdlog)^alpha)/((plnorm(q = x, meanlog = meanlog, sdlog = sdlog)^alpha) + 
        (1 - plnorm(q = x, meanlog = meanlog, sdlog = sdlog))^alpha)))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

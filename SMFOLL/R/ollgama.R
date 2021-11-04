#' The Odd Log-Logistic gamma Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic gamma distribution.
#'
#' @param param Vector with three parameters, first shape parameter, second scale parameter and third odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLgamma
#' @aliases dollgamma
#' @aliases pollgamma
#' @aliases qollgamma
#' @aliases rollgamma
#'
#' @return dollgamma gives the density, pollgamma gives the distribution function, qollgamma quantile function, rollgamma random generation function.
#'
#' @examples
#' pollgamma(c(1,2,3),c(4,0.6))
#'
#' dollgamma(c(1,2,3),c(4,0.6))

dollgamma <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        shape <- param[1]
        scale <- param[2]
        alpha <- param[3]
        
        u = (alpha * dgamma(x = x, shape = shape, scale = scale) * (pgamma(q = x, shape = shape, scale = scale) * (1 - pgamma(q = x, 
            shape = shape, scale = scale)))^(alpha - 1))/(((pgamma(q = x, shape = shape, scale = scale)^alpha) + (1 - pgamma(q = x, 
            shape = shape, scale = scale))^alpha)^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname OLLgamma
#' @export

pollgamma <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        shape <- param[1]
        scale <- param[2]
        alpha <- param[3]
        
        x = (pgamma(q = q, shape = shape, scale = scale)^alpha)/((pgamma(q = q, shape = shape, scale = scale)^alpha) + (1 - pgamma(q = q, 
            shape = shape, scale = scale))^alpha)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLgamma
#' @export

qollgamma <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        shape <- param[1]
        scale <- param[2]
        alpha <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = qgamma(((p^(1/alpha))/(((1 - p)^(1/alpha)) + (p^(1/alpha)))), shape = shape, scale = scale)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLgamma
#' @export

rollgamma <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        shape <- param[1]
        scale <- param[2]
        alpha <- param[3]
        
        u = runif(n)
        x = qgamma(((u^(1/alpha))/(((1 - u)^(1/alpha)) + (u^(1/alpha)))), shape = shape, scale = scale)
        return(x)
    }
    
}

vollgamma <- function(param, x) {
    shape <- param[1]
    scale <- param[2]
    alpha <- param[3]
    
    if (any(shape < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(scale < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * dgamma(x = x, shape = shape, scale = scale) * (pgamma(q = x, shape = shape, scale = scale) * (1 - pgamma(q = x, 
        shape = shape, scale = scale)))^(alpha - 1))/(((pgamma(q = x, shape = shape, scale = scale)^alpha) + (1 - pgamma(q = x, 
        shape = shape, scale = scale))^alpha)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollgamma <- function(param, x, cens) {
    shape <- param[1]
    scale <- param[2]
    alpha <- param[3]
    
    if (any(shape < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(scale < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- ((alpha * dgamma(x = x, shape = shape, scale = scale) * (pgamma(q = x, shape = shape, scale = scale) * (1 - pgamma(q = x, 
        shape = shape, scale = scale)))^(alpha - 1))/(((pgamma(q = x, shape = shape, scale = scale)^alpha) + (1 - pgamma(q = x, 
        shape = shape, scale = scale))^alpha)^2))
    
    S <- (1 - ((pgamma(q = x, shape = shape, scale = scale)^alpha)/((pgamma(q = x, shape = shape, scale = scale)^alpha) + (1 - pgamma(q = x, 
        shape = shape, scale = scale))^alpha)))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

#' The Odd Log-Logistic Weibull Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic Weibull distribution, with three parameters.
#'
#' @param param Vector with three parameters, first the scale parameter, second the shape parameter and the last one the odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLWeibull
#' @aliases pollwei
#' @aliases dollwei
#' @aliases qollwei
#' @aliases rollwei
#'
#' @return  dollwei gives the density, pollwei gives the distribution function, qollwei quantile function, rollwei random generation function.
#'
#' @examples
#' pollwei(c(1,2,3),c(1.2,0.6))
#'
#' dollwei(c(1,2,3),c(1.1,0.6))
#'
#' x <- rollwei(c(1,1,1),1e+3)
#' plot(ecdf(x), lwd = 3)
#' curve(pollwei(c(1,1,1),x), add = TRUE, col = 2, from = 0, lwd = 2)
#' legend('right', legend = c('Aleatória', 'Teórica'), lty = 1, col = 1:2, bty = 'n', lwd = 2)

dollwei <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        u = (alpha * gamma * (x^(gamma - 1)) * ((exp(-(x/lambda)^gamma))^alpha) * (1 - exp(-((x/lambda)^gamma)))^(alpha - 1))/((lambda^(gamma)) * 
            (((1 - exp(-((x/lambda)^gamma)))^alpha) + (exp(-((x/lambda)^gamma)))^alpha)^2)
        if (log) {
            u = log(u)
        }
        return(u)
    }
}

#' @rdname OLLWeibull
#' @export

pollwei <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        x = ((1 - exp(-((q/lambda)^gamma)))^alpha)/(((1 - exp(-((q/lambda)^gamma)))^alpha) + (exp(-((q/lambda)^gamma))^alpha))
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLWeibull
#' @export

qollwei <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = lambda * (log(((p^(1/alpha) + (1 - p)^(1/alpha))/((1 - p)^(1/alpha)))))^(1/gamma)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLWeibull
#' @export

rollwei <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        gamma <- param[2]
        alpha <- param[3]
        
        u = runif(n)
        x = lambda * (log(((u^(1/alpha) + (1 - u)^(1/alpha))/((1 - u)^(1/alpha)))))^(1/gamma)
        return(x)
    }
    
}

vollwei <- function(param, x) {
    lambda <- param[1]
    gamma <- param[2]
    alpha <- param[3]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * gamma * (x^(gamma - 1)) * ((exp(-(x/lambda)^gamma))^alpha) * (1 - exp(-((x/lambda)^gamma)))^(alpha - 1))/((lambda^(gamma)) * 
        (((1 - exp(-((x/lambda)^gamma)))^alpha) + (exp(-((x/lambda)^gamma)))^alpha)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollwei <- function(param, x, cens) {
    lambda <- param[1]
    gamma <- param[2]
    alpha <- param[3]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(gamma < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- ((alpha * gamma * (x^(gamma - 1)) * ((exp(-(x/lambda)^gamma))^alpha) * (1 - exp(-((x/lambda)^gamma)))^(alpha - 1))/((lambda^(gamma)) * 
        (((1 - exp(-((x/lambda)^gamma)))^alpha) + (exp(-((x/lambda)^gamma)))^alpha)^2))
    
    S <- (1 - ((1 - exp(-((x/lambda)^gamma)))^alpha)/(((1 - exp(-((x/lambda)^gamma)))^alpha) + (exp(-((x/lambda)^gamma))^alpha)))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

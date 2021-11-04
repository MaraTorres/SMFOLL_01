#' The Odd Log-Logistic Half-Normal Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic Half-Normal distribution.
#'
#' @param param Vector with three parameters, first shape parameter, second scale parameter and third odd parameter.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name OLLHN
#' @aliases dollhn
#' @aliases pollhn
#' @aliases qollhn
#' @aliases rollhn
#'
#' @return dollhn gives the density, pollhn gives the distribution function, qollhn quantile function, rollhn random generation function.
#'
#' @examples
#' pollhn(c(1,2,3),c(4,0.6))
#'
#' dollhn(c(1,2,3),c(4,0.6))

dollhn <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        alpha <- param[3]
        
        u = (alpha * (exp(-0.5 * (x/theta)^(2 * lambda)) * sqrt(2/pi) * (lambda/x) * (x/theta)^lambda) * ((2 * pnorm(q = (x/theta)^lambda, 
            mean = 0, sd = 1) - 1) * (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1)))^(alpha - 1))/((((2 * pnorm(q = (x/theta)^lambda, 
            mean = 0, sd = 1) - 1)^alpha) + (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1))^alpha)^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname OLLHN
#' @export

pollhn <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        alpha <- param[3]
        
        x = ((2 * pnorm(q = (q/theta)^lambda, mean = 0, sd = 1) - 1)^alpha)/(((2 * pnorm(q = (q/theta)^lambda, mean = 0, sd = 1) - 
            1)^alpha) + (1 - (2 * pnorm(q = (q/theta)^lambda, mean = 0, sd = 1) - 1))^alpha)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLHN
#' @export

qollhn <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        alpha <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        y = (p^(1/alpha))/(((1 - p)^(1/alpha)) + (p^(1/alpha)))
        x = theta * (qnorm(p = ((y + 1) * 0.5), mean = 0, sd = 1))^(1/lambda)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLHN
#' @export

rollhn <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        lambda <- param[1]
        theta <- param[2]
        alpha <- param[3]
        
        
        u = runif(n)
        y = (u^(1/alpha))/(((1 - u)^(1/alpha)) + (u^(1/alpha)))
        x = theta * (qnorm(p = ((y + 1) * 0.5), mean = 0, sd = 1))^(1/lambda)
        return(x)
    }
    
}

vollhn <- function(param, x) {
    lambda <- param[1]
    theta <- param[2]
    alpha <- param[3]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(theta < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (alpha * (exp(-0.5 * (x/theta)^(2 * lambda)) * sqrt(2/pi) * (lambda/x) * (x/theta)^lambda) * ((2 * pnorm(q = (x/theta)^lambda, 
        mean = 0, sd = 1) - 1) * (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1)))^(alpha - 1))/((((2 * pnorm(q = (x/theta)^lambda, 
        mean = 0, sd = 1) - 1)^alpha) + (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1))^alpha)^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollhn <- function(param, x, cens) {
    lambda <- param[1]
    theta <- param[2]
    alpha <- param[3]
    
    if (any(lambda < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(theta < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(alpha < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- ((alpha * (exp(-0.5 * (x/theta)^(2 * lambda)) * sqrt(2/pi) * (lambda/x) * (x/theta)^lambda) * ((2 * pnorm(q = (x/theta)^lambda, 
        mean = 0, sd = 1) - 1) * (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1)))^(alpha - 1))/((((2 * pnorm(q = (x/theta)^lambda, 
        mean = 0, sd = 1) - 1)^alpha) + (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1))^alpha)^2))
    
    S <- (1 - (((2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1)^alpha)/(((2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 
        1)^alpha) + (1 - (2 * pnorm(q = (x/theta)^lambda, mean = 0, sd = 1) - 1))^alpha)))
    
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

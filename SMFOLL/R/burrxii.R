#' The Burr XII Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Burr XII distribution.
#'
#' @param param Vector with three parameters, first the scale parameter and then the two shape parameters.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
#'
#' @name BXII
#' @aliases dbxii
#' @aliases pbxii
#' @aliases qbxii
#' @aliases rbxii
#'
#' @return dbxii gives the density, pbxii gives the distribution function, qbxii quantile function, rbxii random generation function.
#'
#' @examples
#' pbxii(c(4,2.1,3),5)
#'
#' dbxii(c(2,3,3),c(2.1,1.2,0.6))
#'
#' x <- rbxii(c(1,1,1),1e+3)
#' plot(ecdf(x), lwd = 3)
#' curve(pbxii(c(1,1,1),x), add = TRUE, col = 2, from = 0, lwd = 2)
#' legend('right', legend = c('Random', 'Theoretical'), lty = 1, col = 1:2, bty = 'n', lwd = 2)
#'

dbxii <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        
        u = c * k * (s^(-c)) * (1 + ((x/s)^c))^(-k - 1) * x^(c - 1)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

#' @rdname BXII
#' @export

pbxii <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        
        x = 1 - (1 + ((q/s)^c))^(-k)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname BXII
#' @export

qbxii <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = s * ((1 - p)^(-1/k) - 1)^(1/c)
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname BXII
#' @export

rbxii <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        
        u = runif(n)
        x = s * ((1 - u)^(-1/k) - 1)^(1/c)
        return(x)
    }
    
}

vbxii <- function(param, x) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    
    if (any(s < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- c * k * (s^(-c)) * (1 + ((x/s)^c))^(-k - 1) * x^(c - 1)
    
    lv <- log(f)
    
    sum(-lv)
    
}

sbxii <- function(param, x, cens) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    
    if (any(s < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    lv <- cens * log(c * k * (s^(-c)) * (1 + ((x/s)^c))^(-k - 1) * x^(c - 1)) + (1 - cens) * log((1 + (x/s)^c)^(-k))
    
    sum(-lv)
    
}

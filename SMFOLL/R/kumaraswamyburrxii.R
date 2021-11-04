#' The Kumaraswamy Burr XII Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Kumaraswamy Burr XII distribution.
#'
#' @param param Vector with five parameters, first scale parameter, the other ones are shape parameters.
#' @param x Vector of quantiles.
#'
#' @name KwBXII
#' @aliases dkwbxii
#' @aliases pkwbxii
#' @aliases qkwbxii
#' @aliases rkwbxii
#'
#' @return dkwbxii gives the density, pkwbxii gives the distribution function, qkwbxii quantile function, rkwbxii random generation function.
#'
#' @examples
#' pkwbxii(c(1,1,1,1,1),5)
#'
#' dkwbxii(c(2,2,2,2,2),c(2.1,1.2,0.6))

dkwbxii <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {
        
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]
        
        u = a * b * c * k * s^(-c) * x^(c - 1) * (1 + (x/s)^c)^(-k - 1) * (1 - (1 + (x/s)^c)^(-k))^(a - 1) * (1 - (1 - (1 + (x/s)^c)^(-k))^a)^(b - 
            1)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
}

#' @rdname KwBXII
#' @export

pkwbxii <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {
        
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]
        
        x = 1 - (1 - (1 - (1 + (q/s)^c)^(-k))^a)^b
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname KwBXII
#' @export

qkwbxii <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = s * ((1 - (1 - (1 - p)^(1/b))^(1/a))^(-1/k) - 1)^(1/c)
        
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname KwBXII
#' @export

rkwbxii <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0 | param[5] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        b <- param[5]
        
        u = runif(n)
        x = s * ((1 - (1 - (1 - u)^(1/b))^(1/a))^(-1/k) - 1)^(1/c)
        
        return(x)
    }
    
}

vkwbxii <- function(param, x) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    a <- param[4]
    b <- param[5]
    
    if (any(s < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(a < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(b < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- a * b * c * k * s^(-c) * x^(c - 1) * (1 + (x/s)^c)^(-k - 1) * (1 - (1 + (x/s)^c)^(-k))^(a - 1) * (1 - (1 - (1 + (x/s)^c)^(-k))^a)^(b - 
        1)
    
    lv <- log(f)
    sum(-lv)
    
}

skwbxii <- function(param, x, cens) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    a <- param[4]
    b <- param[5]
    
    if (any(s < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(a < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(b < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (a * b * c * k * s^(-c) * x^(c - 1) * (1 + (x/s)^c)^(-k - 1) * (1 - (1 + (x/s)^c)^(-k))^(a - 1) * (1 - (1 - (1 + (x/s)^c)^(-k))^a)^(b - 
        1))
    
    S <- (1 - (1 - (1 + (x/s)^c)^(-k))^a)^b
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

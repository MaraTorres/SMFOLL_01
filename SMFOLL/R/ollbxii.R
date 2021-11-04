#' The Odd Log-Logistic Burr XII Distribution.
#'
#' Provides density, distribution function, quantile function and random generation for the Odd Log-Logistic Burr XII distribution.
#'
#' @param param Vector with four parameters, first the scale parameter, then the two shape parameters and the last one the odd parameter.
#' @param x Vector of quantiles.
#'
#' @name OLLBXII
#' @aliases dollbxii
#' @aliases pollbxii
#' @aliases qollbxii
#' @aliases rollbxii
#'
#' @return dollbxii gives the density, pollbxii gives the distribution function, qollbxii quantile function, rollbxii random generation function.
#'
#'
#' @examples
#' pollbxii(c(4,2.1,3,4),c(4,5))
#'
#' dollbxii(c(2,3,3,3),c(2.1,1.2,0.6))

dollbxii <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0) {
        return(NaN)
    } else {
        
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        
        u = (a * c * k * x^(c - 1) * ((1 - (1 + (x/s)^c)^(-k)) * ((1 + (x/s)^c)^(-k)))^(a - 1))/((s^c) * ((1 + ((x/s)^c))^(k + 1)) * 
            (((1 - (1 + (x/s)^c)^(-k))^a) + (1 + (x/s)^c)^(-k * a))^2)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
}

#' @rdname OLLBXII
#' @export

pollbxii <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0) {
        return(NaN)
    } else {
        
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        
        x = ((1 - (1 + (q/s)^c)^(-k))^a)/(((1 - (1 + (q/s)^c)^(-k))^a) + (1 + (q/s)^c)^(-k * a))
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLBXII
#' @export

qollbxii <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        
        x = s * ((((1 - ((p^(1/a))/((1 - p)^(1/a) + p^(1/a))))^(-1/k)) - 1)^(1/c))
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

#' @rdname OLLBXII
#' @export

rollbxii <- function(param, n) {
    if (param[1] <= 0 | param[2] <= 0 | param[3] <= 0 | param[4] <= 0) {
        return(NaN)
    } else {
        s <- param[1]
        k <- param[2]
        c <- param[3]
        a <- param[4]
        
        u = runif(n)
        x = s * ((((1 - ((u^(1/a))/((1 - u)^(1/a) + u^(1/a))))^(-1/k)) - 1)^(1/c))
        return(x)
    }
    
}

vollbxii <- function(param, x) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    a <- param[4]
    
    if (any(s < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(a < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (a * c * k * x^(c - 1) * ((1 - (1 + (x/s)^c)^(-k)) * ((1 + (x/s)^c)^(-k)))^(a - 1))/((s^c) * ((1 + ((x/s)^c))^(k + 1)) * 
        (((1 - (1 + (x/s)^c)^(-k))^a) + (1 + (x/s)^c)^(-k * a))^2)
    
    lv <- log(f)
    
    sum(-lv)
}

sollbxii <- function(param, x, cens) {
    s <- param[1]
    k <- param[2]
    c <- param[3]
    a <- param[4]
    
    if (any(s < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(k < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(c < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(a < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (a * c * k * x^(c - 1) * ((1 - (1 + (x/s)^c)^(-k)) * ((1 + (x/s)^c)^(-k)))^(a - 1))/((s^c) * ((1 + ((x/s)^c))^(k + 1)) * 
        (((1 - (1 + (x/s)^c)^(-k))^a) + (1 + (x/s)^c)^(-k * a))^2)
    
    S <- ((1 + (x/s)^c)^(-k * a))/(((1 - (1 + (x/s)^c)^(-k))^a) + (1 + (x/s)^c)^(-k * a))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

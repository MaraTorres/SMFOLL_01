dexpo <- function(param, x, log = FALSE) {
    if (param <= 0) {
        return(NaN)
    } else {
        
        u = param * exp(-param * x)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

pexpo <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param <= 0) {
        return(NaN)
    } else {
        
        x = 1 - exp(-param * q)
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

qexpo <- function(param, p, lower.tail = TRUE, log = FALSE) {
    if (param <= 0) {
        return(NaN)
    } else {
        
        if (lower.tail == FALSE) {
            p = 1 - p
        }
        x = -(1/param) * log(1 - p)
        
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

rexpo <- function(param, n) {
    if (param <= 0) {
        return(NaN)
    } else {
        u = runif(n)
        x = -(1/param) * log(1 - u)
        return(x)
    }
    
}

vexpo <- function(param, x) {
    
    if (any(param < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- param * exp(-param * x)
    
    lv <- log(f)
    
    sum(-lv)
}

sexpo <- function(param, x, cens) {
    
    if (any(param < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- (param * exp(-param * x))
    
    S <- exp(-param * x)
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

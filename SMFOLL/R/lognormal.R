dln <- function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        mu <- param[1]
        s <- param[2]
        
        u = exp(-0.5 * ((log(x) - mu)/s)^2)/(sqrt(2 * pi) * x * s)
        
        if (log) {
            u = log(u)
        }
        return(u)
    }
    
}

pln <- function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        mu <- param[1]
        s <- param[2]
        
        x = 1 - pnorm(q = ((mu - log(q))/s), mean = 0, sd = 1)
        
        if (lower.tail == FALSE) {
            x = 1 - x
        }
        if (log) {
            x = log(x)
        }
        return(x)
    }
    
}

vln <- function(param, x) {
    mean <- param[1]
    sd <- param[2]
    
    if (any(mean < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(sd < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- dlnorm(x = x, meanlog = mean, sdlog = sd)
    
    lv <- log(f)
    
    sum(-lv)
}

sln <- function(param, x, cens) {
    mean <- param[1]
    sd <- param[2]
    
    if (any(mean < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    if (any(sd < 1e-20)) 
        return(.Machine$double.xmax^0.5)
    
    f <- dlnorm(x = x, meanlog = mean, sdlog = sd)
    
    S <- (1 - plnorm(q = x, meanlog = mean, sdlog = sd))
    
    lv <- cens * log(f) + (1 - cens) * log(S)
    
    sum(-lv)
}

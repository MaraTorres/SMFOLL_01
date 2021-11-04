dnormal = function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        mu <- param[1]
        sigma <- param[2]

        u <- dnorm(x = x, mean = mu, sd = sigma)

        if (log) {
            u = log(u)
        }
        return(u)
    }
}

pnormal = function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        mu <- param[1]
        sigma <- param[2]

        x <- pnorm(q = q, mean = mu, sd = sigma)

        if (lower.tail == FALSE) {
            x = 1 - x
        }

        if (log) {
            x = log(x)
        }
        return(x)
    }
}

vnorm <- function(param, x) {
    mu <- param[1]
    sigma <- param[2]

    if (any(mu < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- dnorm(x = x, mean = mu, sd = sigma)

    lv <- log(f)

    sum(-lv)
}

snorm <- function(param, x, cens) {
    mu <- param[1]
    sigma <- param[2]

    if (any(mu < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- dnorm(x = x, mean = mu, sd = sigma)

    S <- (1 - pnorm(q = x, mean = mu, sd = sigma))

    lv <- cens * log(f) + (1 - cens) * log(S)

    sum(-lv)
}

dg = function(param, x, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        shape <- param[1]
        scale <- param[2]

        u <- dgamma(x = q, shape = shape, scale = scale)

        if (log) {
            u = log(u)
        }
        return(u)
    }
}

pg = function(param, q, lower.tail = TRUE, log = FALSE) {
    if (param[1] <= 0 | param[2] <= 0) {
        return(NaN)
    } else {
        shape <- param[1]
        scale <- param[2]

        x <- pgamma(q = q, shape = shape, scale = scale)

        if (lower.tail == FALSE) {
            x = 1 - x
        }

        if (log) {
            x = log(x)
        }
        return(x)
    }
}

vgamma <- function(param, x) {
    shape <- param[1]
    scale <- param[2]

    if (any(shape < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(scale < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- dgamma(x = x, shape = shape, scale = scale)

    lv <- log(f)

    sum(-lv)
}

sgamma <- function(param, x, cens) {
    shape <- param[1]
    scale <- param[2]

    if (any(shape < 1e-20))
        return(.Machine$double.xmax^0.5)
    if (any(scale < 1e-20))
        return(.Machine$double.xmax^0.5)

    f <- dgamma(x = x, shape = shape, scale = scale)

    S <- (1 - pgamma(q = x, shape = shape, scale = scale))

    lv <- cens * log(f) + (1 - cens) * log(S)

    sum(-lv)
}

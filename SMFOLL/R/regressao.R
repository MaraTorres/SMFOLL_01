#' Regression to Odd Log-Logistic distributions
#'
#' The function oll.reg shows some regression results
#'
#' @param dist The name of the studied distributions may be 'll', 'bxii', 'ollbxii' or 'ollghn'.
#' @param starts A vector with initial gess.
#' @param y,X,delta The data set.
#' @param ...
#'
#' @name OllReg
#' @aliases oll.reg
#'
#' @return Fit data of a distribution.
#' @seealso \url{http://www.r-project.org}
#'
#' @examples
#' dados <- read.table("http://users.stat.ufl.edu/~winner/data/alc_elim.dat", sep = "",header = FALSE)
#' head(dados)
#' y <- (with(dados, V4))
#' delta <- rep(1, length(y))
#' x0 <- rep(1, length(y))
#' x1 <- (with(dados, V3))
#' x2 <- (with(dados, V2))
#' X = as.matrix(cbind(x0,x1,x2))
#'
#' opt1 = oll.reg(dist = 'll',starts = c(1,1,1,1),y = y,X = X,delta = delta, plot = TRUE); opt1
#' opt2 = oll.reg(dist = 'bxii',starts = c(1,1,1,1,1),y = y,X = X,delta = delta, plot = TRUE); opt2
#' opt3 = oll.reg('ollbxii',c(opt2$par, .9),y,X,delta, plot = TRUE); opt3

oll.reg = function(dist, starts, y, X, delta, domain = c(0, Inf), plot = T) {
  # Geral
  criterios = function(value, d, n) {
    # value -valor do optim d - numero de parametros n - numero da dados
    l = 2 * value
    AIC = l + 2 * d
    BIC = l + d * log(n)
    CAIC = AIC + (2 * (d + 2) * (d + 3))/(n - d - 3)
    AICc = AIC + 2 * (d * (d + 1))/(n - d - 1)
    resul = cbind(l, AIC, CAIC, BIC, AICc)
    return(resul)
  }
  IC = function(parametros, hessiana, n) {
    np = length(parametros)
    Var = solve(hessiana)
    EP = diag(sqrt(Var))
    tvalue = parametros/EP  # cálculo da estatística z
    pval <- 2 * (1 - pt(abs(tvalue), n - np))  # cálculo do valor p teste de Wald
    LI = parametros - qt(0.975, n) * EP
    LS = parametros + qt(0.975, n) * EP
    valorp = pt(tvalue, n, lower.tail = FALSE)
    resul = cbind(parametros, EP, tvalue, pval, valorp, LI, LS)
    return(resul)
  }

  if (dist == "ll"){
    # Regressão LL
    logs <- function(param, y, X, delta) {
      n <- nrow(X)
      j <- ncol(X)

      beta <- param[1:j]
      sigma <- param[j + 1]
      if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)

      xbeta <- X %*% beta
      lv <- rep(0, n)
      for (i in 1:n) {
        z <- (y[i] - xbeta[i])/sigma
        G <- 1 + exp(z)
        f <- (1/sigma) * exp(z) * G^(-2)
        S <- G^(-1)
        if (delta[i] == 1) {
          lv[i] <- log(f)
        } else {
          lv[i] <- log(S)
        }
      }
      soma = -sum(lv)
      return(soma)
    }
    res = function(parametros, y, X, delta, B = 20) {
      # B número de simulações do envelope
      ps = length(parametros)
      n = length(y)
      par = parametros
      D <- matrix(0, nrow = nrow(X), ncol = B)
      random <- function(n) {
        u = runif(n)
        z = log((1/(1 - u)) - 1)

        return(z)
      }  # LL


      mu <- X %*% parametros[1:(ps - 1)] # LL

      z = 1 + exp((y - mu)/parametros[ps]) # LL

      Fll = (1 - z^(-1)) # LL

      Martingale = delta + log(1 - Fll) # LL

      sinal = ifelse(Martingale < 0, -1, 1)

      Deviance = sinal * (-2 * (Martingale + delta * log(delta - Martingale)))^(0.5)
      D[, 1] = Deviance

      yhat = mu + parametros[ps] # LL

      resul = cbind(y, Martingale = Martingale, Deviance = Deviance, yhat)
      colnames(resul) = c("y", "M", "D", "yhat")

      plot(Deviance, xlab = "Index", ylab = "Deviance residual Modificada")
      plot(yhat, Deviance, xlab = "Valores ajustados", ylab = "Deviance residual Modificada")

      for (i in 2:B) {
        ynovo = mu + parametros[ps] * random(n)  # LL

        estimativa <- optim(par = parametros, y = ynovo, X = X, delta = delta,
                            logs, method = "BFGS", hessian = TRUE) # LL
        par = estimativa$par

        mu.novo <- X %*% par[1:ncol(X)]

        z = 1 + exp((y - mu.novo)/par[ps])

        Fll = (1 - z^(-1))  # LL

        Martingale = delta + log(1 - Fll)  # LL

        sinal = ifelse(Martingale < 0, -1, 1)

        D[, i] = sinal * (-2 * (Martingale + delta * log(delta - Martingale)))^(0.5)

      }

      Dordenado <- (apply(D, 2, sort))
      Z <- qnorm(((1:n) - 3/8)/(n + 1/4))
      Zmeio <- qnorm((n + (1:n) + 1/2)/(2 * n + 9/8))
      Dm <- apply(Dordenado, 1, mean)
      Dmin <- apply(Dordenado, 1, min)
      Dmax <- apply(Dordenado, 1, max)
      Dplot <- cbind(Z, Zmeio, Dmin, Dm, Dmax)

      # Envelope
      par(mai = c(1.2, 1.2, 0.5, 0.1))
      plot(Dplot[, 1], Dordenado[, 1])
      lines(Dplot[, 1], Dplot[, 4], col = 2, lwd = 2, lty = 2)
      lines(Dplot[, 1], Dplot[, 3], lty = 2)
      lines(Dplot[, 1], Dplot[, 5], lty = 2)


      return(as.data.frame(resul))
    }
  }
  if (dist == "bxii"){
    # Regressão LBXII
    logs <- function(param, y, X, delta) {
      n <- nrow(X)
      j <- ncol(X)

      beta <- param[1:j]
      sigma <- param[j + 1]
      k <- param[j + 2]
      if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)
      if (any(k < 1e-20))
        return(.Machine$double.xmax^0.5)

      xbeta <- X %*% beta
      lv <- rep(0, n)
      for (i in 1:n) {
        z <- (y[i] - xbeta[i])/sigma
        G <- 1 + exp(z)
        f <- (k/sigma) * exp(z) * G^(-(k + 1))
        S <- G^(-k)
        if (delta[i] == 1) {
          lv[i] <- log(f)
        } else {
          lv[i] <- log(S)
        }
      }
      soma = -sum(lv)
      return(soma)
    }
    res = function(parametros, y, X, delta, B = 20) {
      # B número de simulações do envelope
      ps = length(parametros)
      n = length(y)
      par = parametros
      D <- matrix(0, nrow = nrow(X), ncol = B)
      random <- function(k, n) {
        u = runif(n)
        z = log((1 - u)^(-(1/k)) - 1)

        return(z)
      }  # LBXII


      mu <- X %*% parametros[1:(ps - 2)] # LBXII

      z = 1 + exp((y - mu)/parametros[(ps-1)]) # LBXII

      k = parametros[ps] # LBXII

      Flbxii = (1 - z^(-k)) # LBXII

      Martingale = delta + log(1 - Flbxii) # LBXII

      sinal = ifelse(Martingale < 0, -1, 1)

      Deviance = sinal * (-2 * (Martingale + delta * log(delta - Martingale)))^(0.5)
      D[, 1] = Deviance

      yhat = mu + parametros[(ps-1)] + k  # LBXII

      resul = cbind(y, Martingale = Martingale, Deviance = Deviance, yhat)
      colnames(resul) = c("y", "M", "D", "yhat")

      plot(Deviance, xlab = "Index", ylab = "Deviance residual Modificada")
      plot(yhat, Deviance, xlab = "Valores ajustados", ylab = "Deviance residual Modificada")

      for (i in 2:B) {
        ynovo = mu + parametros[ps - 1] * random(k, n) # LBXII

        estimativa <- optim(par = parametros, y = ynovo, X = X, delta = delta,
                            logs, method = "BFGS", hessian = TRUE) # LBXII
        par = estimativa$par

        mu.novo <- X %*% par[1:ncol(X)]

        z = 1 + exp((y - mu.novo)/par[ps - 1])

        Flbxii = (1 - z^(-par[ps]))  # LBXII

        Martingale = delta + log(1 - Flbxii)  # LBXII

        sinal = ifelse(Martingale < 0, -1, 1)

        D[, i] = sinal * (-2 * (Martingale + delta * log(delta - Martingale)))^(0.5)


      }

      Dordenado <- (apply(D, 2, sort))
      Z <- qnorm(((1:n) - 3/8)/(n + 1/4))
      Zmeio <- qnorm((n + (1:n) + 1/2)/(2 * n + 9/8))
      Dm <- apply(Dordenado, 1, mean)
      Dmin <- apply(Dordenado, 1, min)
      Dmax <- apply(Dordenado, 1, max)
      Dplot <- cbind(Z, Zmeio, Dmin, Dm, Dmax)

      # Envelope
      par(mai = c(1.2, 1.2, 0.5, 0.1))
      plot(Dplot[, 1], Dordenado[, 1])
      lines(Dplot[, 1], Dplot[, 4], col = 2, lwd = 2, lty = 2)
      lines(Dplot[, 1], Dplot[, 3], lty = 2)
      lines(Dplot[, 1], Dplot[, 5], lty = 2)




      return(as.data.frame(resul))
    }
  }
  if (dist == "ollbxii"){
    # Regressão LOLLBXII
    logs <- function(param, y, X, delta) {
      n <- nrow(X)
      j <- ncol(X)

      beta <- param[1:j]
      sigma <- param[j + 1]
      k <- param[j + 2]
      alpha <- param[j + 3]
      if (any(sigma < 1e-20))
        return(.Machine$double.xmax^0.5)
      if (any(k < 1e-20))
        return(.Machine$double.xmax^0.5)
      if (any(alpha < 1e-20))
        return(.Machine$double.xmax^0.5)

      xbeta <- X %*% beta
      lv <- rep(0, n)
      for (i in 1:n) {
        z <- (y[i] - xbeta[i])/sigma
        G <- 1 + exp(z)
        f <- ((alpha * k)/sigma) * exp(z) * ((G^(-(k + 1)) * ((1 - G^(-k)) * G^(-k))^(alpha - 1))/(((1 - G^(-k))^(alpha) + G^(-k * alpha))^2))
        S <- (G^(-k * alpha))/((1 - G^(-k * alpha))^alpha + (G^(-k * alpha))^(-alpha * k))
        if (delta[i] == 1) {
          lv[i] <- log(f)
        } else {
          lv[i] <- log(S)
        }
      }
      soma = -sum(lv)
      return(soma)
    }
    res = function(parametros, y, X, delta, B = 20) {
      # B número de simulações do envelope
      ps = length(parametros)
      n = length(y)
      par = parametros
      D <- matrix(0, nrow = nrow(X), ncol = B)
      random <- function(k, a, n) {
        u = runif(n)
        z = log((((-u/(u - 1))^(1/a) + 1)^(1/k)) - 1)

        return(z)
      }  # LOLLBXII


      mu <- X %*% parametros[1:ncol(X)] # LOLLBXII

      z = 1 + exp((y - mu)/parametros[ps-2]) # LOLLBXII

      k = parametros[ps - 1] # LOLLBXII

      a = parametros[ps] # LOLLBXII

      Flollbxii = ((1 - z^(-k))^a)/((1 - z^(-k))^a + z^(-k * a)) # LOLLBXII

      Martingale = delta + log(1 - Flollbxii) # LOLLBXII

      sinal = ifelse(Martingale < 0, -1, 1)

      Deviance = sinal * (-2 * (Martingale + delta * log(delta - Martingale)))^(0.5)
      D[, 1] = Deviance

      yhat = mu + parametros[ps-2] + k + a # LOLLBXII

      resul = cbind(y, Martingale = Martingale, Deviance = Deviance, yhat)
      colnames(resul) = c("y", "M", "D", "yhat")

      plot(Deviance, xlab = "Index", ylab = "Deviance residual Modificada", las = 1)
      plot(yhat, Deviance, xlab = "Valores ajustados", ylab = "Deviance residual Modificada", las = 1)


      for (i in 2:B) {
        ynovo = mu + parametros[ps - 2] * random(k, a, n)  # LOLLBXII

        estimativa <- optim(par = parametros, y = ynovo, X = X, delta = delta,
                            fn = logs, method = "BFGS", hessian = TRUE) # LOLLBXII
        par = estimativa$par
        mu.novo <- X %*% par[1:ncol(X)]

        z = 1 + exp((y - mu.novo)/par[ncol(X) + 1])

        k.novo = par[ps - 1]  # LOLLBXII
        a.novo = par[ps]  # LOLLBXII
        Flollbxii = (1 - z^(-k.novo))^a.novo/((1 - z^(-k.novo))^a.novo + z^(-k.novo * a.novo)) # LOLLBXII

        Martingale = delta + log(1 - Flollbxii)  # LOLLBXII

        sinal = ifelse(Martingale < 0, -1, 1)

        D[, i] = sinal * (-2 * (Martingale + delta * log(delta - Martingale)))^(0.5)

      }

      Dordenado <- (apply(D, 2, sort))
      Z <- qnorm(((1:n) - 3/8)/(n + 1/4))
      Zmeio <- qnorm((n + (1:n) + 1/2)/(2 * n + 9/8))
      Dm <- apply(Dordenado, 1, mean)
      Dmin <- apply(Dordenado, 1, min)
      Dmax <- apply(Dordenado, 1, max)
      Dplot <- cbind(Z, Zmeio, Dmin, Dm, Dmax)

      # Envelope
      #par(mai = c(1.2, 1.2, 0.5, 0.1))
      plot(Dplot[, 1], Dordenado[, 1], xlab = "Percentil", ylab = "Resíduo Componente do Desvio", las = 1)
      lines(Dplot[, 1], Dplot[, 4], col = 2, lwd = 2, lty = 2)
      lines(Dplot[, 1], Dplot[, 3], col = "darkgreen", lty = 2)
      lines(Dplot[, 1], Dplot[, 5], col = "darkgreen", lty = 2)

      return(as.data.frame(resul))
    }
  }

  if (missingArg(y) == TRUE)
    stop("y missing!")
  if (missingArg(X) == TRUE)
    stop("X missing!")
  if (missingArg(y) == TRUE)
    stop("delta missing!")
  if (TRUE %in% is.nan(y) == TRUE)
    warning("The data have missing information!")
  if (TRUE %in% is.nan(X) == TRUE)
    warning("The data have missing information!")
  if (TRUE %in% is.nan(delta) == TRUE)
    warning("The data have missing information!")
  if (length(domain) != 2)
    stop("The domain must have two arguments!")

  n <- length(delta)

  opt <- optim(par = starts, y = y, X = X, delta = delta,
                logs, method = "BFGS", hessian = TRUE)
  opt$'confidence interval' = IC(opt$par, opt$hessian, n)
  opt$'criteria' = criterios(opt$value, length(opt$par), n)

  if (plot == T)
    res = res(opt$par, y, X, delta)


  return(opt)
}

#' Hypothesis test.
#'
#' Provides the comparison between two fits
#'
#' @param x Vector of quantiles.
#' @param dist1 First distribution.
#' @param result1 First result.
#' @param dist2 Second distribution.
#' @param result2 Second result.
#' @param status Data set with or without censure.
#'
#' @return Statistics, IC and Hypothesis test.
#'
#' @examples
#' data(atuaria)
#' result1=oll.fit(dist = 'bxii',starts = c(1,1,1),data = atuaria$x)
#' test(atuaria$x,'bxii',result1)
#' result2=oll.fit(dist = 'ollbxii',starts = c(result1$par,1),data = atuaria$x)
#' test(atuaria$x,'bxii',result1,'ollbxii',result2)

test = function(x, dist1, result1, dist2 = NULL, result2 = NULL, status = F) {
    if (TRUE %in% is.nan(x) == TRUE) 
        warning("The data have missing information!")
    if (is.null(dist2) == F & is.null(result2) == T) 
        stop("Missing result2. The list must be informed.")
    if (is.null(dist2) == T & is.null(result2) == F) 
        stop("Missing dist2!")
    
    
    if (dist1 == "fr") {
        pdf1 = dfr
        cdf1 = pfr
    }
    if (dist1 == "bxii") {
        pdf1 = dbxii
        cdf1 = pbxii
    }
    if (dist1 == "wei") {
        pdf1 = dwei
        cdf1 = pwei
    }
    if (dist1 == "ll") {
        pdf1 = dll
        cdf1 = pll
    }
    if (dist1 == "expo") {
        pdf1 = dexpo
        cdf1 = pexpo
    }
    if (dist1 == "norm") {
        pdf1 = dnormal
        cdf1 = pnormal
    }
    if (dist1 == "gamma") {
        pdf1 = dg
        cdf1 = pg
    }
    if (dist1 == "ln") {
        pdf1 = dln
        cdf1 = pln
    }
    if (dist1 == "hn") {
        pdf1 = dhn
        cdf1 = phn
    }
    if (dist1 == "ollfr") {
        pdf1 = dollfr
        cdf1 = pollfr
    }
    if (dist1 == "ollbxii") {
        pdf1 = dollbxii
        cdf1 = pollbxii
    }
    if (dist1 == "ollwei") {
        pdf1 = dollwei
        cdf1 = pollwei
    }
    if (dist1 == "ollll") {
        pdf1 = dollll
        cdf1 = pollll
    }
    if (dist1 == "ollexpo") {
        pdf1 = dollexpo
        cdf1 = pollexpo
    }
    if (dist1 == "ollnorm") {
        pdf1 = dollnorm
        cdf1 = pollnorm
    }
    if (dist1 == "ollgamma") {
        pdf1 = dollgamma
        cdf1 = pollgamma
    }
    if (dist1 == "ollln") {
        pdf1 = dollln
        cdf1 = pollln
    }
    if (dist1 == "ollhn") {
        pdf1 = dollhn
        cdf1 = pollhn
    }
    if (dist1 == "bbxii") {
        pdf1 = dbbxii
        cdf1 = pbbxii
    }
    if (dist1 == "kwbxii") {
        pdf1 = dkwbxii
        cdf1 = pkwbxii
    }
    
    if (is.null(dist2) == F) {
        if (dist2 == "fr") {
            pdf2 = dfr
            cdf2 = pfr
        }
        if (dist2 == "bxii") {
            pdf2 = dbxii
            cdf2 = pbxii
        }
        if (dist2 == "wei") {
            pdf2 = dwei
            cdf2 = pwei
        }
        if (dist2 == "ll") {
            pdf2 = dll
            cdf2 = pll
        }
        if (dist2 == "expo") {
            pdf2 = dexpo
            cdf2 = pexpo
        }
        if (dist2 == "norm") {
            pdf2 = dnormal
            cdf2 = pnormal
        }
        if (dist2 == "gamma") {
            pdf2 = dg
            cdf2 = pg
        }
        if (dist2 == "ln") {
            pdf2 = dln
            cdf2 = pln
        }
        if (dist2 == "hn") {
            pdf2 = dhn
            cdf2 = phn
        }
        if (dist2 == "ollfr") {
            pdf2 = dollfr
            cdf2 = pollfr
        }
        if (dist2 == "ollbxii") {
            pdf2 = dollbxii
            cdf2 = pollbxii
        }
        if (dist2 == "ollwei") {
            pdf2 = dollwei
            cdf2 = pollwei
        }
        if (dist2 == "ollll") {
            pdf2 = dollll
            cdf2 = pollll
        }
        if (dist2 == "ollexpo") {
            pdf2 = dollexpo
            cdf2 = pollexpo
        }
        if (dist2 == "ollnorm") {
            pdf2 = dollnorm
            cdf2 = pollnorm
        }
        if (dist2 == "ollgamma") {
            pdf2 = dollgamma
            cdf2 = pollgamma
        }
        if (dist2 == "ollln") {
            pdf2 = dollln
            cdf2 = pollln
        }
        if (dist2 == "ollhn") {
            pdf2 = dollhn
            cdf2 = pollhn
        }
        if (dist2 == "bbxii") {
            pdf2 = dbbxii
            cdf2 = pbbxii
        }
        if (dist2 == "kwbxii") {
            pdf2 = dkwbxii
            cdf2 = pkwbxii
        }
    }
    u = function(x) {
        pnorm((qnorm(x) - mean(qnorm(x)))/sqrt(var(qnorm(x))))
    }
    
    A = function(x) {
        (-n - (1/n) * sum((2 * i - 1) * (log(x)) + (2 * n + 1 - 2 * i) * (log(1 - x)))) * (1 + (0.75/n) + (2.25/(n^(2))))
    }
    
    W = function(x) {
        (sum((x - ((2 * i - 1)/(2 * n)))^2) + (1/(12 * n))) * (1 + (0.5/n))
    }
    
    criterios = function(value, d, n) {
        # value -valor do optim d - numero de parametros n - numero da dados
        l = 2 * value
        AIC = l + 2 * d
        BIC = l + d * log(n)
        CAIC = AIC + (2 * (d + 2) * (d + 3))/(n - d - 3)
        resul = cbind(l, AIC, CAIC, BIC)
        return(resul)
    }
    
    IC = function(parametros, hessiana, n) {
        Var = solve(hessiana)
        EP = diag(sqrt(Var))
        tvalue = parametros/EP
        LI = parametros - qt(0.975, n) * EP
        LS = parametros + qt(0.975, n) * EP
        valorp = pt(tvalue, n, lower.tail = FALSE)
        resul = cbind(parametros, EP, tvalue, valorp, LI, LS)
        return(resul)
    }
    
    x = sort(x)
    n = length(x)
    i = seq(1, n)
    tests = list()
    
    # dist1
    f1 <- pdf1(param = result1[[1]], x = x)
    lv1 <- sum(log(f1))
    F1 <- cdf1(result1[[1]], x)
    
    
    M1 = matrix(c(A(u(F1)), W(u(F1))), nrow = 1)
    colnames(M1) = c("A", "W")
    
    if (status == FALSE) {
        tests$statistics1 = cbind(M1, criterios(result1[[2]], length(result1[[1]]), n))
    } else {
        tests$statistics1 = cbind(criterios(result1[[2]], length(result1[[1]]), n))
    }
    tests$IC1 = IC(result1[[1]], result1[[6]], n = n)
    
    if (is.null(dist2) == F & is.null(result2) == F) {
        
        f2 <- pdf2(param = result2[[1]], x = x)
        lv2 <- sum(log(f2))
        F2 <- cdf2(result2[[1]], x)
        
        M2 = matrix(c(A(u(F2)), W(u(F2))), nrow = 1)
        colnames(M2) = c("A", "W")
        if (status == FALSE) {
            tests$statistics2 = cbind(M2, criterios(result2[[2]], length(result2[[1]]), n))
        } else {
            tests$statistics2 = cbind(criterios(result2[[2]], length(result2[[1]]), n))
        }
        tests$IC2 = IC(result2[[1]], result2[[6]], n = n)
        
        pvalor = pchisq(2 * abs(lv2 - lv1), 1, ncp = 0, lower.tail = F)
        LR = matrix(c(2 * abs(lv2 - lv1), pvalor), nrow = 1)
        colnames(LR) = c("w", "pvalue")
        tests$LR = LR
        
    }
    
    return(tests)
}

#' descritiva - Calculation of some descriptive statistics.
#'
#' The function descritiva shows some descriptive statistics of a numeric vector
#'
#' @param x,cens A numeric vector.
#' @param colb Color box-plot.
#' @param colh Color histogram.
#' @param xlab Axis label.
#' @param breaks Values of the breaks.
#'
#' @return Some descriptive statistics of a numeric vector
#'
#' @seealso \url{http://www.r-project.org}
#'
#' @examples
#'
#' data(atuaria)
#' descritiva(atuaria$x)
#'

descritiva <- function(x, colb = "blue", colh = "green", xlab = "values of the dataset", breaks = 40,
    main = "", cens = FALSE, b = NULL, ...) {
    if (cens == F) {
        if (is.numeric(x) == FALSE) {
            stop("Database must be numeric!")
        }
        my_variable = x
        layout(mat = matrix(c(1, 2), 2, 1, byrow = TRUE), height = c(1, 8))
        par(mar = c(0, 3.1, 1.1, 2.1))
        boxplot(my_variable, horizontal = TRUE, xaxt = "n", col = colb, frame = F)
        par(mar = c(4, 3.1, 1.1, 2.1))
        hist(my_variable, breaks = breaks, col = colh, border = F, main = main, xlab = xlab)
        box()

        desc = list()
        desc$Mean = mean(my_variable)
        desc$Median = median(my_variable)
        desc$Mode = Mode(my_variable)
        desc$Variance = var(my_variable)
        desc$sd = sd(my_variable)
        desc$Minimum = min(my_variable)
        desc$Maximum = max(my_variable)
        desc$n = NROW(my_variable)

        return(desc)

    } else {
        if (is.data.frame(x) == FALSE) {
            stop("Database must be a data frame!")
        }
        if (is.numeric(x[, 1]) == FALSE && is.numeric(x[, 2]) == FALSE) {
            stop("Database must be numeric!")
        }
        if (is.numeric(b) == FALSE) {
            stop("Database must be numeric!")
        }


        fit <- survfit(Surv(x[, 1], x[, 2]) ~ b)

        plot(fit, conf.int = FALSE, xlab = "x", ylab = "S(x)", lwd = 2, ...)

        return(fit)
    }
}

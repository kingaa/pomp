reulermultinom <- function (n = 1, size, rate, dt) {
    tryCatch(
        .Call(R_Euler_Multinom,n,size,rate,dt),
        error = function (e) {
            stop("in ",sQuote("reulermultinom"),": ",
                 conditionMessage(e),call.=FALSE)
        }
    )
}

deulermultinom <- function (x, size, rate, dt, log = FALSE) {
    tryCatch(
        .Call(D_Euler_Multinom,as.matrix(x),size,rate,dt,log),
        error = function (e) {
            stop("in ",sQuote("deulermultinom"),": ",
                 conditionMessage(e),call.=FALSE)
        }
    )
}

rgammawn <- function (n = 1, sigma, dt) {
    tryCatch(
        .Call(R_GammaWN,n,sigma,dt),
        error = function (e) {
            stop("in ",sQuote("rgammwn"),": ",
                 conditionMessage(e),call.=FALSE)
        }
    )
}

library(nimble)

# custom c14 calibration density functions for use in Nimble c14 calibration

dcal <- nimbleFunction(
    run = function(x = double(0),
        c14_err = double(0, default = 20),
        caldate = double(0),
        calcurve = double(2),
        log = integer(0, default = 0)) {
            returnType(double(0))
            a <- which(calcurve[, 1] <= caldate)[1]
            b <- max(which(calcurve[, 1] >= caldate))
            if(a == b){
                logProb <- dnorm(x = calcurve[a, 2],
                                mean = x,
                                sd = sqrt(c14_err^2 + calcurve[a, 3]^2),
                                log = 1)
            }else{
                curve_mean <- ( calcurve[a, 2] * (calcurve[b, 1] - caldate) + 
                                calcurve[b, 2] * (caldate - calcurve[a, 1]) ) / 
                                    ( calcurve[b, 1] - calcurve[a, 1] )
                curve_error <- ( calcurve[a, 3] * (calcurve[b, 1] - caldate) + 
                                calcurve[b, 3] * (caldate - calcurve[a, 1]) ) / 
                                    ( calcurve[b, 1] - calcurve[a, 1] )
                logProb <- dnorm(x = curve_mean,
                                mean = x,
                                sd = sqrt(c14_err^2 + curve_error^2),
                                log = 1)
            }
            if(log) return(logProb)
            else return(exp(logProb)) 
    })

rcal <- nimbleFunction(
    run = function(n = integer(0),
        c14_err = double(0, default = 20),
        caldate = double(0),
        calcurve = double(2)) {
            returnType(double(0))
            if(n != 1) print("rcal only allows n = 1; using n = 1.")
            a <- which(calcurve[, 1] <= caldate)[1]
            b <- max(which(calcurve[, 1] >= caldate))
            if(a == b){
                c14 <- rnorm(n = 1,
                            mean = calcurve[a, 2],
                            sd = sqrt(c14_err^2 + calcurve[a, 3]^2))
            }else{
                curve_mean <- ( calcurve[a, 2] * (calcurve[b, 1] - caldate) + 
                                calcurve[b, 2] * (caldate - calcurve[a, 1]) ) / 
                                    ( calcurve[b, 1] - calcurve[a, 1] )
                curve_error <- ( calcurve[a, 3] * (calcurve[b, 1] - caldate) + 
                                calcurve[b, 3] * (caldate - calcurve[a, 1]) ) / 
                                    ( calcurve[b, 1] - calcurve[a, 1] )
                c14 <- rnorm(n = 1,
                            mean = curve_mean,
                            sd = sqrt(c14_err^2 + curve_error^2))
            }
            return(c14)
    })
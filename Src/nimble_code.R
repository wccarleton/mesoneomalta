library(nimble)

intcal20 <- read.csv("https://intcal.org/curves/intcal20.14c",
                    skip = 10)

names(intcal20) <- c("calbp",
                    "c14bp",
                    "c14_sigma",
                    "f14c",
                    "f14c_sigma")

# custom c14 calibration density functions for use in Nimble c14 calibration

dcal <- nimbleFunction(
    run = function(c14_mean = double(0), c14_err = double(0), caldate = double(0), calcurve = double(1),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- log(rate) - x*rate
        if(log) return(logProb)
        else return(exp(logProb)) 
    })

#to-do if necessary later
#rcal <- nimbleFunction(
#    run = function(n = integer(0), rate = double(0, default = 1)) {
#        returnType(double(0))
#        if(n != 1) print("rmyexp only allows n = 1; using n = 1.")
#        dev <- runif(1, 0, 1)
#        return(-log(1-dev) / rate)
#    })
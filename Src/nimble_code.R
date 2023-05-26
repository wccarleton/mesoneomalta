library(nimble)
library(IntCal)
library(ggplot2)

intcal20 <- read.csv("https://intcal.org/curves/intcal20.14c",
                    skip = 10)

names(intcal20) <- c("calbp",
                    "c14bp",
                    "c14_sigma",
                    "f14c",
                    "f14c_sigma")

# custom c14 calibration density functions for use in Nimble c14 calibration

dcal <- nimbleFunction(
    run = function(c14_mean = double(0),
        c14_err = double(0),
        caldate = double(0),
        calcurve = double(2),
        log = integer(0, default = 0)) {
            returnType(double(0))
            a <- which(calcurve[, 1] <= caldate)[1]
            b <- max(which(calcurve[, 1] >= caldate))
            if(a == b){
                logProb <- dnorm(x = calcurve[a, 2],
                                mean = c14_mean,
                                sd = sqrt(c14_err^2 + calcurve[a, 3]^2),
                                log = T)
            }else{
                curve_mean <- ( calcurve[a, 2] * (calcurve[b, 1] - caldate) + 
                                calcurve[b, 2] * (caldate - calcurve[a, 1]) ) / 
                                    ( calcurve[b, 1] - calcurve[a, 1] )
                curve_error <- ( calcurve[a, 3] * (calcurve[b, 1] - caldate) + 
                                calcurve[b, 3] * (caldate - calcurve[a, 1]) ) / 
                                    ( calcurve[b, 1] - calcurve[a, 1] )
                logProb <- dnorm(x = curve_mean,
                                mean = c14_mean,
                                sd = sqrt(c14_err^2 + curve_error^2),
                                log = T)
            }
            if(log) return(logProb)
            else return(exp(logProb)) 
    })

# testing
cal_density <- data.frame(caldate = seq(5600, 5900, 0.25), caldensity = NA)

for(j in cal_density$caldate){
    jj <- which(cal_density == j)
    cal_density[jj, "caldensity"] <- dcal(c14_mean = 5000, 
                                            c14_err = 20, 
                                            caldate = j, 
                                            calcurve = intcal20)
}

ggplot(cal_density) +
    geom_area(mapping = aes(x = caldate, y = caldensity),
            alpha = 0.75,
            fill = "steelblue") +
    scale_x_reverse() +
    xlim(c(5900, 5600)) +
    theme(text = element_text(size=20))

#to-do if necessary later
#rcal <- nimbleFunction(
#    run = function(n = integer(0),
#        c14_mean = double(0),
#        c14_err = double(0),
#        calcurve = double(2),
#        log = integer(0, default = 0)) {
#            returnType(double(0))
#            a <- which(calcurve[, 1] <= caldate)[1]
#            b <- max(which(calcurve[, 1] >= caldate))
#            if(a == b){
#                logProb <- dnorm(x = calcurve[a, 2],
#                                mean = c14_mean,
#                                sd = sqrt(c14_err^2 + calcurve[a, 3]^2),
#                                log = T)
#            }else{
#                curve_mean <- ( calcurve[a, 2] * (calcurve[b, 1] - caldate) + 
#                                calcurve[b, 2] * (caldate - calcurve[a, 1]) ) / 
#                                    ( calcurve[b, 1] - calcurve[a, 1] )
#                curve_error <- ( calcurve[a, 3] * (calcurve[b, 1] - caldate) + 
#                                calcurve[b, 3] * (caldate - calcurve[a, 1]) ) / 
#                                    ( calcurve[b, 1] - calcurve[a, 1] )
#                logProb <- dnorm(x = curve_mean,
#                                mean = c14_mean,
#                                sd = sqrt(c14_err^2 + curve_error^2),
#                                log = T)
#            }
#            if(log) return(logProb)
#            else return(exp(logProb)) 
#    })


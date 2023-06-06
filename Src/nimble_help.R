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
                                log = T)
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
                                log = T)
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

# I would like to register the distribution explicitly to specify the allowable 
# domain limits, which must be [0, 55000] because it's defined by the calibration 
# curve/data, but I can't get the model to work yet (see below)
#registerDistributions(list(
#    dcal = list(
#        BUGSdist = "dcal(c14_err, caldate, calcurve)",
#        pqAvail = FALSE, 
#        range = c(0, 55000),
#        types = c('calcurve=double(2)')
#        )
#    ))

# testing calibration outside nimble
true_date <- 5000
c14_samples <- rep(NA, 10000)

for(j in 1:10000){
    c14_samples[j] <- rcal(n = 1, 
                    c14_err = 20, 
                    caldate = true_date, 
                    calcurve = intcal20)
}

mean_c14 <- mean(c14_samples)

cal_density <- data.frame(caldate = seq(4900, 5300, 0.25), 
                        caldensity = NA)

for(j in cal_density$caldate){
    jj <- which(cal_density == j)
    cal_density[jj, "caldensity"] <- dcal(x = mean_c14, 
                                            c14_err = 20, 
                                            caldate = j, 
                                            calcurve = intcal20)
}

ggplot(cal_density) +
    geom_area(mapping = aes(x = caldate, y = caldensity),
            alpha = 0.75,
            fill = "steelblue") +
    scale_x_reverse() +
    xlim(c(5300, 4900)) +
    theme(text = element_text(size=20))

# now in a nimble model
# set up a simple test nimble model
gaussian_phase_model <- nimbleCode({
    mu ~ dnorm(mean = 2000, sd = 100)
    for(n in 1:N){
        y[n, 1] ~ dcal(c14_err = y[n, 2], 
                    caldate = mu,
                    calcurve = intcal20[1:9501, 1:3])
    }
})

# simulate some data
true_date <- 5000
N <- 10
c14_samples <- rep(NA, N)

for(j in 1:N){
    c14_samples[j] <- rcal(n = 1, 
                    c14_err = 20, 
                    caldate = true_date, 
                    calcurve = intcal20)
}

c14_data <- cbind(c14_samples, 20)

nimble_data <- list(y = c14_data)

# had to cast intcal20 here as.matrix b/c Nimble won't work with
# the dataframe directly
nimble_const <- list(N = N,
                    intcal20 = as.matrix(intcal20[1:9501, 1:3]))

nimble_inits <- list(mu = 2000)

mcmc_out <- nimbleMCMC(code = gaussian_phase_model,
                    constants = nimble_const,
                    data = nimble_data,
                    inits = nimble_inits,
                    niter = 10000)
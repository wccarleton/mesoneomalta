library(nimble)
library(ggplot2)
library(IntCal)
# get calibration data
intcal20 <- read.csv("https://intcal.org/curves/intcal20.14c",
                    skip = 10)

names(intcal20) <- c("calbp",
                    "c14bp",
                    "c14_sigma",
                    "f14c",
                    "f14c_sigma")

# nimble model

# include the custom nimble distributions
source("./Src/dcal.R")

# set up a test nimble model
cal_model <- nimbleCode({
    mu ~ dunif(0, 55000)#dunif(tt1, tt2)
    sigma ~ dunif(1e-7, 1000)
    for(n in 1:N){
        ybp[n] ~ dnorm(mu, sd = sigma)
        y[n, 1] ~ dcal(c14_err = y[n, 2], 
                    caldate = ybp[n],
                    calcurve = intcal20[1:9501, 1:3])
    }
})

# simulate some data
true_mean_date <- runif(n = 1, 
                        min = 1000, 
                        max = 10000)
N <- 100
true_dates <- rnorm(n = N, 
                    mean = true_mean_date, 
                    sd = 100)
c14_samples <- rep(NA, N)

for(j in 1:N){
    c14_samples[j] <- rcal(n = 1, 
                    c14_err = 20, 
                    caldate = true_dates[j], 
                    calcurve = intcal20)
}

c14_data <- cbind(c14_samples, 20)

#initial cal date guesses...
calibrate_simple <- function(c14,
                    calcurve){
    caldate <- 1:50000
    curve_c14_mean <- approx(x = calcurve[, 1],
                            y = calcurve[, 2],
                            xout = caldate)$y
    curve_err <- approx(x = calcurve[, 1],
                        y = calcurve[, 3],
                        xout = caldate)$y
    d <- dnorm(x = curve_c14_mean,
                mean = c14[1],
                sd = sqrt(c14[2]^2 + curve_err^2))
    return(d)
}

initial_cal_dates <- apply(c14_data, 1, calibrate_simple, calcurve = intcal20)

get_cal_range <- function(d, caldates, epsilon){
    len_caldates <- length(caldates)
    l_idx <- max(which(cumsum(d) <= epsilon))
    h_idx <- len_caldates - max(which(cumsum(rev(d)) <= epsilon))
    l <- caldates[l_idx]
    h <- caldates[h_idx]
    return(c(l, h))
}

cal_ranges <- t(apply(initial_cal_dates, 
                    2, 
                    get_cal_range, 
                    caldates = 1:50000, 
                    epsilon = 1e-15))

widest_date_range <- range(cal_ranges)

nimble_data <- list(y = c14_data,
                    t1 = cal_ranges[, 1],
                    t2 = cal_ranges[, 2],
                    tt1 = widest_date_range[1],
                    tt2 = widest_date_range[2])

mcmc_conf <- configureMCMC(compiled_model)
mcmc_conf$setMonitors(params_to_track)

cal_mcmc <- buildMCMC(mcmc_conf)

mcmc_compiled <- compileNimble(cal_mcmc, 
                            project = compiled_model)

nimble_const <- list(N = N,
                    intcal20 = as.matrix(intcal20[1:9501, 1:3]))

mu_start <- runif(n = 1, 
                widest_date_range[1], 
                widest_date_range[2])

nimble_inits <- list(ybp = rowMeans(cal_ranges),
                    mu = mu_start,#mean(widest_date_range),
                    sigma = 200)

params_to_track <- c("mu", "sigma", "ybp")

#custom mcmc invocation
nimble_model <- nimbleModel(code = cal_model, 
                            name = "cal_model",
                            data = nimble_data, 
                            inits = nimble_inits,
                            constants = nimble_const)

compiled_model <- compileNimble(nimble_model)

mcmc_conf <- configureMCMC(compiled_model)
mcmc_conf$setMonitors(params_to_track)

cal_mcmc <- buildMCMC(mcmc_conf)

mcmc_compiled <- compileNimble(cal_mcmc, 
                            project = compiled_model)

mcmc_out <- runMCMC(mcmc_compiled, 
                    niter = 10000, 
                    nburnin = 2000)

# check some results
plot(mcmc_out[, "ybp[1]"], type= "l", cex.axis = 2)
hist(mcmc_out[, "ybp[1]"], breaks = 100, cex.axis = 2)
calibrate(c14_data[1,1],20, cex.axis = 2)

hist(apply(mcmc_out,1,mean), cex.axis = 2)

plot(mcmc_out[, "mu"], type = "l", cex.axis = 2)
abline(h = true_mean_date, col = "#fcdcdc", lwd = 3)

# one-line mcmc invocation....
#mcmc_out <- nimbleMCMC(code = cal_model,
#                    constants = nimble_const,
#                    data = nimble_data,
#                    inits = nimble_inits,
#                    nburnin = 1000,
#                    niter = 10000,
#                    monitors = params_to_track)


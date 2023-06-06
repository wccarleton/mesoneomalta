library(nimble)
library(ggplot2)
library(IntCal)
library(tidyverse)
library(readxl)

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
    sigma ~ dunif(1e-7, 1000)
    for(j in 1:J){
        b[j] ~ dnorm(0, sd = 1000)
    }
    for(n in 1:N){
        mu[n] <- inprod(b[1:J], x[n, 1:J])
        ybp[n] ~ dnorm(mu[n], sd = sigma)
        y[n, 1] ~ dcal(c14_err = y[n, 2], 
                    caldate = ybp[n],
                    calcurve = intcal20[1:9501, 1:3])
    }
    # for interpolation
    if(!is.na(K)){
        for(k in 1:K){
            mu2[k] <- inprod(b[1:J], x2[k, 1:J])
            ybp2[k] ~ dnorm(mu2[k], sd = sigma)
        }
    }
})

# simulate some data
N <- 10
x <- matrix(c(rep(1, N), 0:(N-1)), ncol = 2)
b <- c(5000, 500)
mu <- x %*% b

true_dates <- rnorm(n = N, 
                    mean = mu, 
                    sd = 500)

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

initial_cal_dates <- apply(c14_data, 
                            1, 
                            calibrate_simple, 
                            calcurve = intcal20)

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

cal_range_means <- rowMeans(cal_ranges)

# predicting ages from undated depths
# probably easiest to maniuplate as a df
df <- as.data.frame(cbind(c14_data, x))
names(df) <- c("c14", "err", "intercept", "depth")

# recast as matrices
df <- as.matrix(df)

# determine positions to sample/predict at
x2 <- seq(0, (N - 1), 0.5)
K <- length(x2)
x2 <- matrix(c(rep(1, K), x2), ncol = 2)

ybp_inits <- approx(x = x[, 2], 
            y = cal_range_means, 
            xout = x2[, 2])$y

nimble_data <- list(y = df[, c(1, 2)],
                    x = df[, c(3, 4)],
                    x2 = x2)

nimble_const <- list(N = N,
                    intcal20 = as.matrix(intcal20[1:9501, 1:3]),
                    J = dim(x)[2],
                    K = K)

nimble_inits <- list(ybp = cal_range_means,
                    ybp2 = ybp_inits,
                    b = c(0, 0),
                    sigma = 200)

params_to_track <- c("b", "sigma", "ybp", "ybp2")

nimble_model <- nimbleModel(code = cal_model, 
                            name = "cal_model",
                            data = nimble_data, 
                            inits = nimble_inits,
                            constants = nimble_const)

compiled_model <- compileNimble(nimble_model)

#customize MCMC if needed
mcmc_conf <- configureMCMC(compiled_model)
mcmc_conf$setMonitors(params_to_track)

cal_mcmc <- buildMCMC(mcmc_conf)

mcmc_compiled <- compileNimble(cal_mcmc, 
                            project = compiled_model)

mcmc_out <- runMCMC(mcmc_compiled, 
                    niter = 20000, 
                    nburnin = 2000)

# check some results
plot(mcmc_out[, "ybp[1]"], type= "l", cex.axis = 2)
hist(mcmc_out[, "ybp[1]"], breaks = 100, cex.axis = 2)
calibrate(df[1,1], 20, cex.axis = 2)

plot(mcmc_out[, "b[2]"], type = "l", cex.axis = 2)
abline(h = b[2], col = "#fcdcdc", lwd = 3)

# extract predicted values
ybp2_idx <- grep("ybp2", colnames(mcmc_out))

predicted_ages <- apply(mcmc_out[, ybp2_idx], 
                    2, 
                    quantile, 
                    probs = c(0.05, 0.95)) 

predicted_ages <- as.data.frame(cbind(x2[, 2], t(predicted_ages)))

names(predicted_ages) <- c("depth", "lower", "upper")

ybp_idx <- grep("ybp\\[", colnames(mcmc_out))
samples_df <- as.data.frame(mcmc_out[, ybp_idx])


scaled_density <- function(x, ...){
    d <- density(x, ...)
    d$y <- d$y / max(d$y)
    return(d)
}

stack_densities <- function(x){
    d <- scaled_density(x)
    return(cbind(d$x, d$y))
}

ybp_densities <- apply(samples_df, 2, stack_densities)

n_density_samples <- dim(ybp_densities)[1] / 2

ybp_densities_long <- data.frame(c14_sample = rep(1:N, each = n_density_samples), 
            ybp = as.vector(ybp_densities[1:n_density_samples, ]),
            d = as.vector(ybp_densities[(n_density_samples + 1):(n_density_samples * 2), ]),
            depth = rep(x[, 2], each = n_density_samples))

samples_df_long <- pivot_longer(samples_df, 
                            cols = everything(), 
                            names_to = "param", 
                            values_to = "sample")

ggplot() +
    geom_ribbon(data = predicted_ages,
                mapping = aes(y = depth, 
                            xmin = lower, 
                            xmax = upper),
                alpha = 0.5,
                fill = "#333375") +
    geom_ribbon(data = ybp_densities_long,
                mapping = aes(x = ybp, 
                            ymin = depth,
                            ymax = -d + depth,
                            group = c14_sample),
                colour = "white",
                alpha = 0.5) +
    scale_y_reverse() +
    scale_x_reverse() +
    labs(title = "Simulated Age-Depth Model",
        x = "cal ybp") +
    theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 30),
            plot.title = element_text(size = 40))

# now the real data

data_path <- "./Data/Salina\ deep\ core.xlsx"
sheets <- excel_sheets(data_path)
sheets

# pull the variables sheet ("Measures")
core_data <- read_excel(data_path, 
                    sheet = sheets[1])

agedepth_df <- data.frame(c14 = core_data[, 3],
                        error = core_data[, 4],
                        depth = core_data[, 1],
                        outlier = core_data[, 7],
                        sample_id = core_data[, 2])

names(agedepth_df) <- c("c14", 
                        "err", 
                        "depth", 
                        "outlier", 
                        "lab_id")

N <- dim(agedepth_df)[1]

x <- matrix(c(rep(1, N), agedepth_df$depth), 
            ncol = 2)

# initial calibration
initial_cal_dates <- apply(agedepth_df[, c(1, 2)], 
                            1, 
                            calibrate_simple, 
                            calcurve = intcal20)

cal_ranges <- t(apply(initial_cal_dates, 
                    2, 
                    get_cal_range, 
                    caldates = 1:50000, 
                    epsilon = 1e-15))

widest_date_range <- range(cal_ranges)

cal_range_means <- rowMeans(cal_ranges)

# select depths to predict ages for

depth_range <- range(x[, 2])
delta <- diff(depth_range) / 50
x2 <- seq(depth_range[1], depth_range[2], by = delta)
K <- length(x2)
x2 <- matrix(c(rep(1, K), x2), ncol = 2)

ybp2_inits <- approx(x = x[, 2], 
            y = cal_range_means, 
            xout = x2[, 2])$y

# collate nimble data, inits, etc

nimble_data <- list(y = agedepth_df[, c(1, 2)],
                    x = x,
                    x2 = x2)

nimble_const <- list(N = N,
                    intcal20 = as.matrix(intcal20[1:9501, 1:3]),
                    J = dim(x)[2],
                    K = K)

nimble_inits <- list(ybp = cal_range_means,
                    ybp2 = ybp2_inits,
                    b = c(0, 0),
                    sigma = 200)

params_to_track <- c("b", "sigma", "ybp", "ybp2")

# build Nimble model

nimble_model <- nimbleModel(code = cal_model, 
                            name = "cal_model",
                            data = nimble_data, 
                            inits = nimble_inits,
                            constants = nimble_const)

compiled_model <- compileNimble(nimble_model)

#customize MCMC if needed
mcmc_conf <- configureMCMC(compiled_model)
mcmc_conf$setMonitors(params_to_track)

cal_mcmc <- buildMCMC(mcmc_conf)

mcmc_compiled <- compileNimble(cal_mcmc, 
                            project = compiled_model)

mcmc_out <- runMCMC(mcmc_compiled, 
                    niter = 20000, 
                    nburnin = 2000)

# check some results
plot(mcmc_out[, "ybp[1]"], type= "l", cex.axis = 2)
hist(mcmc_out[, "ybp[1]"], breaks = 100, cex.axis = 2)
calibrate(agedepth_df[1,1], agedepth_df[1, 2], cex.axis = 2)

plot(mcmc_out[, "b[2]"], type = "l", cex.axis = 2)

# plot age-depth relationship

# extract predicted values
ybp2_idx <- grep("ybp2", colnames(mcmc_out))
ybp_idx <- grep("ybp\\[", colnames(mcmc_out))

predicted_ages <- apply(mcmc_out[, ybp2_idx], 
                    2, 
                    quantile, 
                    probs = c(0.05, 0.95)) 

predicted_ages <- as.data.frame(cbind(x2[, 2], t(predicted_ages)))

names(predicted_ages) <- c("depth", "lower", "upper")

# extract ybp for actual dates

samples_df <- as.data.frame(mcmc_out[, ybp_idx])

ybp_densities <- apply(samples_df, 2, stack_densities)

n_density_samples <- dim(ybp_densities)[1] / 2

# because the depth is in cm, the densities need to be scaled
# again (not just to one) in order to see them on the plot. so
# note the multiplication is scaling for visualization only

ybp_densities_long <- data.frame(c14_sample = rep(1:N, each = n_density_samples), 
            ybp = as.vector(ybp_densities[1:n_density_samples, ]),
            d = as.vector(ybp_densities[(n_density_samples + 1):(n_density_samples * 2), ]) * 100,
            depth = rep(x[, 2], each = n_density_samples))

samples_df_long <- pivot_longer(samples_df, 
                            cols = everything(), 
                            names_to = "param", 
                            values_to = "sample")

ggplot() +
    geom_ribbon(data = predicted_ages,
                mapping = aes(y = depth, 
                            xmin = lower, 
                            xmax = upper),
                alpha = 0.5,
                fill = "#333375") +
    geom_ribbon(data = ybp_densities_long,
                mapping = aes(x = ybp, 
                            ymin = depth,
                            ymax = -d + depth,
                            group = c14_sample),
                colour = "white",
                alpha = 0.5) +
    scale_y_reverse() +
    scale_x_reverse() +
    labs(title = "Salina Deep Age-Depth Model",
        x = "cal ybp") +
    theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 30),
            plot.title = element_text(size = 40))
ggsave(filename = "Output/agedepth_custom.pdf")

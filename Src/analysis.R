### LIBRARIES ##################################################################
library(nimble)
library(ggplot2)
library(IntCal)
library(tidyverse)
library(readxl)
library(ggplot2)
library(oxcAAR)
library(gtools)
quickSetupOxcal()

### AGE DEPTH MODEL ############################################################

#### PRE-ANALYSIS PREP #########################################################

# get calibration data
intcal20 <- read.csv("https://intcal.org/curves/intcal20.14c",
                    skip = 10)

names(intcal20) <- c("calbp",
                    "c14bp",
                    "c14_sigma",
                    "f14c",
                    "f14c_sigma")

# nimble model

# include the custom nimble distribution for radiocarbon dates
source("./Src/dcal.R")

##### AGE DEPTH REGRESSION MODEL ###############################################
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

##### SIMULATED CHECKS #########################################################

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

#### SALINA DEEP AGE DEPTH ANALYSIS ############################################
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
ggsave(filename = "Output/agedepth_custom.pdf", device = "pdf")
ggsave(filename = "Output/agedepth_custom.png", device = "png")

#### ALTERNATIVE AGE DEPTH MODEL ###############################################

## Age-depth model
data_path <- "./Data/Salina\ deep\ core.xlsx"
sheets <- excel_sheets(data_path)
sheets

# pull the variables sheet ("Measures")
core_data <- read_excel(data_path, 
                    sheet = sheets[1])

agedepth_df <- data.frame(c14 = core_data[, 3],
                        error = core_data[, 4],
                        depth = core_data[, 1],
                        thickness = 1,
                        curve = "intcal20",
                        outlier = core_data[, 7],
                        sample_id = core_data[, 2])

names(agedepth_df) <- c("c14_ybp", 
                    "error", 
                    "depth_cm", 
                    "thickness", 
                    "curve", 
                    "outlier", 
                    "ids")

predict_depths <- seq(range(core_data[, 1])[1],
                        range(core_data[, 1])[2],
                        5)

## Run Bchronology
agedepth_model <- Bchronology(agedepth_df$c14_ybp,
                        agedepth_df$error,
                        agedepth_df$depth_cm,
                        agedepth_df$thickness,
                        agedepth_df$curve,
                        ids = agedepth_df$ids,
                        predictPositions = predict_depths,
                        iterations = 50000,
                        thin = 1)

## Summarize the model
agedepth_mean <- apply(agedepth_model$thetaPredict,
                        2,
                        mean)

agedepth_quant <- t(apply(agedepth_model$thetaPredict,
                            2,
                            quantile,
                            probs = c(0.05, 0.95)))

agedepth_summary <- data.frame(Depth = predict_depths,
                                    Mean = agedepth_mean,
                                    L05 = agedepth_quant[, 1],
                                    U95 = agedepth_quant[, 2])

## Write out results
write.table(agedepth_summary,
            file = "./Output/agedepth_Bchron.csv",
            sep = ",",
            row.names = F)

## Plot the results
ggplot(agedepth_summary) +
    geom_ribbon(mapping = aes(xmin = L05, xmax = U95, y = Depth),
                fill = "steelblue") +
    geom_line(mapping = aes(y = Depth, x = Mean)) +
    labs(x = "Year BP",
        y = "Depth (cm)") +
    scale_y_reverse() +
    scale_x_reverse() +
    theme_minimal()
ggsave(device = "pdf",
    file = "./Output/agedepth_Bchron.pdf")

## Bchron autoplot
plot(agedepth_model)
ggsave(device = "pdf",
    file = "./Output/agedepth_Bchron_auto.pdf")

# Bchron model 2---outliers out

ind <- which(agedepth_df$outlier == 0)

agedepth_model_2 <- Bchronology(agedepth_df$c14_ybp[ind],
                        agedepth_df$error[ind],
                        agedepth_df$depth_cm[ind],
                        agedepth_df$thickness[ind],
                        agedepth_df$curve[ind],
                        ids = agedepth_df$ids[ind],
                        predictPositions = predict_depths,
                        iterations = 100000,
                        thin = 10)

plot(agedepth_model_2)
ggsave(device = "pdf",
    file = "./Output/agedepth_Bchron_auto_2.pdf")

### OXCAL STUFF ################################################################

# This and the next anlaysis involve the use of OxCal. Rather than using it
# and then pulling data manually from OxCal into R, we decided to use an 
# R package with wrappers to make this process more easily repeatable

#### OxCal utility functions ###################################################

# extract oxcal list element names, paramater names, and types 
oxcal_params_table <- function(nested_list, bcad = TRUE) {
  # Find all elements that match the pattern 'ocd[' using grep
  ocd_indices <- grep("^ocd\\[", names(nested_list))

  # Initialize vectors to store results
  ocd_names <- rep(NA, length(ocd_indices))
  ocd_types <- rep(NA, length(ocd_indices))
  ocd_agreement <- rep(NA, length(ocd_indices))
  ocd_overall_agreement <- rep(NA, length(ocd_indices))
  ocd_model_agreement <- rep(NA, length(ocd_indices))
  ocd_matrices <- vector("list", length(ocd_indices))  # Initialize a list to store matrices
  ocd_matrices_lik <- vector("list", length(ocd_indices))  # Initialize a list to store matrices
  expected_age <- rep(NA, length(ocd_indices))
  hdi_95 <- rep(NA, length(ocd_indices))
  lo_95_modelled <- rep(NA, length(ocd_indices))
  up_95_modelled <- rep(NA, length(ocd_indices))
  lo_95_unmodelled <- rep(NA, length(ocd_indices))
  up_95_unmodelled <- rep(NA, length(ocd_indices))

  # Loop over indices explicitly to extract 'name', 'type', and matrices of ages and posterior probabilities
  for (i in seq_along(ocd_indices)) {
    index <- ocd_indices[i]
    if (!is.null(nested_list[[index]]$name)) {
      ocd_names[i] <- nested_list[[index]]$name
    }
    if (index == 1){
      ocd_types[i] <- "model"
    }else{
      if (!is.null(nested_list[[index]]$type)) {
        ocd_types[i] <- nested_list[[index]]$type
      }
    }
    
    # Extract values and create matrix
    if (!is.null(nested_list[[index]]$posterior$prob)) {
      posterior_prob <- nested_list[[index]]$posterior$prob
      plot_prob <- posterior_prob / max(posterior_prob)  # Normalize to max 1 for plotting
      start_age_bp <- nested_list[[index]]$posterior$start
      age_resolution <- nested_list[[index]]$posterior$resolution
      ages <- seq(from = start_age_bp, by = age_resolution, length.out = length(posterior_prob))
      ocd_matrices[[i]] <- cbind(ages, plot_prob)  # Combine ages and probabilities into a matrix for plotting
      
      # Normalize posterior_prob to integrate to 1 for calculations
      posterior_prob <- posterior_prob / sum(posterior_prob * age_resolution)
      
      # Calculate expected age (weighted mean)
      expected_age[i] <- sum(ages * posterior_prob) / sum(posterior_prob)
      
      # Calculate the 95% HDR for modelled ages
      age_prob_matrix <- cbind(ages, posterior_prob)
      age_prob_matrix <- age_prob_matrix[order(-age_prob_matrix[, 2]), ]  # Sort by posterior_prob descending
      cumsum_probs <- cumsum(age_prob_matrix[, 2] * age_resolution)
      hdr_indices <- which(cumsum_probs <= 0.95)
      lo_95_modelled[i] <- min(age_prob_matrix[hdr_indices, 1])
      up_95_modelled[i] <- max(age_prob_matrix[hdr_indices, 1])
      hdi_95[i] <- abs(up_95_modelled[i] - lo_95_modelled[i])
    }

    if (!is.null(nested_list[[index]]$likelihood$prob)) {
      likelihood_prob <- nested_list[[index]]$likelihood$prob
      lik_prob <- likelihood_prob / max(likelihood_prob)  # Normalize to max 1 for plotting
      lik_start_age_bp <- nested_list[[index]]$likelihood$start
      lik_age_resolution <- nested_list[[index]]$likelihood$resolution
      lik_ages <- seq(from = lik_start_age_bp, by = lik_age_resolution, length.out = length(lik_prob))
      ocd_matrices_lik[[i]] <- cbind(lik_ages, lik_prob)

      # Normalize unmodelled densities to integrate to 1 for calculations
      lik_prob <- lik_prob / sum(lik_prob * lik_age_resolution)

      # Calculate the 95% HDR for unmodelled ages
      age_prob_matrix <- cbind(lik_ages, lik_prob)
      age_prob_matrix <- age_prob_matrix[order(-age_prob_matrix[, 2]), ]  # Sort by posterior_prob descending
      cumsum_probs <- cumsum(age_prob_matrix[, 2] * lik_age_resolution)
      hdr_indices <- which(cumsum_probs <= 0.95)
      lo_95_unmodelled[i] <- min(age_prob_matrix[hdr_indices, 1])
      up_95_unmodelled[i] <- max(age_prob_matrix[hdr_indices, 1])
      #hdi_95_unmodelled[i] <- abs(upper_95_hdr[i] - lower_95_hdr[i])
    }

    # extract agreement indeces
    if(!is.null(nested_list[[index]]$posterior$overallAgreement)){
      ocd_overall_agreement[i] <- nested_list[[index]]$posterior$overallAgreement
    }
    if(!is.null(nested_list[[index]]$posterior$modelAgreement)){
      ocd_model_agreement[i] <- nested_list[[index]]$posterior$modelAgreement
    }
    if(!is.null(nested_list[[index]]$posterior$agreement)){
      ocd_agreement[i] <- nested_list[[index]]$posterior$agreement
    }
  }

  if(!bcad){
    l <- length(expected_age)
    expected_age[5:l] <- expected_age[5:l] - 1950
    lo_95_unmodelled[5:l] <- lo_95_unmodelled[5:l] - 1950
    up_95_unmodelled[5:l] <- up_95_unmodelled[5:l] - 1950
    lo_95_modelled[5:l] <- lo_95_modelled[5:l] - 1950
    up_95_modelled[5:l] <- up_95_modelled[5:l] - 1950
  }

  # Create the tibble with extracted values
  result_table <- tibble(
    ocd = names(nested_list)[ocd_indices],  # ocd[n]
    index = ocd_indices,                    # List indices for reference
    name = ocd_names,                       # Extracted 'name'
    type = ocd_types,                       # Extracted 'type'
    expected_age = round(expected_age, 2),
    hdi_95_modelled = hdi_95,
    lo_95_unmodelled = lo_95_unmodelled,
    up_95_unmodelled = up_95_unmodelled,
    lo_95_modelled = lo_95_modelled,
    up_95_modelled = up_95_modelled,
    agreement = ocd_agreement,
    overall_agreement = ocd_overall_agreement,
    model_agreement = ocd_model_agreement,
    likelihood = ocd_matrices_lik, # combined age and likelihood (unmodelled cal dates)
    postprob = ocd_matrices # Combined matrix of ages and posterior probabilities
  )
  return(result_table)
}

# example
# oxcal_table <- oxcal_params_table(result_data)
# access postprob like this:
# oxcal_table[8,]$postprob[[1]]

### REGIONAL MESO-NEO MODEL ####################################################
#### DATA WRANGLING ############################################################
# pull in the data from Excel sheets and CSVs as needed
# get sheet names
data_path <- "./Data/radiocarbon_database_regional.xlsx"
sheets <- excel_sheets(data_path)
sheets

# pull the variables sheet ("Measures")
central_med <- read_excel(data_path, 
                    sheet = sheets[1])

# create a new dataframe for input into OxCal. Will have to include an arbitrary
# unique key/name column for the dates because not all of the dates in the
# spreadsheet have lab numbers.

c14_dates <- data.frame(Name = central_med$ProjCode,#1:dim(central_med)[1],
        Date = central_med$date,
        Uncertainty = central_med$error,
        Culture = central_med$Culture,
        Region = central_med$Region,
        Marine = central_med$'Marine?',
        DeltaR = central_med$'ΔR',
        DeltaR_Err = central_med$'ΔR error')

# break up the data by regions and save as seperate CSVs for easier input into
# oxcal scripts

# meso first
c14_dates_meso <- subset(c14_dates, Culture == "Mesolithic")
for(j in unique(c14_dates$Region)){
    write.csv(c14_dates_meso[grep(j, c14_dates_meso$Region), 1:3],
        file = paste("./Data/meso_",j,".csv",sep=""),
        row.names = FALSE)
}

# next neo
c14_dates_neo <- subset(c14_dates, Culture == "Neolithic")
for(j in unique(c14_dates$Region)){
    write.csv(c14_dates_neo[grep(j, c14_dates_neo$Region), 1:3],
        file = paste("./Data/neo_",j,".csv",sep=""),
        row.names = FALSE)
}

#### OXCAL MODELLING ###########################################################
# OxCal scripts were written separately (see /Src) and so here we just call
# them with oxcAAR functions and process the results for plotting

# grab all files in /Src with meso|neo in the filename
filenames <- list.files("Src/")
subset_filenames <- grep("meso|neo", filenames, value = TRUE)

oxcal_tables <- list()

for (j in seq_along(subset_filenames)) {
  # for each file, pass the OxCal code to oxcal, grab and process the result
  filename <- subset_filenames[j]

  modelname <- tools::file_path_sans_ext(filename)

  file_lines <- readLines(paste("Src/", filename, sep = ""), warn = FALSE)

  # Collapse all lines into a single string
  oxcal_code <- paste(file_lines, collapse = "\n")

  oxcal_results_file <- executeOxcalScript(oxcal_code)

  oxcal_results_txt <- readOxcalOutput(oxcal_results_file)

  oxcal_results_list <- parseFullOxcalOutput(oxcal_results_txt)

  oxcal_table <- oxcal_params_table(oxcal_results_list)

  oxcal_table$modelname <- modelname

  oxcal_tables[[j]] <- oxcal_table
}

oxcal_table_final <- bind_rows(oxcal_tables)

report_variables <- c("name",
                        "type",
                        "expected_age",
                        "hdi_95",
                        "lower_95_hdr",
                        "upper_95_hdr",
                        "agreement",
                        "overall_agreement",
                        "model_agreement")

report_variables_idx <- grep(paste(report_variables, collapse = "|"), 
                                names(oxcal_table_final))

# save main results
write.table(oxcal_table_final[, report_variables_idx], 
          file = "Output/meso_neo_model_results.csv",
          row.names = FALSE)

# we're only interested in plotting the boundaries for the end of the mesolithic
# and start of the neolithic
meso_end_idx <- grep("meso_end", oxcal_table_final$name)
neo_start_idx <- grep("neo_start", oxcal_table_final$name)

# Create a long format table for plotting
plot_data <- data.frame(name = character(), 
                        period = character(),
                        region = character(),
                        ages = numeric(), 
                        posterior_prob = numeric())

# Populate the long format table by extracting data from each row of oxcal_table
for (i in meso_end_idx) {
  current_name <- oxcal_table_final[i, 'name']
  current_matrix <- oxcal_table_final$postprob[[i]]
  
  # Convert matrix to data frame and add other variables
  current_df <- data.frame(ages = current_matrix[, 1], posterior_prob = current_matrix[, 2])
  current_df$period <- "Mesolithic End"
  current_df$region <- gsub("_meso_end", "", current_name)
    
  # Append to the plot_data data frame
  plot_data <- rbind(plot_data, current_df)
}

# Populate the long format table by extracting data from each row of oxcal_table
for (i in neo_start_idx) {
  current_name <- oxcal_table_final[i, 'name']
  current_matrix <- oxcal_table_final$postprob[[i]]
  
  # Convert matrix to data frame and add other variables
  current_df <- data.frame(ages = current_matrix[, 1], posterior_prob = current_matrix[, 2])
  current_df$period <- "Neolithic Start"
  current_df$region <- gsub("_neo_start", "", current_name)
    
  # Append to the plot_data data frame
  plot_data <- rbind(plot_data, current_df)
}

# reorder regions for nice plotting
regions_order <- c("N_Italy",
                  "C_Italy",
                  "S_Italy",
                  "Sicily",
                  "Corsica",
                  "Sardinia",
                  "Malta")
plot_data$region <- factor(plot_data$region, 
                          levels = regions_order)

plt <- ggplot(plot_data, mapping = aes(x = ages - 1950, y = posterior_prob, fill = period)) +
        geom_area(position = "identity") +
        facet_grid(region ~ .) +
        labs(x = "year cal. BP") +
        theme_minimal()

plt

ggsave(filename = "./Output/modelled_meso_neo_boundaries.pdf",
        device = "pdf")

ggsave(filename = "./Output/modelled_meso_neo_boundaries.png",
        device = "png")

### LATNIJA PHASE MODEL ########################################################
#### OXCAL MODELLING ###########################################################
# read in the OxCal model
# Read the file contents as lines
file_lines <- readLines("Src/LAT_main.oxcal", warn = FALSE)

# Collapse all lines into a single string
oxcal_code <- paste(file_lines, collapse = "\n")

# look at the code
cat(oxcal_code)

# pass it to oxcal
oxcal_results_file <- executeOxcalScript(oxcal_code)

# bring in OxCal results for plotting
oxcal_results_txt <- readOxcalOutput(oxcal_results_file)

# save results just in case
writeLines(oxcal_results_txt, con = "Output/oxcal_latnija_model_results.txt")

# this can of course be read back in
oxcal_results_txt <- readLines("Output/oxcal_latnija_model_results.txt")

# parse the oxcal text into an R list
result_data <- parseFullOxcalOutput(oxcal_results_txt)

#### EXTRACT RESULTS AND PLOT ##################################################
# table relating oxcal data list elements to model parameter names and 
# indeces
oxcal_table_full <- oxcal_params_table(result_data, bcad = FALSE)

report_variables <- c("name",
                        "type",
                        "expected_age",
                        "lo_95_unmodelled",
                        "up_95_unmodelled",
                        "lo_95_modelled",
                        "up_95_modelled",
                        "hdi_95_modelled",
                        "agreement",
                        "overall_agreement",
                        "model_agreement")

report_variables_idx <- grep(paste(report_variables, collapse = "|"), 
                        names(oxcal_table_full))

write.table(oxcal_table_full[, report_variables_idx], 
            file = "Output/latnija_oxcal_results_table.csv", 
            row.names = FALSE)

oxcal_table <- oxcal_table_full[-which(oxcal_table_full$type == "model" | is.na(oxcal_table_full$type)),]

# Create a long format table for plotting
plot_data <- data.frame(name = character(), 
                        ages = numeric(), 
                        density = numeric(),
                        density_type = character(),
                        offset = numeric())

# Populate the long format table by extracting data from each row of oxcal_table
for (i in 1:nrow(oxcal_table)) {
  current_name <- oxcal_table$name[i]
  postprob_matrix <- oxcal_table$postprob[[i]]
  lik_matrix <- oxcal_table$likelihood[[i]]
  
  # Convert matrix to data frame and add the 'name' column
  if (!is.null(postprob_matrix)) {
    # first for the posterior
    postprob_df <- data.frame(ages = postprob_matrix[, 1], 
                              density = postprob_matrix[, 2])
    postprob_df$name <- current_name
    postprob_df$density_type <- "posterior"
    postprob_df$offset <- i  # Add an offset for vertical separation in the plot
    # again for the likelihood (unmodelled cal date)
    if(!is.null(lik_matrix)) {
      lik_df <- data.frame(ages = lik_matrix[, 1], 
                            density = lik_matrix[, 2])
      lik_df$name <- current_name
      lik_df$density_type <- "likelihood"
      lik_df$offset <- i  # Add an offset for vertical separation in the plot
    }
    # Append to the plot_data data frame
    plot_data <- rbind(plot_data, postprob_df, lik_df)
  }
}

# plot_data contains posteriors for model parameters we don't want to plot
plot_data <- plot_data[-which(plot_data$name == "General" | plot_data$name == ""),]

# Read the updated CSV with the order and type columns
oxcal_sample_map <- read.csv("Data/sample_model_map_ordered.csv")

# Merge the updated sample map with the plot data
merged_plot_data <- merge(plot_data, oxcal_sample_map, by = "name", all.x = TRUE)

# Calculate expected years for each density
expected_years <- merged_plot_data %>%
  group_by(name) %>%
  summarise(expected_year = sum(ages * density) / sum(density))

# Merge expected years into the sample map
oxcal_sample_map <- merge(oxcal_sample_map, expected_years, by = "name", all.x = TRUE) %>%
    arrange(order)

# For each sequence and phase, sort the events by expected year to improve the plot
phasing <- c("6a", "6b", "5a", "5b", "4", "3a", "3b")

for(j in 1:length(phasing)){
  event_rows <- as.numeric(row.names(subset(oxcal_sample_map, phase_num == phasing[j] &
                                            type == "event")))
  oxcal_sample_map[event_rows,] <- oxcal_sample_map[event_rows,] %>% 
                                    arrange(expected_year)
}

oxcal_sample_map$order <- 1: dim(oxcal_sample_map)[1]

# Merge the updated sample map with the plot data
merged_plot_data <- merge(plot_data, oxcal_sample_map, by = "name", all.x = TRUE)

# Set a threshold for density
threshold <- 0.2  # Adjust this value as needed for visibility

# Determine the position for labels by isolating ages where densities exceed the threshold
label_positions <- subset(merged_plot_data, density_type == "posterior") %>%
  group_by(name) %>%
  filter(density > threshold) %>%  
  filter(ages == max(ages)) %>%           
  mutate(label_x = ages - 1950, label_y = order)

# Plot using ggplot2 with vertical offsets
ggplot(merged_plot_data, aes(x = ages - 1950, group = name, fill = density_type)) +
  # Plot posterior as a filled ribbon
  geom_ribbon(data = subset(merged_plot_data, density_type == "posterior"),
              aes(ymin = order, ymax = density + order, fill = type), alpha = 0.5) +
  
  # Plot likelihood as a gray line
  geom_line(data = subset(merged_plot_data, density_type == "likelihood"),
            aes(y = density + order), color = "black", size = 0.1) +

  # Add labels
  geom_text(data = label_positions, aes(x = label_x + 100, y = label_y + 0.2, label = name), 
            hjust = 0, nudge_x = 10, size = 2) +
  
  # Additional plot settings
  xlim(-12000, -5500) +  
  theme_minimal() + 
  labs(x = "Ages (BP)", y = "Posterior Probability (Normalized)", 
       title = "OxCal Phase Model Posteriors", 
       fill = "Density Type") +
  
  # Customize the appearance of facets and legend
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = ""
  ) +
  
  # Specify fill colors for posterior based on type
  scale_fill_manual(values = c("date" = "steelblue", "boundary" = "lightgreen"), 
                    name = "Posterior Type") +
  
  # Specify color for likelihood lines based on density_type
  scale_color_manual(values = c("likelihood" = "black"), 
                     name = "Density Type", 
                     labels = c("Likelihood")) +
  
  # Guide adjustment to combine legends
  guides(fill = guide_legend(override.aes = list(color = NA)),
         color = guide_legend(override.aes = list(fill = NA, size = 1)))

ggsave(filename = "Output/posterior_probabilities_plot.pdf",
    width = 90,
    height = 170,
    units = "mm",
    dpi = 300)

ggsave(filename = "Output/posterior_probabilities_plot.png",
    width = 90,
    height = 170,
    units = "mm",
    dpi = 300)

# Plot using ggplot2 with vertical offsets, focus on only the key boundaries
boundaries_only <- subset(merged_plot_data, type == "boundary")

# Filter out rows where name is "=IV_st"
boundaries_only <- subset(boundaries_only, name != "=IV_st")

boundaries_only <- subset(boundaries_only, ages >= -(9000 - 1950) & ages <= -(7000 - 1950))

# Sort the dataframe by 'order' in ascending order
boundaries_only <- boundaries_only[order(boundaries_only$order), ]

# Get the unique values of 'order' and create a mapping to consecutive integers
unique_order <- sort(unique(boundaries_only$order))
order_mapping <- setNames(seq_along(unique_order), unique_order)

# Replace 'order' in the dataframe with the relative order based on the mapping
boundaries_only$order <- order_mapping[as.character(boundaries_only$order)]

# Determine the position for labels by isolating ages where densities exceed the threshold
label_positions_bounds <- boundaries_only %>%
  group_by(name) %>%
  filter(density == max(density)) %>%  
  filter(ages == mean(ages)) %>%           
  mutate(label_x = ages - 1950, label_y = 1)

# one label needs shifting to avoid overlap and since it's only one
# a pragmatic solution is to change its y coordinates manually
label_positions_bounds[which(label_positions_bounds$name == "IIIA_st"), "label_y"] <- 1.05

ggplot(boundaries_only, aes(x = ages - 1950, group = name)) +
  # Plot posterior as a filled ribbon
  geom_ribbon(aes(ymin = 0, ymax = density, fill = name), alpha = 0.5) +

  # Add labels
  geom_text(data = label_positions_bounds, aes(x = label_x, y = label_y, label = name), 
          size = 2) +
  
  # Additional plot settings
  theme_minimal() + 
  labs(x = "Ages (BP)", y = "Posterior Probability (Normalized)", 
       title = "Key Boundaries") +
  
  # Customize the appearance of facets and legend
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = ""
  )

ggsave(filename = "Output/posterior_probabilities_bounds_plot.pdf",
    width = 90,
    height = 40,
    units = "mm",
    scale = 1.5,
    dpi = 300)

ggsave(filename = "Output/posterior_probabilities_bounds_plot.png",
    width = 90,
    height = 40,
    units = "mm",
    scale = 1.5,
    dpi = 300)
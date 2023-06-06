library(nimble)
library(IntCal)
library(ggplot2)

source(file = "dcal.R")

intcal20 <- read.csv("https://intcal.org/curves/intcal20.14c",
                    skip = 10)

names(intcal20) <- c("calbp",
                    "c14bp",
                    "c14_sigma",
                    "f14c",
                    "f14c_sigma")

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
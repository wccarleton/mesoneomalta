# load libraries
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(readxl)
library(Bchron)

# data wrangling

# pull in the data from Excel sheets and CSVs as needed

# get sheet names
data_path <- "./Data/FINALRadiocarbonDatabase.xlsx"
sheets <- excel_sheets(data_path)
sheets

# pull the variables sheet ("Measures")
central_med <- read_excel(data_path, 
                    sheet = sheets[1])

# create a new dataframe for input into OxCal. Will have to include an arbitrary
# unique key/name column for the dates because not all of the dates in the
# spreadsheet have lab numbers.

c14_dates <- data.frame(Name = 1:dim(central_med)[1],
        Date = central_med$date,
        Uncertainty = central_med$error,
        Culture = central_med$Culture,
        Region = central_med$Region,
        Marine = central_med$'Marine?',
        DeltaR = central_med$'ΔR',
        DeltaR_Err = central_med$'ΔR error')

# break up the data by regions and save as seperate CSVs for easier input into
# oxcal

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

# extracting coords for plotting samples in QGIS
# coords in the spreadsheet appear to have the wrong headers, so corrected
# here
c14_coords <- data.frame(idx = c14_dates$Name,
                        Lat = central_med$Longitude,
                        Long = central_med$Latitude)

write.csv(c14_coords,
        file="./Data/c14dates_coords.csv", 
        row.names = F)

# modelling occurs in OxCal

# bring in OxCal results for plotting

oxcal_results_path <- "OxCal/"

files <- list.files(oxcal_results_path)

posteriors <- grep("prior", files)

posts <- list()

for(j in 1:length(posteriors)){
        f <- files[posteriors[j]]
        p <- paste(oxcal_results_path, f, sep = "")
        df <- as.data.frame(read.table(file = p, header = F))
        names(df) <- c("date", "density")
        r <- regexpr("meso|neo", files[posteriors[j]])
        period <- substr(files[posteriors[j]], r[1], 
                        r[1] + attr(r, "match.length") - 1)
        df$period <- period
        idx <- regexpr("meso|neo", 
                files[posteriors[j]])[1] - 2
        rgn <- substr(files[posteriors[j]], 1, idx)
        df$region <- rgn
        posts[[j]] <- df
}

posts <- do.call(rbind, posts)

regions_ordered <- c("N_Italy",
                "C_Italy",
                "S_Italy",
                "Sicily",
                "Corsica",
                "Sardinia",
                "Malta")

posts$region <- factor(posts$region, levels = regions_ordered)
posts$date <- abs(posts$date)

plt <- ggplot(posts, mapping = aes(x = 1950 + date, y = density, fill = period)) +
        geom_area(position = "identity") +
        facet_grid(region ~ .) +
        labs(x = "years BP") +
        scale_x_reverse() +
        theme_minimal()

plt

ggsave(filename = "./Output/modelled_meso_neo_boundaries.pdf",
        device = "pdf")

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

names(agedepth_df) <- c("c14_ybp", "error", "depth_cm", "thickness", "curve", "outlier", "ids")

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

# try with Bacon

bacon_agedepth_df <- data.frame(lab_ID = core_data[, 2],
                                Age = core_data[, 3],
                                Error = core_data[, 4],
                                Depth = core_data[, 1],
                                cc = 1)

write.csv(bacon_agedepth_df, 
        file = "Bacon/SalinaDeep.csv", 
        row.names = F)

Bacon(core = "SalinaDeep", 
        coredir = "Bacon/",
        thick = 50,
        ssize = 50000,
        burnin = 30000,
        d.by = 50,
        acc.mean = 2,
        acc.shape = 1.5)
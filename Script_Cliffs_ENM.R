rm(list=ls())
set.seed(10)
setwd("~")
# Libraries ####
pacman::p_load(predicts, sdm, tidyverse, mapview, devtools, countrycode, rgbif, CoordinateCleaner,
               ggsci, usethis, usdm, spatialEco, colorRamps, terra, randomForest, leaflet, spdep, MLeval,
               sp, ggpubr, ggspatial, plotly, ggnewscale, readr, maps, viridis, kernelshap, shapviz,
               rmapshaper, ggthemes, tidyterra, ggsignif, pdp, pROC, caret, ranger, tuneRanger)

# Occurrences and Raster loading ####
Species_data <- vect("~/Cliff_Occurrences_ESTR25831.shp")
Predictive_rasters <- rast("~/Predictive_Europe.tif")
Europe_Limits <- vect("~/Coastline_Europe_25381.shp")
Balearic_Limits <- vect("~/Coastline_BalearicIslands_25831.shp")

# Pseudo-Absence calculation ####
absence_multiplier = 2
buffers <- buffer(Species_data, 1000)

Background_area <- erase(Limites_Europe, aggregate(buffers))

Absence_points <- terra::spatSample(Background_area, 
                                    size = absence_multiplier*nrow(Species_data),
                                    method = "random")

Absence_points$Presence <- 0
Species_data$Presence <- 1

Occ_cliffs <- terra::union(Absence_points, Species_data)

# Extract data from rasters for each points ####
occs <- Occ_cliffs[,"Presence"]

extracted_predictors_raw <- as.data.frame(terra::extract(Predictive_rasters, occs, xy=T, bind=T))

extracted_predictors_raw2 <- extracted_predictors_raw %>% 
  select(-"Presence",-"x",-"y", -"Geo", -"Land_Use")
extracted_predictors <- extracted_predictors_raw %>% 
  select("Presence","x","y","Geo","Land_Use")

# Fill missing data
extracted_NA_kNN <- caret::preProcess(extracted_predictors_raw2,
                                      method="bagImpute")
missing_points <- predict(extracted_NA_kNN, extracted_predictors_raw2)

cliff_data <- cbind(extracted_predictors, missing_points)

# Predictive variables selection ####
Not_Numeric_spg <- cliff_data %>%
  select_if(is.factor)

loop_val <- 100
while (loop_val >= 10) {
  Numeric_spg <- cliff_data %>%
    select(-"Presence",-"x",-"y") %>%
    select_if(.,is.numeric)
  
  VIF_numeric <- usdm::vif(Numeric_spg[,6:length(Numeric_spg)]) %>%
    arrange(VIF)
  VIF_numeric 
  
  Selected_variables <- VIF_numeric[-c(length(VIF_numeric$VIF)),]
  
  cliff_selected_data <- cliff_data %>%
                  select("Presence","x","y",names(Numeric_spg[,1:5]),
                         Selected_variables$Variables,names(Not_Numeric_spg))
  
  loop_val <- VIF_numeric[c(length(VIF_numeric$VIF)),2]
  print(paste("max VIF at", loop_val))
}

# Random Forest model ####
cliff_selected_data$Presence <- as.factor(cliff_selected_data$Presence)
levels(cliff_selected_data$Presence) <- c("Ausence","Presence")

rf.task <-  makeClassifTask(data = cliff_selected_data, target = "Presence")
res <- tuneRanger(rf.task, measure = list(multiclass.brier), num.trees = 500,
                  num.threads = 20, iters = 100, save.file.path = NULL)

cliff_RF_model <- caret::train(
  cliff_selected_data[,-1],     # This way the factor rasters are keep together
  as.factor(cliff_selected_data$Presence),
  method = "ranger",
  trControl = trainControl(method = "repeatedcv",   # k-fold cross-validation
                           number = 40,             # 40 folds
                           repeats = 40,            # repeated 40 times
                           allowParallel = T,
                           classProbs=T,
                           returnData = T,
                           savePredictions = "final"),
  tuneGrid = expand.grid(mtry = res$recommended.pars[,1],
                         min.node.size = res$recommended.pars[,2],
                         splitrule = "gini"),
  num.trees = 500,
  num.threads = 20,
  importance = 'impurity'
)
cliff_RF_model
cliff_RF_model$results

conf_RF <- confusionMatrix(cliff_RF_model)$table
ERROR_positive=(1-round(conf_RF[3]/conf_RF[1],4))*100
ERROR_negative=(1-round(conf_RF[2]/conf_RF[4],4))*100

resample_stats <- thresholder(ranger_model,
                              threshold = seq(0.2, 0.6, by = 0.01),
                              final = TRUE)
threshold_value <- as.numeric(resample_stats$prob_threshold[min(which(resample_stats$Precision > 0.9))])

# Current prediction for Balearic Islands ####
Predictive_rasters_Balearic <- crop(Predictive_rasters, Balearic_Limits)

ENM_Balearic_current <- terra::predict(Predictive_rasters_Balearic, cliff_RF_model, 
                                      type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# Future prediction for Balearic Islands ####
# CMCC 2021-2040 ####
Files_bioc_2040=list.files("~/Bioclim/Future predictions/2040/CMCC", ".tif", full.names = TRUE)

# SPP 126
CMCC_2040_ssp126 <- rast(Files_bioc_2040[1])
CMCC_pred_2040_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2040_ssp126)
ENM_CMCC_2040_ssp126 <- terra::predict(CMCC_pred_2040_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
CMCC_2040_ssp245 <- rast(Files_bioc_2040[2])
CMCC_pred_2040_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2040_ssp245)
ENM_CMCC_2040_ssp245 <- terra::predict(CMCC_pred_2040_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
CMCC_2040_ssp370 <- rast(Files_bioc_2040[3])
CMCC_pred_2040_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2040_ssp370)
ENM_CMCC_2040_ssp370 <- terra::predict(CMCC_pred_2040_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
CMCC_2040_ssp585 <- rast(Files_bioc_2040[4])
CMCC_pred_2040_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2040_ssp585)
ENM_CMCC_2040_ssp585 <- terra::predict(CMCC_pred_2040_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# CMCC 2041-2060 ####
Files_bioc_2060=list.files("~/Bioclim/Future predictions/2060/CMCC", ".tif", full.names = TRUE)

# SPP 126
CMCC_2060_ssp126 <- rast(Files_bioc_2060[1])
CMCC_pred_2060_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2060_ssp126)
ENM_CMCC_2060_ssp126 <- terra::predict(CMCC_pred_2060_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
CMCC_2060_ssp245 <- rast(Files_bioc_2060[2])
CMCC_pred_2060_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2060_ssp245)
ENM_CMCC_2060_ssp245 <- terra::predict(CMCC_pred_2060_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
CMCC_2060_ssp370 <- rast(Files_bioc_2060[3])
CMCC_pred_2060_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2060_ssp370)
ENM_CMCC_2060_ssp370 <- terra::predict(CMCC_pred_2060_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
CMCC_2060_ssp585 <- rast(Files_bioc_2060[4])
CMCC_pred_2060_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2060_ssp585)
ENM_CMCC_2060_ssp585 <- terra::predict(CMCC_pred_2060_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# CMCC 2061-2080 ####
Files_bioc_2080=list.files("~/Bioclim/Future predictions/2080/CMCC", ".tif", full.names = TRUE)

# SPP 126
CMCC_2080_ssp126 <- rast(Files_bioc_2080[1])
CMCC_pred_2080_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2080_ssp126)
ENM_CMCC_2080_ssp126 <- terra::predict(CMCC_pred_2080_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
CMCC_2080_ssp245 <- rast(Files_bioc_2080[2])
CMCC_pred_2080_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2080_ssp245)
ENM_CMCC_2080_ssp245 <- terra::predict(CMCC_pred_2080_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
CMCC_2080_ssp370 <- rast(Files_bioc_2080[3])
CMCC_pred_2080_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2080_ssp370)
ENM_CMCC_2080_ssp370 <- terra::predict(CMCC_pred_2080_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
CMCC_2080_ssp585 <- rast(Files_bioc_2080[4])
CMCC_pred_2080_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2080_ssp585)
ENM_CMCC_2080_ssp585 <- terra::predict(CMCC_pred_2080_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# CMCC 2081-2100 ####
Files_bioc_2100=list.files("~/Bioclim/Future predictions/2100/CMCC", ".tif", full.names = TRUE)

# SPP 126
CMCC_2100_ssp126 <- rast(Files_bioc_2100[1])
CMCC_pred_2100_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2100_ssp126)
ENM_CMCC_2100_ssp126 <- terra::predict(CMCC_pred_2100_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
CMCC_2100_ssp245 <- rast(Files_bioc_2100[2])
CMCC_pred_2100_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2100_ssp245)
ENM_CMCC_2100_ssp245 <- terra::predict(CMCC_pred_2100_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
CMCC_2100_ssp370 <- rast(Files_bioc_2100[3])
CMCC_pred_2100_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2100_ssp370)
ENM_CMCC_2100_ssp370 <- terra::predict(CMCC_pred_2100_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
CMCC_2100_ssp585 <- rast(Files_bioc_2100[4])
CMCC_pred_2100_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           CMCC_2100_ssp585)
ENM_CMCC_2100_ssp585 <- terra::predict(CMCC_pred_2100_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# MPI 2021-2040 ####
Files_bioc_2040=list.files("~/Bioclim/Future predictions/2040/MPI", ".tif", full.names = TRUE)

# SPP 126
MPI_2040_ssp126 <- rast(Files_bioc_2040[1])
MPI_pred_2040_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2040_ssp126)
ENM_MPI_2040_ssp126 <- terra::predict(MPI_pred_2040_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
MPI_2040_ssp245 <- rast(Files_bioc_2040[2])
MPI_pred_2040_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2040_ssp245)
ENM_MPI_2040_ssp245 <- terra::predict(MPI_pred_2040_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
MPI_2040_ssp370 <- rast(Files_bioc_2040[3])
MPI_pred_2040_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2040_ssp370)
ENM_MPI_2040_ssp370 <- terra::predict(MPI_pred_2040_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
MPI_2040_ssp585 <- rast(Files_bioc_2040[4])
MPI_pred_2040_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2040_ssp585)
ENM_MPI_2040_ssp585 <- terra::predict(MPI_pred_2040_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# MPI 2041-2060 ####
Files_bioc_2060=list.files("~/Bioclim/Future predictions/2060/MPI", ".tif", full.names = TRUE)

# SPP 126
MPI_2060_ssp126 <- rast(Files_bioc_2060[1])
MPI_pred_2060_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2060_ssp126)
ENM_MPI_2060_ssp126 <- terra::predict(MPI_pred_2060_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
MPI_2060_ssp245 <- rast(Files_bioc_2060[2])
MPI_pred_2060_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2060_ssp245)
ENM_MPI_2060_ssp245 <- terra::predict(MPI_pred_2060_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
MPI_2060_ssp370 <- rast(Files_bioc_2060[3])
MPI_pred_2060_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2060_ssp370)
ENM_MPI_2060_ssp370 <- terra::predict(MPI_pred_2060_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
MPI_2060_ssp585 <- rast(Files_bioc_2060[4])
MPI_pred_2060_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2060_ssp585)
ENM_MPI_2060_ssp585 <- terra::predict(MPI_pred_2060_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# MPI 2061-2080 ####
Files_bioc_2080=list.files("~/Bioclim/Future predictions/2080/MPI", ".tif", full.names = TRUE)

# SPP 126
MPI_2080_ssp126 <- rast(Files_bioc_2080[1])
MPI_pred_2080_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2080_ssp126)
ENM_MPI_2080_ssp126 <- terra::predict(MPI_pred_2080_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
MPI_2080_ssp245 <- rast(Files_bioc_2080[2])
MPI_pred_2080_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2080_ssp245)
ENM_MPI_2080_ssp245 <- terra::predict(MPI_pred_2080_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
MPI_2080_ssp370 <- rast(Files_bioc_2080[3])
MPI_pred_2080_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2080_ssp370)
ENM_MPI_2080_ssp370 <- terra::predict(MPI_pred_2080_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
MPI_2080_ssp585 <- rast(Files_bioc_2080[4])
MPI_pred_2080_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2080_ssp585)
ENM_MPI_2080_ssp585 <- terra::predict(MPI_pred_2080_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
# MPI 2081-2100 ####
Files_bioc_2100=list.files("~/Bioclim/Future predictions/2100/MPI", ".tif", full.names = TRUE)

# SPP 126
MPI_2100_ssp126 <- rast(Files_bioc_2100[1])
MPI_pred_2100_ssp126 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2100_ssp126)
ENM_MPI_2100_ssp126 <- terra::predict(MPI_pred_2100_ssp126, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 245
MPI_2100_ssp245 <- rast(Files_bioc_2100[2])
MPI_pred_2100_ssp245 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2100_ssp245)
ENM_MPI_2100_ssp245 <- terra::predict(MPI_pred_2100_ssp245, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 370
MPI_2100_ssp370 <- rast(Files_bioc_2100[3])
MPI_pred_2100_ssp370 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2100_ssp370)
ENM_MPI_2100_ssp370 <- terra::predict(MPI_pred_2100_ssp370, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)

# SPP 585
MPI_2100_ssp585 <- rast(Files_bioc_2100[4])
MPI_pred_2100_ssp585 <- c(Predictive_rasters_Balearic[[1:6]],
                           MPI_2100_ssp585)
ENM_MPI_2100_ssp585 <- terra::predict(MPI_pred_2100_ssp585, cliff_RF_model,
                                       type="prob", cores=10, na.rm=TRUE, se.fit=TRUE)
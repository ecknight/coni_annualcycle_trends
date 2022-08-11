library(tidyverse)
library(rgee)
library(sf)
library(readr)
library(lubridate)

#1. Read in region areas & points----
region <- read_sf("Data/PopulationPolygons.shp")

pts <- read.csv("Data/PopulationPoints.csv") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326)

#2. Initialize rgee----
ee_Initialize(gcs=TRUE)
ee_check()

#3. Create GCS bucket----
# project_id <- ee_get_earthengine_path() %>%
#   list.files(., "\\.json$", full.names = TRUE) %>%
#   jsonlite::read_json() %>%
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("coni_trend", projectId = project_id)

#4. Send polygons to GEE---
poly <- sf_as_ee(region)
point <- sf_as_ee(pts)

#5. Drought data----

#5a. Read in images---
drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date('2000-01-01', '2019-12-31'))$select('pdsi')$toBands()

#5b. Get mean within region polygons---
drght.mn <- drght$reduceRegions(reducer=ee$Reducer$mean(),
                                         collection=poly,
                                         scale=4638)

task_vector <- ee_table_to_gcs(collection=drght.mn,
                               bucket="coni_trend",
                               fileFormat = "CSV")
task_vector$start()
ee_monitoring(task_vector, max_attempts=1000)
ee_gcs_to_local(task = task_vector, dsn="Data/Covariates/PopulationCovariates_Polygon_drought.csv")

#5c. Get point values----
# drght.pt <- ee_extract(
#   x = drght,
#   y = point,
#   scale = 4638,
#   sf = FALSE
# )
# 
# write.csv(drght.pt, "Data/Covariates/PopulationCovariates_Point_drought.csv")

#6. ALAN data----

#6a. Get data & scale from 0 to 1----
#Older dataset, scale from 0 to 1 (dmsp is 0 to 63)
scaledsmp <- function(feature){
  feature$divide(63)
}
alan.dsmp <- ee$ImageCollection("NOAA/DMSP-OLS/NIGHTTIME_LIGHTS")$select('stable_lights')$map(scaledsmp)$toBands()

#Newer dataset, scale from 0 to 1 (viirs is -1.5 to	193565)
scaleviirs <- function(feature){
  feature$subtract(-1.5)$divide(193565.2)
}
alan.viirs <- ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG")$select('avg_rad')$map(scaleviirs)$toBands()

#6b. Stack & make an image----
alan <- alan.dsmp$addBands(alan.viirs)

#6c. Get mean within region polygons---
alan.mn <- alan$reduceRegions(reducer=ee$Reducer$mean(),
                                collection=poly,
                                scale=500)

task_vector <- ee_table_to_gcs(collection=alan.mn,
                               bucket="coni_trend",
                               fileFormat = "CSV")
task_vector$start()
ee_monitoring(task_vector, max_attempts=1000)
ee_gcs_to_local(task = task_vector, dsn="Data/Covariates/PopulationCovariates_Polygon_alan.csv")

# #6d. Get point values----
# alan.pt <- ee_extract(
#   x = alan,
#   y = point,
#   scale = 4638,
#   sf = FALSE
# )
# 
# write.csv(alan.pt, "Data/Covariates/PopulationCovariates_Point_alan.csv")

#7. Tree regrowth----

#7a. Set up year loop----
years <- c(2001:2019)
for(i in 1:(length(years)-1)){
  
  year1 <- years[i]
  year2 <- years[i+1]
  
  #7b. Get images for t and t-1----
  modis1 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(year1, '-01-01'), paste0(year1,'-12-31')))$select('LC_Type1')$mean()
  modis2 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(year2, '-01-01'), paste0(year2,'-12-31')))$select('LC_Type1')$mean()
  
  #7c. Remap to treed and untreed habitat---
  modis.nottree <- modis1$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,2,2,2,2,2,2,2,0,2,0,2,2))
  modis.tree<- modis2$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0))
  
  #7d. Take difference between years and remp----
  modis.treeconv <- modis.nottree$subtract(modis.tree)$remap(c(-1,0,1,2), c(0,0,1,0))

  #7e. Get mean within region polygons---
  tree.mn <- modis.treeconv$reduceRegions(reducer=ee$Reducer$mean(),
                                collection=poly,
                                scale=500)
  
  task_vector <- ee_table_to_gcs(collection=tree.mn,
                                 bucket="coni_trend",
                                 fileFormat = "CSV")
  task_vector$start()
  ee_monitoring(task_vector, max_attempts=1000)
  ee_gcs_to_local(task = task_vector, dsn=paste0("Data/Covariates/PopulationCovariates_Polygon_tree_", year2, ".csv"))
  
  #7f. Get point values----
  # tree.pt <- ee_extract(
  #   x = modis.treeconv,
  #   y = point,
  #   scale = 500,
  #   sf = FALSE
  # )
  # 
  # write.csv(tree.pt, paste0("Data/Covariates/PopulationCovariates_Point_tree_", year2, ".csv"))
  
}

#8. Crop conversion----

#8a. Set up year loop----
years <- c(2001:2019)
for(i in 1:(length(years)-1)){
  
  year1 <- years[i]
  year2 <- years[i+1]
  
  #8b. Get images for t and t-1----
  modis1 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(year1, '-01-01'), paste0(year1,'-12-31')))$select('LC_Type1')$mean()
  modis2 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(year2, '-01-01'), paste0(year2,'-12-31')))$select('LC_Type1')$mean()
  
  #8c. Remap to cropd and uncropd habitat---
  modis.notcrop <- modis1$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(2,2,2,2,2,2,2,2,2,2,2,0,0,2,0,2,2))
  modis.crop <- modis2$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  
  #8d. Take difference between years and remp----
  modis.cropconv <- modis.notcrop$subtract(modis.crop)$remap(c(-1,0,1,2), c(0,0,1,0))
  
  #8e. Get mean within region polygons---
  crop.mn <- modis.cropconv$reduceRegions(reducer=ee$Reducer$mean(),
                                         collection=poly,
                                         scale=500)
  
  task_vector <- ee_table_to_gcs(collection=crop.mn,
                                 bucket="coni_trend",
                                 fileFormat = "CSV")
  task_vector$start()
  ee_monitoring(task_vector, max_attempts=1000)
  ee_gcs_to_local(task = task_vector, dsn=paste0("Data/Covariates/PopulationCovariates_Polygon_crop_", year2, ".csv"))
  
  #8f. Get point values----
  # crop.pt <- ee_extract(
  #   x = modis.cropconv,
  #   y = point,
  #   scale = 500,
  #   sf = FALSE
  # )
  # 
  # write.csv(crop.pt, paste0("Data/Covariates/PopulationCovariates_Point_crop_", year2, ".csv"))
  
}

#9. Put everything together----
files <- data.frame(file = list.files("Data/Covariates")) %>% 
  separate(file, into=c("descrip", "type", "cov", "year", "filetype"), remove=FALSE) %>% 
  mutate(year = ifelse(year=="csv", NA, year)) %>% 
  dplyr::select(-descrip, -filetype)

#9a. Drought----
cov.drght <- read.csv("Data/Covariates/PopulationCovariates_Polygon_drought.csv") %>% 
  dplyr::select(-'.geo', -system.index) %>% 
  pivot_longer(cols=X200001_pdsi:X201912_pdsi, names_to="band", values_to="pdsi") %>% 
  mutate(year = as.numeric(str_sub(band, 2, 5)),
         month = as.numeric(str_sub(band, 6, 7)),
         year = ifelse(stage=="winter" & month <=4, year-1, year),
         use = case_when(stage=="breed" & month %in% c(4:8) ~ 1,
                         stage=="fall" & month %in% c(8:12) ~ 1,
                         stage=="winter" & month %in% c(11:4) ~ 1,
                         stage=="spring" & month %in% c(3:5) ~ 1)) %>% 
  dplyr::filter(use==1) %>% 
  group_by(area, pop, stage, type, year) %>% 
  summarize(pdsi.mn = mean(pdsi),
            pdsi.cv = sd(pdsi)/mean(pdsi)) %>% 
  ungroup()

#9b. ALAN----
cov.alan <- read.csv("Data/Covariates/PopulationCovariates_Polygon_alan.csv") %>% 
  dplyr::select(-'.geo', -system.index) %>% 
  pivot_longer(cols=X20140101_avg_rad:F182013_stable_lights, names_to="band", values_to="alan") %>% 
  mutate(source = ifelse(str_sub(band,1, 1)=="X", "dsmp", "viirs"),
         year = case_when(source=="dsmp" ~ str_sub(band, 2, 5),
                          source=="viirs" ~ str_sub(band, 4, 7)),
         year = as.numeric(year)) %>% 
  group_by(area, pop, stage, type, year) %>% 
  summarize(alan = mean(alan)) %>% 
  ungroup()

#9c. Tree & crop----
files.tree <- files %>% 
  dplyr::filter(cov %in% c("tree", "crop")) %>% 
  mutate(file = paste0("Data/Covariates/", file))

cov.treecrop <- read_csv(files.tree$file, id="file") %>% 
  separate(file, into=c("data", "covariates", "descrip", "shape", "cov", "year", "filetype")) %>%
  dplyr::select(area, pop, stage, type, year, cov, mean) %>% 
  pivot_wider(names_from=cov, values_from=mean) %>% 
  mutate(year = as.numeric(year))

#9d. Put together---
cov <- cov.drght %>% 
  left_join(cov.alan) %>% 
  left_join(cov.treecrop) %>% 
  dplyr::filter(year %in% c(2002:2019))

cov.year <- cov %>% 
  mutate(year = ifelse(stage %in% c("breed", "fall", "winter"), year+1, year)) %>% 
  rename(cropt1 = crop, treet1 = tree) %>% 
  left_join(cov %>% 
            dplyr::filter(stage=="breed") %>% 
            rename(cropt2 = crop, treet2 = tree) %>% 
            dplyr::select(area, pop, stage, type, year, cropt2, treet2))

#10. Join to monitoring data----

#10a. Read in data & filter----
dat <- read.csv("Data/MonitoringData_Offsets.csv") %>% 
  mutate(datetime = ymd_hms(datetime),
         year = year(datetime)) %>% 
  dplyr::filter(pop %in% unique(cov$pop))
  
#10b. Join to covariates----
dat.cov <- dat %>% 
  inner_join(cov.year)

write.csv(dat.cov, "Data/MonitoringData_Covariates.csv")
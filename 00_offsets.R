library(tidyverse)
library(ggmap)
library(lubridate)
library(suncalc)
library(detect)
library(lutz)
library(naturecounts)

options(scipen=9999)

map.theme <- theme_nothing() +
  theme(text=element_text(size=12, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text = element_blank())

#TO DO: FIX METHODOLOGY MATRIX ISSUE####
#TO DO: ADD IN RECOGNIZER DATA####

#A. CONI DATABASE DATA#####

#1. Read in CONI DB data dictionaries----
datatype <- read.csv("Offsets/DD_DataType.csv")
rectype <- read.csv("Offsets/DD_Recognizer.csv")
durtype <- read.csv("Offsets/DD_MethodDuration.csv")
durdesc <- read.csv("Offsets/DD_DescripDuration.csv")
disttype <- read.csv("Offsets/DD_MethodDistance.csv")
distdesc <- read.csv("Offsets/DD_DescripDistance.csv")

#2. Read in CONI DB data----
locs <- read.csv("Offsets/Location_SS.csv")
samp <- read.csv("Offsets/Sample_PKEY.csv") 
proj <- read.csv("Offsets/Project_PCODE.csv")
bird <- read.csv("Offsets/Bird_PC.csv") %>% 
  mutate(SPECIES = "CONI")
det <- read.csv("Offsets/Detection.csv")

#3. Identify survey type for each location----
type <- locs %>% 
  rename(PCODE = PCODE_fkey) %>% 
  left_join(proj %>% 
              rename(PCODE = PCODE_pkey) %>% 
              dplyr::select(PCODE, DATATYPE, RECOGNIZER, DURMETH, DISTMETH)) %>% 
  left_join(datatype %>% 
              rename(DATATYPE = DATAID)) %>% 
  left_join(rectype %>% 
              rename(RECOGNIZER = RECOGID)) %>% 
  left_join(durtype %>% 
              rename(DURMETH = DURMETHID)) %>% 
  mutate(type = paste0(DATA_DESC, "-", RECOG_DESC))

#4. Filter to usable surveys----
#must be time structured or recognizer-processed
#Can't be broadcast surveys
#Can't be highest score filtered (Sam's data)
#Remove a bunch with nonsensical duration data

CONI.use <- type %>% 
  dplyr::filter(PCODE %in% c(3, 5, 7, 11, 18, 28, 59)) %>% 
  left_join(samp %>% 
              rename(LOCATION_pkey = LOCATION_fkey) %>% 
              dplyr::select(-COMMENTS)) %>% 
  inner_join(bird %>% 
               rename(SAMPLE_pkey = SAMPLE_fkey) %>% 
               dplyr::select(-COMMENTS)) %>% 
  inner_join(durdesc %>% 
               rename(DURATION = DURID) %>% 
               mutate(DURATION = as.character(DURATION))) %>% 
  dplyr::filter(!(Dur_Start==0 & DUR_end==3 & DURMETH=="CC"),
                !(Dur_Start==3 & DUR_end==6 & DURMETH=="CC"),
                !(Dur_Start==8 & DUR_end==10 & DURMETH=="CC"),
                !Dur_Start=="unk",
                !DUR_end=="unk",
                !(Dur_Start==0 & DUR_end==10 & DURMETH=="G"))

#5. Select first detection of each individual and tidy----
CONI.tidy <-  CONI.use %>% 
  mutate(HR = ifelse(MIN==60, HR+1, HR),
         MIN = ifelse(MIN==60, 0, MIN),
         DD = ifelse(MM==15, 15, DD),
         MM = ifelse(MM==15, 6, MM),
         date = ymd_hm(paste(YYYY, MM, DD, HR, MIN, sep="-")),
         Max_Duration = case_when(DURMETH=="CC" ~ 6,
                                  DURMETH=="G" ~ 10,
                                  DURMETH=="II" ~ 6),
         Max_Distance = Inf) %>% 
  rename(Sample_ID = SAMPLE_pkey, Latitude = Y, Longitude = X, Time_Method = DURMETH, Start_Duration = Dur_Start, End_Duration = DUR_end, Abundance = ABUND) %>% 
  mutate(Time_Level = case_when(Time_Method=="CC" ~ as.numeric(End_Duration),
                                Time_Method=="G" & Start_Duration==0 ~ 1,
                                Time_Method=="G" & Start_Duration==3 ~ 2,
                                Time_Method=="G" & Start_Duration==5 ~ 3,
                                Time_Method=="II" & Start_Duration==0 ~ 1,
                                Time_Method=="II" & Start_Duration==3 ~ 2)) %>% 
  arrange(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Max_Duration, Abundance, BIRDID, Time_Level) %>% 
  group_by(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Max_Duration, Abundance, BIRDID) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup() %>% 
  dplyr::select(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Start_Duration, End_Duration, Max_Duration, Time_Level, Abundance)

#6. Calculate time since sunset & sunrise----
#Filter out surveys with no survey time, get local timezone
CONI.tz <- CONI.tidy %>% 
  mutate(tz=tz_lookup_coords(Latitude, Longitude, method="accurate"))

#Loop through timezones to calculate time since sunrise
tzs <- unique(CONI.tz$tz)

CONI.list <- list()
for(i in 1:length(tzs)){
  
  CONI.i <- CONI.tz %>% 
    dplyr::filter(tz==tzs[i]) %>% 
    rename(lat = Latitude, lon = Longitude) %>% 
    mutate(date = as.Date(date))
  
  CONI.i$sunset <- getSunlightTimes(data=CONI.i, keep="sunset", tz=tzs[i])$sunset
  CONI.i$sunrise <- getSunlightTimes(data=CONI.i, keep="sunrise", tz=tzs[i])$sunrise
  CONI.i$sunset <- as.POSIXct(as.character(CONI.i$sunset), tz=tzs[i])
  CONI.i$sunrise <- as.POSIXct(as.character(CONI.i$sunrise), tz=tzs[i])
  CONI.i$date <- as.POSIXct(as.character(CONI.i$date), tz=tzs[i])
  CONI.i$tsss <- as.numeric(difftime(CONI.i$date, CONI.i$sunset), units="hours")
  CONI.i$tssr <- as.numeric(difftime(CONI.i$date, CONI.i$sunrise), units="hours")
  
  CONI.list[[i]] <- CONI.i
}

CONI.sun <- do.call(rbind, CONI.list) %>% 
  mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
         tsss = ifelse(tsss < -12, tsss+24, tsss),
         tssr = ifelse(tssr > 12, tssr-24, tssr),
         tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
  dplyr::select(-tz, -sunrise, -sunset) %>% 
  rename(Latitude = lat, Longitude = lon)

hist(CONI.sun$tsss)
hist(CONI.sun$tssr) #looks good

#B. BAM DATA####

#1. Read in BAM data dictionaries----
time <- read.csv("offsets/time_lookup.csv") %>% 
  rename(Time_Method = Method, Time_Level = Level, Time_Range = Range, Time_Description = Description)

dist <- read.csv("offsets/distance_lookup.csv") %>% 
  rename(Distance_Method = Method, Distance_Level = Level, Distance_Range = Range, Distance_Description = Description)
  
#2. Read in diurnal BAM data----
counts <- read.csv("offsets/coni_counts.csv") %>% 
  left_join(time) %>% 
  left_join(dist)

#3. Filter BAM to time structured data----
BAM.use <- counts %>% 
  dplyr::filter(!Time_Description %in% c("0-5 min",
                                         "0-6 min",
                                         "1 min intervals to 10 mins"),
                !is.na(Time_Description)) %>% 
  mutate(date = ymd_hms(UTC, tz="GMT")) %>% 
  dplyr::select(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Start_Duration, End_Duration, Max_Duration, Time_Level, Abundance)
#Took out 1 min intervals to 10 minutes because there is only one obs of this time type

#4. Calculate BAM time since sunset & sunrise----
BAM.sun <- getSunlightTimes(data=BAM.use %>% 
                              rename(lat=Latitude, lon=Longitude) %>% 
                              mutate(date=as.Date(date)),
                            tz="GMT") %>% 
  dplyr::select(sunrise, sunset) %>% 
  cbind(BAM.use) %>% 
  mutate(tsss = as.numeric(difftime(date, sunset), units="hours"),
         tssr = as.numeric(difftime(date, sunrise), units="hours"),
         tsss = ifelse(tsss > 12, tsss-24, tsss),
         tsss = ifelse(tsss < -12, tsss+24, tsss),
         tssr = ifelse(tssr > 12, tssr-24, tssr),
         tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
  dplyr::select(-sunrise, -sunset)

hist(BAM.sun$tsss)
hist(BAM.sun$tssr) #looks good

#C. NIGHTJAR SURVEY NETWORK DATA####

#1. Read in data----
nsn.raw <- read.csv("offsets/NSN_CONI.csv")

#2. Select first detection of each individual and tidy----
nsn.use <- nsn.raw %>% 
  mutate(start_time = mdy_hm(start_time),
         end_time = mdy_hm(end_time),
         date = case_when(stop_number==1 ~ start_time,
                          stop_number==10 ~ end_time),
         Duration = (end_time - start_time)/9,
         date = start_time + Duration*(stop_number-1),
         Sample_ID=paste0("NSN-", survey_stop_id),
         Max_Distance = Inf,
         Time_Method="CC",
         Max_Duration=6) %>% 
  pivot_longer(cols=c(heard_minute_1:heard_minute_6),
               names_to="Time_Level",
               values_to="Abundance") %>% 
  mutate(Time_Level = as.numeric(str_sub(Time_Level, -1, -1)),
         Start_Duration = Time_Level-1,
         End_Duration = Time_Level) %>% 
  rename(Latitude = latitude, Longitude = longitude) %>% 
  dplyr::filter(Abundance > 0) %>% 
  arrange(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Max_Duration, birdacct_id, Time_Level) %>% 
  dplyr::group_by(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Max_Duration, birdacct_id) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup() %>% 
  dplyr::select(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Start_Duration, End_Duration, Max_Duration, Time_Level, Abundance) %>% 
  dplyr::filter(!is.na(date),
                Latitude!="blank") %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>% 
  dplyr::filter(Latitude > 0, Latitude < 49)

#3. Calculate time since sunset & sunrise----
#Filter out surveys with no survey time, get local timezone
nsn.tz <- nsn.use %>% 
  mutate(tz=tz_lookup_coords(Latitude, Longitude, method="accurate")) %>% 
  dplyr::filter(tz!="Asia/Shanghai; Asia/Urumqi")

#Loop through timezones to calculate time since sunrise
tzs <- unique(nsn.tz$tz)

nsn.list <- list()
for(i in 1:length(tzs)){
  
  nsn.i <- nsn.tz %>% 
    dplyr::filter(tz==tzs[i]) %>% 
    rename(lat = Latitude, lon = Longitude) %>% 
    mutate(date = as.Date(date))
  
  nsn.i$sunset <- getSunlightTimes(data=nsn.i, keep="sunset", tz=tzs[i])$sunset
  nsn.i$sunrise <- getSunlightTimes(data=nsn.i, keep="sunrise", tz=tzs[i])$sunrise
  nsn.i$sunset <- as.POSIXct(as.character(nsn.i$sunset), tz=tzs[i])
  nsn.i$sunrise <- as.POSIXct(as.character(nsn.i$sunrise), tz=tzs[i])
  nsn.i$date <- as.POSIXct(as.character(nsn.i$date), tz=tzs[i])
  nsn.i$tsss <- as.numeric(difftime(nsn.i$date, nsn.i$sunset), units="hours")
  nsn.i$tssr <- as.numeric(difftime(nsn.i$date, nsn.i$sunrise), units="hours")
  
  nsn.list[[i]] <- nsn.i
}

nsn.sun <- do.call(rbind, nsn.list) %>% 
  mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
         tsss = ifelse(tsss < -12, tsss+24, tsss),
         tssr = ifelse(tssr > 12, tssr-24, tssr),
         tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
  dplyr::select(-tz, -sunrise, -sunset) %>% 
  rename(Latitude = lat, Longitude = lon)

hist(nsn.sun$tsss)
hist(nsn.sun$tssr) #looks good

#D. CANADIAN NIGHTJAR SURVEY DATA####

#1. Get data----
#Datasets I have access to
#dbs <- nc_count(username = "ecknight")
#password: CONICOPO2010

#raw <- nc_data_dl(request_id=206938, username = "ecknight", fields_set = "extended", info = "CONI data for offsets")

#write.csv(raw, "Offsets/naturecounts.csv", row.names = FALSE)
raw <- read.csv("Offsets/naturecounts.csv")

#2. Wrangle data----
cns.use <- raw %>% 
  dplyr::filter(ScientificName=="Chordeiles minor",
                collection=="NIGHTJAR") %>% 
  rename(Latitude=DecimalLatitude, Longitude=DecimalLongitude) %>% 
  mutate(date = ymd(paste0(YearCollected,"-", MonthCollected, "-", DayCollected)),
         hour = floor(as.numeric(TimeCollected)),
         minute = round((as.numeric(TimeCollected)-hour)*60),
         date = ymd_hm(paste(date, " ", hour,":",minute))) %>% 
  dplyr::select(record_id, Latitude, Longitude, date, ObservationCount2, ObservationCount3, ObservationCount4, ObservationCount5, ObservationCount6, ObservationCount7) %>% 
  pivot_longer(cols=ObservationCount2:ObservationCount7,
               names_to="interval",
               values_to="observation") %>% 
  dplyr::filter(observation%in%c("C", "U", "V", "W"),
                !is.na(Latitude),
                !is.na(Longitude)) %>% 
  mutate(Sample_ID=paste0("CNS-", record_id),
         Max_Distance=Inf,
         Time_Method="CC",
         Start_Duration = as.numeric(str_sub(interval, -1, -1))-2,
         End_Duration = Start_Duration + 1,
         Max_Duration = 6,
         Time_Level = End_Duration,
         Abundance = 1)
  
#3. Calculate time since sunset & sunrise----
#Filter out surveys with no survey time, get local timezone
cns.tz <- cns.use %>% 
  mutate(tz=tz_lookup_coords(Latitude, Longitude, method="accurate"))

#Loop through timezones to calculate time since sunrise
tzs <- unique(cns.tz$tz)

cns.list <- list()
for(i in 1:length(tzs)){
  
  cns.i <- cns.tz %>% 
    dplyr::filter(tz==tzs[i]) %>% 
    rename(lat = Latitude, lon = Longitude) %>% 
    mutate(date = as.Date(date))
  
  cns.i$sunset <- getSunlightTimes(data=cns.i, keep="sunset", tz=tzs[i])$sunset
  cns.i$sunrise <- getSunlightTimes(data=cns.i, keep="sunrise", tz=tzs[i])$sunrise
  cns.i$sunset <- as.POSIXct(as.character(cns.i$sunset), tz=tzs[i])
  cns.i$sunrise <- as.POSIXct(as.character(cns.i$sunrise), tz=tzs[i])
  cns.i$date <- as.POSIXct(as.character(cns.i$date), tz=tzs[i])
  cns.i$tsss <- as.numeric(difftime(cns.i$date, cns.i$sunset), units="hours")
  cns.i$tssr <- as.numeric(difftime(cns.i$date, cns.i$sunrise), units="hours")
  
  cns.list[[i]] <- cns.i
}

cns.sun <- do.call(rbind, cns.list) %>% 
  mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
         tsss = ifelse(tsss < -12, tsss+24, tsss),
         tssr = ifelse(tssr > 12, tssr-24, tssr),
         tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
  dplyr::select(-tz, -sunrise, -sunset) %>% 
  rename(Latitude = lat, Longitude = lon) %>% 
  dplyr::select(colnames(nsn.sun))

hist(cns.sun$tsss)
hist(cns.sun$tssr) #looks good

#E. COMBINE DATASETS & FORMAT####

#1. Put datasets together----
all <- rbind(BAM.sun, CONI.sun, nsn.sun, cns.sun) %>% 
  unique() %>% 
  mutate(doy = yday(date)) %>% 
  dplyr::filter(doy < 230)

#2. Plot geographic distribution----
plot.dat <- all %>% 
  dplyr::select(Latitude, Longitude) %>% 
  mutate(type="human") %>% 
  rbind(CONI.use %>% 
          dplyr::select(X, Y, type) %>% 
          rename(Longitude = X, Latitude = Y)) %>% 
  mutate(method = case_when(type %in% c("human", "Human - passive-NA") ~ "human",
                            !is.na(type) ~ "ARU"))

can <- map_data("world") %>% 
  dplyr::filter(region %in% c("Canada", "USA"))

plot.distribution <- ggplot(plot.dat) +
  geom_polygon(data = can, aes(x=long, y = lat, group = group), alpha = 0.7) + 
  geom_point(aes(x=Longitude, y=Latitude, colour=method), alpha=0.5) +
  map.theme +
  theme(legend.position = "right") +
  xlab("") +
  ylab("") + 
  xlim(c(-180, -50))
plot.distribution

ggsave(plot.distribution, filename="Figs/LocationsPlot.jpeg", width=16, height=8)

#3. Plot temporal coverage----
plot.time <- ggplot(all) +
  geom_hex(aes(x=doy, y=tsss))
plot.time

#4. Format observation data----
all.m <- all %>% 
  group_by(Sample_ID, Time_Method, Time_Level) %>% 
  summarize(Abundance = sum(Abundance)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=Time_Level,
              values_from=Abundance,
              names_prefix = "t",
              values_fill=0,
              names_sort=TRUE) %>% 
  dplyr::filter(!is.na(Sample_ID))

#4. Format design data----
all.time <- all %>%  
  dplyr::select(Time_Method, Max_Duration, End_Duration, Time_Level) %>% 
  unique()
time.d <- reshape2::dcast(all.time, Time_Method + Max_Duration ~ Time_Level, value.var = "End_Duration")

all.d <- all.m %>% 
  dplyr::select(Sample_ID, Time_Method) %>% 
  left_join(time.d) 

#5. Change zeros to NAs in the observation matrix----
cols <- 1:6
for (i in cols){
  
  indices <- which(all.m[,i+2] == 0)
  
  all.m[indices, i+2] <- ifelse(is.na(all.d[indices, i+3]), NA, 0)
}

#6. Change to matrices---
m <- all.m %>% 
  dplyr::select(-Sample_ID, -Time_Method) %>% 
  as.matrix()
row.names(m) <- all.m$Sample_ID

d <- all.d %>% 
  dplyr::select(-Sample_ID, -Time_Method, -Max_Duration) %>% 
  as.matrix()
row.names(d) <- all.d$Sample_ID

#7. Format covariates----
c <- all.m %>% 
  dplyr::select(Sample_ID, Time_Method) %>% 
  left_join(all) %>% 
  mutate(doy = yday(date),
         hour = as.numeric(hour(date))) %>% 
  rename(lat = Latitude, lon = Longitude) %>% 
  dplyr::select(Sample_ID, lat, lon, doy, hour, tsss, tssr) %>% 
  mutate(doy2 = doy^2,
         tsss2 = tsss^2,
         tssr2 = tssr^2,
         ts = (hour +2)/8,
         sin = sin(ts*2*pi),
         cos = cos(ts*2*pi)) %>% 
  unique()

tsss <- c$tsss
ts <- c$ts
tssr <- c$tssr
lat <- c$lat
sin <- c$sin
cos <- c$cos

#F. MODEL####

#1. Select how to model time of day----
m1 = cmulti(m | d ~ 1 + poly(tsss, 2), type = "rem")
m2 = cmulti(m | d ~ 1 + sin, type = "rem")
m3 = cmulti(m | d ~ 1 + cos, type = "rem")
m4 = cmulti(m | d ~ 1 + sin + cos, type = "rem")
sapply(list(m1, m2, m3, m4), AIC) #tsss2



m1 = cmulti(m | d ~ 1, type="rem")
m2 = cmulti(m | d ~ 1 + c$lat, type = "rem")
m3 = cmulti(m | d ~ 1 + c$doy, type="rem")
m4 = cmulti(m | d ~ 1 + c$doy + c$doy2, type = "rem")
m5 = cmulti(m | d ~ 1 + c$tssr, type="rem")
m6 = cmulti(m | d ~ 1 + c$tssr + c$tssr2, type="rem")
m7 = cmulti(m | d ~ 1 + c$sin, type="rem")
m8 = cmulti(m | d ~ 1 + c$cos, type="rem")
m9 = cmulti(m | d ~ 1 + c$sin + c$cos, type="rem")

removal.list <- list(m1, m3, m4, m5, m6, m7, m8, m9, m10)

#2. Compare with AIC----
removal.aic <- data.frame(df=sapply(removal.list, function(z) length(coef(z))),
                          AIC=sapply(removal.list, AIC),
                          loglik = sapply(removal.list, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(c)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight)
removal.aic

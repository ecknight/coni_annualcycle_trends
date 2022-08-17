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

#TO DO: CLEAN UP MODEL SELECTION####

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
                !(Dur_Start==0 & DUR_end==10 & DURMETH=="G"),
                !HR %in% c(12, 13))

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
  dplyr::select(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Start_Duration, End_Duration, Max_Duration, Time_Level, Abundance) %>% 
  mutate(source = "National CONI database")

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
  # mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
  #        tsss = ifelse(tsss < -12, tsss+24, tsss),
  #        tssr = ifelse(tssr > 12, tssr-24, tssr),
  #        tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
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
  mutate(date = ymd_hms(UTC, tz="UTC")) %>% 
  dplyr::select(Sample_ID, Latitude, Longitude, date, Max_Distance, Time_Method, Start_Duration, End_Duration, Max_Duration, Time_Level, Abundance) %>% 
  mutate(source = "Boreal Avian Modelling project")
#Took out 1 min intervals to 10 minutes because there is only one obs of this time type

#4. Calculate BAM time since sunset & sunrise----
BAM.sun <- getSunlightTimes(data=BAM.use %>% 
                              rename(lat=Latitude, lon=Longitude) %>% 
                              mutate(date=as.Date(date)),
                            tz="GMT") %>% 
  dplyr::select(sunrise, sunset) %>% 
  cbind(BAM.use) %>% 
  mutate(tsss = as.numeric(difftime(date, sunset), units="hours"),
         tssr = as.numeric(difftime(date, sunrise), units="hours")) %>% 
  # mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
  #        tsss = ifelse(tsss < -12, tsss+24, tsss),
  #        tssr = ifelse(tssr > 12, tssr-24, tssr),
  #        tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
  dplyr::select(-sunrise, -sunset)

hist(BAM.sun$tsss)
hist(BAM.sun$tssr) #looks good

#C. NIGHTJAR SURVEY NETWORK DATA####

#1. Read in data----
nsn.raw <- read.csv("offsets/NSN_CONI.csv")
nsn.route <- read.csv("offsets/NSN_all.csv")

#2. Calculate # of stops per route----
nsn.n <- nsn.route %>% 
  dplyr::select(route_id, Year, survey_id, stop_number) %>% 
  unique() %>% 
  group_by(route_id, Year, survey_id) %>% 
  summarize(n=n()) %>% 
  ungroup()
table(nsn.n$n)
#Ok so they're all 10 stops! Weird

#2. Select first detection of each individual and tidy----
nsn.use <- nsn.raw %>% 
  mutate(start.time1 = mdy_hm(start_time),
         end.time1 = mdy_hm(end_time),
         end.time2 = case_when((hour(start.time1) > hour(end.time1)) & 
                                 (yday(start.time1)==yday(end.time1)) ~ end.time1 + days(1),
                              !is.na(end.time1) ~ end.time1),
         end.time3 = case_when(hour(start.time1) %in% c(5:12) ~ end.time2 + hours(12),
                              !is.na(end.time2) ~ end.time2),
         start.time2 = case_when(hour(start.time1) %in% c(5:12) ~ start.time1 + hours(12),
                                !is.na(start.time1) ~ start.time1),
         datetime = case_when(stop_number==1 ~ start.time2,
                          stop_number==10 ~ end.time3),
         Duration = (end.time3 - start.time2)/9,
         datetime = start.time2 + Duration*(stop_number-1),
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
  arrange(Sample_ID, Latitude, Longitude, datetime, Max_Distance, Time_Method, Max_Duration, birdacct_id, Time_Level) %>% 
  dplyr::group_by(Sample_ID, Latitude, Longitude, datetime, Max_Distance, Time_Method, Max_Duration, birdacct_id) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup() %>% 
  dplyr::filter(as.numeric(Duration) > 300,
                as.numeric(Duration) < 1500) %>% 
  dplyr::select(Sample_ID, Latitude, Longitude, datetime, Max_Distance, Time_Method, Start_Duration, End_Duration, Max_Duration, Time_Level, Abundance) %>% 
  dplyr::filter(!is.na(datetime),
                Latitude!="blank") %>% 
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>% 
  dplyr::filter(Latitude > 0, Latitude < 49) %>% 
  mutate(source = "Nightjar Survey Network")
  # mutate(starthour = hour(start.time2), 
  #        endhour = hour(end.time3),
  #        hour = hour(date),
  #        hours = abs(starthour - endhour))

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
    mutate(date = as.Date(datetime))
  
  nsn.i$sunset <- getSunlightTimes(data=nsn.i, keep="sunset", tz=tzs[i])$sunset
  nsn.i$sunrise <- getSunlightTimes(data=nsn.i, keep="sunrise", tz=tzs[i])$sunrise
  nsn.i$sunset <- as.POSIXct(as.character(nsn.i$sunset), tz=tzs[i])
  nsn.i$sunrise <- as.POSIXct(as.character(nsn.i$sunrise), tz=tzs[i])
  nsn.i$date <- as.POSIXct(as.character(nsn.i$datetime), tz=tzs[i])
  nsn.i$tsss <- as.numeric(difftime(nsn.i$date, nsn.i$sunset), units="hours")
  nsn.i$tssr <- as.numeric(difftime(nsn.i$date, nsn.i$sunrise), units="hours")
  
  nsn.list[[i]] <- nsn.i
}

nsn.sun <- do.call(rbind, nsn.list) %>% 
  # mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
  #        tsss = ifelse(tsss < -12, tsss+24, tsss),
  #        tssr = ifelse(tssr > 12, tssr-24, tssr),
  #        tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
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
         datetime = ymd_hm(paste(date, " ", hour,":",minute))) %>% 
  dplyr::select(record_id, Latitude, Longitude, datetime, ObservationCount2, ObservationCount3, ObservationCount4, ObservationCount5, ObservationCount6, ObservationCount7) %>% 
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
         Abundance = 1) %>% 
  mutate(source = "Canadian Nightjar Survey")
  
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
    mutate(date = as.Date(datetime))
  
  cns.i$sunset <- getSunlightTimes(data=cns.i, keep="sunset", tz=tzs[i])$sunset
  cns.i$sunrise <- getSunlightTimes(data=cns.i, keep="sunrise", tz=tzs[i])$sunrise
  cns.i$sunset <- as.POSIXct(as.character(cns.i$sunset), tz=tzs[i])
  cns.i$sunrise <- as.POSIXct(as.character(cns.i$sunrise), tz=tzs[i])
  cns.i$date <- as.POSIXct(as.character(cns.i$datetime), tz=tzs[i])
  cns.i$tsss <- as.numeric(difftime(cns.i$date, cns.i$sunset), units="hours")
  cns.i$tssr <- as.numeric(difftime(cns.i$date, cns.i$sunrise), units="hours")
  
  cns.list[[i]] <- cns.i
}

cns.sun <- do.call(rbind, cns.list) %>% 
  # mutate(tsss = ifelse(tsss > 12, tsss-24, tsss),
  #        tsss = ifelse(tsss < -12, tsss+24, tsss),
  #        tssr = ifelse(tssr > 12, tssr-24, tssr),
  #        tssr = ifelse(tssr < -12, tssr+24, tssr)) %>% 
  dplyr::select(-tz, -sunrise, -sunset) %>% 
  rename(Latitude = lat, Longitude = lon) %>% 
  dplyr::select(colnames(nsn.sun))

hist(cns.sun$tsss)
hist(cns.sun$tssr) #looks good

#E. COMBINE DATASETS & VISUALIZE####

#1. Put datasets together----
all <- rbind(BAM.sun, CONI.sun) %>% 
  rbind(nsn.sun %>% 
          dplyr::select(-datetime)) %>% 
  rbind(cns.sun %>% 
          dplyr::select(-datetime) %>% 
          dplyr::filter(tsss > -5)) %>% 
  unique() %>% 
  mutate(doy = yday(date),
         Start_Duration = as.numeric(Start_Duration),
         End_Duration = as.numeric(End_Duration)) %>% 
  dplyr::filter(doy < 230)

write.csv(all, "Data/OffsetData.csv", row.names = FALSE)

#2. Adjust tsss & tssr----
ggplot(all) +
  geom_point(aes(x=tsss, y=tssr, colour=Time_Method))
 
all$tsss <- ifelse(all$tsss < -5, all$tsss+24, all$tsss)
all$tssr <- ifelse(all$tssr >10, all$tssr-24, all$tssr)
 
 ggplot(all) +
   geom_point(aes(x=tsss, y=tssr, colour=Time_Method))


#3. Plot geographic distribution----
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

#4. Plot temporal coverage----
plot.tsss <- ggplot(all) +
  geom_hex(aes(x=doy, y=tsss))
plot.tsss
plot.tssr <- ggplot(all) +
  geom_hex(aes(x=doy, y=tssr))
plot.tssr

#5. Explore temporal vals----
ggplot(all) +
  geom_hex(aes(x=tsss, y=Start_Duration))
ggplot(all) +
  geom_hex(aes(x=tssr, y=Start_Duration))
ggplot(all) +
  geom_hex(aes(x=doy, y=Start_Duration))

#6. Determine peaks for tsss & tsssr----
mon <- read.csv("Data/MonitoringData.csv") %>% 
  mutate(pres = ifelse(count > 0, 1, 0),
         latr = round(lat, -1),
         latr = case_when(latr==20 ~ 30,
                          latr==70 ~ 60,
                          !is.na(latr) ~ latr))

ggplot(mon) +
  geom_smooth(aes(x=tsss, y=pres, colour=factor(latr))) +
  geom_jitter(aes(x=tsss, y=pres, colour=factor(latr))) +
  facet_wrap(~latr)

ggplot(mon) +
  geom_smooth(aes(x=tssr, y=pres, colour=factor(latr))) +
  facet_wrap(~latr)

#F. FORMAT FOR MODELLING####

#1. Format observation data----
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

#2. Format design data----
all.time <- all %>%  
  dplyr::select(Time_Method, Max_Duration, End_Duration, Time_Level) %>% 
  unique()
time.d <- reshape2::dcast(all.time, Time_Method + Max_Duration ~ Time_Level, value.var = "End_Duration")

all.d <- all.m %>% 
  dplyr::select(Sample_ID, Time_Method) %>% 
  left_join(time.d) 

#3. Change zeros to NAs in the observation matrix----
cols <- 1:6
for (i in cols){
  
  indices <- which(all.m[,i+2] == 0)
  
  all.m[indices, i+2] <- ifelse(is.na(all.d[indices, i+3]), NA, 0)
}

#4. Change to matrices---
m <- all.m %>% 
  dplyr::select(-Sample_ID, -Time_Method) %>% 
  as.matrix()
row.names(m) <- all.m$Sample_ID

d <- all.d %>% 
  dplyr::select(-Sample_ID, -Time_Method, -Max_Duration) %>% 
  as.matrix()
row.names(d) <- all.d$Sample_ID

#5. Format covariates----
c <- all.m %>% 
  dplyr::select(Sample_ID, Time_Method) %>% 
  left_join(all) %>% 
  rename(lat = Latitude, lon = Longitude) %>% 
  dplyr::select(Sample_ID, lat, lon, date, tsss, tssr) %>% 
  mutate(doy = yday(date),
         start = hour(date) + minute(date)/60) %>% 
  unique()

tsss <- c$tsss
tssr <- c$tssr
ds <- c$doy/365
lats <- (c$lat)/(64)
lons <- (c$lon+139)/(139-61)

#G. MODEL####

#1. Day of year----
md1 = cmulti(m | d ~ 1 + ds, type="rem")
md2 = cmulti(m | d ~ 1 + poly(ds, 2), type = "rem")
md3 = cmulti(m | d ~ 1 + poly(ds, 2)+ poly(ds, 2):lats, type = "rem")
md4 = cmulti(m | d ~ 1 + poly(ds, 2) + poly(ds, 2):lons, type = "rem")
md5 = cmulti(m | d ~ 1 + poly(ds, 2) + poly(ds, 2):lons + poly(ds, 2):lats, type = "rem")
md6 = cmulti(m | d ~ 1 + + poly(ds, 2) + poly(ds, 2):lons:lats, type = "rem")
sapply(list(md1, md2, md3, md4, md5, md6), AIC)
#2nd order plus latitude interaction

#2. Find best combination of sunset & sunrise peaks----
loops <- expand.grid(sins.off = c(-1:2),
                     sinr.off = c(-1:2))

removal.list.sr <- list()
for(i in 1:nrow(loops)){
  
  sins.i <- sin((tsss-loops$sins.off[i])/24*2*pi)
  sinr.i <- sin((tssr-loops$sinr.off[i])/24*2*pi)
  
  removal.list.sr[[i]]  <- cmulti(m | d ~ 1 + sins.i + sins.i:lats + sinr.i + sinr.i:lats, type="rem")
  
  print(paste0("Finished model ", i, " of ", nrow(loops)))
  
}

removal.aic <- data.frame(df=sapply(removal.list.sr, function(z) length(coef(z))),
                          AIC=sapply(removal.list.sr, AIC),
                          loglik = sapply(removal.list.sr, function(z) logLik(z))) %>%
  mutate(AICc = AIC + (2*df^2+2*df) / (length(c)-df-1),
         delta = AICc - min(AICc),
         rellike = exp(-.5*delta),
         weight = rellike/sum(rellike)) %>%
  mutate_at(c("loglik", "AICc", "delta", "weight"), ~round(., 2)) %>%
  dplyr::select(df, loglik, AICc, delta, weight) %>% 
  cbind(loops)
removal.aic

sins = sin((tsss+2)/24*2*pi)
sinr = sin((tssr-0.5)/24*2*pi)

plot(tsss, sins)
plot(tssr, sinr)

#3. Time of day----
mt1 = cmulti(m | d ~ 1, type = "rem")
mt2 = cmulti(m | d ~ 1 + poly(tsss, 2), type="rem")
mt3 = cmulti(m | d ~ 1 + poly(tsss, 2) + + poly(tsss, 2):lats, type="rem")
mt4 = cmulti(m | d ~ 1 + sins, type="rem")
mt5 = cmulti(m | d ~ 1 + sins + sins:lats, type="rem")
mt6 = cmulti(m | d ~ 1 + sinr, type="rem")
mt7 = cmulti(m | d ~ 1 + sinr + sinr:lats, type="rem")
mt8 = cmulti(m | d ~ 1 + sinr + sins, type="rem")
mt9 = cmulti(m | d ~ 1 + sinr + sinr:lats + sins + sins:lats, type="rem")

removal.list <- list(mt1, mt2, mt3, mt4, mt5, mt6, mt7, mt8, mt9)
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

#4. Full model list----
m1 = cmulti(m | d ~ 1, type="rem")
m2 = cmulti(m | d ~ 1 + poly(ds, 2) + poly(ds, 2):lats, type = "rem")
m3 = cmulti(m | d ~ 1 + sinr + sinr:lats + sins + sins:lats, type="rem")
m4 = cmulti(m | d ~ 1 + sinr + sinr:lats + sins + sins:lats + poly(ds, 2) + poly(ds, 2):lats, type="rem")

removal.list <- list(m1, m2, m3, m4, m5)
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

#H. VISUALIZE####

#1. Data frame of new data----
dat <- expand.grid(tsss = seq(-12, 12, 1),
                   tssr = seq(-12, 12, 1),
                   ds=seq(0.27, 0.59, 0.01),
                   lat = c(25, 45, 65)) %>% 
  mutate(day = ds*365,
         duration=3,
         lats = (lat)/(64),
         sins = sin((tsss+2)/24*2*pi),
         sinr = sin((tssr+0.5)/24*2*pi))

#2. Predict for survey data----
dat$pr <- predict(m4, newdata = dat, type="response")

#Turn singing rates into probabilities requires total # minutes
dat$p <- 1-exp(-dat$duration*dat$pr)

#3. Visualize----
dat$lat <- factor(dat$lat, levels=c(25, 45, 65), labels=c("25 degrees N", "45 degrees N", "65 degrees N"))

ggplot(dat) +
  geom_raster(aes(x=day, y=tsss, fill=p)) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(~lat) +
  scale_fill_viridis_c(name="Availability\n
                       for detection\n
                       in 6-minute\n
                       survey") +
  xlab("Day of year") +
  ylab("Time since sunset")

write.csv(dat, "Data/OffsetPredictionsForFigure.csv", row.names=FALSE)

#I. SAVE####

#1. Save best model----
saveRDS(m4, "Data/BestOffsetModel.RDS")

#Ok this is still far from perfect, but let's go with it for now.
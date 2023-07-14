library(tidyverse)
library(lubridate)
library(lutz)
library(suncalc)
library(detect)
library(sf)
library(data.table)
library(ggmap)

options(scipen = 99999)

#Thoughts:
#needs to be at route level so abundances are not so zero-inflated and because station-level data unavailable for US
#so how to equalize effort between survey protocols (10-12 * 6 minutes vs 50 * 3 minutes)
#and how to choose TOD to calculate offset? Take mean of all surveys
#and does analysis use repeated measures or are data points unlinked between years? If so, should I be putting EVERYTHING in?
#Don't forget Clark has observer... not sure I need this for CONI TBH

#1. Get BBS data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bbs.bird <- bbs_data$bird %>% 
  dplyr::filter(AOU==4200)

bbs <- bbs_data$route %>% 
  left_join(bbs.bird) %>% 
  mutate(count = ifelse(is.na(SpeciesTotal), 0, SpeciesTotal),
         startdate = ymd_hm(paste0(Year, "-", Month, "-", Day, " ", str_sub(StartTime, -4, -3), ":", str_sub(StartTime, -2, -1))),
         enddate = ymd_hm(paste0(Year, "-", Month, "-", Day, " ", str_sub(EndTime, -4, -3), ":", str_sub(EndTime, -2, -1))),
         interval = interval(startdate, enddate),
         meandate = int_start(interval) + (int_end(interval) - int_start(interval))/2,
         method = "BBS",
         n = 50) %>%
  rename(route = Route, lat = Latitude, lon = Longitude, obs = ObsN) %>% 
  dplyr::select(method, route, lat, lon, obs, startdate, enddate, meandate, count, n) %>% 
  dplyr::filter(!is.na(meandate), 
                !is.na(lat),
                !is.na(lon))

#2. Get CNS data----
cns.raw <- read.csv("Offsets/naturecounts.csv")

cns.bird <- cns.raw %>% 
  dplyr::filter(ScientificName=="Chordeiles minor") %>% 
  dplyr::select(SiteCode, SurveyAreaIdentifier, latitude, longitude, survey_year, survey_month, survey_day, TimeCollected, CollectorNumber, ObservationCount)

cns.stop <- cns.raw %>% 
  dplyr::select(ProtocolCode, SiteCode, SurveyAreaIdentifier, latitude, longitude, survey_year, survey_month, survey_day, TimeCollected, CollectorNumber) %>% 
  unique() %>% 
  dplyr::filter(!is.na(TimeCollected)) %>% 
  left_join(cns.bird) %>% 
  separate(TimeCollected, into=c("hour", "minute"), remove=FALSE) %>% 
  mutate(count = ifelse(is.na(ObservationCount), 0, ObservationCount),
         minute = ifelse(ProtocolCode=="Nightjars", round(as.numeric(paste0("0.",minute))*60), as.numeric(minute)),
         minute = ifelse(is.na(minute), 0, minute),
         date = ymd(paste0(survey_year, "-", survey_month, "-", survey_day)),
         correcteddate = as_date(ifelse(hour < 12, date + days(1), date)),
         datetime = ymd_hm(paste0(correcteddate, " ", hour, ":", minute))) %>% 
  rename(route = SiteCode, lat = latitude, lon = longitude, obs=CollectorNumber) %>% 
  dplyr::filter(!is.na(datetime), 
                !is.na(lat),
                !is.na(lon)) %>% 
  group_by(route, lat, lon, obs, date, datetime) %>% 
  summarize(count = sum(count)) %>% 
  ungroup()

cns <- cns.stop %>% 
  group_by(route, obs, date) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon),
            startdate = min(datetime),
            enddate = max(datetime),
            count = sum(count),
            n=n()) %>% 
  ungroup() %>% 
  mutate(interval = interval(startdate, enddate),
         meandate = int_start(interval) + (int_end(interval) - int_start(interval))/2,
         method="CNS") %>% 
  dplyr::select(method, route, lat, lon, obs, startdate, enddate, meandate, count, n)

#3. Get NSN data----
nsn.raw <- read.csv("offsets/NSN_all.csv") %>% 
  mutate(start.time1 = mdy_hm(start_time),
         end.time1 = mdy_hm(end_time),
         end.time2 = case_when((hour(start.time1) > hour(end.time1)) & 
                                 (yday(start.time1)==yday(end.time1)) ~ end.time1 + days(1),
                               !is.na(end.time1) ~ end.time1),
         end.time3 = case_when(hour(start.time1) %in% c(5:12) ~ end.time2 + hours(12),
                               !is.na(end.time2) ~ end.time2),
         start.time2 = case_when(hour(start.time1) %in% c(5:12) ~ start.time1 + hours(12),
                                 !is.na(start.time1) ~ start.time1),
         start_time = start.time2,
         end_time = end.time3) 

nsn.bird <- nsn.raw %>% 
  dplyr::filter(common_name=="Common Nighthawk") %>% 
  group_by(route_id, start_time, end_time, latitude, longitude, user_id) %>% 
  summarize(count = n()) %>% 
  ungroup()

nsn.stop <- nsn.raw %>% 
  dplyr::select(route_id, start_time, end_time, latitude, longitude, user_id) %>% 
  unique() %>% 
  left_join(nsn.bird) %>% 
  mutate(interval = interval(start_time, end_time),
         meandate = int_start(interval) + (int_end(interval) - int_start(interval))/2,
         latitude = as.numeric(latitude),
         longitude = as.numeric(longitude),
         count = ifelse(is.na(count), 0, count)) %>% 
  rename(route = route_id, obs = user_id, lat = latitude, lon = longitude, startdate = start_time, enddate = end_time) %>% 
  dplyr::filter(!is.na(meandate), 
                !is.na(lat),
                !is.na(lon))

nsn <- nsn.stop %>% 
  group_by(route, obs, meandate, startdate, enddate) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon),
            count = sum(count),
            n=n()) %>% 
  ungroup() %>% 
  mutate(method="NSN") %>% 
  dplyr::select(method, route, lat, lon, obs, startdate, enddate, meandate, count, n)

#4. Combine----
all <- rbind(bbs, cns, nsn) %>% 
  dplyr::filter(lat > 0,
                lat < 100,
                lon < 100) %>% 
  mutate(lon = ifelse(lon > 0, -lon, lon))

#5. Calculate time since sunset----

#Get local timezone
all.tz <- all %>% 
  mutate(tz=tz_lookup_coords(lat, lon, method="accurate"))

#Loop through timezones to calculate time since sunrise
tzs <- unique(all.tz$tz)

all.list <- list()
for(i in 1:length(tzs)){
  
  all.i <- all.tz %>% 
    dplyr::filter(tz==tzs[i]) %>% 
    mutate(date = as.Date(startdate))
  
  all.i$sunset <- getSunlightTimes(data=all.i, keep="sunset", tz=tzs[i])$sunset
  all.i$sunrise <- getSunlightTimes(data=all.i, keep="sunrise", tz=tzs[i])$sunrise
  all.i$sunset <- as.POSIXct(as.character(all.i$sunset), tz=tzs[i])
  all.i$sunrise <- as.POSIXct(as.character(all.i$sunrise), tz=tzs[i])
  all.i$date <- as.POSIXct(as.character(all.i$startdate), tz=tzs[i])
  all.i$tsss <- as.numeric(difftime(all.i$date, all.i$sunset), units="hours")
  all.i$tssr <- as.numeric(difftime(all.i$date, all.i$sunrise), units="hours")
  
  all.list[[i]] <- all.i
}

all.sun <- do.call(rbind, all.list) %>% 
  mutate(tsss = ifelse(tsss < -5, tsss+24, tsss),
         tssr = ifelse(tssr > 10, tssr-24, tssr)) %>% 
  dplyr::select(-tz, -sunrise, -sunset) %>% 
  dplyr::filter(!(method=="CNS" & tsss > 10),
                !(method=="NSN" & tsss > 12),
                !(method=="NSN" & tsss < -2))

ggplot(all.sun) +
  geom_histogram(aes(x=tsss, fill=method)) +
  facet_wrap(~method, scales="free_y")


#6. Create line for each survey from route-level data----
all.stop.list <- list()
for(i in 1:nrow(all.sun)){
  
  n.i <- all.sun[i,"n"]
  
  all.stop.list[[i]] <- all.sun %>% 
    slice(rep(i,n.i))
  
}

all.stop <- rbindlist(all.stop.list) %>% 
  group_by(method, route, lat, lon, obs, startdate, enddate, tsss) %>% 
  mutate(stop = row_number()) %>% 
  ungroup() %>% 
  mutate(duration = (enddate - startdate)/n,
         datetime = startdate + duration*(stop-1),
         tsss = as.numeric(tsss + duration*(stop-1)),
         tssr = as.numeric(tssr + duration*(stop-1)))

#7. Calculate offsets----

mod <- readRDS("Data/BestOffsetModel.RDS")

all.off <- all.stop %>% 
  mutate(doy = yday(datetime),
         ds = doy/365,
         length = ifelse(method=="BBS", n*3, n*6),
         sins = sin((tsss+2)/24*2*pi),
         sinr = sin((tssr+0.5)/24*2*pi),
         lats = (lat)/(64),
         pres = ifelse(count > 0, 1, 0)) %>% 
  dplyr::filter(!is.na(tsss),
                !is.na(ds))

all.off$pr <- predict(mod, newdata = all.off, type="response")

all.off$p <- 1-exp(-all.off$length*all.off$pr)

#8. Calculate mean p per route----
all.out <- all.off %>% 
  group_by(method, route, lat, lon, obs, startdate, enddate, meandate, count, n) %>% 
  summarize(pr = mean(pr),
            p = mean(p),
            tsss = mean(tsss),
            tssr = mean(tssr)) %>% 
  ungroup

write.csv(all.out, "Data/MonitoringData.csv")

#9. Visualize offsets----
ggplot(all.out) +
  geom_histogram(aes(x=duration, fill=method)) +
  facet_wrap(~method, scales="free_y")

ggplot(all.out) +
  geom_histogram(aes(x=pr, fill=method)) +
  facet_wrap(~method, scales="free_y")

ggplot(all.out) +
  geom_histogram(aes(x=p, fill=method)) +
  facet_wrap(~method, scales="free_y")

ggplot(all.out) +
#  geom_point(aes(x=count, y=p, colour=method)) +
  geom_smooth(aes(x=p, y=pres), method="lm") +
  facet_wrap(~method, scales="free")

all.out %>% 
  mutate(count.mn = count/duration) %>% 
  group_by(method) %>% 
  summarize(mean = mean(count.mn),
            sd = sd(count.mn))

#8. Assign to population----
pops <- read_sf("Data/NaturalPops.shp") %>% 
  st_make_valid()

all.sf <- all.out %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_intersection(pops) 

all.pop <- all.sf %>% 
  rename(pop = id) %>% 
  dplyr::select(method, route, obs, startdate, enddate, pop) %>% 
  data.frame() %>% 
  left_join(all.off)

#9. Save----
write.csv(all.out, "Data/MonitoringData_Offsets.csv", row.names = FALSE)

#10. Plot----
raw <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv")

bird <- raw %>% 
  dplyr::select(PinpointID, BandLat, BandLong) %>% 
  unique() %>% 
  st_as_sf(coords=c("BandLong", "BandLat"), crs=4326) %>% 
  st_intersection(pops) %>% 
  rename(pop=id) %>% 
  data.frame() %>% 
  dplyr::select(PinpointID, pop) %>% 
  rbind(data.frame(PinpointID=c(437, 442, 448, 455, 456, 457, 458, 489, 491),
                   pop = c(1, 2, 2, 6, 6, 6, 6, 6, 6)))

all.pop.use <- all.pop %>% 
  dplyr::filter(pop %in% bird$pop) %>% 
  group_by(pop, method, lon, lat) %>% 
  summarize(years = n()) %>% 
  ungroup() 

table(all.pop.use$method, all.pop.use$pop)

can <- map_data("world") %>% 
  dplyr::filter(region %in% c("Canada", "USA"))

map.theme <- theme_nothing() +
  theme(text=element_text(size=12, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text = element_blank())

all.pop.use$method <- factor(all.pop.use$method, levels=c("BBS", "CNS", "NSN"), labels=c("Breeding Bird Survey", "Canadian Nightjar Survey", "Nightjar Survey Network"))

plot.data <- ggplot(all.pop.use) +
  geom_polygon(data = can, aes(x=long, y = lat, group = group), alpha = 0.7) + 
  geom_point(aes(x=lon, y=lat, colour=factor(pop)), alpha=0.5) +
  scale_colour_viridis_d(name="Deployment\npopulation") +
  map.theme +
  theme(legend.position = "bottom") +
  xlab("") +
  ylab("") + 
  xlim(c(-180, -50)) +
  facet_wrap(~method)
plot.data

ggsave(plot.data, filename="Figs/MonitoringDataPlot.jpeg", width=24, height=8)

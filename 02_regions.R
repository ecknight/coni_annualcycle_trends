library(tidyverse)
library(sf)
library(adehabitatHR)

options(scipen=99999)

#1. Read in tracking data----
raw <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv")

#2. Get pop id for each bird----
pops <- read_sf("Data/NaturalPops.shp") %>% 
  st_make_valid()

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

#add points not within polygons manually: Yukon = region2, NWT = region1, SK = region6

#3. Pull points for each season----

#3a. Breeding----
pt.breed <- raw %>% 
  dplyr::filter(Season2 %in% c("Breed1", "Breed2")) %>% 
  dplyr::select(PinpointID, Season2, Type, Lat, Long) %>% 
  left_join(bird) %>% 
  mutate(stage = "breed")

#3b. Winter----
pt.winter <- raw %>% 
  dplyr::filter(Season2 %in% c("Winter", "Winter1", "Winter2")) %>% 
  dplyr::select(PinpointID, Season2, Type, Lat, Long) %>% 
  left_join(bird) %>% 
  mutate(stage = "winter")

#3c. Spring stopover----
pt.spring <- raw %>% 
  dplyr::filter(Season2=="SpringMig",
                Lat > 0,
                Lat < 12) %>% 
  dplyr::select(PinpointID, Season2, Type, Lat, Long) %>% 
  left_join(bird) %>% 
  mutate(stage = "spring")

#3d. Fall stopover----
pt.fall <- raw %>% 
  dplyr::filter(Season2=="FallMig",
                doy >= 290,
                doy <= 307) %>% 
  dplyr::select(PinpointID, Season2, Type, Lat, Long) %>% 
  left_join(bird) %>% 
  mutate(stage = "fall")

#3e. Put together----
pt <- rbind(pt.breed, pt.winter, pt.spring, pt.fall)

#4. Create regions for each nonbreeding season----
table(pt$pop, pt$stage)

ggplot(pt) +
  geom_point(aes(x=Long, y=Lat, colour=factor(pop))) +
  facet_wrap(~stage, scales="free")

write.csv(pt, "Data/PopulationPoints.csv", row.names=FALSE)

#4a. MCP for stages with at least 5 locations and at least 3 individuals----
pt.n <- pt %>% 
  group_by(pop, stage) %>% 
  summarize(pts = n()) %>% 
  ungroup() %>% 
  full_join(pt %>% 
              dplyr::select(pop, stage, PinpointID) %>% 
              unique() %>% 
              group_by(pop, stage) %>% 
              summarize(inds = n()) %>% 
              ungroup())

pt.mcp <- pt.n %>% 
  dplyr::filter((pts >= 5 & inds >= 3),
                stage!="breed") %>% 
  left_join(pt) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  mutate(ID = paste0(pop, "-", stage)) %>% 
  dplyr::select(ID) %>% 
  as_Spatial()

mcp.95 <- mcp(pt.mcp, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  separate(id, into=c("pop", "stage")) %>% 
  mutate(type="mcp")

mcp.95 %>% 
  data.frame() %>% 
  group_by(stage) %>% 
  summarize(mean = mean(area),
            sd = sd(area)) %>% 
  mutate(radius = sqrt(mean/pi))

table(mcp.95$pop, mcp.95$stage)

#4b. Buffer with mean area radius for other pts----
pt.buff <- pt.n %>% 
  dplyr::filter((pts < 5 | inds < 3),
                stage!="breed") %>% 
  left_join(pt) %>% 
  group_by(pop, stage) %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat)) %>% 
  ungroup() %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=350000) %>% 
  mutate(type="buffer", 
         area = pi*350^2) %>% 
  dplyr::select(colnames(mcp.95))

#5. Put polygons together with breeding region natural populations----
pt.area <- pops %>% 
  st_area() %>% 
  cbind(pops) %>% 
  dplyr::filter(id %in% unique(pt$pop)) %>% 
  rename(pop=id, area='.') %>% 
  mutate(stage = "breed", 
         area = as.numeric(area)/1000000,
         type = "natpop") %>% 
  dplyr::select(colnames(mcp.95)) %>% 
  st_transform(3857) %>% 
  rbind(mcp.95, pt.buff)

ggplot(pt.area) +
  geom_sf(aes(fill=pop)) +
  facet_wrap(~stage)

write_sf(pt.area, "Data/PopulationPolygons.shp")

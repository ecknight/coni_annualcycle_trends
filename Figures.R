library(tidyverse)
library(sf)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))


#Availability offset figure----

dat <- read.csv("Data/OffsetPredictionsForFigure.csv")

ggplot(dat) +
  geom_raster(aes(x=day, y=tsss, fill=p)) +
  geom_hline(aes(yintercept=0)) +
  ylim(c(-2, 10)) +
  facet_wrap(~lat) +
  scale_fill_viridis_c(name="Availability\nfor detection in\n3 minute survey") +
  xlab("Day of year") +
  ylab("Time since sunset") +
  my.theme +
  theme(legend.position = "right")

ggsave("Figs/Availability.jpeg", width = 10, height = 6)

#Removal modelling data figure----
all <- read.csv("Data/OffsetData.csv")

plot.dat <- all %>% 
  dplyr::select(Latitude, Longitude, source) %>% 
  unique()

can <- map_data("world") %>% 
  dplyr::filter(region %in% c("Canada", "USA"))

plot.distribution <- ggplot(plot.dat) +
  geom_polygon(data = can, aes(x=long, y = lat, group = group), alpha = 0.5) + 
  geom_point(aes(x=Longitude, y=Latitude, colour=source), alpha=0.5) +
  scale_colour_viridis_d(option="plasma", name="Database") +
  map.theme +
  theme(legend.position = "left")  +
  xlab("") +
  ylab("") + 
  xlim(c(-170, -50)) +
  ylim(c(25, 80))
plot.distribution

ggsave(plot.distribution, filename="Figs/RemovalModelData.jpeg", width=7, height=4)

#Population trend drivers----
coefs <- read.csv("MLFullModelResults_Summary.csv") %>% 
  mutate(method = factor(method, levels=c("BBS", "All"), labels=c("BBS only", "Integrated")),
         cov = factor(cov, levels=c("pdsis", "trees", "crops"),
                      labels = c("Drought index", "Forest growth", "Crop conversion")),
         season = factor(season, levels=c("breed", "fall", "winter", "spring"),
                         labels=c("Breeding", "Postbreeding\nmigration", "Nonbreeding", "Prebreeding\nmigration")),
         est = ifelse(cov=="Drought index", -est, est))

plot.coefs <- ggplot(coefs) +
  geom_raster(aes(x=season, y=cov, fill=sig*100)) +
  facet_wrap(~method) +
  scale_fill_gradient2(low="red", mid="white", high = "blue", midpoint=0, name="Direction of\neffect &\n% of\npopulations\nwith p < 0.05") +
  my.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=45, vjust = 1, hjust = 1),
        axis.text.y=element_text(angle=45, vjust = 0.5, hjust = 1))
plot.coefs

ggsave(plot.coefs, filename="Figs/EffectResults.jpeg", width=8, height=4)

#Population trend drivers - spatial----
regions <- read_sf("Data/PopulationPolygons_Results.shp") %>% 
  st_transform(crs=4326) %>% 
  dplyr::filter(method=="All") %>% 
  mutate(method = factor(method, levels=c("BBS", "All"), labels=c("BBS only", "Integrated")),
                  cov = factor(cov, levels=c("pdsis", "trees", "crops"),
                               labels = c("Drought index", "Forest growth", "Crop conversion")),
                  stage = factor(stage, levels=c("breed", "fall", "winter", "spring"),
                                  labels=c("Breeding", "Postbreeding\nmigration", "Nonbreeding", "Prebreeding\nmigration")),
                  est = ifelse(cov=="Drought index", -est, est),
                  direction = ifelse(est < 0, -1, 1),
                  plot = direction * sig) 

world <- map_data("world")

plot.region <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), alpha = 0.5) + 
  geom_sf(data=regions, aes(fill=est, alpha = factor(sig))) +
#  geom_sf(data=dplyr::filter(regions, sig==1), aes(fill=est)) +
  facet_grid(stage~cov, switch="y") +
  scale_fill_gradient2(low="red", mid="white", high = "blue", midpoint=0, name="Effect on\nabundance") +
  scale_alpha_manual(values = c(0.3, 1.0), labels=c("Not significant", "Significant"), name="") +
  my.theme +
  theme(legend.position = "bottom",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())  +
  xlab("") +
  ylab("") +
  xlim(c(-170, -30)) +
  ylim(c(-60, 80)) 
plot.region

ggsave(plot.region, filename="Figs/EffectResults - Spatial.jpeg", width = 12, height = 16)


#Population trend results----
regions <- read_sf("Data/PopulationPolygons_Results_Breed.shp") %>% 
  st_transform(crs=4326) %>% 
  mutate(method = factor(method, levels=c("BBS", "All"), labels=c("BBS only", "Integrated")),
         sig = factor(sig, levels=c(1,0), labels=c("Significant", "Nonsignificant")))

can <- map_data("world") %>% 
  dplyr::filter(region %in% c("Canada", "USA"))

plot.trend <- ggplot() +
  geom_polygon(data = can, aes(x=long, y = lat, group = group), alpha = 0.5) + 
  geom_sf(data=regions, aes(fill=est)) +
  facet_grid(sig~method, switch="y") +
  scale_fill_gradient2(low="red", mid="white", high = "blue", midpoint=0, name="Population\ntrend") +
  my.theme +
  theme(legend.position = "bottom",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())  +
  xlab("") +
  ylab("") + 
  xlim(c(-170, -50)) +
  ylim(c(25, 80)) 

ggsave(plot.trend, filename="Figs/TrendResults.jpeg", width=6, height=6)

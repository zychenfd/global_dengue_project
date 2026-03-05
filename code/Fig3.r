#==load package==
library(RevGadgets)
library(tracerer)
library(stringr)
library(tidyverse)
library(readxl)
library(Biostrings)
library(rworldmap)
library(lubridate)
library(ggpubr)
library(ape)
library(ggtree)
library(phangorn)
library(phytools)
library(treeio)
library(zoo)
library(ggsci)
library(scales)
library(sf)
library(rnaturalearth)

theme_plot <- theme(panel.grid = element_blank(),
                    text = element_text(size = 7),
                    panel.border = element_rect(linewidth = 0.2, fill = "transparent"),
                    plot.background = element_rect(linewidth = 0.2, fill = "transparent"),
                    axis.line.x = element_line(color = "black", linewidth = 0.2),
                    axis.line.y = element_line(color = "black", linewidth = 0.2),
                    axis.ticks = element_line(color = "black", linewidth = 0.2))

#==1. define color==
color <- c("#8DD3C7","#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
color1 <- c(pal_aaas("default", alpha = 1)(10))
col_value <- c("Support" = color1[5],
               "Strong support" = color1[1],
               "Very strong support" = color1[4],
               "Decisive support" = color1[6])
col_value1 <- c("Florida, US" = "#8DD3C7",
                "West Europe" = "#CCEBC5",
                "South Europe" = "#00A6FF",
                "Endemic America" = "#FB8072",
                "Uruguay/North Argentina" = "#FCCDE5",
                "Endemic Africa" = "#FDB462",
                "Central China" = "#AEA200",
                "South China" =  "#B3DE69",
                "Endemic Asia" = "#BC80BD",
                "Queensland, Australia" = "#80B1D3",
                "Japan" =  "#FFFFB3")

#==2. read data==(Data are not provided in GitHub)
state <- as.data.frame(readRDS("../data/map_data/map_fig2.rds"))[,1:2] %>% distinct() #get sub-location and region name
MJ_geo1 <- readRDS("../data/genome/final_data_202512/Markov_jump/DENV_Jumps_geo1_per_period.rds") %>% #sub-location level
  mutate(period = str_replace_all(period,"P","epoch")) %>%
  left_join(state, by = c("startLocation" = "geo1")) %>%
  left_join(state, by = c("endLocation" = "geo1"))

MJ_geo3 <- readRDS("../data/genome/final_data_202512/Markov_jump/DENV_Jumps_geo3_per_period.rds") %>% #region level
  mutate(period = str_replace_all(period,"P","epoch"))

#==3. Markov Jump comparsion between sub-location and region==
MJ_1 <- MJ_geo1 %>%
  group_by(serotype, scheme, region_final.x, region_final.y, period) %>%
  summarise(jump_mean_per_year_geo1 = sum(jump_mean_per_year)) %>%
  filter(region_final.x != region_final.y) %>%
  rename("region_final.x" = "startLocation", "region_final.y" = "endLocation") %>%
  mutate(startLocation = sapply(str_split(str_remove_all(startLocation, " |/|-"),","), function(x) x[1])) %>% 
  mutate(startLocation = ifelse(startLocation %in% "UruguayNorthArgentina", "NorthArgentinaUruguay", startLocation)) %>%
  mutate(endLocation = sapply(str_split(str_remove_all(endLocation, " |/|-"),","), function(x) x[1])) %>% 
  mutate(endLocation = ifelse(endLocation %in% "UruguayNorthArgentina", "NorthArgentinaUruguay", endLocation)) %>%
  left_join(MJ_geo3[,c(1:4,7:8)]) %>%
  filter(scheme == "phylogeo_sampling_infection")

MJ_2 <- MJ_1 %>%
  group_by(scheme, startLocation, endLocation) %>%
  summarise(jump_mean_per_year_geo1 = sum(jump_mean_per_year_geo1),
            jump_mean_per_year_geo3 = sum(jump_mean_per_year))

range(MJ_2$jump_mean_per_year_geo1)
range(MJ_2$jump_mean_per_year_geo3)
cut(MJ_2$jump_mean_per_year_geo1, breaks = c(-0.1, 5, 10, 25, 50, 200),right = T,
    labels = c("[0, 5]" ,"(5, 10]", "(10, 25]","(25, 50]","> 50")) -> MJ_2$jump_mean_per_year_geo1_2
cut(MJ_2$jump_mean_per_year_geo3, breaks = c(-0.1,  5, 10, 25, 50, 200),right = T,
    labels = c("[0, 5]", "(5, 10]", "(10, 25]","(25, 50]","> 50")) -> MJ_2$jump_mean_per_year_geo3_2
factor(MJ_2$startLocation, levels = rev(c("Florida","WestEurope","SouthEurope","CentralChina","SouthChina","Japan",
                                      "EndemicAmerica","EndemicAfrica","EndemicAsia",
                                      "NorthArgentinaUruguay","Queensland" ))) -> MJ_2$startLocation
factor(MJ_2$endLocation, levels = rev(c("Florida","WestEurope","SouthEurope","CentralChina","SouthChina","Japan",
                                      "EndemicAmerica","EndemicAfrica","EndemicAsia",
                                      "NorthArgentinaUruguay","Queensland" ))) -> MJ_2$endLocation
text <- rev(c("Florida, US","West Europe","South Europe","Central China","South China","Japan",
              "Endemic America","Endemic Africa","Endemic Asia",
              "Uruguay/North Arg","Queensland, Aus"))
text1 <- round(cor.test(MJ_1$jump_mean_per_year, MJ_1$jump_mean_per_year_geo1)[4][[1]], 3)

ggplot(MJ_1) +
  geom_smooth(mapping = aes(x = jump_mean_per_year, y = jump_mean_per_year_geo1), method = "lm",color =  "grey", linewidth = 0.3)+
  geom_point(aes(x = jump_mean_per_year, y = jump_mean_per_year_geo1, color = serotype),alpha = 0.6)+
  annotate("text", x = 1, y = max(MJ_1$jump_mean_per_year_geo1, MJ_1$jump_mean_per_year)*0.9, label = paste0("ρ = ",text1), hjust = 0, size = 2)+
  scale_x_continuous(breaks = seq(0,25,5)) +
  scale_y_continuous(breaks = seq(0,25,5)) +
  theme_bw()+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(face = "bold", margin = margin(b = 2, t = 0)),
        legend.text = element_text(margin = margin(r = 0, l = 0)),
        legend.key.size = unit(0.3,"cm"),
        legend.title.position = "top",
        plot.tag = element_text(size = 8, face = "bold"))+
  theme_plot+
  labs(x = "Intensity of Markov jump inferred\nunder regional levels",
       y = "Intensity of Markov jump inferred\nunder sub-location levels",
       tag = "c")+
  scale_color_manual("Serotypes", values = color1) -> fig3

ggplot()+
  geom_tile(data = MJ_2 %>% filter(scheme == "phylogeo_sampling_infection"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year_geo1))+
  scale_fill_gradientn("Intensity of\nMarkov jump",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5")),
                       na.value="white",limits = c(0, 102))+
  theme_bw()+
  theme_plot+
  scale_x_discrete(label = text)+
  scale_y_discrete(label = text)+
  theme(plot.margin = margin(0,0.15,0,0, "cm"),
        text = element_text(size = 7),
        panel.grid = element_blank(),
        legend.title = element_text(face = "bold", margin = margin(b = 2, t = 0)),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.key.size = unit(0.3,"cm"),
        legend.title.position = "top",
        legend.frame = element_rect(linewidth = 0.2, fill = "transparent"),
        legend.ticks = element_line(linewidth = 0.2),
        plot.tag = element_text(size = 8, face = "bold"))+
  labs(x = "To",
       y = "From",
       subtitle = "Inferred under sub-location levels",
       tag = "a")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black"))-> fig1
ggplot()+
  geom_tile(data = MJ_2 %>% filter(scheme == "phylogeo_sampling_infection"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year_geo3)) +
  scale_fill_gradientn("Intensity of\nMarkov jump",colors = rev(c("#E64B35FF","#0061A3", "#CDDDE5")),
                       na.value="white",limits = c(0, 102))+
  theme_bw()+
  theme_plot+
  scale_x_discrete(label = text)+
  scale_y_discrete(label = text)+
  theme(plot.margin = margin(0,0.15,0,0, "cm"),
        text = element_text(size = 7),
        panel.grid = element_blank(),
        legend.title = element_text(face = "bold", margin = margin(b = 2, t = 0)),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title.position = "top",
        legend.frame = element_rect(linewidth = 0.2, fill = "transparent"),
        legend.ticks = element_line(linewidth = 0.2),
        plot.tag = element_text(size = 8, face = "bold"))+
  labs(x = "To",
       y = "From",
       subtitle = "Inferred under regional levels",
       tag = "b")+
  guides(fill = guide_colorbar(frame.colour = "black",
                               ticks.colour = "black")) -> fig2

#==4. Read more data at regional level==
BF <- read.csv("../data/genome/final_data_202512/Bayes_factors/DENV_BF_geo3_per_period.csv")
dat <- left_join(readRDS("../data/genome/final_data_202512/Markov_jump/DENV_Jumps_geo3_per_period.rds") %>%
                   mutate(period = str_replace_all(period,"P","epoch")) %>%
                   rename("startLocation" = "FROM",
                          "endLocation" = "TO",
                          "period" = "epoch"),
                 BF) %>%
  mutate(jump_mean_per_year1 = ifelse(BAYES_FACTOR < 3, NA, jump_mean_per_year))

#==5. map==
#global map
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  mutate(adm0_a3 = ifelse(adm0_a3 == "SDS", "SSD" ,adm0_a3)) %>%
  filter(!sovereignt %in% c("China", "Taiwan")) %>% 
  filter(admin != "France") %>%
  filter(admin != "Antarctica")
world <- st_transform(world, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

#regional map
map <- readRDS("../data/map_data/map_fig2.rds")
map_island <- map %>% filter(geo1 == "OceaniaIslands") %>% filter(admin != "Papua New Guinea")
map <- map %>%
  filter(!(geo1 == "OceaniaIslands" &  admin != "Papua New Guinea")) %>%
  group_by(region_final) %>%
  summarise()
map1 <- rbind(map, map_island[,c(1,3)])
map1 <- st_transform(map1, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

#China map
china <- st_read("../../../../data/ChinaAdminDivisonSHP-master/1. Country/country.shp")
china <- st_transform(china, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
nine <- st_read("../../../../../data/Geographic/nine.shp")
nine <- st_transform(nine, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

#--central point--
st_centroid(map) -> cen_point
cen_point$point <- as.character(cen_point$geometry)
cen_point$lon <-  as.numeric(str_trim(sapply(str_split(str_remove_all(cen_point$point,"\\(|\\)|c"),","),function(x) x[1])))
cen_point$lat <-  as.numeric(str_trim(sapply(str_split(str_remove_all(cen_point$point,"\\(|\\)|c"),","),function(x) x[2])))

cen_point$lat[cen_point$region_final == "South China"] <- cen_point$lat[cen_point$region_final == "South China"] - 1.5
cen_point$lat[cen_point$region_final == "Central China"] <- cen_point$lat[cen_point$region_final == "Central China"] + 2
cen_point$lon[cen_point$region_final == "Central China"] <- cen_point$lon[cen_point$region_final == "Central China"] + 6
cen_point$lat[cen_point$region_final == "Endemic Africa"] <- cen_point$lat[cen_point$region_final == "Endemic Africa"] - 5
cen_point$lon[cen_point$region_final == "Endemic Africa"] <- cen_point$lon[cen_point$region_final == "Endemic Africa"] + 4
cen_point$lat[cen_point$region_final == "South Europe"] <- cen_point$lat[cen_point$region_final == "South Europe"] + 1.5
cen_point$lon[cen_point$region_final == "Florida, US"] <- cen_point$lon[cen_point$region_final == "Florida, US"] + 0.9
cen_point$lon[cen_point$region_final == "Japan"] <- cen_point$lon[cen_point$region_final == "Japan"] + 2.5
cen_point$lat[cen_point$region_final == "Endemic Asia"] <- 0
cen_point$lon[cen_point$region_final == "Endemic Asia"] <- 114
cen_point <- as.data.frame(cen_point)
cen_point <- cen_point[,c(1,4,5)]

#==6. Final data==
# Significant diffusion rates were identified using Bayes factors (BF): ≥1,000 was deemed as decisive support, 
# 100 ≤ BF < 1000 as very strong support, 10 ≤ BF < 100 as strong support and 3 ≤ BF < 10 as supported. 
# Ref: https://pubmed.ncbi.nlm.nih.gov/31843889/

factor(cen_point$region_final, levels = unique(cen_point$region_final)[c(5,11,9,1,8,6,3,2,4,10,7)]) -> cen_point$region_final
dat_final <- left_join(dat, 
                       cen_point %>% mutate(region_final = sapply(str_split(str_remove_all(region_final, " |/|-"),","), function(x) x[1])) %>% 
                         mutate(region_final = ifelse(region_final %in% "UruguayNorthArgentina", "NorthArgentinaUruguay", region_final)) ,
                       by = c("FROM" = "region_final") ) %>%
  left_join(cen_point %>% mutate(region_final = sapply(str_split(str_remove_all(region_final, " |/|-"),","), function(x) x[1])) %>% 
              mutate(region_final = ifelse(region_final %in% "UruguayNorthArgentina", "NorthArgentinaUruguay", region_final)) ,
            by = c("TO" = "region_final") ) %>%
  mutate(BAYES_FACTOR1 = cut(BAYES_FACTOR, breaks = c(-Inf,3,10,100,1000,Inf),
                             right = T, labels = c("No support", "Support", "Strong support", 
                                                   "Very strong support", "Decisive support"))) 
min_v <- 0.5
max_v <- max(dat_final$jump_mean_per_year[dat_final$BAYES_FACTOR >= 10 & dat_final$scheme == "phylogeo_sampling_infection"])

deme1 <- dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1")
unique(deme1$FROM)
deme2 <- dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2")
unique(deme2$FROM)
deme3 <- dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3")
unique(deme3$FROM)
deme4 <- dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4")
unique(deme4$FROM)

#==7. figure plot==
#DENV1
ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(epoch == "epoch1") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4, order = 1,keyheight = unit(0, "cm")),
         color = guide_legend(ncol = 1, order = 2,keyheight = unit(0, "cm")),
         linewidth = guide_legend(ncol = 1, order = 3,keyheight = unit(0, "cm")))+
  theme_void()+
  geom_point(data = cen_point,
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump                                         ", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  labs(subtitle = "Pre-pandemic", y = "DENV1")+
  scale_color_manual("Bayes factors (BF)", values = col_value)-> p1

ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(epoch == "epoch2") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point,
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  labs(subtitle = "Pandemic")+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p2


ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(epoch == "epoch3") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4, order = 1,keyheight = unit(0.1, "cm")),
         color = guide_legend(ncol = 1, order = 2,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1, order = 3,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point,
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  labs(subtitle = "Post-pandemic")+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p3

ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China", "Uruguay/North Argentina")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2") %>% filter(epoch == "epoch1") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina")),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  labs(y = "DENV2")+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p4

ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China", "Uruguay/North Argentina")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2") %>% filter(epoch == "epoch2") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina", "Uruguay/North Argentina")), 
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p5


ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China", "Uruguay/North Argentina")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2") %>% filter(epoch == "epoch3") %>%
               filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4, order = 1,keyheight = unit(0, "cm")),
         color = guide_legend(ncol = 1, order = 2,keyheight = unit(0, "cm")),
         linewidth = guide_legend(ncol = 1, order = 3,keyheight = unit(0, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina")),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump                                         ", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p6


ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China","Uruguay/North Argentina","West Europe")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3") %>% filter(epoch == "epoch1") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina","West Europe")),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  labs(y = "DENV3")+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p7

ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China", "Uruguay/North Argentina","West Europe")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3") %>% filter(epoch == "epoch2") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina","West Europe")),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p8


ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China", "Uruguay/North Argentina","West Europe")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3") %>% filter(epoch == "epoch3") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina","West Europe")),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p9


ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China","Uruguay/North Argentina","West Europe", "Central China", "Endemic Africa", "South Europe")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4") %>% filter(epoch == "epoch1") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina","West Europe", "Central China", "Endemic Africa", "South Europe" )),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  labs(y = "DENV4")+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p10

ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China","Uruguay/North Argentina","West Europe", "Central China", "Endemic Africa", "South Europe")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4") %>% filter(epoch == "epoch2") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina","West Europe", "Central China", "Endemic Africa", "South Europe" )),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p11


ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China","Uruguay/North Argentina","West Europe", "Central China", "Endemic Africa", "South Europe")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4") %>% filter(epoch == "epoch3") %>% filter(BAYES_FACTOR >= 3) %>% filter(jump_mean_per_year >= 0.5),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4,keyheight = unit(0.4, "cm")),
         color = guide_legend(ncol = 1,keyheight = unit(0.4, "cm")),
         linewidth = guide_legend(ncol = 1,keyheight = unit(0.01, "cm")))+
  theme_void()+
  geom_point(data = cen_point %>% filter(!region_final %in% c("Uruguay/North Argentina","West Europe", "Central China", "Endemic Africa", "South Europe" )),
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p12


#legend
main_data <- data.frame(
  x1 = c(100, 50),  
  y1 = c(100, 90),   
  x2 = c(50, 100),  
  y2 = c(100, 90))
ggplot(main_data) +
  geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    curvature = 0.3, 
    size = 0.1,
    arrow = arrow(length = unit(0.1, "cm"), type = "closed"))+
  scale_y_continuous(limits = c(0,200))+
  coord_cartesian(clip = "off")+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  annotate("text", x = 50 ,y = 150, label = "Dispersion direction", hjust = 0, size = 2.4) -> p0

ggplot(main_data) +
  scale_x_continuous(limits = c(0,100))+
  coord_cartesian(clip = "off")+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  annotate("text", x = 50 ,y = 150, label = "Intensity of Markov jump inferred\nunder regional levels", hjust = 0.5, size = 2.3) -> p0_tmp

ggplot() + 
  annotate("rect", xmin = -180, xmax = 180, ymin = -23.27, ymax = 23.27,size = 0.1, alpha = 0.3, fill = "lightblue")+
  annotate("rect", xmin = -180, xmax = 180, ymin = -29, ymax = -50,size = 0.1, alpha = 0.1, fill = "red")+
  annotate("rect", xmin = -180, xmax = 180, ymin = 29, ymax = 50,size = 0.1, alpha = 0.1, fill = "red")+
  geom_sf(data = world, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(!region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = china, color = "black", fill = "white", linewidth = 0.01)+
  geom_sf(data = map1 %>% filter(region_final %in% c("Central China","South China")),
          fill = "grey" ,color = "black", size = 0.02, lwd = 0.01)+
  geom_sf(data = nine, color="black",size = 0.01, lwd = 0.005)+
  scale_x_continuous(expand = c(0.01,0),limits=c(-180,180))+
  scale_y_continuous(expand = c(0.01,0))+
  geom_curve(data = dat_final %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(epoch == "epoch1"),
             aes(x = as.double(lon.x), 
                 y = as.double(lat.x), 
                 xend = as.double(lon.y), 
                 yend = as.double(lat.y),
                 linewidth = jump_mean_per_year,
                 color = BAYES_FACTOR1),
             # arrow = arrow(type = "open", length = unit(1.5, "mm"), angle = 30),
             alpha = 1,
             curvature = 0.35)+
  scale_fill_manual("Geographic regions", values = col_value1)+
  guides(fill = guide_legend(nrow = 4, order = 1,keyheight = unit(0, "cm")),
         color = guide_legend(ncol = 1, order = 2,keyheight = unit(0, "cm")),
         linewidth = guide_legend(ncol = 1, order = 3,keyheight = unit(0, "cm")))+
  theme_void()+
  geom_point(data = cen_point,
             aes(x = lon, y = lat, fill = region_final), shape = 21, size = 1.2, stroke = 0.02)+
  theme(plot.margin = margin(0,0,0,0, "cm"),
        text = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title.position = "top",
        legend.position = "bottom",
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.spacing.x = unit(0.05, "cm"),
        # axis.title.y = element_text(angle = 90),
        plot.tag = element_text(size = 8, face = "bold"))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Intensity of Markov jump                                         ", range = c(0.15,1.25), limits = c(min_v, max_v), breaks = c(5,10,15,20))+
  scale_color_manual("Bayes factors (BF)", values = col_value) -> p_legend
legend1 <- get_legend(p_legend) 
p_legend <- ggplot() + theme_void() + cowplot::draw_grob(legend1)

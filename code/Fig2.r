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
library(patchwork)

#==1. define color and geographic level==
color <- rev(c("#FB8072", "#FDB462", "#BC80BD", "#80B1D3"))
color1 <- c(pal_aaas("default", alpha = 1)(10))
col_value <- c("[0, 0.05]" = "grey92",
               "(0.05, 0.1]" = "#CDDDE5",
               "(0.1, 0.5]" = "#7BABCA",
               "(0.5, 1]" = "#2879B0",
               "(1, 5]" =  "#2E5C8C",
               "> 5" =  "#FDB462")

map <- as.data.frame(readRDS("../data/map_data/map_fig2.rds"))[,1:2] %>% distinct()
levels <- rev(c(map$geo1[map$region_final == "Endemic America"],map$geo1[map$region_final == "Endemic Africa"],unique(map$geo1[map$region_final == "Endemic Asia"]),
                "Florida", "NorthArgentinaUruguay","SouthEurope", map$geo1[map$region_final == "South China"], "Queensland",
                "WestEurope", "CentralChina", "Japan"))

#==2. read data==(Data are not provided in GitHub)
MJ <- readRDS("../data/genome/final_data_202512/Markov_jump/DENV_Jumps_geo1_per_period.rds") %>%
  filter(!startLocation %in% c("WestEurope", "CentralChina", "Japan")) %>%
  filter(!endLocation %in% c("WestEurope", "CentralChina", "Japan")) %>%
  mutate(period = str_replace_all(period,"P","epoch")) %>%
  left_join(map, by = c("startLocation" = "geo1"))  %>%
  left_join(map, by = c("endLocation" = "geo1"))

MJ_emp <- MJ %>% 
  filter(region_final.x == region_final.y) %>%
  group_by(scheme, period, region_final.x) %>%
  summarise(mean = mean(jump_mean_per_year),
            n = n())

factor(MJ$startLocation, levels = levels) -> MJ$startLocation
factor(MJ$endLocation, levels = levels) -> MJ$endLocation
range(MJ$jump_mean_per_year, na.rm = T)
cut(MJ$jump_mean_per_year, breaks = c(-0.1, 0.05, 0.1, 0.5, 1, 5, 25),right = T,
    labels = c("[0, 0.05]","(0.05, 0.1]", "(0.1, 0.5]" ,"(0.5, 1]", "(1, 5]","> 5")) -> MJ$jump_mean_per_year1

#==3. Collate data==
deme1 <- data.frame(deme = levels) %>% filter(deme %in% unique(MJ$startLocation[MJ$scheme == "phylogeo_sampling_infection" & MJ$serotype == "DENV1"])) %>%
  left_join(map, by = c("deme" = "geo1")) %>%
  mutate(region_final1 = ifelse(deme %in% c("Florida", "NorthArgentinaUruguay","SouthEurope", map$geo1[map$region_final == "South China"], "Queensland"),
                                "Fringe  \nareas   ", region_final)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Asia"), "Endemic\nAsia    ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Africa"), "Endemic\nAfrica   ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic America"), "Endemic\nAmerica", region_final1)) %>%
  mutate(position = 1:n()) %>%
  group_by(region_final1) %>%
  mutate(num = n()) %>%
  ungroup() %>%
  group_by(region_final1) %>%
  mutate(mid_num = (min(position)+max(position))/2) %>%
  mutate(position1 = ifelse( (position <= mid_num + (num-3)/2) & (position >= mid_num - (num-3)/2), NA, position))
deme1_pos_min <- sort((deme1 %>% group_by(region_final1) %>% summarise(min(position1, na.rm = T)) %>% select(2))[[1]])
deme1_pos_max <- sort((deme1 %>% group_by(region_final1) %>% summarise(max(position1, na.rm = T)) %>% select(2))[[1]])

deme2 <- data.frame(deme = levels) %>% filter(deme %in% unique(MJ$startLocation[MJ$scheme == "phylogeo_sampling_infection" & MJ$serotype == "DENV2"])) %>%
  left_join(map, by = c("deme" = "geo1")) %>%
  mutate(region_final1 = ifelse(deme %in% c("Florida", "NorthArgentinaUruguay","SouthEurope", map$geo1[map$region_final == "South China"], "Queensland"),
                                "Fringe  \nareas   ", region_final)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Asia"), "Endemic\nAsia    ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Africa"), "Endemic\nAfrica   ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic America"), "Endemic\nAmerica", region_final1)) %>%
  mutate(position = 1:n()) %>%
  group_by(region_final1) %>%
  mutate(num = n()) %>%
  ungroup() %>%
  group_by(region_final1) %>%
  mutate(mid_num = (min(position)+max(position))/2) %>%
  mutate(position1 = ifelse( (position <= mid_num + (num-3)/2) & (position >= mid_num - (num-3)/2), NA, position))
deme2_pos_min <- sort((deme2 %>% group_by(region_final1) %>% summarise(min(position1, na.rm = T)) %>% select(2))[[1]])
deme2_pos_max <- sort((deme2 %>% group_by(region_final1) %>% summarise(max(position1, na.rm = T)) %>% select(2))[[1]])

deme3 <- data.frame(deme = levels) %>% filter(deme %in% unique(MJ$startLocation[MJ$scheme == "phylogeo_sampling_infection" & MJ$serotype == "DENV3"])) %>%
  left_join(map, by = c("deme" = "geo1")) %>%
  mutate(region_final1 = ifelse(deme %in% c("Florida", "NorthArgentinaUruguay","SouthEurope", map$geo1[map$region_final == "South China"], "Queensland"),
                                "Fringe  \nareas   ", region_final)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Asia"), "Endemic\nAsia    ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Africa"), "Endemic\nAfrica   ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic America"), "Endemic\nAmerica", region_final1)) %>%
  mutate(position = 1:n()) %>%
  group_by(region_final1) %>%
  mutate(num = n()) %>%
  ungroup() %>%
  group_by(region_final1) %>%
  mutate(mid_num = (min(position)+max(position))/2) %>%
  mutate(position1 = ifelse( (position <= mid_num + (num-3)/2) & (position >= mid_num - (num-3)/2), NA, position))
deme3_pos_min <- sort((deme3 %>% group_by(region_final1) %>% summarise(min(position1, na.rm = T)) %>% select(2))[[1]])
deme3_pos_max <- sort((deme3 %>% group_by(region_final1) %>% summarise(max(position1, na.rm = T)) %>% select(2))[[1]])

deme4 <- data.frame(deme = levels) %>% filter(deme %in% unique(MJ$startLocation[MJ$scheme == "phylogeo_sampling_infection" & MJ$serotype == "DENV4"])) %>%
  left_join(map, by = c("deme" = "geo1")) %>%
  mutate(region_final1 = ifelse(deme %in% c("Florida", "NorthArgentinaUruguay","SouthEurope", map$geo1[map$region_final == "South China"], "Queensland"),
                                "Fringe  \nareas   ", region_final)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Asia"), "Endemic\nAsia    ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic Africa"), "Endemic\nAfrica   ", region_final1)) %>%
  mutate(region_final1 = ifelse(region_final1 %in% c("Endemic America"), "Endemic\nAmerica", region_final1)) %>%
  mutate(position = 1:n()) %>%
  group_by(region_final1) %>%
  mutate(num = n()) %>%
  ungroup() %>%
  group_by(region_final1) %>%
  mutate(mid_num = (min(position)+max(position))/2) %>%
  mutate(position1 = ifelse( (position <= mid_num + (num-3)/2) & (position >= mid_num - (num-3)/2), NA, position))
deme4_pos_min <- sort((deme4 %>% group_by(region_final1) %>% summarise(min(position1, na.rm = T)) %>% select(2))[[1]])
deme4_pos_max <- sort((deme4 %>% group_by(region_final1) %>% summarise(max(position1, na.rm = T)) %>% select(2))[[1]])

#==4. plot==
ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(period == "epoch1"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme1$deme[!is.na(deme1$position1)], labels = NULL, expand = c(0.02,0))+
  scale_y_discrete("DENV1\n\nFrom\n\n\n\n", breaks = deme1$deme[!is.na(deme1$position1)], labels = NULL, expand = c(0.02,0))+
  coord_cartesian(clip = "off")+
  annotate("text", x= -Inf, y = unique(deme1$mid_num), label = unique(deme1$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme1$mid_num), y = -Inf, label = unique(deme1$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme1_pos_min,
           yend = deme1_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme1_pos_min,
           xend = deme1_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme1_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme1_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme1_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme1_pos_max, color = color, size = 0.6)+
  theme_bw()+
  # guides(fill = guide_colorbar(frame.colour = "black",
  #                              ticks.colour = "black"))+
  labs(title = "Pre-pandemic period")+
  geom_vline(xintercept = deme1_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme1_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    axis.title = element_text(face = "bold"),
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p1

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(period == "epoch2"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme1$deme[!is.na(deme1$position1)], labels = NULL, expand = c(0.02,0))+
  scale_y_discrete("DENV1\n\nFrom\n\n\n\n", breaks = deme1$deme[!is.na(deme1$position1)], labels = NULL, expand = c(0.02,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme1$mid_num), label = unique(deme1$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme1$mid_num), y = -Inf, label = unique(deme1$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme1_pos_min,
           yend = deme1_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme1_pos_min,
           xend = deme1_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme1_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme1_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme1_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme1_pos_max, color = color, size = 0.6)+
  theme_bw()+
  labs(title = "Pandemic period")+
  geom_vline(xintercept = deme1_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme1_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p2

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV1") %>% filter(period == "epoch3"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme1$deme[!is.na(deme1$position1)], labels = NULL, expand = c(0.02,0))+
  scale_y_discrete("DENV1\n\nFrom\n\n\n\n", breaks = deme1$deme[!is.na(deme1$position1)], labels = NULL, expand = c(0.02,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme1$mid_num), label = unique(deme1$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme1$mid_num), y = -Inf, label = unique(deme1$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme1_pos_min,
           yend = deme1_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme1_pos_min,
           xend = deme1_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme1_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme1_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme1_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme1_pos_max, color = color, size = 0.6)+
  theme_bw()+
  labs(title = "Post-pandemic period")+
  geom_vline(xintercept = deme1_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme1_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p3

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2") %>% filter(period == "epoch1"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme2$deme[!is.na(deme2$position1)], labels = NULL, expand = c(0.02,0))+
  scale_y_discrete("DENV2\n\nFrom\n\n\n\n", breaks = deme2$deme[!is.na(deme2$position1)], labels = NULL, expand = c(0.02,0))+
  coord_cartesian(clip = "off")+
  annotate("text", x= -Inf, y = unique(deme2$mid_num), label = unique(deme2$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme2$mid_num), y = -Inf, label = unique(deme2$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme2_pos_min,
           yend = deme2_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme2_pos_min,
           xend = deme2_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme2_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme2_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme2_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme2_pos_max, color = color, size = 0.6)+
  theme_bw()+
  # guides(fill = guide_colorbar(frame.colour = "black",
  #                              ticks.colour = "black"))+
  geom_vline(xintercept = deme2_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme2_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    axis.title = element_text(face = "bold"),
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p4

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2") %>% filter(period == "epoch2"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme2$deme[!is.na(deme2$position1)], labels = NULL, expand = c(0.02,0))+
  scale_y_discrete("DENV2\n\nFrom\n\n\n\n", breaks = deme2$deme[!is.na(deme2$position1)], labels = NULL, expand = c(0.02,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme2$mid_num), label = unique(deme2$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme2$mid_num), y = -Inf, label = unique(deme2$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme2_pos_min,
           yend = deme2_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme2_pos_min,
           xend = deme2_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme2_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme2_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme2_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme2_pos_max, color = color, size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme2_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme2_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p5

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV2") %>% filter(period == "epoch3"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Intensity of Markov jump per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme2$deme[!is.na(deme2$position1)], labels = NULL, expand = c(0.02,0))+
  scale_y_discrete("DENV2\n\nFrom\n\n\n\n", breaks = deme2$deme[!is.na(deme2$position1)], labels = NULL, expand = c(0.02,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme2$mid_num), label = unique(deme2$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme2$mid_num), y = -Inf, label = unique(deme2$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme2_pos_min,
           yend = deme2_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme2_pos_min,
           xend = deme2_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme2_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme2_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme2_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme2_pos_max, color = color, size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme2_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme2_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  guides(fill = guide_legend(nrow = 1)) +
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p6

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3") %>% filter(period == "epoch1"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme3$deme[!is.na(deme3$position1)], labels = NULL, expand = c(0.025,0))+
  scale_y_discrete("DENV3\n\nFrom\n\n\n\n", breaks = deme3$deme[!is.na(deme3$position1)], labels = NULL, expand = c(0.025,0))+
  coord_cartesian(clip = "off")+
  annotate("text", x= -Inf, y = unique(deme3$mid_num), label = unique(deme3$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme3$mid_num), y = -Inf, label = unique(deme3$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme3_pos_min,
           yend = deme3_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme3_pos_min,
           xend = deme3_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme3_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme3_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme3_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme3_pos_max, color = color, size = 0.6)+
  theme_bw()+
  # guides(fill = guide_colorbar(frame.colour = "black",
  #                              ticks.colour = "black"))+
  geom_vline(xintercept = deme3_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme3_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    axis.title = element_text(face = "bold"),
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p7

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3") %>% filter(period == "epoch2"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme3$deme[!is.na(deme3$position1)], labels = NULL, expand = c(0.025,0))+
  scale_y_discrete("DENV3\n\nFrom\n\n\n\n", breaks = deme3$deme[!is.na(deme3$position1)], labels = NULL, expand = c(0.025,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme3$mid_num), label = unique(deme3$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme3$mid_num), y = -Inf, label = unique(deme3$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme3_pos_min,
           yend = deme3_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme3_pos_min,
           xend = deme3_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme3_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme3_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme3_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme3_pos_max, color = color, size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme3_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme3_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p8

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV3") %>% filter(period == "epoch3"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\n\n\nTo", breaks = deme3$deme[!is.na(deme3$position1)], labels = NULL, expand = c(0.025,0))+
  scale_y_discrete("DENV3\n\nFrom\n\n\n\n", breaks = deme3$deme[!is.na(deme3$position1)], labels = NULL, expand = c(0.025,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme3$mid_num), label = unique(deme3$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  # annotate("text", x= unique(deme3$mid_num), y = -Inf, label = unique(deme3$region_final1), hjust = 1, vjust = 1, angle = 45, size = 2, color = color)+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme3_pos_min,
           yend = deme3_pos_max,
           linewidth = 0.3,
           color = color)+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme3_pos_min,
           xend = deme3_pos_max,
           linewidth = 0.3,
           color = color)+
  # annotate("point", x = -Inf, y =  deme3_pos_min, color = color, size = 0.6)+
  # annotate("point", x = -Inf, y =  deme3_pos_max, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme3_pos_min, color = color, size = 0.6)+
  # annotate("point", y = -Inf, x =  deme3_pos_max, color = color, size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme3_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme3_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p9

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4") %>% filter(period == "epoch1"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\nTo", breaks = deme4$deme[!is.na(deme4$position1)], labels = NULL, expand = c(0.035,0))+
  scale_y_discrete("DENV4\n\nFrom\n\n\n\n", breaks = deme4$deme[!is.na(deme4$position1)], labels = NULL, expand = c(0.035,0))+
  coord_cartesian(clip = "off")+
  annotate("text", x= -Inf, y = unique(deme4$mid_num), label = unique(deme4$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color[c(1,2,4)])+
  annotate("text", x= unique(deme4$mid_num), y = -Inf, label = str_remove_all(unique(deme4$region_final1)," "), hjust = 0.5, vjust = 1.3, size = 2, color = color[c(1,2,4)])+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme4_pos_min,
           yend = deme4_pos_max,
           linewidth = 0.3,
           color = color[c(1,2,4)])+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme4_pos_min,
           xend = deme4_pos_max,
           linewidth = 0.3,
           color = color[c(1,2,4)])+
  # annotate("point", x = -Inf, y =  deme4_pos_min, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", x = -Inf, y =  deme4_pos_max, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", y = -Inf, x =  deme4_pos_min, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", y = -Inf, x =  deme4_pos_max, color = color[c(1,2,4)], size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme4_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme4_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    axis.title = element_text(face = "bold"),
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p10

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4") %>% filter(period == "epoch2"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\nTo", breaks = deme4$deme[!is.na(deme4$position1)], labels = NULL, expand = c(0.035,0))+
  scale_y_discrete("DENV4\n\nFrom\n\n\n\n", breaks = deme4$deme[!is.na(deme4$position1)], labels = NULL, expand = c(0.035,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme4$mid_num), label = unique(deme4$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  annotate("text", x= unique(deme4$mid_num), y = -Inf, label = str_remove_all(unique(deme4$region_final1)," "), hjust = 0.5, vjust = 1.3, size = 2, color = color[c(1,2,4)])+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme4_pos_min,
           yend = deme4_pos_max,
           linewidth = 0.3,
           color = color[c(1,2,4)])+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme4_pos_min,
           xend = deme4_pos_max,
           linewidth = 0.3,
           color = color[c(1,2,4)])+
  # annotate("point", x = -Inf, y =  deme4_pos_min, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", x = -Inf, y =  deme4_pos_max, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", y = -Inf, x =  deme4_pos_min, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", y = -Inf, x =  deme4_pos_max, color = color[c(1,2,4)], size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme4_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme4_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    axis.title = element_text(face = "bold"),
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold"))-> p11

ggplot()+
  geom_tile(data = MJ %>% filter(scheme == "phylogeo_sampling_infection") %>% filter(serotype == "DENV4") %>% filter(period == "epoch3"),
            aes(x = endLocation, y = startLocation, fill = jump_mean_per_year1))+
  scale_fill_manual("Number of transition event per year",values = col_value,
                    na.value="white")+
  scale_x_discrete("\n\nTo", breaks = deme4$deme[!is.na(deme4$position1)], labels = NULL, expand = c(0.035,0))+
  scale_y_discrete("DENV4\n\nFrom\n\n\n\n", breaks = deme4$deme[!is.na(deme4$position1)], labels = NULL, expand = c(0.035,0))+
  coord_cartesian(clip = "off")+
  # annotate("text", x= -Inf, y = unique(deme4$mid_num), label = unique(deme4$region_final1), hjust = 1.2, vjust = 0.5, size = 2, color = color)+
  annotate("text", x= unique(deme4$mid_num), y = -Inf, label = str_remove_all(unique(deme4$region_final1)," "), hjust = 0.5, vjust = 1.3, size = 2, color = color[c(1,2,4)])+
  annotate("segment", x = -Inf, xend = -Inf, 
           y = deme4_pos_min,
           yend = deme4_pos_max,
           linewidth = 0.3,
           color = color[c(1,2,4)])+
  annotate("segment", y = -Inf, yend = -Inf, 
           x = deme4_pos_min,
           xend = deme4_pos_max,
           linewidth = 0.3,
           color = color[c(1,2,4)])+
  # annotate("point", x = -Inf, y =  deme4_pos_min, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", x = -Inf, y =  deme4_pos_max, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", y = -Inf, x =  deme4_pos_min, color = color[c(1,2,4)], size = 0.6)+
  # annotate("point", y = -Inf, x =  deme4_pos_max, color = color[c(1,2,4)], size = 0.6)+
  theme_bw()+
  geom_vline(xintercept = deme4_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  geom_hline(yintercept = deme4_pos_min[-1] - 0.5, linetype = 2, color = "grey",linewidth = 0.2)+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  theme(
    axis.title = element_text(face = "bold"),
    legend.key.size = unit(0.3,"cm"),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 7),
    text = element_text(size = 7),
    panel.border = element_rect(fill = "transparent", color = "transparent"),
    axis.ticks = element_line(linewidth = 0.3),
    plot.tag = element_text(face = "bold")) -> p12

#==load package==
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
library(coda)
library(sf)
library(RevGadgets)

#==define settings==
color1 <- c(pal_aaas("default", alpha = 1)(10))
color2 <- c(pal_nejm("default", alpha = 1)(8))
Sys.setlocale('LC_TIME', 'C')

date_match <- data.frame(date_week = seq(as.Date("2000-01-06"),as.Date("2024-12-26"),7))
date_match$ISO_YEAR <- isoyear(date_match$date)
date_match$ISO_WEEK <- isoweek(date_match$date)

#==HPD function==
summarise_hpd_lower <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[1])
}

summarise_hpd_upper <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[2])
}

#===========================
#===========DENV1===========
#===========================
#==read data== (Data are not provided in GitHub)
files1 = c(Sys.glob("../data/genome/final_data_202512/china_sampling_all/TL/output_single/*.tsv"))
files2 = c(Sys.glob("../data/genome/final_data_202512/china_sampling_all/TL/output_all/*.tsv"))

#based on single MCC tree
DENV_TL_duration <- c()
DENV_TL_size <- c()
for (i in files1) {
  # i = files1[1]
  DENV_tmp <- read.delim(i) %>%
    mutate(duration = last_seen - first_seen) %>%
    mutate(duration1 = cut(duration, breaks = c(-Inf, 7/365, 60/365, 180/365, 1, Inf),
                           right = T, labels = c("(0, 7]", "(7, 60]", "(60, 180]", "(180, 365]", "> 365"))) %>%
    mutate(ntaxa1 = cut(ntaxa, breaks = c(0, 1, 5, 25, 100, Inf), 
                        right = T, labels = c("Singleton", "(1, 5]", "(5, 25]", "(25, 100]", "> 100"))) %>%
    mutate(import_date = ((tmrca - ptmrca)/2) + ptmrca) %>%
    mutate(import_date1 = decimal2Date(import_date)) %>%
    mutate(detection_lag = first_seen - tmrca) %>%
    mutate(lag_define = str_remove_all(str_split(i, "\\/")[[1]][8], "DENV1_fixroot_tl_keep_height_|DENV2_fixroot_tl_keep_height_|DENV1_fixroot_tl_mean_height_|DENV2_fixroot_tl_mean_height_|year.tsv")) %>%
    mutate(serotype = substr(str_split(i, "\\/")[[1]][8], 1, 5 )) %>%
    mutate(tree_type = substr(str_split(i, "\\/")[[1]][8], 18, 28 ) )
  
  DENV_tmp_duration <- DENV_tmp %>%
    group_by(serotype,tree_type, lag_define, duration1) %>%
    summarise(n = n()) %>%
    mutate(total = sum(n)) %>%
    mutate(prop = n/total)
  DENV_tmp_size <- DENV_tmp %>%
    group_by(serotype, tree_type, lag_define, ntaxa1) %>%
    summarise(n = n()) %>%
    mutate(total = sum(n)) %>%
    mutate(prop = n/total)
  
  DENV_TL_duration <- rbind(DENV_TL_duration, DENV_tmp_duration)
  DENV_TL_size <- rbind(DENV_TL_size, DENV_tmp_size)
}

#based on all posterior trees
DENV_TL_duration_all <- c()
DENV_TL_size_all <- c()
for (i in files2) {
  # i = files2[1]
  DENV_tmp <- read.delim(i) %>%
    mutate(duration = last_seen - first_seen) %>%
    mutate(duration1 = cut(duration, breaks = c(-Inf, 7/365, 60/365, 180/365, 1, Inf),
                           right = T, labels = c("(0, 7]", "(7, 60]", "(60, 180]", "(180, 365]", "> 365"))) %>%
    mutate(ntaxa1 = cut(ntaxa, breaks = c(0, 1, 5, 25, 100, Inf), 
                        right = T, labels = c("Singleton", "(1, 5]", "(5, 25]", "(25, 100]", "> 100"))) %>%
    mutate(import_date = ((tmrca - ptmrca)/2) + ptmrca) %>%
    mutate(import_date1 = decimal2Date(import_date)) %>%
    mutate(detection_lag = first_seen - tmrca) %>%
    mutate(lag_define = str_remove_all(str_split(i, "\\/")[[1]][8], "DENV1_fixroot_tl_all_|DENV2_fixroot_tl_all_|year.tsv")) %>%
    mutate(serotype = substr(str_split(i, "\\/")[[1]][8], 1, 5))
  
  tmp_duration <- data.frame(tree = rep(unique(DENV_tmp$tree), length(unique(DENV_tmp$duration1))),
                        duration1 = rep(unique(DENV_tmp$duration1), each = length(unique(DENV_tmp$tree)))) %>%
    mutate(lag_define = str_remove_all(str_split(i, "\\/")[[1]][8], "DENV1_fixroot_tl_all_|DENV2_fixroot_tl_all_|year.tsv"),
           serotype = substr(str_split(i, "\\/")[[1]][8], 1, 5))
  tmp_size <- data.frame(tree = rep(unique(DENV_tmp$tree), length(unique(DENV_tmp$ntaxa1))),
                         ntaxa1 = rep(unique(DENV_tmp$ntaxa1), each = length(unique(DENV_tmp$tree)))) %>%
    mutate(lag_define = str_remove_all(str_split(i, "\\/")[[1]][8], "DENV1_fixroot_tl_all_|DENV2_fixroot_tl_all_|year.tsv"),
           serotype = substr(str_split(i, "\\/")[[1]][8], 1, 5))
  
  DENV_tmp_duration1 <- DENV_tmp %>%
    group_by(serotype, tree, lag_define, duration1) %>%
    summarise(n = n()) %>%
    right_join(tmp_duration) %>%
    mutate(n = ifelse(is.na(n), 0,  n))
  DENV_tmp_size1 <- DENV_tmp %>%
    group_by(serotype, tree, lag_define, ntaxa1) %>%
    summarise(n = n()) %>%
    right_join(tmp_size) %>%
    mutate(n = ifelse(is.na(n), 0,  n))
  
  DENV_tmp_duration2 <- DENV_tmp_duration1 %>%
    group_by(serotype, lag_define, duration1) %>%
    summarise(duration_mean = mean(n),
              duration_median = median(n),
              duration_lower = summarise_hpd_lower(n),
              duration_upper = summarise_hpd_upper(n),
              duration_min = min(n),
              duration_max = max(n))  
  DENV_tmp_size2 <- DENV_tmp_size1 %>%
    group_by(serotype, lag_define, ntaxa1) %>%
    summarise(size_mean = mean(n),
              size_median = median(n),
              size_lower = summarise_hpd_lower(n),
              size_upper = summarise_hpd_upper(n),
              size_min = min(n),
              size_max = max(n))
  
  DENV_TL_duration_all <- rbind(DENV_TL_duration_all, DENV_tmp_duration2)
  DENV_TL_size_all <- rbind(DENV_TL_size_all, DENV_tmp_size2)
}

tmp1 <- DENV_TL_size_all %>%
  group_by(serotype, lag_define) %>%
  mutate(total = sum(size_mean))
tmp2 <- DENV_TL_duration_all %>%
  group_by(serotype, lag_define) %>%
  mutate(total = sum(duration_mean))

#==plot (a-d)==
ggplot(DENV_TL_duration_all %>% filter(serotype == "DENV1") %>% filter(duration1 == "> 365")) + 
  annotate("rect", xmin = 9/12, xmax = 15/12, ymin = -0.5, ymax = 18,size = 0.1, alpha = 0.3, fill = "lightblue")+
  geom_line(aes(x = as.numeric(lag_define), y = duration_mean), color = "black", width = 0.5, linewidth = 0.1)+
  geom_errorbar(aes(x = as.numeric(lag_define), ymin = duration_lower, ymax = duration_upper), width = 0.05, linewidth = 0.1)+
  geom_point(aes(x = as.numeric(lag_define), y = duration_mean), color = "black", width = 0.5, linewidth = 0.1)+
  scale_y_continuous(breaks = seq(0,18,3), limits = c(-0.5,18), expand = c(0,0))+
  scale_x_continuous(breaks = seq(0.25,2,0.25),
                     labels = seq(3,24,3))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0.1,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV1", x = "Detection lags (months)", y = "No. of TLs with duration > 1 year", tag = "b") -> p1

ggplot(DENV_TL_size_all %>% filter(serotype == "DENV1") %>% filter(ntaxa1 == "> 100")) + 
  annotate("rect", xmin = 9/12, xmax = 15/12, ymin = -0.5, ymax = 18,size = 0.1, alpha = 0.3, fill = "lightblue")+
  geom_line(aes(x = as.numeric(lag_define), y = size_mean), color = "black", width = 0.5, linewidth = 0.1)+
  geom_errorbar(aes(x = as.numeric(lag_define), ymin = size_lower, ymax = size_upper), width = 0.05, linewidth = 0.1)+
  geom_point(aes(x = as.numeric(lag_define), y = size_mean), color = "black", width = 0.5, linewidth = 0.1)+
  scale_y_continuous(breaks = seq(0,18,3), limits = c(-0.5,18), expand = c(0,0))+
  scale_x_continuous(breaks = seq(0.25,2,0.25),
                     labels = seq(3,24,3))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0.1,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV1", x = "Detection lags (months)", y = "No. of TLs with size > 100", tag = "c")-> p2

ggplot(DENV_TL_duration_all %>% filter(serotype == "DENV2") %>% filter(duration1 == "> 365")) + 
  annotate("rect", xmin = 9/12, xmax = 15/12, ymin = -0.5, ymax = 18,size = 0.1, alpha = 0.3, fill = "lightblue")+
  geom_line(aes(x = as.numeric(lag_define), y = duration_mean), color = "black", width = 0.5, linewidth = 0.1)+
  geom_errorbar(aes(x = as.numeric(lag_define), ymin = duration_lower, ymax = duration_upper), width = 0.05, linewidth = 0.1)+
  geom_point(aes(x = as.numeric(lag_define), y = duration_mean), color = "black", width = 0.5, linewidth = 0.1)+
  scale_y_continuous(breaks = seq(0,18,3), limits = c(-0.5,18), expand = c(0,0))+
  scale_x_continuous(breaks = seq(0.25,2,0.25),
                     labels = seq(3,24,3))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0.1,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV2", x = "Detection lags (months)", y = "No. of TLs with duration > 1 year", tag = "d") -> p3

ggplot(DENV_TL_size_all %>% filter(serotype == "DENV2") %>% filter(ntaxa1 == "> 100")) + 
  annotate("rect", xmin = 9/12, xmax = 15/12, ymin = -0.5, ymax = 18,size = 0.1, alpha = 0.3, fill = "lightblue")+
  geom_line(aes(x = as.numeric(lag_define), y = size_mean), color = "black", width = 0.5, linewidth = 0.1)+
  geom_errorbar(aes(x = as.numeric(lag_define), ymin = size_lower, ymax = size_upper), width = 0.05, linewidth = 0.1)+
  geom_point(aes(x = as.numeric(lag_define), y = size_mean), color = "black", width = 0.5, linewidth = 0.1)+
  scale_y_continuous(breaks = seq(0,18,3), limits = c(-0.5,18), expand = c(0,0))+
  scale_x_continuous(breaks = seq(0.25,2,0.25),
                     labels = seq(3,24,3))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0.1,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV2", x = "Detection lags (months)", y = "No. of TLs with size > 100", tag = "e")-> p4

#==plot (e-h)==
ggplot(DENV_TL_duration_all %>% filter(serotype == "DENV1") %>% filter(lag_define %in% c(0.75,1,1.25))) + 
  geom_bar(aes(x = duration1, y = duration_mean, fill = lag_define), stat = "identity", position = "dodge", color = "black", width = 0.6, linewidth = 0.1)+
  geom_errorbar(aes(x =  duration1, ymin = duration_lower, ymax = duration_upper, group = lag_define), width = 0.1, linewidth = 0.1, position = position_dodge(0.6)) +
  scale_y_continuous(limits = c(-10,1200), expand = c(0,0))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.2,"cm"),
        legend.position = c(0.7,0.8),
        legend.title = element_text(margin = margin(b = 2, t = 0)),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV1", x = "Duration of TLs (days)", y = "No. of TLs", tag = "f")+
  scale_fill_manual("Detection lags\n(months)", values = c("#8DD3C7","#FFFFB3", "#BEBADA"),
                    labels = c("0.75" = "9",
                               "1" = "12",
                               "1.25" = "15")) -> p5

ggplot(DENV_TL_size_all %>% filter(serotype == "DENV1") %>% filter(lag_define %in% c(0.75,1,1.25))) + 
  geom_bar(aes(x = ntaxa1, y = size_mean, fill = lag_define), stat = "identity", position = "dodge", color = "black", width = 0.6, linewidth = 0.1)+
  geom_errorbar(aes(x =  ntaxa1, ymin = size_lower, ymax = size_upper, group = lag_define), width = 0.1, linewidth = 0.1, position = position_dodge(0.6)) +
  scale_y_continuous(limits = c(-10,1200), expand = c(0,0))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.2,"cm"),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV1", x = "Size of TLs", y = "No. of TLs", tag = "g")+
  guides(fill = F) +
  scale_fill_manual("Detection lags (months)", values = c("#8DD3C7","#FFFFB3", "#BEBADA"))-> p6

ggplot(DENV_TL_duration_all %>% filter(serotype == "DENV2") %>% filter(lag_define %in% c(0.75,1,1.25))) + 
  geom_bar(aes(x = duration1, y = duration_mean, fill = lag_define), stat = "identity", position = "dodge", color = "black", width = 0.6, linewidth = 0.1)+
  geom_errorbar(aes(x =  duration1, ymin = duration_lower, ymax = duration_upper, group = lag_define), width = 0.1, linewidth = 0.1, position = position_dodge(0.6)) +
  scale_y_continuous(limits = c(-10,1200), expand = c(0,0))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.2,"cm"),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV2", x = "Duration of TLs (days)", y = "No. of TLs", tag = "h")+
  guides(fill = F) +
  scale_fill_manual("Detection lags (months)", values = c("#8DD3C7","#FFFFB3", "#BEBADA")) -> p7

ggplot(DENV_TL_size_all %>% filter(serotype == "DENV2") %>% filter(lag_define %in% c(0.75,1,1.25))) + 
  geom_bar(aes(x = ntaxa1, y = size_mean, fill = lag_define), stat = "identity", position = "dodge", color = "black", width = 0.6, linewidth = 0.1)+
  geom_errorbar(aes(x =  ntaxa1, ymin = size_lower, ymax = size_upper, group = lag_define), width = 0.1, linewidth = 0.1, position = position_dodge(0.6)) +
  scale_y_continuous(limits = c(-10,1200), expand = c(0,0))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.2,"cm"),
        text = element_text(size = 7),
        plot.margin = margin(0,0,0,0, "cm"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))+
  labs(subtitle = "DENV2", x = "Size of TLs", y = "No. of TLs", tag = "i")+
  guides(fill = F) +
  scale_fill_manual("Detection lags (months)", values = c("#8DD3C7","#FFFFB3", "#BEBADA")) -> p8

#==plot (i-j)==
#DENV1
DENV1_TL <- read.delim(files1[15]) %>% #15 mean height MCC tree
  mutate(duration = last_seen - first_seen) %>%
  mutate(duration1 = cut(duration, breaks = c(-Inf, 7/365, 60/365, 180/365, 1, Inf),
                         right = T, labels = c("(0, 7]", "(7, 60]", "(60, 180]", "(180, 365]", "> 365"))) %>%
  mutate(ntaxa1 = cut(ntaxa, breaks = c(0, 1, 5, 25, 100, Inf), 
                      right = T, labels = c("Singleton", "(1, 5]", "(5, 25]", "(25, 100]", "> 100"))) %>%
  mutate(import_date = ((tmrca - ptmrca)/2) + ptmrca) %>%
  mutate(import_date1 = decimal2Date(import_date)) %>%
  mutate(detection_lag = first_seen - tmrca)

DENV1_meta <- read.csv("../data/genome/final_data_202512/china_sampling_all/DENV1/DENV1_china_v3.csv") %>%
  filter(geography == "China") %>%
  mutate(lineage = NA)
for (i in 1:nrow(DENV1_meta)) {
  seq <- DENV1_meta$new_name1[i]
  DENV1_meta$lineage[i] <- DENV1_TL$lineage[which(grepl(seq, DENV1_TL$taxa, fixed = TRUE))]
}

DENV1_TL1 <- DENV1_TL %>% 
  filter(import_date >= 2010) %>% 
  filter(ntaxa > 5) %>%
  filter(duration > 7/365) %>%
  mutate(first_seen = first_seen - 0.5/365,
         last_seen = last_seen + 0.5/365) %>%
  mutate(duration2 = cut(duration, breaks = c(-Inf, 1, Inf),
                         right = T, labels = c("Less than 1 year", "More than 1 year")))

DENV1_TL1  <- DENV1_TL1[order(DENV1_TL1$import_date),]
factor(DENV1_TL1$lineage, levels = DENV1_TL1$lineage) -> DENV1_TL1$lineage1
DENV1_meta <- DENV1_meta %>% filter(lineage %in% unique(DENV1_TL1$lineage))
factor(DENV1_meta$lineage, levels = DENV1_TL1$lineage) -> DENV1_meta$lineage
DENV1_meta <- left_join(DENV1_meta, DENV1_TL1[,16:17] %>% distinct(), by = c("lineage" = "lineage1") )

ggplot(DENV1_TL1)+
  annotate("rect", xmin = 0, xmax = nrow(DENV1_TL1), ymin = 2020, ymax = 2022,size = 0.1, alpha = 0.1, fill = "blue")+
  geom_errorbar(aes(x = lineage1, ymin = ptmrca, ymax = tmrca),linewidth = 1, width = 0, alpha = 0.6, color = "grey80")+
  geom_errorbar(aes(x = lineage1, ymin = import_date, ymax = first_seen),linewidth = 0.2, width = 0, linetype = "dotted", color = "grey50")+
  geom_errorbar(aes(x = lineage1, ymin = first_seen, ymax = last_seen, color = duration2), width = 0)+
  geom_point(DENV1_meta, mapping = aes(x = lineage, y = decimal_date(as.Date(Collection.date1)), color = duration2), alpha = 0.2) +
  coord_flip()+
  theme_bw()+
  scale_y_continuous("Date", breaks = seq(2011,2025,2), labels = seq(2011,2025,2), expand = c(0.02,0), limits = c(2011,2025))+
  scale_x_discrete(expand = c(0.015,0))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = "transparent", color =  "black",linewidth = 0.1),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.border =  element_blank(),
        panel.grid  = element_blank(),
        # panel.grid.minor.y  = element_blank(),
        legend.position = c(0.14,0.6),
        legend.key.size = unit(0.3,"cm"),
        plot.margin = margin(0,0.1,0,0.1, "cm"),
        legend.title = element_text(size = 7),
        legend.background = element_rect(fill = "transparent", color =  "transparent"),
        legend.key = element_rect(fill = "transparent", color =  "transparent"),
        axis.line.x = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        text = element_text(size = 7))+
  scale_color_manual("Duration of TLs", values = color2[c(2:1)])+
  annotate("text", y = 2011, x = nrow(DENV1_TL1)-1, label = "DENV1, detection lag set to 1 year", hjust = 0, size = 2.5)+
  # scale_linewidth_discrete("Size of TLs", range = c(0.15,1))+
  labs(tag = "j ")-> p9

#DENV2
DENV2_TL <- read.delim(files1[31]) %>% #31 mean height MCC tree
  mutate(duration = last_seen - first_seen) %>%
  mutate(duration1 = cut(duration, breaks = c(-Inf, 7/365, 60/365, 180/365, 1, Inf),
                         right = T, labels = c("(0, 7]", "(7, 60]", "(60, 180]", "(180, 365]", "> 365"))) %>%
  mutate(ntaxa1 = cut(ntaxa, breaks = c(0, 1, 5, 25, 100, Inf), 
                      right = T, labels = c("Singleton", "(1, 5]", "(5, 25]", "(25, 100]", "> 100"))) %>%
  mutate(import_date = ((tmrca - ptmrca)/2) + ptmrca) %>%
  mutate(import_date1 = decimal2Date(import_date)) %>%
  mutate(detection_lag = first_seen - tmrca)

DENV2_meta <- read.csv("../data/genome/final_data_202512/china_sampling_all/DENV2/DENV2_china_v3.csv") %>%
  filter(geography == "China") %>%
  mutate(lineage = NA)
for (i in 1:nrow(DENV2_meta)) {
  seq <- DENV2_meta$new_name1[i]
  DENV2_meta$lineage[i] <- DENV2_TL$lineage[which(grepl(seq, DENV2_TL$taxa, fixed = TRUE))]
}

DENV2_TL1 <- DENV2_TL %>% 
  filter(import_date >= 2010) %>% 
  filter(ntaxa > 5) %>%
  filter(duration > 7/365) %>%
  mutate(first_seen = first_seen - 0.5/365,
         last_seen = last_seen + 0.5/365) %>%
  mutate(duration2 = cut(duration, breaks = c(-Inf, 1, Inf),
                         right = T, labels = c("Less than 1 year", "More than 1 year")))

DENV2_TL1  <- DENV2_TL1[order(DENV2_TL1$import_date),]
factor(DENV2_TL1$lineage, levels = DENV2_TL1$lineage) -> DENV2_TL1$lineage1
DENV2_meta <- DENV2_meta %>% filter(lineage %in% unique(DENV2_TL1$lineage))
factor(DENV2_meta$lineage, levels = DENV2_TL1$lineage) -> DENV2_meta$lineage
DENV2_meta <- left_join(DENV2_meta, DENV2_TL1[,16:17] %>% distinct(), by = c("lineage" = "lineage1") )

ggplot(DENV2_TL1)+
  annotate("rect", xmin = 0, xmax = nrow(DENV2_TL1), ymin = 2020, ymax = 2022,size = 0.1, alpha = 0.1, fill = "blue")+
  geom_errorbar(aes(x = lineage1, ymin = ptmrca, ymax = tmrca),linewidth = 1, width = 0, alpha = 0.6, color = "grey80")+
  geom_errorbar(aes(x = lineage1, ymin = import_date, ymax = first_seen),linewidth = 0.2, width = 0, linetype = "dotted", color = "grey50")+
  geom_errorbar(aes(x = lineage1, ymin = first_seen, ymax = last_seen, color = duration2), width = 0)+
  geom_point(DENV2_meta, mapping = aes(x = lineage, y = decimal_date(as.Date(Collection.date1)), color = duration2), alpha = 0.2) +
  coord_flip()+
  theme_bw()+
  scale_y_continuous("Date", breaks = seq(2011,2025,2), labels = seq(2011,2025,2), expand = c(0.02,0), limits = c(2011,2025))+
  scale_x_discrete(expand = c(0.015,0))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = "transparent", color =  "black",linewidth = 0.1),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.border =  element_blank(),
        panel.grid  = element_blank(),
        # panel.grid.minor.y  = element_blank(),
        legend.position = c(0.12,0.7),
        legend.key.size = unit(0.3,"cm"),
        plot.margin = margin(0,0.1,0,0.1, "cm"),
        legend.title = element_text(size = 7),
        legend.background = element_rect(fill = "transparent", color =  "transparent"),
        axis.line.x = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        text = element_text(size = 7))+
  scale_color_manual("Duration of TLs", values = color2[c(2:1)])+
  guides(color = F)+
  annotate("text", y = 2011, x = nrow(DENV2_TL1)-0.5, label = "DENV2, detection lag set to 1 year", hjust = 0, size = 2.5)+
  # scale_linewidth_discrete("Size of TLs", range = c(0.15,1))+
  labs(tag = "k ")-> p10

#legend
ggplot() +
  geom_errorbar(aes(x = 1, ymin = 1, ymax = 2 ),linewidth = 2.5,  width = 0, color = "grey",alpha = 0.6,)+
  geom_errorbar(aes(x = 1, ymin = 1.5, ymax = 2.5),linewidth = 0.2,linetype = "dotted",  width = 0, color = "black")+
  geom_errorbar(aes(x = 1, ymin = 2.5, ymax = 5),linewidth = 2.5,  width = 0, color = color2[2])+
  geom_point(aes(x = 1, y = c(2.5,2.6,2.7,4,4.5,4.7,5)), size = 4, color = color2[2], alpha = 0.5)+
  coord_flip(clip = "off")+
  scale_y_continuous(limits = c(1,5))+
  scale_x_continuous(limits = c(0.95,1.05))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 5),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  annotate("text", x = 1.027 ,y = 1, label = "TPMRCA", hjust = 0.5, size = 1.5)+
  annotate("text", x = 1.027 ,y = 2, label = "TMRCA", hjust = 0.5, size = 1.5)+
  annotate("text", x = 1.027 ,y = 3.75, label = "Duration of TLs", hjust = 0.5, size = 1.5)+
  annotate("text", x = 0.973 ,y = 1.5, label = "Importation date", hjust = 0.5, size = 1.5)+
  annotate("text", x = 0.973 ,y = 5, label = "Case of TLs", hjust = 0.5, size = 1.5)+
  geom_segment(aes(x= 1.02, y= 1, xend= 1.007 , yend= 1), arrow = arrow(length=unit(0.1, 'cm'), type = "closed"),lwd= 0.28)+
  geom_segment(aes(x= 1.02, y= 2, xend= 1.007 , yend= 2), arrow = arrow(length=unit(0.1, 'cm'), type = "closed"),lwd= 0.28)+
  geom_segment(aes(x= 1.02, y= 3.75, xend= 1.007 , yend= 3.75), arrow = arrow(length=unit(0.1, 'cm'), type = "closed"),lwd= 0.28)+
  geom_segment(aes(x= 0.98, y= 1.5, xend= 0.993 , yend= 1.5), arrow = arrow(length=unit(0.1, 'cm'), type = "closed"),lwd= 0.28) +
  geom_segment(aes(x= 0.98, y= 5, xend= 0.993 , yend= 5), arrow = arrow(length=unit(0.1, 'cm'), type = "closed"),lwd= 0.28)-> p0

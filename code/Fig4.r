#==load package==
library(stringr)
library(tidyverse)
library(readxl)
library(lubridate)
library(ggpubr)
library(ape)
library(ggtree)
library(phangorn)
library(phytools)
library(treeio)
library(zoo)
library(ggraph)
library(igraph)
library(ggsci)
library(patchwork)
library(scales)
library(showtext)
library(RevGadgets)
library(tracerer)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf"); font.files()
showtext_auto()

#==1. define color==
color <- c("#8DD3C7","#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
color1 <- c(pal_aaas("default", alpha = 1)(10))
arrow_colors <- c("Strong support" ="grey80", "Very strong support" =color1[4], "Decisive support" = color1[6])
col_value1 <- c("Florida" = "#8DD3C7",
                "SouthEurope" = "#00A6FF",
                "EndemicRegion" = "#FB8072",
                "NorthArgentinaUruguay" = "#FCCDE5",
                "SouthChina" =  "#B3DE69",
                "Queensland" = "#80B1D3")

#==2. read data==(Data are not provided in GitHub)
#Markov jump
MJ <- readRDS("../data/genome/final_data_202512/Markov_jump/DENV_Jumps_temspa_no_indices.rds") %>%
  rename("FROM" = "startLocation",
         "TO" = "endLocation") %>%
  select(c(1,2,3,6))

#Bayes factors
dat_DENV <- read.csv("../data/genome/final_data_202512/Bayes_factors/BF_temspa_no_indices.csv") %>%
  select(c(1,3:5)) %>%
  filter(!str_detect(FROM, "20002009")) %>%
  filter(!str_detect(TO, "20002009"))

#==3. Plot setting==
#DENV1
nodes_DENV1 <- data.frame(name = unique(dat_DENV$FROM[dat_DENV$serotype == "DENV1"])) %>%
  mutate(location = gsub("[^A-Za-z]", "", name),
         year = as.integer(gsub("[^0-9]", "", name))) %>%
  left_join(data.frame(year = 2010:2024, x = 1:15)) %>%
  left_join(data.frame(location = c("NorthArgentinaUruguay", "Queensland", "EndemicRegion", "Florida", "SouthEurope", "SouthChina"), y = c(0.8, 1.2, 2, 2.8, 3.2, 3.6))) 

dat_DENV1_2 <- dat_DENV %>%
  filter(serotype == "DENV1") %>%
  mutate(year_from = as.integer(gsub("[^0-9]", "", FROM))) %>%
  mutate(year_to = as.integer(gsub("[^0-9]", "", TO)))%>%
  filter(year_from <= year_to) %>%
  filter(year_to <= (year_from + 2)) %>%
  left_join(MJ %>% filter(serotype == "DENV1")) %>%
  filter(BAYES_FACTOR >= 100) %>%
  mutate(BAYES_FACTOR1 = cut(BAYES_FACTOR, breaks = c(-Inf,3,10,100,1000,Inf),
                             right = T, labels = c("No support", "Support", "Strong support", "Very strong support", "Decisive support")))

edges_DENV1 <- data.frame(from = dat_DENV1_2$FROM,
                          to = dat_DENV1_2$TO,
                          intensity = dat_DENV1_2$BAYES_FACTOR1,
                          weight = dat_DENV1_2$jump_mean) %>%
  mutate(weight1 = as.numeric(as.character(cut(weight, breaks = c(0,5,10,25,50,Inf),
                                               right = T, labels = c(1:5)))))

graph_DENV1 <- graph_from_data_frame(edges_DENV1, directed = TRUE, vertices = nodes_DENV1)
layout_DENV1 <- create_layout(graph_DENV1, layout = "manual", x = nodes_DENV1$x, y = nodes_DENV1$y)

#DENV2
nodes_DENV2 <- data.frame(name = unique(dat_DENV$FROM[dat_DENV$serotype == "DENV2"])) %>%
  mutate(location = gsub("[^A-Za-z]", "", name),
         year = as.integer(gsub("[^0-9]", "", name))) %>%
  left_join(data.frame(year = 2010:2024, x = 1:15)) %>%
  left_join(data.frame(location = c("NorthArgentinaUruguay", "Queensland", "EndemicRegion", "Florida", "SouthEurope", "SouthChina"), y = c(0.8, 1.2, 2, 2.8, 3.2, 3.6))) 

dat_DENV2_2 <- dat_DENV %>%
  filter(serotype == "DENV2") %>%
  mutate(year_from = as.integer(gsub("[^0-9]", "", FROM))) %>%
  mutate(year_to = as.integer(gsub("[^0-9]", "", TO)))%>%
  filter(year_from <= year_to) %>%
  filter(year_to <= (year_from + 2)) %>%
  left_join(MJ %>% filter(serotype == "DENV2")) %>%
  filter(BAYES_FACTOR >= 100) %>%
  mutate(BAYES_FACTOR1 = cut(BAYES_FACTOR, breaks = c(-Inf,3,10,100,1000,Inf),
                             right = T, labels = c("No support", "Support", "Strong support", "Very strong support", "Decisive support")))

edges_DENV2 <- data.frame(from = dat_DENV2_2$FROM,
                          to = dat_DENV2_2$TO,
                          intensity = dat_DENV2_2$BAYES_FACTOR1,
                          weight = dat_DENV2_2$jump_mean) %>%
  mutate(weight1 = as.numeric(as.character(cut(weight, breaks = c(0,5,10,25,50,Inf),
                                               right = T, labels = c(1:5)))))

graph_DENV2 <- graph_from_data_frame(edges_DENV2, directed = TRUE, vertices = nodes_DENV2)

layout_DENV2 <- create_layout(graph_DENV2, layout = "manual", x = nodes_DENV2$x, y = nodes_DENV2$y)

#DENV3
nodes_DENV3 <- data.frame(name = unique(dat_DENV$FROM[dat_DENV$serotype == "DENV3"])) %>%
  mutate(location = gsub("[^A-Za-z]", "", name),
         year = as.integer(gsub("[^0-9]", "", name))) %>%
  left_join(data.frame(year = 2010:2024, x = 1:15)) %>%
  left_join(data.frame(location = c("NorthArgentinaUruguay", "Queensland", "EndemicRegion", "Florida", "SouthEurope", "SouthChina"), y = c(0.8, 1.2, 2, 2.8, 3.2, 3.6))) 

dat_DENV3_2 <- dat_DENV %>%
  filter(serotype == "DENV3") %>%
  mutate(year_from = as.integer(gsub("[^0-9]", "", FROM))) %>%
  mutate(year_to = as.integer(gsub("[^0-9]", "", TO)))%>%
  filter(year_from <= year_to) %>%
  filter(year_to <= (year_from + 2)) %>%
  left_join(MJ %>% filter(serotype == "DENV3")) %>%
  filter(BAYES_FACTOR >= 100) %>%
  mutate(BAYES_FACTOR1 = cut(BAYES_FACTOR, breaks = c(-Inf,3,10,100,1000,Inf),
                             right = T, labels = c("No support", "Support", "Strong support", "Very strong support", "Decisive support")))

edges_DENV3 <- data.frame(from = dat_DENV3_2$FROM,
                          to = dat_DENV3_2$TO,
                          intensity = dat_DENV3_2$BAYES_FACTOR1,
                          weight = dat_DENV3_2$jump_mean) %>%
  mutate(weight1 = as.numeric(as.character(cut(weight, breaks = c(0,5,10,25,50,Inf),
                                               right = T, labels = c(1:5)))))

graph_DENV3 <- graph_from_data_frame(edges_DENV3, directed = TRUE, vertices = nodes_DENV3)

layout_DENV3 <- create_layout(graph_DENV3, layout = "manual", x = nodes_DENV3$x, y = nodes_DENV3$y)

#DENV4
nodes_DENV4 <- data.frame(name = unique(dat_DENV$FROM[dat_DENV$serotype == "DENV4"])) %>%
  mutate(location = gsub("[^A-Za-z]", "", name),
         year = as.integer(gsub("[^0-9]", "", name))) %>%
  left_join(data.frame(year = 2010:2024, x = 1:15)) %>%
  left_join(data.frame(location = c("NorthArgentinaUruguay", "Queensland", "EndemicRegion", "Florida", "SouthEurope", "SouthChina"), y = c(0.8, 1.2, 2, 2.8, 3.2, 3.6))) 

dat_DENV4_2 <- dat_DENV %>%
  filter(serotype == "DENV4") %>%
  mutate(year_from = as.integer(gsub("[^0-9]", "", FROM))) %>%
  mutate(year_to = as.integer(gsub("[^0-9]", "", TO)))%>%
  filter(year_from <= year_to) %>%
  filter(year_to <= (year_from + 2)) %>%
  left_join(MJ %>% filter(serotype == "DENV4")) %>%
  filter(BAYES_FACTOR >= 100) %>%
  mutate(BAYES_FACTOR1 = cut(BAYES_FACTOR, breaks = c(-Inf,3,10,100,1000,Inf),
                             right = T, labels = c("No support", "Support", "Strong support", "Very strong support", "Decisive support")))

edges_DENV4 <- data.frame(from = dat_DENV4_2$FROM,
                          to = dat_DENV4_2$TO,
                          intensity = dat_DENV4_2$BAYES_FACTOR1,
                          weight = dat_DENV4_2$jump_mean) %>%
  mutate(weight1 = as.numeric(as.character(cut(weight, breaks = c(0,5,10,25,50,Inf),
                                               right = T, labels = c(1:5)))))

graph_DENV4 <- graph_from_data_frame(edges_DENV4, directed = TRUE, vertices = nodes_DENV4)

layout_DENV4 <- create_layout(graph_DENV4, layout = "manual", x = nodes_DENV4$x, y = nodes_DENV4$y)

#==4. plot==
min_v <- min(edges_DENV1$weight1,edges_DENV2$weight1,edges_DENV3$weight1,edges_DENV4$weight1)
max_v <- max(edges_DENV1$weight1,edges_DENV2$weight1,edges_DENV3$weight1,edges_DENV4$weight1)

ggraph(layout_DENV1) +
  # geom_hline(yintercept = 1:4, color = "gray90", linewidth = 0.3) +
  # geom_vline(xintercept = 1:15, color = "gray90", linewidth = 0.3) +
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 2.6, ymax = 3.8, fill = "red", alpha = 0.1)+
  annotate("text", x = 12, y = 3, label = "Northern\nHemisphere", size = 3, hjust = 0.5)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 1.5, ymax = 2.5, fill = "lightblue", alpha = 0.3)+
  annotate("text", x = 12, y = 1.65, label = "Tropical areas", size = 3, hjust = 0.5)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 0.6, ymax = 1.4, fill = "red", alpha = 0.1)+
  annotate("text", x = 12, y = 1, label = "Southern\nHemisphere", size = 3, hjust = 0.5)+
  geom_node_point(aes(fill = location), size = 4.5, color = "black", shape = 21)+
  geom_edge_arc2(aes(color = intensity, width = weight1),
                 arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 30),
                 end_cap = circle(2, "mm"),
                 start_cap = circle(2, "mm"),
                 alpha = 0.8,
                 # show.legend = TRUE,
                 strength = 0.25)+
  # geom_node_text(aes(label = gsub("_.*", "", name)), 
  #                size = 3,
  #                nudge_y = -0.15,
  #                color = "grey", 
  #                check_overlap = TRUE)+
  theme_graph()+
  scale_edge_color_manual(name = "Bayes factors", values = arrow_colors, breaks = c("Very strong support", "Decisive support")) +
  scale_fill_manual(name = "Locations",
                    values = col_value1,
                    labels = c("Endemic Regions","Florida, US","Uruguay/North Argentina","Queensland, Australia","South China","South Europe"))+
  scale_x_continuous(breaks = 1:15,
                     labels = 2010:2024,
                     name = "Year",
                     limits = c(0.5, 15.5),
                     expand = expansion(0)) +
  scale_y_continuous(labels = NULL,
                     limits = c(0.6, 3.8),
                     expand = expansion(0))+
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        legend.title.position  = "top",
        plot.margin = margin(0, 0, 0, 5),
        text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial"),
        axis.text = element_text(size = 9))+
  scale_edge_width_continuous("Intensity of Markov jump", range = c(0.01,1.5), limits = c(min_v, max_v), breaks = c(1:5),
                              label = c("(0, 5]", "(5, 10]", "(10, 25]", "(25, 50]", "> 50"))+
  guides(fill = guide_legend(nrow = 2, order = 1),
         edge_color = guide_legend(nrow = 2, order = 2),
         edge_width = guide_legend(nrow = 2, order = 3))+
  labs(y = "DENV1") -> p1

ggraph(layout_DENV2) +
  # geom_hline(yintercept = 1:4, color = "gray90", linewidth = 0.3) +
  # geom_vline(xintercept = 1:15, color = "gray90", linewidth = 0.3) +
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 2.6, ymax = 3.8, fill = "red", alpha = 0.1)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 1.5, ymax = 2.5, fill = "lightblue", alpha = 0.3)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 0.6, ymax = 1.4, fill = "red", alpha = 0.1)+
  geom_node_point(aes(fill = location), size = 4.5, color = "black", shape = 21)+
  geom_edge_arc2(aes(color = intensity, width = weight1),
                 arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 30),
                 end_cap = circle(2, "mm"),
                 start_cap = circle(2, "mm"),
                 alpha = 0.8,
                 # show.legend = TRUE,
                 strength = 0.25) +
  # geom_node_text(aes(label = gsub("_.*", "", name)), 
  #                size = 3,
  #                nudge_y = -0.15,
  #                color = "grey", 
  #                check_overlap = TRUE)+
  theme_graph()+
  scale_edge_color_manual(name = "Bayes factors", values = arrow_colors, breaks = c("Very strong support", "Decisive support")) +
  scale_fill_manual(name = "Locations",
                    values = col_value1,
                    labels = c("Endemic Regions","Florida, US","Queensland, Australia","South China","South Europe"))+
  scale_x_continuous(breaks = 1:15,
                     labels = 2010:2024,
                     name = "Year",
                     limits = c(0.5, 15.5),
                     expand = expansion(0)) +
  scale_y_continuous(labels = NULL,
                     limits = c(0.6, 3.8),
                     expand = expansion(0))+
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        legend.title.position  = "top",
        plot.margin = margin(0, 0, 0, 5),
        text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial"),
        axis.text = element_text(size = 9)) +
  scale_edge_width_continuous("Intensity of Markov jump", range = c(0.01,1.5), limits = c(min_v, max_v), breaks = c(1:5),
                              label = c("(0, 5]", "(5, 10]", "(10, 25]", "(25, 50]", "> 50"))+
  guides(fill = guide_legend(nrow = 2, order = 1),
         edge_color = guide_legend(nrow = 2, order = 2),
         edge_width = guide_legend(nrow = 2, order = 3))+
  labs(y = "DENV2")-> p2

ggraph(layout_DENV3) +
  # geom_hline(yintercept = 1:4, color = "gray90", linewidth = 0.3) +
  # geom_vline(xintercept = 1:15, color = "gray90", linewidth = 0.3) +
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 2.6, ymax = 3.8, fill = "red", alpha = 0.1)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 1.5, ymax = 2.5, fill = "lightblue", alpha = 0.3)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 0.6, ymax = 1.4, fill = "red", alpha = 0.1)+
  geom_node_point(aes(fill = location), size = 4.5, color = "black", shape = 21)+
  geom_edge_arc2(aes(color = intensity, width = weight1),
                 arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 30),
                 end_cap = circle(2, "mm"),
                 start_cap = circle(2, "mm"),
                 alpha = 0.8,
                 # show.legend = TRUE,
                 strength = 0.25) +
  # geom_node_text(aes(label = gsub("_.*", "", name)), 
  #                size = 3,
  #                nudge_y = -0.15,
  #                color = "grey", 
  #                check_overlap = TRUE)+
  theme_graph()+
  scale_edge_color_manual(name = "Bayes factors", values = arrow_colors, breaks = c("Very strong support", "Decisive support"))+
  scale_fill_manual(name = "Locations",
                    values = col_value1,
                    labels = c("Endemic Regions","Florida, US","Queensland, Australia","South China","South Europe"))+
  scale_x_continuous(breaks = 1:15,
                     labels = 2010:2024,
                     name = "Year",
                     limits = c(0.5, 15.5),
                     expand = expansion(0)) +
  scale_y_continuous(labels = NULL,
                     limits = c(0.6, 3.8),
                     expand = expansion(0))+
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 5),
        legend.title.position  = "top",
        text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial"),
        axis.text = element_text(size = 9)) +
  scale_edge_width_continuous("Intensity of Markov jump", range = c(0.01,1.5), limits = c(min_v, max_v), breaks = c(1:5),
                              label = c("(0, 5]", "(5, 10]", "(10, 25]", "(25, 50]", "> 50"))+
  guides(fill = guide_legend(nrow = 2, order = 1),
         edge_color = guide_legend(nrow = 2, order = 2),
         edge_width = guide_legend(nrow = 2, order = 3))+
  labs(y = "DENV3")-> p3

ggraph(layout_DENV4) +
  # geom_hline(yintercept = 1:4, color = "gray90", linewidth = 0.3) +
  # geom_vline(xintercept = 1:15, color = "gray90", linewidth = 0.3) +
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 2.6, ymax = 3.8, fill = "red", alpha = 0.1)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 1.5, ymax = 2.5, fill = "lightblue", alpha = 0.3)+
  annotate("rect", xmin = 0.5, xmax = 15.5, ymin = 0.6, ymax = 1.4, fill = "red", alpha = 0.1)+
  geom_node_point(aes(fill = location), size = 4.5, color = "black", shape = 21)+
  geom_edge_arc2(aes(color = intensity, width = weight1),
                 arrow = arrow(type = "closed", length = unit(2, "mm"), angle = 30),
                 end_cap = circle(2, "mm"),
                 start_cap = circle(2, "mm"),
                 alpha = 0.8,
                 # show.legend = TRUE,
                 strength = 0.25) +
  # geom_node_text(aes(label = gsub("_.*", "", name)), 
  #                size = 3,
  #                nudge_y = -0.15,
  #                color = "grey", 
  #                check_overlap = TRUE)+
  theme_graph()+
  scale_edge_color_manual(name = "Bayes factors",
                          values = arrow_colors, breaks = c("Very strong support", "Decisive support")) +
  scale_fill_manual(name = "Locations",
                    values = col_value1,
                    labels = c("Endemic Regions","Florida, US","South China"))+
  scale_x_continuous(breaks = 1:15,
                     labels = 2010:2024,
                     name = "Year",
                     limits = c(0.5, 15.5),
                     expand = expansion(0)) +
  scale_y_continuous(labels = NULL,
                     limits = c(0.6, 3.8),
                     expand = expansion(0))+
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 0, 0, 5),
        legend.title.position  = "top",
        text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"),
        legend.title = element_text(family = "Arial"),
        axis.text = element_text(size = 9)) +
  scale_edge_width_continuous("Intensity of Markov jump", range = c(0.01,1.5), limits = c(min_v, max_v), breaks = c(1:5),
                              label = c("(0, 5]", "(5, 10]", "(10, 25]", "(25, 50]", "> 50"))+
  guides(fill = guide_legend(nrow = 2, order = 1),
         edge_color = guide_legend(nrow = 2, order = 2),
         edge_width = guide_legend(nrow = 2, order = 3))+
  labs(y = "DENV4") -> p4

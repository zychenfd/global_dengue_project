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
library(coda)
library(reshape2)

#==0. function==
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

#==1. define color==
color <- c("#8DD3C7","#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
color1 <- c(pal_aaas("default", alpha = 1)(10))
point_colors <- c("No support" = "white",
                  "Support" = color1[5],
                  "Strong support" = color1[1],
                  "Very strong support" = color1[4],
                  "Decisive support" = color1[6])

#==2. read data from BEAST output== (Data are not provided in GitHub)
log_DENV1 <- rbind(readTrace("D:/dengue202512/country_even/RUN_1/DENV1_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]],
                   readTrace("D:/dengue202512/country_even/RUN_2/DENV1_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]])
log_DENV2 <- rbind(readTrace("D:/dengue202512/country_even/RUN_1/DENV2_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]],
                   readTrace("D:/dengue202512/country_even/RUN_2/DENV2_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]])
log_DENV3 <- rbind(readTrace("D:/dengue202512/country_even/RUN_1/DENV3_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]],
                   readTrace("D:/dengue202512/country_even/RUN_2/DENV3_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]])
log_DENV4 <- rbind(readTrace("D:/dengue202512/country_even/RUN_1/DENV4_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]],
                   readTrace("D:/dengue202512/country_even/RUN_2/DENV4_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]],
                   readTrace("D:/dengue202512/country_even/RUN_3/DENV4_country_even_GLM1.geo.rates.log",burnin = 0.1)[[1]])

indicators_DENV1 <- log_DENV1[, grep("geo.coefIndicators", colnames(log_DENV1))]
indicators_DENV2 <- log_DENV2[, grep("geo.coefIndicators", colnames(log_DENV2))]
indicators_DENV3 <- log_DENV3[, grep("geo.coefIndicators", colnames(log_DENV3))]
indicators_DENV4 <- log_DENV4[, grep("geo.coefIndicators", colnames(log_DENV4))]

product_DENV1 <- log_DENV1[, grep("Times", colnames(log_DENV1))]
product_DENV2 <- log_DENV2[, grep("Times", colnames(log_DENV2))]
product_DENV3 <- log_DENV3[, grep("Times", colnames(log_DENV3))]
product_DENV4 <- log_DENV4[, grep("Times", colnames(log_DENV4))]

pip_DENV1 <- colMeans(indicators_DENV1)
pip_DENV2 <- colMeans(indicators_DENV2)
pip_DENV3 <- colMeans(indicators_DENV3)
pip_DENV4 <- colMeans(indicators_DENV4)

#==3. Bayes factors==
prior_odds <- 0.06106908933829369/(1-0.06106908933829369)
bf_DENV1 <- (pip_DENV1 / (1 - pip_DENV1))/prior_odds
bf_DENV2 <- (pip_DENV2 / (1 - pip_DENV2))/prior_odds
bf_DENV3 <- (pip_DENV3 / (1 - pip_DENV3))/prior_odds
bf_DENV4 <- (pip_DENV4 / (1 - pip_DENV4))/prior_odds

# BF Table
results <- rbind(data.frame(Predictor = paste("geo.coefficients", 1:11, sep=""),
                            PIP = pip_DENV1, BF = bf_DENV1, serotype = "DENV1"),
                 data.frame(Predictor = paste("geo.coefficients", 1:11, sep=""),
                            PIP = pip_DENV2, BF = bf_DENV2, serotype = "DENV2"),
                 data.frame(Predictor = paste("geo.coefficients", 1:11, sep=""),
                            PIP = pip_DENV3, BF = bf_DENV3, serotype = "DENV3"),
                 data.frame(Predictor = paste("geo.coefficients", 1:11, sep=""),
                            PIP = pip_DENV4, BF = bf_DENV4, serotype = "DENV4"))

# Product Table
log_table_DENV1 <- c()
for (i in 1:11) {
  tmp <- data.frame(Predictor = paste0("geo.coefficients",i),
                    mean = mean(product_DENV1[,i]),
                    serotype = "DENV1",
                    HPD_low =  summarise_hpd_lower(product_DENV1[,i]),
                    HPD_upp = summarise_hpd_upper(product_DENV1[,i]))
  log_table_DENV1 <- rbind(log_table_DENV1, tmp)
}

log_table_DENV2 <- c()
for (i in 1:11) {
  tmp <- data.frame(Predictor = paste0("geo.coefficients",i),
                    mean = mean(product_DENV2[,i]),
                    serotype = "DENV2",
                    HPD_low =  summarise_hpd_lower(product_DENV2[,i]),
                    HPD_upp = summarise_hpd_upper(product_DENV2[,i]))
  log_table_DENV2 <- rbind(log_table_DENV2, tmp)
}

log_table_DENV3 <- c()
for (i in 1:11) {
  tmp <- data.frame(Predictor = paste0("geo.coefficients",i),
                    mean = mean(product_DENV3[,i]),
                    serotype = "DENV3",
                    HPD_low =  summarise_hpd_lower(product_DENV3[,i]),
                    HPD_upp = summarise_hpd_upper(product_DENV3[,i]))
  log_table_DENV3 <- rbind(log_table_DENV3, tmp)
}

log_table_DENV4 <- c()
for (i in 1:11) {
  tmp <- data.frame(Predictor = paste0("geo.coefficients",i),
                    mean = mean(product_DENV4[,i]),
                    serotype = "DENV4",
                    HPD_low =  summarise_hpd_lower(product_DENV4[,i]),
                    HPD_upp = summarise_hpd_upper(product_DENV4[,i]))
  log_table_DENV4 <- rbind(log_table_DENV4, tmp)
}

#==4. Plot==
label <- c("Air traffic",
           "Distance",
           "Index P (O)",
           "Population size (O)",
           "Population size (D)",
           "Maximum precipitation (O)",
           "Maximum precipitation (D)",
           "Same region",
           "Sharing a border",
           "Mosquito occurrence (O)",
           "Mosquito occurrence (D)")

label_order <- rev(c("Population size (O)",
                     "Population size (D)",
                     "Distance",
                     "Same region",
                     "Sharing a border",
                     "Air traffic",
                     "Mosquito occurrence (O)",
                     "Mosquito occurrence (D)",
                     "Index P (O)",
                     "Maximum precipitation (O)",
                     "Maximum precipitation (D)"))

log_table <- rbind(log_table_DENV1,log_table_DENV2,log_table_DENV3,log_table_DENV4) %>%
  left_join(results) %>%
  mutate(BAYES_FACTOR = cut(BF, breaks = c(-Inf,3,10,100,1000,Inf),
                            right = T, labels = c("No support", "Support", "Strong support", 
                                                  "Very strong support", "Decisive support"))) %>%
  left_join(data.frame(Predictor = paste0("geo.coefficients",seq(1,11,1)),
                       variables = label))

log_table1 <- log_table %>% filter(serotype == "DENV1"); factor(log_table1$variables, levels = label_order) -> log_table1$variables
log_table2 <- log_table %>% filter(serotype == "DENV2"); factor(log_table2$variables, levels = label_order) -> log_table2$variables
log_table3 <- log_table %>% filter(serotype == "DENV3"); factor(log_table3$variables, levels = label_order) -> log_table3$variables
log_table4 <- log_table %>% filter(serotype == "DENV4"); factor(log_table4$variables, levels = label_order) -> log_table4$variables

product_DENV1_long <- melt(product_DENV1) %>% 
  left_join(data.frame(variable = paste0("geo.coefficientsTimesIndicators",seq(1,11,1)), variables = label)) %>%
  filter(variables %in% unique(log_table1$variables))
product_DENV2_long <- melt(product_DENV2) %>% 
  left_join(data.frame(variable = paste0("geoy.coefficientsTimesIndicators",seq(1,11,1)), variables = label)) %>%
  filter(variables %in% unique(log_table2$variables))
product_DENV3_long <- melt(product_DENV3) %>% 
  left_join(data.frame(variable = paste0("geo.coefficientsTimesIndicators",seq(1,11,1)), variables = label)) %>%
  filter(variables %in% unique(log_table3$variables))
product_DENV4_long <- melt(product_DENV4) %>% 
  left_join(data.frame(variable = paste0("geo.coefficientsTimesIndicators",seq(1,11,1)), variables = label)) %>%
  filter(variables %in% unique(log_table4$variables))

min_value <- min(product_DENV1_long$value,product_DENV2_long$value,product_DENV3_long$value,product_DENV4_long$value)
max_value <- max(product_DENV1_long$value,product_DENV2_long$value,product_DENV3_long$value,product_DENV4_long$value)

ggplot(log_table1)+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", linewidth = 0.2)+
  geom_errorbar(aes(x = variables, ymin = HPD_low, ymax = HPD_upp), width = 0.2, linewidth =0.1)+
  geom_point(aes(x = variables, y = mean, fill = BAYES_FACTOR), shape =21, size = 2, stroke = 0.01)+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-max_value, max_value))+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.05),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        text = element_text(size = 7))+
  scale_fill_manual("Bayes factors", values = point_colors)+
  labs( x = "Predictor", y = "Coefficient × Inclusion probability")-> p1

ggplot(log_table2)+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", linewidth = 0.2)+
  geom_errorbar(aes(x = variables, ymin = HPD_low, ymax = HPD_upp), width = 0.2, linewidth =0.1)+
  geom_point(aes(x = variables, y = mean, fill = BAYES_FACTOR), shape =21, size = 2, stroke = 0.01)+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-max_value, max_value))+
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.05),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 7))+
  scale_fill_manual("Bayes factors", values = point_colors)+
  labs(x = "Predictor", y = "Coefficient × Inclusion probability") -> p2

ggplot(log_table3)+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", linewidth = 0.2)+
  geom_errorbar(aes(x = variables, ymin = HPD_low, ymax = HPD_upp), width = 0.2, linewidth =0.1)+
  geom_point(aes(x = variables, y = mean, fill = BAYES_FACTOR), shape =21, size = 2, stroke = 0.01)+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-max_value, max_value))+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.05),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        text = element_text(size = 7))+
  scale_fill_manual("Bayes factors", values = point_colors)+
  labs(x = "Predictor", y = "Coefficient × Inclusion probability")-> p3

ggplot(log_table4 )+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", linewidth = 0.2)+
  geom_errorbar(aes(x = variables, ymin = HPD_low, ymax = HPD_upp), width = 0.2, linewidth =0.1)+
  geom_point(aes(x = variables, y = mean, fill = BAYES_FACTOR), shape =21, size = 2, stroke = 0.01)+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-max_value, max_value))+
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        text = element_text(size = 7))+
  scale_fill_manual("Bayes factors", values = point_colors)+
  labs( x = "Predictor", y = "Coefficient × Inclusion probability") -> p4

#Legend plot
ggplot(log_table)+
  geom_point(aes(x = 0, y = 0, fill = BAYES_FACTOR), shape =21, size = 2, stroke = 0.01)+
  scale_fill_manual("Bayes factors", values = point_colors)+
  scale_x_continuous(limits = c(-0.05,-0.025))+
  theme(legend.position = c(0.5,0.5))+
  theme(panel.background = element_rect(fill = "transparent", color =  "transparent"),
        plot.background = element_rect(fill = "transparent", color =  "transparent"),
        legend.background = element_rect(fill = "transparent", color =  "transparent"),
        panel.grid = element_blank(),
        text = element_text(size = 7),
        legend.title.position = "left",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  guides(fill = guide_legend(nrow = 1)) -> p0

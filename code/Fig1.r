#==load packages==
library(stringr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(ggsci)
library(scales)
library(patchwork)
library(rgdal)
library(cowplot) 

#=========================
#==1. define environment==
#=========================
Sys.setlocale('LC_TIME', 'C')
colors <- c(pal_aaas("default", alpha = 1)(10)); show_col(colors)
colors1 <- c(pal_nejm("default", alpha = 1)(10)); show_col(colors1)
colors2 <- c("#8DD3C7","#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
col_value <- c("Florida, US" = "#8DD3C7",
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

#=======================
#==2. define geography==
#=======================
study_area <- read.csv("../data/geo_data/iso3_included.csv") %>% 
  mutate(iso3 = substr(alpha.3,1,3)) %>%
  mutate(region_final = ifelse(str_detect(region_final,"Endemic"), region_final, geo1))
world_region <- read.csv("../data/geo_data/list_country.csv")

#========================
#======3. Epi data======= (Data are not provided in GitHub due to data sharing agreement)
#========================
#==0). Europe data (Italy and France)==
EU1 <- read.csv("../data/epi_data/Europe/ECDC_surveillance_data_Dengue_Local_20251226.csv") %>% 
  filter(RegionName %in% c("France","Italy","Spain")) %>% 
  mutate(NumValue = as.integer(NumValue)) %>%
  group_by(RegionName, Time, Population) %>%
  summarise(cases = sum(NumValue)) %>%
  as.data.frame() 

EU2 <- read.csv("../data/epi_data/Europe/ECDC_surveillance_data_Dengue_Travel_20251226.csv") %>% 
  filter(RegionName %in% c("France","Italy","Spain")) %>% 
  mutate(NumValue = as.integer(NumValue)) %>%
  group_by(RegionName, Time, Population) %>%
  summarise(cases = sum(NumValue)) %>%
  as.data.frame() 

EU_data <- rbind(EU1, EU2) %>%
  filter(!str_detect(Time, "2008|2009")) %>%
  mutate(date_lab = paste(format(as.Date(paste0(Time,"-01")), "%b"), substr(Time, 1, 4))) %>%
  mutate(Population = str_remove_all(Population, " cases")) %>%
  mutate(Population = str_replace_all(Population, "-", " ")) %>%
  rename("TravelStatus" = "Population", "country" = "RegionName") %>%
  left_join(world_region[,c(1,3)], by = c("country" = "name")) %>%
  rename("iso3" = "alpha.3") %>%
  mutate(year_loc = paste0(substr(date_lab, 4,8), iso3)) %>%
  select(-"Time")
  
#==1). US data (includes U.S. Territories and Freely Associated States, by Travel status)==
ArboNET <- read_xlsx("../data/epi_data/US/Dengue cases - ArboNET.xlsx") %>%
  filter(Year <= 2024); sum(ArboNET$DengueCases[ArboNET$Year <= 2024])
ArboNET1 <- ArboNET %>% 
  filter(Jurisdictions == 1) %>% 
  filter(State %in% c("FL")) %>%
  group_by(Year, Month, TravelStatus) %>%
  summarise(DengueCases = sum(DengueCases)) %>%
  mutate(State = "US")

US_data <- rbind(ArboNET1, ArboNET %>% filter(Jurisdictions %in% 2:3)) %>%
  left_join(world_region[,1:3], by = c("State" = "alpha.2")) %>%
  mutate(date_lab = paste(Month, Year),
         year_loc = paste0(Year,alpha.3)) %>%
  rename("country" = "name",
         "iso3" = "alpha.3",
         "cases" = "DengueCases") %>%
  dplyr::select(c("date_lab","year_loc","country","iso3","cases","TravelStatus")) 

#==2). China data (includes local and imported cases)==
#Zhejiang
ZJ_tmp <- read_xlsx("../data/epi_data/China/Zhejiang_2005-2024_1226_lwz.xlsx")
ZJ <- rbind(read_xlsx("../data/epi_data/China/Zhejiang_2005-2020.xlsx"),
            read_xlsx("../data/epi_data/China/Zhejiang_2005-2024.xlsx")[2885:3548,]) %>%
  select(c(12,13,17)) %>%
  mutate(local_case_check = ZJ_tmp$local_case,
         local_case1 = ZJ_tmp$case)
table(ZJ$local_case == ZJ$local_case_check)
ZJ <- ZJ[,c(1,2,5)] #imported cases from other provinces are removed
names(ZJ) <- c("date_onset", "date_diagnosis", "TravelStatus")
ZJ <- ZJ %>% 
  filter(!is.na(TravelStatus)) %>%
  filter(TravelStatus != 2) %>%
  mutate(date_onset = as.Date(date_onset)) %>%
  filter(date_onset >= as.Date("2010-01-01")) %>%
  mutate(date_lab = paste(format(date_onset, "%b"), format(date_onset, "%Y"))) %>%
  group_by(date_lab, TravelStatus) %>%
  summarise(cases = n()) %>%
  mutate(TravelStatus = ifelse(TravelStatus == 0, "Travel associated", "Locally acquired"))

#Guangdong
GD1 <- rbind(read_xlsx("../data/epi_data/China/Guangdong_2005_202507.xlsx", sheet = 2)) %>%
  gather(key = "Year", value = "cases",2:16) %>%
  mutate(cases = ifelse(is.na(cases),0, cases)) %>% 
  mutate(TravelStatus = "Locally acquired")
GD2 <- rbind(read_xlsx("../data/epi_data/China/Guangdong_2005_202507.xlsx", sheet = 3)) %>%
  gather(key = "Year", value = "cases",2:16) %>%
  mutate(cases = ifelse(is.na(cases),0, cases)) %>% 
  mutate(TravelStatus = "Travel associated")
GD <- rbind(GD1, GD2) %>%
  mutate(date_lab = paste(format(as.Date(paste0(Year,"-",Month,"-01")), "%b"),
                          format(as.Date(paste0(Year,"-",Month,"-01")), "%Y"))) %>%
  select(c(3:5))

#Taiwan
TW <- read_excel("../data/epi_data/China/dengue-data_TWN.xlsx")[,1:7] %>%
  mutate(local_cases = all_cases - import_cases) %>%
  filter(!str_detect(date, "2020|2021|2022")) %>%
  filter(!str_detect(date, "2016")) %>%
  select(c(2,6,7)) %>%
  gather(key = "TravelStatus", value = "cases", "local_cases", "import_cases") %>%
  mutate(TravelStatus = ifelse(TravelStatus == "import_cases", "Travel associated", "Locally acquired"))

#Fujian
FJ <- read_excel("../data/epi_data/China/Fujuan_by_onset_date.xlsx") %>%
  filter(Year %in% 2010:2024) %>%
  mutate(date_lab = paste(format(as.Date(paste0(Year,"-",Month,"-01")), "%b"),
                          format(as.Date(paste0(Year,"-",Month,"-01")), "%Y"))) %>%
  mutate(Local_or_import = ifelse(Local_or_import == "Import", "Travel associated", "Locally acquired")) %>%
  select(c(4:6)) 
names(FJ)[1:2] <- c("TravelStatus","cases")

#Yunnan
YN <- read_xlsx("../data/epi_data/China/Yunnan_2010_2025.xlsx") %>%
  mutate(date_lab = paste(format(as.Date(paste0(year,"-",month,"-01")), "%b"),
                          format(as.Date(paste0(year,"-",month,"-01")), "%Y"))) %>%
  select(c(4,5,7)) %>%
  gather(key = "TravelStatus", value = "cases", "local", "import") %>%
  mutate(TravelStatus = ifelse(TravelStatus == "import", "Travel associated", "Locally acquired"))

#South China as a whole
CHN_data <- rbind(ZJ %>% mutate(iso3 = "CHNZhejiang"),
                  GD %>% mutate(iso3 = "CHNGuangdong"),
                  TW %>% mutate(iso3 = "TWN"),
                  FJ %>% mutate(iso3 = "CHNFujian"),
                  YN %>% mutate(iso3 = "CHNYunnan")) %>%
  group_by(iso3, date_lab, TravelStatus) %>%
  summarise(cases = sum(cases)) %>%
  right_join(data.frame(date_lab = format(rep(seq(as.Date("2010-01-01"), as.Date("2024-12-01"), by = "1 month"), 2),"%b %Y"),
                        TravelStatus = rep(c("Travel associated", "Locally acquired"), each = 180)))%>%
  mutate(cases = ifelse(is.na(cases),0,cases)) %>%
  mutate(year_loc = paste0(substr(date_lab, 5,8), "CHN"),
         # iso3 = "CHN",
         country = "China")

#==3). Australia data (cannot distinguish imported cases from local cases, but have that in yearly resolution)==
AUS_data <- read_xlsx("../data/epi_data/Australia/Denguecase-AUSQLD.xlsx") %>% 
  filter(iso3 == "AUS") %>%
  mutate(year_loc = paste0(substr(date_lab,5,8), iso3)) %>% 
  dplyr::select(c("date_lab", "country", "iso3", "cases", "year_loc")) %>%
  mutate(TravelStatus = "Mix")
AUS_data$cases[is.na(AUS_data$cases)] <- 0

AUS_per_year <- read.csv("../data/epi_data/Australia_local_case.csv", header = F) %>%
  mutate(V1 = round(V1, 0)) %>%
  mutate(local_case = round(V2, 0)) %>%
  add_row(V1 = 2021:2022,
          local_case = 0) %>%
  mutate(import_case =  c(1510,1602,1533,1648,
                          2207,1120,928,1452,
                          220,11,528)) %>%
  select(c(1,3,4)) %>%
  gather(key = "TravelStatus", value = "cases", "local_case", "import_case") %>%
  group_by(V1) %>%
  mutate(total_case = sum(cases))

ggplot(AUS_per_year) +
  geom_bar(aes(x = V1, y = cases*100/total_case, fill = TravelStatus), stat = "identity",width = 0.5)+
  theme_bw()+
  scale_x_continuous("Year",breaks = seq(2011,2022,1))+
  scale_y_continuous(breaks = c(0,25,50,75,100))+
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0.03,0.3,0.03,0.3, "cm"),
        text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5))+
  labs(y = "Prop. of DENV cases (%)", tag = "j", subtitle = "Australia")+
  scale_fill_manual(values = colors1[c(3,6)]) -> fig_AUS

#==4). Global cases from WHO (mainly are local cases except for Australia, Réunion, China)==
#NOTE: Negative values were observed in Americas; Data in China was unavailable before 2024
data <- read_xlsx("../data/epi_data/dengue-global-data-2025-12-26.xlsx") %>% 
  mutate(year_loc = paste0(substr(date_lab,5,8), iso3)) %>%
  filter(!iso3 %in% c(unique(US_data$iso3), "CHN", "AUS", "REU", "FRA", "ITA", "ESP")) %>% # without imported data in China in 2024; Australia and Réunion do not distinguish between imported and locally acquired cases
  filter(country != "Autonomous Region of Madeira") %>%
  filter(!is.na(cases)) %>%
  dplyr::select(c("date_lab","country","iso3","cases","year_loc")) # Data reported as of 09 December 2025; local cases in Uruguay 
tmp_data <- data %>% filter(cases <= -10)
data1 <- data %>% filter(!year_loc %in% unique(tmp_data$year_loc)) 
data1$cases[data1$cases < 0] <- 0
data1 <- data1 %>% filter(year_loc != "2015TTO") #Trinidad and Tobago 2015 removed due to potential miss-reporting

#==5). Final data==
data2 <- rbind(data1 %>% mutate(TravelStatus = "Locally acquired"),
               EU_data,
               US_data[,3:8],
               CHN_data,
               AUS_data) %>%
  mutate(date = as.Date(paste(substr(date_lab, 5,8),substr(date_lab, 1,3), "01"), format = c("%Y %b %d"))) %>%
  filter(date < as.Date("2025-01-01")) %>%
  mutate(month = substr(date_lab, 1, 3)) %>%
  filter(iso3 %in% study_area$alpha.3) %>%
  left_join(study_area[(4:6)], by = c("iso3"="alpha.3")) %>%
  mutate(region_final = ifelse(str_detect(iso3, "CHN|TWN"), "South China", region_final))

data2 <- data2 %>%
  group_by(region_final, date_lab, date, month, TravelStatus, iso3) %>%
  summarise(cases = sum(cases)) %>%
  mutate(year_loc = paste0(substr(date,1,4),iso3))

data2 <- data2 %>% filter(year_loc != "2024TWN") #potential miss/bias-reporting in Taiwan 2024

#5.1 local cases
data2_local <- data2 %>% filter(TravelStatus %in% c("Locally acquired")) 
tmp <- data.frame(table(data2_local$year_loc)); table(tmp$Freq)
data3_local <- data2_local %>% 
  filter(year_loc %in% tmp$Var1[tmp$Freq == 12]) %>%
  group_by(year_loc) %>%
  mutate(total_cases = sum(cases)) %>%
  mutate(ratio = cases*100/total_cases) %>%
  filter(total_cases >= 15) #30
tmp1 <- data.frame(table(data3_local$year_loc[data3_local$cases > 0]))
data4_local <- data3_local %>%
  filter(year_loc %in% tmp1$Var1[tmp1$Freq >= 3]) %>% #4
  filter(!str_detect(year_loc, "2020")) %>%
  filter(!str_detect(year_loc, "2021")) %>% 
  filter(!str_detect(year_loc, "2022")) 
factor(data4_local$month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")) -> data4_local$month
data4_local$month1 <- as.numeric(data4_local$month)

data5_local <- data4_local %>% 
  group_by(region_final, month) %>%
  summarise(average_ratio = mean(ratio),
            sd_ratio = sd(ratio),
            median_ratio = median(ratio),
            range_low = quantile(ratio, 0.25),
            range_upp = quantile(ratio, 0.75))

#5.2 imported cases
data2_import <- data2 %>% filter(TravelStatus %in% c("Travel associated", "Mix")) %>%
  filter(region_final %in% c("Florida","NorthArgentinaUruguay", "Queensland", "South China", "SouthEurope"))
tmp <- data.frame(table(data2_import$year_loc)); table(tmp$Freq)
data3_import <- data2_import %>% 
  filter(year_loc %in% tmp$Var1[tmp$Freq == 12]) %>%
  group_by(year_loc) %>%
  mutate(total_cases = sum(cases)) %>%
  mutate(ratio = cases*100/total_cases) %>%
  filter(total_cases >= 15)
tmp1 <- data.frame(table(data3_import$year_loc[data3_import$cases > 0]))
data4_import <- data3_import %>%
  filter(year_loc %in% tmp1$Var1[tmp1$Freq >= 3]) %>%
  filter(!str_detect(year_loc, "2020")) %>%
  filter(!str_detect(year_loc, "2021")) %>%
  filter(!str_detect(year_loc, "2022")) 
factor(data4_import$month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")) -> data4_import$month
data4_import$month1 <- as.numeric(data4_import$month)
table(data4_import$iso3[data4_import$region_final == "SouthEurope"])

data5_import <- data4_import %>% 
  group_by(region_final, month) %>%
  summarise(average_ratio = mean(ratio),
            sd_ratio = sd(ratio),
            median_ratio = median(ratio),
            range_low = quantile(ratio, 0.25),
            range_upp = quantile(ratio, 0.75))

#============================
#=========4. Index P=========
#============================
indexP_dat_other <- readRDS("../output/IndexP_by_month_202512.rds")[,c(1:6,9)] %>%
  mutate(month = as.integer(substr(year_month1, 6, 7)))
indexP_dat <- indexP_dat_other %>%
  filter(!geo1 %in% c("Hainan","Guangxi")) %>%
  left_join(unique(study_area[,c(4,6)])) %>%
  mutate(region_final = ifelse(str_detect(alpha.3, "CHN|TWN"), "South China", region_final))
indexP_dat <- indexP_dat %>%
  group_by(region_final, year_month1, month, alpha.3) %>%
  summarise(index_value = sum(index_value * num_grid)/sum(num_grid)) %>%
  mutate(year_loc = paste0(substr(year_month1,1,4),alpha.3)) %>%
  filter(!region_final %in% c("Japan","WestEurope"))

#===========================
#===5. PLOT: fitted curve===
#===========================
ggplot() +
  geom_smooth(data4_import %>% filter(region_final %in% c("Florida")),
              mapping = aes(x = month1, y = ratio, color = "Travel associated", fill = "Travel associated"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(data4_local %>% filter(region_final %in% c("Florida")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("Florida")),
              mapping = aes(x = month, y = index_value * 7.5, color = "Index P", fill = "Index P"), #10.8
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  # annotate("text", x = 8.475, y = 28, label = "Two-month lag", size = 1.8)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/7.5,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-5,53))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = c(0.3,0.7),
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "Florida, US", x = "",  tag = "b",
       y = "Prop. of monthly cases to\nall cases in each year (%)") -> p1

ggplot() +
  geom_smooth(data4_import %>% filter(region_final %in% c("SouthEurope")),
              mapping = aes(x = month1, y = ratio, color = "Travel associated", fill = "Travel associated"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(data4_local %>% filter(region_final %in% c("SouthEurope")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("SouthEurope")),
              mapping = aes(x = month, y = index_value * 122, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  # annotate("text", x = 8.475, y = 28, label = "Two-month lag", size = 1.8)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/122,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-5,53))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "South Europe", x = "",  tag = "c",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p1_eu

ggplot() +
  geom_smooth(data4_import %>% filter(region_final %in% c("South China")),
              mapping = aes(x = month1, y = ratio, color = "Travel associated", fill = "Travel associated"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(data4_local %>% filter(region_final %in% c("South China")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("South China")),
              mapping = aes(x = month, y = index_value * 23, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/23,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-5,53))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "South China", x = "",  tag = "d",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p2

ggplot() +
  geom_smooth(data4_local %>% filter(region_final %in% c("Endemic America")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("Endemic America")) %>%
                filter(alpha.3 %in% unique(data4_local$iso3[data4_local$region_final == "Endemic America"])),
              mapping = aes(x = month, y = index_value * 8.05, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/8.05,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-2,22))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F, color = F)+
  scale_linetype_manual(values = 2)+
  theme(legend.position = c(0.25,0.8),
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.6, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "Endemic America", x = "",  tag = "e",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p3

ggplot() +
  geom_smooth(data4_local %>% filter(region_final %in% c("Endemic Africa")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("Endemic Africa")) %>%
                filter(alpha.3 %in% unique(data4_local$iso3[data4_local$region_final == "Endemic Africa"])),
              mapping = aes(x = month, y = index_value * 6.7, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/6.7,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-2,22))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "Endemic Africa", x = "",  tag = "f",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p4

ggplot() +
  geom_smooth(data4_local %>% filter(region_final %in% c("Endemic Asia")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("Endemic Asia")) %>%
                filter(alpha.3 %in% unique(data4_local$iso3[data4_local$region_final == "Endemic Asia"])),
              mapping = aes(x = month, y = index_value * 7.95, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/7.95,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-2,22))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "Endemic Asia", x = "",  tag = "g",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p5

ggplot() +
  geom_smooth(data4_local %>% filter(region_final %in% c("NorthArgentinaUruguay")),
              mapping = aes(x = month1, y = ratio, color = "Locally acquired", fill = "Locally acquired"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("NorthArgentinaUruguay")),
              mapping = aes(x = month, y = index_value * 41.5, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/41.5,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-5,53))+
  scale_color_manual(values = colors1[c(1,6,3)])+
  scale_fill_manual(values = colors1[c(1,6,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "Uruguay/North Argentina", x = "Month",  tag = "h",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p6

ggplot() +
  geom_smooth(data4_import %>% filter(region_final %in% c("Queensland")),
              mapping = aes(x = month1, y = ratio, color = "Travel associated", fill = "Travel associated"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  geom_smooth(indexP_dat %>% filter(region_final %in% c("Queensland")),
              mapping = aes(x = month, y = index_value * 8, color = "Index P", fill = "Index P"),
              method = "gam",
              formula = y ~ s(x, bs = "cc", k = 12),
              alpha = 0.2,
              linewidth = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb)+
  scale_y_continuous(expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/8,
                                         name="Index P"))+
  coord_cartesian(ylim = c(-5,53))+
  scale_color_manual(values = colors1[c(1,3)])+
  scale_fill_manual(values = colors1[c(1,3)])+
  theme_bw()+
  guides(fill = F)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 7))+
  labs(subtitle = "Queensland, Australia", x = "Month",  tag = "i",
       y = "Prop. of monthly cases to\nall cases in each year (%)")-> p7

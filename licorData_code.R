# install.packages("ggridges") # you only need to do this once
library(Rmisc) # for summarySE
library(tidyverse) # for readr, dplyr, ggplot...
library(lubridate) # for dealing with dates and times (e.g., hms)
library(ggridges) # for ridge plots
library(gridExtra) # for plotting multiple plots in a single frame

setwd("~/Desktop/andreas18ssequencing/")

# Load the file
d0 <- read_csv("newTimes.csv")
head(d0)
summary(d0)

# make depth a factor, not a continuous variable
d0$depth <- as.factor(d0$depth)

# get rid of spaces in site names
d0$siteName <- gsub(" ", "", d0$siteName)
summary(d0)

# plot licor input over time
d0 %>% ggplot(aes(x = localTime, y = licorData, group = siteType, color = siteType)) +
  geom_line()+
  facet_grid(siteName ~ .)+
  theme_bw()

# Get rid of the first 5 measurements of each depth for each site (time during descent)----
d1 <- d0 %>% group_by(siteName,depth) %>% slice(5:n())
head(d1)

# Gid rid of measurements at Bastimentos South at 6 ft after 12:35 (meter was taken out but not turned off)
endBS6 <- hms("12:35:00")

d2 <- d1 %>%
  filter(!(siteName=="BASTIMSOUTH" & depth==6 & hms(localTime)>endBS6))
head(d2)

# Get rid of any measurements before 11 to 2
dstart <- hms("11:00:00")
dend <- hms("14:00:00")

d3 <- d2 %>%
  filter(cloudy=="N") %>% # try analysis with and without this line. Removing "cloudy" points removes 1198 observations
  filter(!(hms(localTime)<dstart), 
         !(hms(localTime)>dend) )


table(as.factor(d2$cloudy))
# N    Y
# 801 1198
table(as.factor(d3$cloudy))
# N 
# 520

# Times look good?
# before filtering ...
summary(hms(d2$localTime))
# min is 955 am , max is 238 pm

# after filtering...
summary(hms(d3$localTime))
# min time is 11H 0M 4S
# max time is 13H 59M 50S
# looks good!

# Subset for values within N (1? 2?) standard deviations of the mean for each depth and site ----
sd_subset <- d3 %>%
  ungroup() %>% # this line is important! the data were previously grouped, so we must ungroup before grouping again :)
  group_by(siteName,depth) %>%
  mutate(mean = mean(licorData), sd = sd(licorData))  %>%
  filter(licorData > mean-(1*sd), licorData < mean+(1*sd) )

head(sd_subset)
summary(sd_subset)

sd_subset %>% ggplot(aes(x = localTime, y = licorData, group = siteType, color = siteType)) +
  geom_line()+
  facet_grid(siteName ~ .)+
  theme_bw()

sd_subset %>% ggplot(aes(x = siteName, y = licorData)) +
  geom_boxplot(aes(x = siteName, y = licorData, color = depth)) +
  geom_jitter(width = 0.4, size = 0.1) +
  theme_bw()

# Ridge plot for all data (not subset by standard deviaton)
plot_all_data <- d3 %>%
  ggplot(aes(x = licorData, y = siteName, fill = paste(siteType,depth))) +
  ggtitle("Plot All Data")+
  geom_density_ridges(jittered_points=TRUE, scale = .95, rel_min_height = .01,
                      point_shape = "|", point_size = 2, size = 0.1,
                      position = position_points_jitter(height = 0)) +
  # scale_fill_manual(values=c("tomato2", "salmon", "royalblue1", "royalblue4", "deepskyblue1"))+
  theme_ridges(center = T)

# Ridge plot for data subset by standard devistion
plot_sd_subset <- sd_subset %>%
   ggplot(aes(x = licorData, y = siteName, fill = paste(siteType,depth))) +
   ggtitle("Plot Data Subset by SD")+
   geom_density_ridges(jittered_points=TRUE, scale = .95, rel_min_height = .01,
                       point_shape = "|", point_size = 2, size = 0.1,
                       position = position_points_jitter(height = 0)) +
   # scale_fill_manual(values=c("tomato2", "salmon", "royalblue1", "royalblue4", "deepskyblue1"))+
   theme_ridges(center = T)

# plot both side by side
grid.arrange(plot_all_data, plot_sd_subset, ncol = 2)
# Based on this, I feel like we shouldn't subset based on standard deviation. It doesn't solve our problem sites and it emphasizes differences in our good sites that maybe shouldn't exist

# so I'm going to finalize "d3" as the "licor_data" for downstream analyses...for now ;)
licor_data <- d3

# plot each site type
licor_data %>% ggplot(aes(x = siteType, y = licorData)) +
  geom_boxplot(aes(x = siteType, y = licorData, fill = depth)) +
  geom_jitter(width = 0.4, size = 0.1) +
  theme_bw()

licor_data %>% ggplot(aes(x = depth, y = licorData)) +
  geom_boxplot(aes(x = depth, y = licorData, fill = siteType)) +
  # geom_jitter(width = 0.4, size = 0.1) +
  scale_fill_manual(values = c("salmon", "royalblue4")) +
  theme_bw()

licor_data %>% ggplot(aes(x = siteType, y = licorData)) +
  geom_boxplot(aes(x = siteType, y = licorData, fill = siteType)) +
  scale_fill_manual(values = c("salmon", "royalblue4")) +
  theme_bw()

# plot mean ± standard deviation for each site and depth
sumstats <- licor_data %>%
  ungroup() %>%
  filter(cloudy=="N") %>%
  summarySE("licorData", c("siteName","depth"))
sumstats
str(sumstats)

ggplot(data = sumstats, 
       aes(x = depth,
           y = licorData, 
           ymin = licorData-sd, 
           ymax = licorData+sd,       
           fill = siteName)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(position = position_dodge(1), colour="black", width = 0.3, size = 0.5) +
  theme_bw()

ggplot(data = sumstats, 
       aes(x = siteName,
           y = licorData, 
           ymin = licorData-sd, 
           ymax = licorData+sd,       
           fill = depth)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(position = position_dodge(1), colour="black", width = 0.3, size = 0.5) +
  theme_bw()


# STOP ------
# test code below

# Subset for depths at each site?
# 
# fiv <- d [ which(d$depth=='5'), ]
# six <- d [ which(d$depth== '6'), ]
# eig <- d [ which(d$depth== '8'), ]
# nin <- d [ which(d$depth== '9'), ]
# ten <- d [ which(d$depth== '10'), ]
# fif <- d [ which(d$depth== '15'), ]
# sev <- d [ which(d$depth== '17'), ]

# Subet for good times at each site?
STRIfiv <- subset(sd_subset, siteName== 'STRI' & depth== '5' & cloudy=='N')

# STRIsev <- subset(d, siteName== 'STRI' & depth== '17')
# PLeig <- subset(d, siteName== 'PUNTA LAUREL' & depth== '8')
# PDten <- subset(d, siteName== 'PUNTA DONATO' & depth== '10')
# BNfiv<- subset(d, siteName== 'BASTIM NORTH' & depth== '5')
# BNfif<- subset(d, siteName== 'BASTIM NORTH' & depth== '15')
# BSsix<- subset(d, siteName== 'BASTIM SOUTH' & depth== '6')
# BSten<- subset(d, siteName== 'BASTIM SOUTH' & depth== '10')
# PIsix<- subset(d, siteName== 'POPA ISLAND' & depth== '6')
# DMeig<- subset(d, siteName== 'DRAGO MAR' & depth== '8')
# DMnin<- subset(d, siteName== 'DRAGO MAR' & depth== '9')
# DMten<- subset(d, siteName== 'DRAGO MAR' & depth== '10')

# plot by depth 

STRIfiv %>% ggplot(aes(x = localTime, y = licorData)) +
  geom_line()+
  theme_bw()

# Set limits for times 
# the first couple of minutes are duing the dive down
# it got cloudy at 12:15

start <- hms("11:30:00")
end <- hms("12:10:00")

head(d)

d %>% 
  filter(hms(localTime) >= start & hms(localTime) < end, siteName == "STRI", depth == "5", cloudy == "N") %>%
  ggplot(aes(x = localTime, y = licorData)) +
  geom_point() +
  geom_line()+
  theme_bw()


# mean and standard deviation
d %>% 
  filter(hms(localTime) >= start & hms(localTime) < end, siteName == "STRI", depth == "5", cloudy == "N") %>%
  summarySE("licorData")

d %>% 
  filter(hms(localTime) >= start & hms(localTime) < end, siteName == "STRI", depth == "5", cloudy == "N") %>%
  c(summarize(mean(licorData)+sd(licorData)), summarize(mean(licorData)-sd(licorData)))


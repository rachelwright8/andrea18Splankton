library(Rmisc) # for summarySE
library(tidyverse) # for readr, dplyr, ggplot...
library(lubridate) # for dealing with dates and times (e.g., hms)

setwd("andreas18ssequencing/")

# Load the file
d0 <- read_csv("newTimes.csv")
head(d0)

# plot lycor input over time
d0 %>% ggplot(aes(x = localTime, y = licorData, group = siteType, color = siteType)) +
  geom_line()+
  facet_grid(siteName ~ .)+
  theme_bw()

# Get rid of the first 5 measurements of each depth for each site (time during descent)----
d <- d0 %>% group_by(siteName,depth) %>% slice(5:n())

# Subset for values within 2 standard deviations of the mean for each depth and site ----
sd_subset <- d %>%
  filter(cloudy=="N") %>%
  group_by(siteName,depth) %>%
  mutate(mean = mean(licorData), sd = sd(licorData)) %>%
  filter(licorData > mean-(2*sd), licorData < mean+(2*sd) )

sd_subset %>% ggplot(aes(x = localTime, y = licorData, group = siteType, color = siteType)) +
  geom_line()+
  facet_grid(siteName ~ .)+
  theme_bw()

str(sd_subset)
sd_subset$depth <- as.factor(sd_subset$depth)

sd_subset %>% ggplot(aes(x = siteName, y = licorData)) +
  geom_boxplot(aes(x = siteName, y = licorData, color = depth)) +
  geom_jitter(width = 0.4, size = 0.1) +
  theme_bw()

# plot each site type
sd_subset %>% ggplot(aes(x = siteType, y = licorData)) +
  geom_boxplot(aes(x = siteType, y = licorData, fill = depth)) +
  # geom_jitter(width = 0.4, size = 0.1) +
  theme_bw()

sd_subset %>% ggplot(aes(x = depth, y = licorData)) +
  geom_boxplot(aes(x = depth, y = licorData, fill = siteType)) +
  scale_fill_manual(values = c("salmon", "royalblue4")) +
  theme_bw()

sd_subset %>% ggplot(aes(x = siteType, y = licorData)) +
  geom_boxplot(aes(x = siteType, y = licorData, fill = siteType)) +
  scale_fill_manual(values = c("salmon", "royalblue4")) +
  theme_bw()

# plot mean ± standard deviation for each site and depth
sumstats <- sd_subset %>% 
  filter(cloudy=="N") %>%
  summarySE("licorData", c("siteName","depth"))
sumstats
str(sumstats)
sumstats$depth <- as.factor(sumstats$depth) # recast depth as factor instead of integer

ggplot(data = sumstats, 
       aes(x = depth,
           y = licorData, 
           ymin = licorData-sd, 
           ymax = licorData+sd,       
           fill = siteName)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(position = position_dodge(), colour="black") +
  theme_bw()

ggplot(data = sumstats, 
       aes(x = siteName,
           y = licorData, 
           ymin = licorData-sd, 
           ymax = licorData+sd,       
           fill = depth)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(position = position_dodge(), colour="black") +
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


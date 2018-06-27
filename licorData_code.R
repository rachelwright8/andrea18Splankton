library(Rmisc) # for summarySE
library(tidyverse) # for readr, dplyr, ggplot...
library(lubridate) # for dealing with dates and times (e.g., hms)

setwd("andreas18ssequencing/")

# Load the file
d <- read_csv("newTimes.csv")
head(d)

# plot lycor input over time
d %>% ggplot(aes(x = localTime, y = licorData, group = siteType, color = siteType)) +
  geom_line()+
  facet_grid(siteName ~ .)+
  theme_bw()



# STOP 6/26/18 -------
# Subset for depths at each site?
# Subet for good times at each site?




# Set limits for times -----
# the first couple of minutes are duing the dive down
# it got cloudy at 12:15

start <- hms("11:30:00")
end <- hms("12:10:00")

d %>% 
  filter(hms(localTime) >= start & hms(localTime) < end) %>%
  ggplot(aes(x = localTime, y = input)) +
  geom_point() +
  geom_line()+
  theme_bw()


# mean and standard deviation
d %>% 
  filter(hms(localTime) >= start & hms(localTime) < end) %>%
  summarySE("input")



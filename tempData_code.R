# setwd("~/Desktop/andreas18ssequencing/")

# Load libraries
library(tidyverse)
library(stringr)

# Load data
temp_data <- read.csv("minmeanmax_tempdata.csv")

# Reformate date
temp_data <- temp_data %>% separate(date, c("month", "day", "year"), "/")
head(temp_data)

summary(temp_data)

# Add site names
temp_data <- temp_data %>% mutate(site_name = ifelse(site_number=="IR1", "PuntaDonato", ifelse(site_number == "IR2", "STRIPoint", ifelse(site_number== "IR3", "Cristobal", ifelse(site_number=="OR3", "PopaIsland", "DragoMar")))))

temp_data$site_name <- as.factor(temp_data$site_name)

summary(temp_data)
  
# Subset data for just June 2015
sub <- temp_data %>% filter(month=="6" & year=="15")

# Plot temp over time
ggplot(sub, aes(x = day, y = Max, group = site_name, color = site_name)) +
  labs(x = "Day in June", y = "Max Temperature °C")+
  geom_line()+
  geom_point()+
  theme_bw()


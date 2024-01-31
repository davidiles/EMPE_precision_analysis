# install/load necessary packages
my.packs <- c(
  
  # Bayesian analysis
  'jagsUI','MCMCvis',
  
  # Data manipulation
  'tidyverse','reshape2','lubridate','readxl',
  
  
  # Spatial analysis
  'rgeos','raster','sp','sf',
  
  # Plotting
  'ggrepel','ggthemes','ggpubr','scales','viridis')

if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/EMPE_Global/analysis")

rm(list=ls())


# ****************************************************************************************************************
# ****************************************************************************************************************
# PART 1: FORMAT DATA AND PREPARE FOR ANALYSIS 
# ****************************************************************************************************************
# ****************************************************************************************************************


# --------------------------------------------------
# Colony attributes
# --------------------------------------------------

year_range <- 2009:2018
n_years <- length(year_range)

# --------------------------------------------------
# Colony attributes
# --------------------------------------------------

colony_attributes <- read.csv("../data/colony_attributes.csv")

# --------------------------------------------------
# Read in satellite data
# --------------------------------------------------

sat <- read_xlsx("../data/empe_satellite_2023-05-25.xlsx") %>%
  mutate(area_m2 = as.numeric(area_m2)) %>%
  subset(!is.na(area_m2) ) # & img_year %in% year_range & site_id %in% colony_attributes$site_id
sat$adult_count <- as.numeric(sat$area_m2)
sat$img_qualit <- as.numeric(sat$img_qualit)
sat <- sat %>% dplyr::rename(year = img_year)
sat$obs_id <- 1:nrow(sat)

# --------------------------------------------------
# Read in aerial data
# --------------------------------------------------

aer <- read_xlsx("../data/empe_aerial_2023-05-25.xlsx") %>%
  subset(count_type == "adults") %>%
  rename(adult_count = penguin_count, adult_accuracy = accuracy)
aer$adult_count <- as.numeric(aer$adult_count)
aer$obs_id <- 1:nrow(aer)

# --------------------------------------------------
# Format survey dates
# --------------------------------------------------

aer$Date <- lubridate::ymd(paste(aer$year,aer$month,aer$day, sep="-"))
sat$Date <- lubridate::ymd(paste(sat$year,sat$img_month,sat$img_day, sep="-"))

aer$yday <- yday(aer$Date) 
sat$yday <- yday(sat$Date) 

# --------------------------------------------------
# Only select aerial observations from suitable year ranges
# --------------------------------------------------

aer_orig <- aer # Dataframe to store full dataset (only a subset of these will be used in analysis)
sat_orig <- sat # Dataframe to store full dataset (only a subset of these will be used in analysis)

# Retain surveys within appropriate date range (2009:2018)
aer1 <- aer %>%
  dplyr::select(-reference) %>%
  subset(!is.na(adult_count) & year %in% year_range & site_id %in% colony_attributes$site_id)

# Retain surveys within appropriate date range (2009:2018)
sat1 <- sat %>%
  subset(!is.na(area_m2) & year %in% year_range & site_id %in% colony_attributes$site_id)

# Numeric year identifier (2009 = 1 ... 2010 = 2 ... 2018 = 10)
aer1$year_number <- match(aer1$year,year_range)
sat1$year_number <- match(sat1$year,year_range)

# Surveys excluded based on year range
aer_excluded_YearRange <- subset(aer_orig, !(obs_id %in% aer1$obs_id))
sat_excluded_YearRange <- subset(sat_orig, !(obs_id %in% sat1$obs_id))

# --------------------------------------------------
# Exclude surveys outside ideal annual "survey window" (Sep 1 to Nov 1 each year)
# --------------------------------------------------

min_date <- lubridate::ymd("2023-09-01") %>% lubridate::yday()
max_date <- lubridate::ymd("2023-11-30") %>% lubridate::yday()

# aerial surveys within appropriate year and month range
aer2 <- subset(aer1, yday >= min_date & yday <= max_date)

# satellite surveys within appropriate year and month range
# Note: also remove any "apparent colony absences" detected from satellite after nov 1 (colony status is too uncertain)
nov1_yday <- lubridate::ymd("2023-11-01") %>% lubridate::yday()

sat2 <- subset(sat1, 
               (yday >= min_date & yday <= max_date) &
                 !(yday >= nov1_yday & area_m2 == 0))  

# Surveys excluded based on day range
aer_excluded_DayRange <- subset(aer1, !(obs_id %in% aer2$obs_id))
sat_excluded_DayRange <- subset(sat1, !(obs_id %in% sat2$obs_id))

aer_excluded_DayRange %>% dim() # 1 - aug to nov
sat_excluded_DayRange %>% dim() # 23 - aug to nov


# --------------------------------------------------
# Remove any additional surveys as necessary (e.g., inadequate colony coverage, extremely poor visual conditions, etc)
# --------------------------------------------------

sat2 <- subset(sat2, is.na(remove))

# --------------------------------------------------
# Create numeric identifiers for each colony (for inclusion in model)
# --------------------------------------------------

# Only include sites that have been surveyed at least once between 2009 and 2018
sites <- c(aer2$site_id,sat2$site_id) %>% unique() %>% sort()
n_sites <- length(sites)

# Order colonies by longitude
colony_attributes <- subset(colony_attributes, site_id %in% sites) %>%
  arrange(lon) %>%
  mutate(site_id_factor = fct_reorder(site_id, lon))
colony_attributes$site_number <- 1:nrow(colony_attributes)

# --------------------------------------------------
# Final, cleaned datasets 
# --------------------------------------------------

aer <- left_join(aer2, colony_attributes, by = c("site_id"))
sat <- left_join(sat2, colony_attributes, by = c("site_id"))

# --------------------------------------------------
# Create "day of season" covariate
# --------------------------------------------------
baseline_day <- mean(c(min_date,max_date))
aer$day_of_season <- aer$yday - baseline_day
sat$day_of_season <- sat$yday - baseline_day

lubridate::ymd("2009-01-01") + 288

# --------------------------------------------------
# Package data for JAGS
# --------------------------------------------------

jags.data <- list( n_years = n_years,
                   n_sites = n_sites,
                   
                   # aerial counts of adults
                   n_obs_aerial = length(aer$adult_count),
                   adult_count = aer$adult_count,
                   aerial_site = aer$site_number,
                   aerial_year = aer$year_number,
                   aerial_month = aer$month,
                   aerial_DoS = aer$day_of_season,
                   
                   # satellite counts
                   n_obs_satellite = length(sat$area_m2),
                   img_qual = sat$img_qualit,
                   satellite = sat$area_m2,
                   satellite_site = sat$site_number,
                   satellite_year = sat$year_number,
                   satellite_month = sat$img_month,
                   satellite_DoS = sat$day_of_season
)

# For calculating log-linear trend (least squares regression line)
XX <-  cbind(rep(1,jags.data$n_years),1:jags.data$n_years)
jags.data$regression_weights <- matrix(c(0,1),1,2) %*% solve(t(XX) %*% XX) %*% t(XX)

# --------------------------------------------------
# Save data for subsequent analysis
# --------------------------------------------------

save.image("output/EMPE_data_formatted.RData")

# ****************************************************************************************************************
# ****************************************************************************************************************
# # PART 2: EXPLORATORY ANALYSIS / DATA VISUALIZATION
# ****************************************************************************************************************
# ****************************************************************************************************************

# -----------------------------------------
# Examine number of colonies surveyed per year
# -----------------------------------------

aer_annual <- aer %>%
  group_by(site_id,year) %>%
  summarize(n_aer = n())

sat_annual <- sat %>%
  group_by(site_id,year) %>%
  summarize(n_sat = n())

colonies_surveyed_annual <- full_join(aer_annual,sat_annual)
colonies_surveyed_annual[is.na(colonies_surveyed_annual)] <- 0
colonies_surveyed_annual <- colonies_surveyed_annual %>%
  mutate(aer_only = (n_aer>0) & (n_sat == 0),
         sat_only = (n_aer == 0) & (n_sat > 0),
         both = (n_aer>0) & (n_sat>0)) %>%
  group_by(year) %>%
  summarize(aer_only = sum(aer_only),
            sat_only = sum(sat_only),
            both = sum(both),
            total = length(unique(site_id)))

colonies_surveyed_annual <- colonies_surveyed_annual %>%
  pivot_longer(cols = aer_only:total, names_to = "survey_type", values_to = "n") %>%
  subset(survey_type != "total")

colonies_surveyed_annual$survey_type <- factor(colonies_surveyed_annual$survey_type,
                                               levels = c("aer_only","sat_only","both"))

plot1 <- ggplot(colonies_surveyed_annual, aes(x = year, y = n, fill = survey_type)) +
  geom_bar(stat = "identity", col = "black")+
  scale_x_continuous(breaks = year_range)+
  ylab("# colonies")+
  xlab("Year")+
  scale_fill_manual(values = c("white","gray80","black"), name = "Survey Type",
                    labels = c("Aerial Only","Satellite Only","Both"))+
  theme_bw()+
  ggtitle("Number of colonies surveyed per year")

png("output/data_viz/plot1_nsurveys.png", units = "in", res = 500, 
    width = 8, height = 5)
print(plot1)
dev.off()

# -----------------------------------------
# Examine phenology of surveys, and how it has changed over time
# -----------------------------------------

aer_colony_months <- aer %>%
  group_by(site_id,year,month) %>%
  summarize(n_aer = n()) %>%
  group_by(year,month) %>%
  summarize(n = sum(n_aer>0))
aer_colony_months$month.name <- month.name[aer_colony_months$month]
aer_colony_months$month.name <- factor(aer_colony_months$month.name, levels = c("August","September","October","November","December"))

p2_aer <- ggplot(aer_colony_months, aes(x = year, y = n, fill = month.name)) +
  geom_bar(stat = "identity", width = 0.2)+
  scale_fill_manual(values = c("darkblue","dodgerblue","orange","red","darkred"), name = "Month", drop = FALSE)+
  theme_bw()+
  ggtitle("Aerial")+
  ylab("Number of colonies surveyed")

sat_colony_months <- sat %>%
  group_by(site_id,year,img_month) %>%
  summarize(n_sat = n()) %>%
  dplyr::rename(month = img_month) %>%
  group_by(year,month) %>%
  summarize(n = sum(n_sat>0))
sat_colony_months$month.name <- month.name[sat_colony_months$month]
sat_colony_months$month.name <- factor(sat_colony_months$month.name, levels = c("August","September","October","November","December"))
p2_sat <- ggplot(sat_colony_months, aes(x = year, y = n, fill = month.name)) +
  geom_bar(stat = "identity", width = 0.2)+
  scale_fill_manual(values = c("darkblue","dodgerblue","orange","red","darkred"), name = "Month", drop = FALSE)+
  theme_bw()+
  ggtitle("Satellite")+
  ylab("Number of colonies surveyed")

# combine into single figure
plot2 <- ggarrange(p2_aer,p2_sat,nrow=2,align="hv")

png("output/data_viz/plot2_survey_phenology.png", units = "in", res = 500, 
    width = 8, height = 8)
print(plot2)
dev.off()

# -----------------------------------------
# Raw observations at each colony
# - point type for satellite or aerial
# - point colour for "day of year"
# - point size for satellite image quality
# - observed absences with a red vertical line
# -----------------------------------------
aer$img_qualit = 3
plot3 <- ggplot() +
  
  # Satellite observations
  geom_vline(data = subset(sat, area_m2 == 0), aes(xintercept = year), col = "gray80", size=2)+
  geom_point(data = sat, aes(x = year, y = area_m2, shape = "Satellite count (area in m2)", col = yday,
                             size = factor(img_qualit)), stroke = 0.5)+
  
  # Aerial observations
  geom_point(data = aer, aes(x = year, y = adult_count, shape = "Aerial count (adult)", col = yday,size = factor(img_qualit)))+
  
  # Scales
  scale_shape_manual(name = 'Survey Type', values =c('Satellite count (area in m2)'=10,'Aerial count (adult)'= 17))+
  scale_color_gradientn(name = 'Date of survey\n(day of year)', colors = magma(10)[1:9])+
  scale_size_manual(values = c(0.5,1,2), name = "Image Quality",
                    labels = c("1 - Poor","2 - Moderate","3 - Good"))+
  scale_x_continuous(limits = range(year_range))+
  
  ylab("Count")+
  xlab("Year")+
  facet_wrap(site_id_factor~., scales = "free",nrow=10,ncol=5)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, col = "transparent")

pdf("output/data_viz/plot3_raw_data.pdf", width = 20, height = 20)
print(plot3)
dev.off()
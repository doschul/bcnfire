#### validate actual burning patterns ####

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

# load all shp files in subfolder of ./data/wildfires/wildfires_shp
library(sf)
library(tidyverse)
library(stringr)
library(mapview)


# load grd shape
grd <- st_read("./data/rdat/grd.gpkg")

# convez hull polygon
ch <- st_convex_hull(st_union(grd)) %>%
  st_as_sf() %>%
  st_make_valid()

# par(mfrow = c(1, 1))
# 
# wf_df <- read.csv("./data/wildfires/Incendis_forestals_a_Catalunya._Anys_2011-2023_20250725.csv", stringsAsFactors = FALSE)
# wf_df$total_burned_ha <- rowSums(wf_df[, c("HAARBRADES",  "HANOARBRAD", "HANOFOREST", "HAFORESTAL")], na.rm = TRUE)
# 
# table(wf_df$total_burned_ha > 0)
# table(wf_df$total_burned_ha > 4)
# 
# hist(wf_df$total_burned_ha[wf_df$total_burned_ha > 0 & wf_df$total_burned_ha < 1])

wildfires <- list.files("./data/wildfires/wildfires_shp", pattern = ".shp$", 
                        recursive = TRUE, full.names = TRUE, include.dirs = TRUE) |>
  lapply(st_read) |>
  bind_rows()

wildfires <- wildfires |>
  filter(!is.na(CODI_FINAL)) |>
  mutate(
    year = case_when(nchar(DATA_INCEN) == 7 ~ as.integer(paste0("20", substr(DATA_INCEN, 6, 7))),
                     nchar(DATA_INCEN) == 8 ~ as.integer(paste0("20", substr(DATA_INCEN, 7, 8))),
                     nchar(DATA_INCEN) >= 10 ~ as.integer(substr(DATA_INCEN, 7, 10)),
                     TRUE ~ NA),
    month = case_when(nchar(DATA_INCEN) == 7 ~ as.integer(substr(DATA_INCEN, 3, 4)),
                     nchar(DATA_INCEN) >= 8 ~ as.integer(substr(DATA_INCEN, 4, 5)),
                     TRUE ~ NA),
    # day is anything before first slash
    day = as.integer(substring(DATA_INCEN, 1, regexpr("/", DATA_INCEN) - 1)),
    date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d"),
    area_ha = as.numeric(st_area(geometry)) / 10000,  # convert from m² to ha
    MUNICIPI = str_to_title(MUNICIPI)  # capitalize first letter of each word
  ) |>
  select(CODI_FINAL, MUNICIPI, date, year, month, area_ha, geometry)


# transform crs
wildfires <- st_transform(wildfires, crs = st_crs(grd)) %>%
  st_make_valid()


# wildfires within study region between 2015 and 2025
wf_aoi <- wildfires |>
  # valid geom only
  filter(st_is_valid(geometry)) |>
  filter(year >= 2015 & year <= 2025)

wf_aoi <- wf_aoi[st_intersects(wf_aoi, ch, sparse = FALSE), ]


mapview(wf_aoi, zcol = "year", legend = TRUE, layer.name = "Wildfires") +
  # add convex hull polygon, no fill, just border
  mapview(ch, color = "black", lwd = 2, layer.name = "Convex Hull")
  

hist(wf_aoi$area_ha, breaks = 15, main = "Distribution of Burned Area (ha)", xlab = "Area (ha)")


View(wf_aoi)







##### Ignition probability #####



wildfire_centroids <- st_centroid(wf_aoi) |>
  st_transform(crs = st_crs(grd))

overlay <- st_join(grd, wildfire_centroids)
overlay$wildfire <- ifelse(is.na(overlay$CODI_FINAL), 0, 1)

overlay_rest <- overlay[overlay$lc_artificial < 0.1, ]
table(overlay_rest$wildfire)

# density of access grouped by wildfire
library(ggplot2)

ggplot(overlay_rest, aes(x = road_distance, fill = as.factor(wildfire))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Road Distance by Wildfire Occurrence",
       x = "Road distance",
       fill = "Wildfire Occurrence [m]") +
  theme_minimal()



# logistic regression, exclude urban areas
m1 <- glm(wildfire ~ ignition_probability_surface, data = overlay_rest, family = binomial(link = "logit"))
summary(m1)






#### test SAL cost assumptions ####
# two random correlated variables with rho = 0.5
library(mvtnorm)
library(tidyverse)

set.seed(123)
n <- 10000

mu <- c(0, 0)
rho <- 0

sigma <- matrix(c(1, rho, rho, 1), nrow = 2)

data <- rmvnorm(n, mean = mu, sigma = sigma)
colnames(data) <- c("slope", "access")

# create a data frame
df <- data.frame(slope = data[, 1], access = data[, 2]) |>
  # scale both to 0-1, making new variables
  mutate(slope_sc = (slope - min(slope)) / (max(slope) - min(slope)),
         access_sc = (access - min(access)) / (max(access) - min(access))) |>
  # add scaled variables
  mutate(slope_access_add = slope_sc + access_sc) |>
  # multiply scaled variables
  mutate(slope_access_mult = slope_sc * access_sc)

# make two plots in one
par(mfrow = c(2, 1))
hist(df$slope_access_add, main = "Sum of Scaled Variables", xlab = "Sum")
hist(df$slope_access_mult, main = "Product of Scaled Variables", xlab = "Product")

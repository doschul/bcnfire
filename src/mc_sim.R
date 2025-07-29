# decision Support 

rm(list=ls())

library(decisionSupport)
library(tidyverse)

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcnfire")

# https://cran.r-project.org/web/packages/decisionSupport/vignettes/wildfire_example.html
#format(read.csv("./src/decisionSupport/wildfire_input_table.csv", sep = ";")[,1:5],scientific=FALSE)


dat <- read.csv("./src/decisionSupport/wildfire_input_table.csv", sep = ";")
write.csv(dat, file = "./src/decisionSupport/wildfire_input.csv")


load("./out/mc_bau_lc_long.RData")

burn.sim <- function(){
  
  # load data from parent environment
  mc_bau_lc_long <- get("mc_bau_lc_long", envir = parent.frame())
  #load("./out/mc_bau_lc_long.RData")
  
  start_cell <- sample(unique(mc_bau_lc_long$ignition_cell), 1)

  urban_burned <- mc_bau_lc_long %>%
    filter(time == time_horizon) %>%
    filter(landcover == "urban") %>%
    filter(ignition_cell == start_cell) %>%
    mutate(house_cost = houses_per_ha * house_price * burned * 4,
           ffight_cost = fire_fight_cost_per_ha * burned * 4)
  
  # costs per scenario
  c.burn.houses.bau <- urban_burned %>%
    filter(scenario == "BAU")
  
  c.burn.houses.sal <- urban_burned %>%
    filter(scenario == "SAL")
  
  #total_area <- 59446 * 4 # total area in ha
  
  c.logging <- salvage_logging_cost_per_ha * total_wildland_area
  
  # sum up costs of salvage logging and benefits of business as usual
  c.tot.sal <- c.logging + c.burn.houses.sal$house_cost + c.burn.houses.sal$ffight_cost
  c.tot.bau <- c.burn.houses.bau$house_cost + c.burn.houses.bau$ffight_cost
  
  # calculate NPV
  NPV.dif <- c.tot.bau - c.tot.sal
  
  # return NPV
  return(NPV.dif)
}



decisionSupport(inputFilePath = "./src/decisionSupport/wildfire_input.csv", #input file with estimates
                outputPath = "./out/decisionSupport", #output folder
                write_table = TRUE,
                welfareFunction = burn.sim,
                numberOfModelRuns = 1000,
                functionSyntax = "plainNames")

e1 <- read.csv("./src/decisionSupport/wildfire_input.csv") %>%
  as.estimate()

# mcres1 <- mcSimulation(estimate = e1,
#                        model_function = burn.sim,
#                        numberOfModelRuns = 1000, 
#                        functionSyntax = "plainNames")

# load MC results
mc_res <- read.csv("./out/decisionSupport/mcSimulationResults.csv")


# plot histogram of NPV

# remove 1% extreme outliers on both ends
mc_res <- mc_res %>%
  filter(output_1 > quantile(output_1, 0.01) & output_1 < quantile(output_1, 0.99))
  


ggplot(data = mc_res, aes(x = output_1)) +
  geom_histogram(fill = "lightblue", bins = 100) +
  # horizontal line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", lwd = 1.5) +
  labs(title = "Histogram of NPV",
       x = "NPV",
       y = "Frequency") +
  theme_minimal()


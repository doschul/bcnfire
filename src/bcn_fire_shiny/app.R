
##### load data #####

#setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/bcn_fire")

# get functions
source("bcn_funs.R")

# libraries
library(tidyverse)
library(shiny)
library(leaflet)
library(dplyr)
library(sf)



# load grid
grd <- st_read("grd_bau.gpkg") %>%
  # transform cs to lon lat
  st_transform(., 4326)

# load neighbors
load("neighbor_idx.RData")

##### Shiny app #####


ui <- fluidPage(
  leafletOutput("map"),
  sliderInput("slider", "Select a Time horizon:", min = 1, max = 24, value = 4),
  verbatimTextOutput("selected_cell")
)

server <- function(input, output) {
  
  selected_cell <- reactiveVal(NULL)
  burned_cells <- reactiveVal(0)
  partially_burned_cells <- reactiveVal(0)
  
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addPolygons(data = grd,
                  layerId = ~id,
                  color = "black",
                  fillColor = "black",
                  weight = 1,
                  fillOpacity = 0.1,
                  fill = T)
  })
  
  observeEvent(input$map_shape_click, {
    clicked_cell <- input$map_shape_click$id
    selected_cell(clicked_cell)
  })
  
  observe({
    req(selected_cell())
    
    # run get burners with ignition cell and t 10
    burners <- get_burners(time_horizon = input$slider,
                           ignition_cell = selected_cell(),
                           cell_size = 200, # in meters
                           time_step = 15,    # in minutes
                           scenario = "bau", # select scenario (bau, sal)
                           full = FALSE, # return full data frame or last column
                           grd = grd,
                           neighbor_idx = neighbor_idx)
    
    # set 0 burner to NA
    burners$burning[burners$burning == 0] <- NA
    
    # count cells completely burned (value of 1)
    burned_cells(sum(burners$burning == 1, na.rm = TRUE))
    partially_burned_cells(sum(burners$burning, na.rm = TRUE))
    
    pal <- colorNumeric(palette = colorRampPalette(c("yellow", "red"))(100), domain = seq(0, 1, length.out = 10))
    
    # Update the polygons based on the clicked cell
    leafletProxy("map") %>%
      addPolygons(data = burners,
                  layerId = ~id,
                  color = "darkgrey",
                  fillColor = ~pal(burning),
                  fillOpacity = 0.8)
  })
  
  output$selected_cell <- renderPrint({
    paste("Selected Cell ID:", selected_cell()) 
    paste("Burned cells:", burned_cells())
    paste("Incluidng partially burned cells:", partially_burned_cells())
  })
}

shinyApp(ui, server)















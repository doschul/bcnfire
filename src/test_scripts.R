

#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyWidgets)

time_labels <- expand.grid(
  hours = 0:23,  # Up to 23 hours
  minutes = seq(0, 55, by = 5) # 5-minute intervals
) %>%
  mutate_at(1:2, ~ifelse(nchar(.) == 1, paste0('0', .), .)) %>%
  mutate(time_label = paste0(hours, ':', minutes)) %>%
  arrange(hours, minutes) %>% # Important for correct ordering
  pull(time_label)



ui <- fluidPage(
  sliderTextInput(
    "mySliderText",
    "Select a time (hours:minutes)",
    choices = time_labels,
    width = "400px" # Adjust width for better display
  ),
  textOutput('selection')
)

server <- function(input, output, session) {
  
  output$selection <- renderText({
    req(input$mySliderText)
    
    my_hours <- as.numeric(substr(input$mySliderText, 1, 2))
    my_minutes <- as.numeric(substr(input$mySliderText, 4, 5))
    
    total_minutes <- 60 * my_hours + my_minutes
    
    out <- paste0(
      'You selected ', my_hours, ' hours, ', my_minutes,
      ' minutes (', scales::comma(total_minutes), ' total minutes)'
    )
    
    out
    
  })
  
}

shinyApp(ui, server)

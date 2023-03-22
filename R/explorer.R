#' Simple graphical interface to explore group response curves
#'
#' This function launches a simple graphical interface to explore fire response
#' curves. Given a selected group, the response to each fire regime component
#' (frequency, severity and time since fire) based on expert-elicited data are
#' shown, together with the resulting overall response distribution(s). Overall
#' responses can be displayed for multiple fire regimes, but at present these
#' can only differ in time since and must share the same frequency and severity
#' values. This restriction will be relaxed in a future version of the package.
#'
#' The \code{shiny} package must be installed to use this function.
#'
#' @importFrom ggplot2 labs theme_bw
#'
#' @export
#
explorer <- function() {
  if (!("shiny" %in% installed.packages()[, "Package"])) {
    stop("You need to install the 'shiny' package to use this function")
  }

  library(shiny)

  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput("groupId", label = "Group", choices = 1:18),

        selectInput("fireFreq", label = "Fire frequency (50 years)",
                    choices = sort(unique(GroupOverallResponse$frequency))),

        selectInput("fireSeverity", label = "Severity of last fire",
                    choices = sort(unique(GroupOverallResponse$severity))),

        sliderInput("fireTSF", label = "Years since last fire",
                    min = 0, max = max(GroupOverallResponse$tsf),
                    value = 0,
                    step = 1,
                    round = TRUE),

        actionButton("addTSF", "Add", class = "btn-sm"),
        actionButton("resetTSF", "Reset", class = "btn-sm")
      ),

      mainPanel(
        plotOutput("componentGraphs"),
        plotOutput("overallGraph")
      )
    )
  )

  server <- function(input, output, session) {
    grp <- reactive({
      as.integer( input$groupId )
    })

    freq <- reactive({
      as.integer( input$fireFreq )
    })

    sev <- reactive({
      as.integer( input$fireSeverity )
    })

    tsfs <- reactiveVal(0)

    observeEvent(input$addTSF, {
      curVals <- tsfs()
      newVal <- as.integer( input$fireTSF )
      if (!(newVal %in% curVals)) tsfs(c(curVals, newVal))
    })

    observeEvent(input$resetTSF, {
      newVal <- as.integer( input$fireTSF )
      tsfs(newVal)
    })

    output$componentGraphs <- renderPlot({
      fireresponse::draw_response_curves(grp()) +
        labs(title = glue::glue("Group {grp()} response curves")) +
        theme_bw()
    })

    output$overallGraph <- renderPlot({
      fireresponse::draw_overall_response(grp(),
                                          frequency = freq(),
                                          severity = sev(),
                                          tsf = tsfs(),
                                          draw_pzero = TRUE) +

        labs(title = glue::glue("Group {grp()} overall response")) +
        theme_bw()
    })
  }

  shinyApp(ui, server)
}

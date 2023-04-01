#' Simple graphical interface to explore group response curves
#'
#' This function launches a simple graphical interface to explore fire response
#' curves. Given a selected group, the response to each fire regime component
#' (frequency, severity and time since fire) based on expert-elicited data are
#' shown, together with the resulting overall response distribution(s). Overall
#' responses can be displayed for multiple fire regimes, but presently only for
#' one group at a time.
#'
#' The \code{shiny} package must be installed to use this function.
#'
#' @importFrom ggplot2 labs theme theme_classic element_text element_rect
#'
#' @export
#
explorer <- function() {
  if (!("shiny" %in% installed.packages()[, "Package"])) {
    stop("You need to install the 'shiny' package to use this function")
  }

  library(shiny)

  ggplot2::theme_set( theme_classic() +
                        theme(strip.text.x = element_text(size = 18),
                              axis.text = element_text(size = 12),
                              axis.title = element_text(size = 16),
                              legend.text = element_text(size = 12),
                              legend.title = element_text(size = 16),
                              plot.background = element_rect(fill = "white")) )

  ui <- fluidPage(
    theme = bslib::bs_theme(bootswatch = "journal"),

    sidebarLayout(
      sidebarPanel(
        selectInput("groupId", label = "Group", choices = 1:18),

        sliderInput("fireFreq", label = "Fire frequency (50 years)",
                    min = 0, max = max(GroupOverallResponse$frequency),
                    value = 0,
                    step = 1,
                    round = TRUE),

        sliderInput("fireSeverity", label = "Severity of last fire",
                    min = 0, max = max(GroupOverallResponse$severity),
                    value = 0,
                    step = 1,
                    round = TRUE),

        sliderInput("fireTSF", label = "Years since last fire",
                    min = 0, max = max(GroupOverallResponse$tsf),
                    value = 0,
                    step = 1,
                    round = TRUE),

        actionButton("addCurve", "Add", class = "btn-sm"),
        actionButton("reset", "Reset", class = "btn-sm")
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

    regimes <- reactiveVal()
    m <- matrix(c(0,0,0), nrow=1)
    regimes(m)

    observeEvent(input$addCurve, {
      cur_regimes <- regimes()

      f <- as.integer( input$fireFreq )
      s <- as.integer( input$fireSeverity )
      t <- as.integer( input$fireTSF )

      same <- apply(cur_regimes, MARGIN = 1, identical, c(f,s,t))
      if (!any(same)) {
        # add the regime as a new matrix row
        cur_regimes <- rbind(cur_regimes, c(f, s, t))
        regimes(cur_regimes)
      }
    })

    observeEvent(input$reset, {
      f <- as.integer( input$fireFreq )
      s <- as.integer( input$fireSeverity )
      t <- as.integer( input$fireTSF )
      regimes( matrix(c(f,s,t), nrow=1) )
    })

    output$componentGraphs <- renderPlot({
      fireresponse::draw_response_curves(grp()) +
        labs(title = glue::glue("Group {grp()} response curves"))
    })

    output$overallGraph <- renderPlot({
      cur_regimes <- regimes()
      frequency <- cur_regimes[,1]
      severity <- cur_regimes[,2]
      tsf <- cur_regimes[,3]

      fireresponse::draw_overall_response(grp(),
                                          frequency = frequency,
                                          severity = severity,
                                          tsf = tsf,
                                          cross = FALSE,
                                          draw_pzero = TRUE) +

        labs(title = glue::glue("Group {grp()} overall response"))
    })
  }

  shinyApp(ui, server)
}

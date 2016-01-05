library(shiny)
library(clusterPower)
library(xtable)
library(ggplot2)

# this function allows text inputs to be on the same line; default has them on new lines
textInputRow<-function (inputId, label, value = "") 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value, class="input-small"))
}

shinyUI(fluidPage(
  withMathJax(),
  # section below allows in-line LaTeX via $ in mathjax.
  tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                <script type='text/javascript' async src='path-to-mathjax/MathJax.js?config=TeX-AMS_CHTML'></script>")),
  # we don't like the default css, so we're going to tinker with it a little bit
  tags$head(
    tags$style(HTML("
                    # changing table borders to gray instead of black, increasing cell padding 
                    
                    table, td, th, tr {
                    border: 1px solid gray;
                    }
                    td, th {
                    padding: 5px;
                    }
                    shiny-table-class {
                    font: sans-serif;
                    }
                    
                    "))
    ), 
  titlePanel("Power for Cluster-Randomized Trials"),
  sidebarLayout(
    sidebarPanel(
      h5("Model:"),
      p('$y_{ij} = \\beta_0 + x_{ij}\\delta + b_{0i} +e_{ij} $'),
      tags$small('where $b_{0i} \\sim N(0, \\sigma_b^2)$, $e_{ij} \\sim N(0, \\sigma^2)$ 
                 and ICC = $\\sigma_b^2 / (\\sigma_b^2 + \\sigma^2)$'),
      br(), br(),
      textInput("N","Enter the number of subjects per cluster (N):",value="100"), # number of subjects per clusters
      textInput("M","Enter the number of clusters (M):",value="3"), # number of clusters
      textInput("d","Enter the difference $\\delta$ in mean between the two arms:",value="0.2"), # difference in mean between the two arms
      # alternative hypothesis (null is 0)
      
      # radio button: ICC or sigma-b
      radioButtons("rhosigmab", "Use either an ICC or $\\sigma_b^2$ value:", list("ICC","$\\sigma_b^2$" = "sigmab"), inline = TRUE),
      conditionalPanel(condition = "input.rhosigmab == 'ICC'",
                       textInputRow("ICC","Enter ICC value:",value="0.001")), # intra-cluster correlation coefficient
      conditionalPanel(condition = "input.rhosigmab == 'sigmab'",
                       textInputRow("sigmab","Enter $\\sigma_b^2$ value:", value="0.001001")),
      textInputRow("sigma","Enter $\\sigma^2$ value:",value="1"),
      
      ## alpha can be one of the three options
      radioButtons("alpha","Enter alpha value:",c("0.025", "0.05","0.10"),selected="0.05",inline=TRUE),
      
      ## options: range of sigma-b, use CV (0 otherwise)
      checkboxGroupInput("options", label="Options:", choices = c("Cluster size is variable" = "useCV", "Plot multiple values for $\\sigma_b^2$"
                                                                  = "multipleICC")),
      conditionalPanel(condition = "input.options == 'useCV' || input.options.length > 1",
                       textInput("CV", label="Enter the coefficient of variation (CV):")),
      conditionalPanel(condition = "input.options == 'multipleICC' || input.options.length > 1",
                       # calling it rho and ICC because I thought it was, and sigma-b is super long
                       # it's not actually rho though
                       textInputRow("rho1", "Minimum:"),
                       tags$head(tags$style(type="text/css", "#rho1 {width: 75px}")),
                       textInputRow("rho2", "Maximum:"),
                       tags$head(tags$style(type="text/css", "#rho2 {width: 75px}")),
                       textInputRow("numval", "Number of values:", value="10"),
                       tags$head(tags$style(type="text/css", "#numval {width: 50px}")))
      ),
    #  actionButton("calc", "Calculate!")
    mainPanel(
      tabsetPanel(
        #         tabPanel("Numeric Output",
        #                  br(),
        #                  h4("Analytic Power:"),
        #                  textOutput("ap"),
        #                  br(),
        #                  h4("Approximate Power:"),
        #                  textOutput("dp"),
        #                  br(),
        #                  br()),
        tabPanel("Table",
                 dataTableOutput("table")),
        tabPanel("Plot",
                 plotOutput("plot")),
        tabPanel("Simulations",
                 textOutput("simulation"))
      )))))
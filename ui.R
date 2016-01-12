library(shiny)
library(shinyjs)
library(DT)

# this function allows text inputs to be on the same line; default has them on new lines
numericInputRow<-function (inputId, label, value) 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value, class="input-small"))
}

shinyUI(fluidPage(
  withMathJax(),
  useShinyjs(),
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
      actionLink("model","Model Equation"),
      hidden(div(id="modeltext", p('$y_{ij} = \\beta_0 + x_{ij}\\delta + b_{0i} +e_{ij} $'),
      tags$small('where $b_{0i} \\sim N(0, \\sigma_b^2)$, $e_{ij} \\sim N(0, \\sigma^2)$ 
                 and ICC = $\\sigma_b^2 / (\\sigma_b^2 + \\sigma^2)$'))),
      br(), br(),
      numericInput("N","Enter the number of subjects per cluster (N):",value=100), # number of subjects per clusters
      numericInput("M","Enter the number of clusters (M):",value=3), # number of clusters
      numericInput("d","Enter the difference $\\delta$ in mean between the two arms:",value=0.2), # difference in mean between the two arms
      # alternative hypothesis (null is 0)
      
      # radio button: ICC or sigma-b
      radioButtons("rhosigmab", "Use either an ICC or $\\sigma_b^2$ value:", list("ICC","$\\sigma_b^2$" = "sigmab"), inline = TRUE),
      conditionalPanel(condition = "input.rhosigmab == 'ICC'",
                       numericInputRow("ICC","Enter ICC value:",value=0.001)), # intra-cluster correlation coefficient
      conditionalPanel(condition = "input.rhosigmab == 'sigmab'",
                       numericInputRow("sigmab","Enter $\\sigma_b^2$ value:", value=0.001001)),
      numericInputRow("sigma","Enter $\\sigma^2$ value:",value=1),
      
      ## alpha can be one of the three options
      radioButtons("alpha","Enter alpha value:",c(0.025, 0.05,0.10),selected=0.05,inline=TRUE),
      
      ## options: range of rho/sigma-b, use CV (0 otherwise)
      checkboxGroupInput("options", label="Options:", choices = c("Cluster size is variable" = "useCV", "Use multiple values for the ICC or $\\sigma_b^2$"
                                                                  = "multipleICC")),
      conditionalPanel(condition = "input.options == 'useCV' || input.options.length > 1",
                       numericInput("CV", label="Enter the coefficient of variation (CV):", value=0)),
      conditionalPanel(condition = "input.options == 'multipleICC' || input.options.length > 1",
                       # calling it rho, but if sigmab is selected, it'll run as sigmab
                       numericInputRow("rho1", "Minimum:", value=0),
                       tags$head(tags$style(type="text/css", "#rho1 {width: 75px}")),
                       numericInputRow("rho2", "Maximum:", value=1),
                       tags$head(tags$style(type="text/css", "#rho2 {width: 75px}")),
                       numericInputRow("numval", "Number of values:", value=10),
                       tags$head(tags$style(type="text/css", "#numval {width: 50px}"))),
      actionButton("calc", "Calculate!")),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Table",
                 br(),
                 textOutput("tablefiller"),
                 DT::dataTableOutput("table"),
                 br(),
                 actionButton("save", "Save this Calculation"),
                 actionButton("clearall", "Clear All Saved Values"),
                 downloadButton("download",label="Download as CSV File"),
                 br(), br(),
                 p("Abbreviations:"),
                 helpText("Sigma-B: $\\sigma_b^2$"),
                 helpText("Sigma: $\\sigma^2$"),
                 helpText("ICC: Intra-cluster Correlation Coefficient"),
                 helpText("CV: Coefficent of Variation")),
        tabPanel("Plot",
                 plotOutput("plot")),
        tabPanel("Simulations",
                 br(),
                 numericInput("nsims","Number of Simulations to Run:", value=100),
                 actionButton("run", "Run"),
                 br(),br(),
                 p(textOutput("text")),
                 br(),
                  textOutput("simulation"))
      )))))
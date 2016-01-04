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

ui <- fluidPage(
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
      ))))

server <- function(input, output){

  
  ICC <- reactive({
    if (input$rhosigmab == 'ICC') as.numeric(input$ICC) else  {
      as.numeric(input$sigmab)/(as.numeric(input$sigma)+as.numeric(input$sigmab))
    }
  })
  
  CV <- reactive({
    if (length(input$options) == 0) 0 else {
    if ('useCV' %in% input$options) as.numeric(input$CV) else 0 }
  })
  
  DP <- reactive ({ # design effect power
    round(as.numeric(power.t.test(n = as.numeric(input$N)*as.numeric(input$M)/(1+(as.numeric(input$N)-1)*ICC()),
                                          delta=as.numeric(input$d), sig.level=as.numeric(input$alpha))$power),2)
    # approximate version of analytic power
    })
  
  AP <- reactive ({
    df <- 2*(as.numeric(input$M)-1) # degrees of freedom
    deff <- 1+((CV()^2+1)*as.numeric(input$N)-1)*ICC() # correction factor with CV
    lambda <- (as.numeric(input$d)/as.numeric(input$sigma)) / sqrt(2*deff/(as.numeric(input$M)*as.numeric(input$N)))
    nullq <- qt(as.numeric(input$alpha)/2, df, ncp=0)  # quantile of the alpha error, under the null, ncp=0
    analyticpower <- pt(nullq, df, lambda, lower.tail = TRUE) + pt(-nullq, df, lambda, lower.tail = FALSE)
    round(analyticpower,2)
  })
  
# if sigma-b isn't input
  SB <- reactive ({
    if (input$rhosigmab == 'ICC') (as.numeric(input$ICC)*as.numeric(input$sigma))/(as.numeric(input$ICC)+1) else {
    as.numeric(input$sigmab) }
  })
  
   
# values are reactive so I can use them wherever
  output$ap <- renderText(AP())
  output$dp <- renderText(DP())
  
# b in sigma-b means random effects
  
  # not just one input any more
  # input min and max
  # output a table of five values between
  # option: output a value or plot a range
  
  ### table output ###
  # difference n.clusters  n.per.cluster sigma.b sigma icc cv approx.power analytic.power
  icctable <- reactive({
  icctab <- data.frame(delta=as.numeric(input$d), M=as.numeric(input$M),N=as.numeric(input$N), 
                         SB=round(SB(),4), sigma=as.numeric(input$sigma), ICC=round(ICC(),4), CV=CV(),
                         DP=DP(), AP=AP())
  colnames(icctab) <- c("Difference", "Number of Clusters", "Number per Cluster", "$\\sigma_b^2$",
                          "$\\sigma^2$", "Intra-cluster Correlation Coefficient", "Coefficent of Variation",
                          "Approximate Power", "Analytic Power")
  icctab
  })
  output$table <- renderDataTable(icctable(),options=list(paging=FALSE,searching=FALSE,
                                                          ordering=0, processing=0, info=0))
  # options for rendertable: include.rownames=FALSE, digits=c(0,2,0,0,4,4,4,2,2,2)
  
  # mixed effects model -> used for cluster-randomized trials
  # add an extra term and some extra notation
  # assume that individuals are not independent in the model
  # regular sigma comes from the normal distribution of e_ij (residual error)
  # sigmab is a parameter of the normally distributed b_i (random effect)
  
  
  ### plot output ###
  # vector of analytic powers from this equation using Vectorize()
  # feed that to the plot function
  # !! approximate power plotted over on top !! - haven't done this yet
  output$plot <- renderPlot({
    rho.values <- seq(as.numeric(input$rho1),as.numeric(input$rho2),length.out=as.numeric(input$numval))
    power <- function(rho){
      df <- 2*(as.numeric(input$M)-1) # degrees of freedom
      ICC <- rho/(as.numeric(input$sigma)+rho)
      deff <- 1+((CV()^2+1)*as.numeric(input$N)-1)*ICC # correction factor with CV
      lambda <- (as.numeric(input$d)/as.numeric(input$sigma)) / sqrt(2*deff/(as.numeric(input$M)*as.numeric(input$N)))
      nullq <- qt(as.numeric(input$alpha)/2, df, ncp=0)  # quantile of the alpha error, under the null, ncp=0
      ap <- pt(nullq, df, lambda, lower.tail = TRUE) + pt(-nullq, df, lambda, lower.tail = FALSE)
      return(ap)
    }
    v.power <- Vectorize(power)
    apr.power <- function(rho){
      ICC <- rho/(as.numeric(input$sigma)+rho)
      power.t.test(n = as.numeric(input$N)*as.numeric(input$M)/(1+(as.numeric(input$N)-1)*ICC),
                              delta=as.numeric(input$d), sig.level=as.numeric(input$alpha))$power
    }
    v.apr.power <- Vectorize(apr.power)
    power.df <- data.frame(rho.values,v.power(rho.values),v.apr.power(rho.values))
    colnames(power.df) <- c("Values", "Power", "Approximate")
    ggplot() + geom_point(data=power.df,aes(x=Values,y=Power),color="blue") + 
      geom_point(data=power.df,aes(x=Values,y=Approximate),color="green")
    
  })
  
  
  ### simulations ###
  output$simulation <- renderText({

    p.try <- power.sim.normal(n.sim=100, effect.size=as.numeric(input$d), alpha=.05,
                              n.clusters=as.numeric(input$N), n.periods=1,
                              cluster.size=as.numeric(input$M),
                              period.effect = .7, period.var = 0,
                              btw.clust.var=as.numeric(input$sigmab), indiv.var=as.numeric(input$sigma),
                              verbose=TRUE,
                              estimation.function=random.effect)

    p.try$power
    binom.test(p.try$power*nsims, nsims)$conf.int[1:2]

  })
  
  }

shinyApp(ui = ui, server = server)




### original code ###
# normanal <- function(N,M,d,rho,sigma,alpha=.05) {
#  df <- 2*(M-1)
#  deff <- 1+(N-1)*rho
#  lambda <- (d/sigma) / sqrt(2*deff/(M*N))
#  nullq <- qt(alpha/2, df, ncp=0)  # quantile of the alpha error, under the null, ncp=0
#  analyticpower <- pt(nullq, df, lambda, lower.tail = TRUE) + pt(-nullq, df, lambda, lower.tail = FALSE)
#  deffpower <- power.t.test(n = N*M/deff, delta=d, sig.level=alpha)$power
#  return(list(ap=analyticpower, dp=deffpower))
#}
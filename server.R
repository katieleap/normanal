library(shiny)
library(clusterPower)
library(ggplot2)

shinyServer(function(input, output, session){
  observeEvent(input$model, {
    toggle("modeltext")
  })
  
  ICC <- eventReactive(input$calc, {
    if (input$rhosigmab == 'ICC') as.numeric(input$ICC) else  {
      as.numeric(input$sigmab)/(as.numeric(input$sigma)+as.numeric(input$sigmab))
    }
  })
  
  CV <- eventReactive(input$calc, {
    if (length(input$options) == 0) 0 else {
      if ('useCV' %in% input$options) as.numeric(input$CV) else 0 }
  })
  
  DP <- eventReactive(input$calc, { # design effect power
    round(as.numeric(power.t.test(n = as.numeric(input$N)*as.numeric(input$M)/(1+(as.numeric(input$N)-1)*ICC()),
                                  delta=as.numeric(input$d), sig.level=as.numeric(input$alpha))$power),2)
    # approximate version of analytic power
  })
  
  AP <- eventReactive(input$calc, {
    df <- 2*(as.numeric(input$M)-1) # degrees of freedom
    deff <- 1+((CV()^2+1)*as.numeric(input$N)-1)*ICC() # correction factor with CV
    lambda <- (as.numeric(input$d)/as.numeric(input$sigma)) / sqrt(2*deff/(as.numeric(input$M)*as.numeric(input$N)))
    nullq <- qt(as.numeric(input$alpha)/2, df, ncp=0)  # quantile of the alpha error, under the null, ncp=0
    analyticpower <- pt(nullq, df, lambda, lower.tail = TRUE) + pt(-nullq, df, lambda, lower.tail = FALSE)
    round(analyticpower,2)
  })
  
  # if sigma-b isn't input
  SB <- eventReactive(input$calc, {
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
  
  # calctable runs when we ask it to calculate
  # it also needs to append to anything that has already been saved, but without becoming permanent
  calctable <- eventReactive(input$calc, {
    calctab <- data.frame(delta=as.numeric(input$d), M=as.numeric(input$M),N=as.numeric(input$N), 
                         SB=round(SB(),4), sigma=as.numeric(input$sigma), ICC=round(ICC(),4), CV=CV(),
                         DP=DP(), AP=AP())
    colnames(calctab) <- c("Difference", "Number of Clusters", "Number per Cluster", "$\\sigma_b^2$",
                          "$\\sigma^2$", "Intra-cluster Correlation Coefficient", "Coefficent of Variation",
                          "Approximate Power", "Analytic Power")
    calctab
  })

  # savetable runs when we ask it to save and should append the new calculation to the old info
  # need a way to have calctable append to savetable if we run a new calculaton
  savetable <- eventReactive(input$save, {
    savetab <- data.frame(savetable(),calctable())
    savetab
  })
  
  # clearall makes an empty table and forgets everything we've done
  clearall <- eventReactive(input$clearall, {
    emptytab <- data.frame()
    colnames(emptytab) <- c("Difference", "Number of Clusters", "Number per Cluster", "$\\sigma_b^2$",
                           "$\\sigma^2$", "Intra-cluster Correlation Coefficient", "Coefficent of Variation",
                           "Approximate Power", "Analytic Power")
    emptytab
  })
  
  output$table <- renderDataTable(calctable(),options=list(paging=FALSE,searching=FALSE,
                                                           ordering=0, processing=0, info=0))

# we probably want a save as csv option, yeah?
  
  
  # options for rendertable: include.rownames=FALSE, digits=c(0,2,0,0,4,4,4,2,2,2)
  
  # mixed effects model -> used for cluster-randomized trials
  # add an extra term and some extra notation
  # assume that individuals are not independent in the model
  # regular sigma comes from the normal distribution of e_ij (residual error)
  # sigmab is a parameter of the normally distributed b_i (random effect)
  
  
  ### plot output ###
  # vector of analytic powers from this equation using Vectorize()
  # feed that to the plot function

  output$plot <- renderPlot({
    validate(
      need(input$options == 'multipleICC' || input$options.length > 1, "Input multiple values to plot!")
    )
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
    validate(
      need(input$nsims,"Enter the number of simulations!")
    )
    simulate()
  })
  simulate <- eventReactive(input$run, {
    withCallingHandlers({
      shinyjs::text("text", "")
    p.try <- suppressWarnings(power.sim.normal(n.sim=100, effect.size=as.numeric(input$d), alpha=.05,
                              n.clusters=as.numeric(input$N), n.periods=1,
                              cluster.size=as.numeric(input$M),
                              period.effect = .7, period.var = 0,
                              btw.clust.var=as.numeric(input$sigmab), indiv.var=as.numeric(input$sigma),
                              verbose=TRUE,
                              estimation.function=random.effect))
    p.try$power
    nsims <- as.numeric(input$nsims)
    binom.test(p.try$power*nsims, nsims)$conf.int[1:2]
    },
    message = function(m) {
      shinyjs::text(id = "text", text = m$message)
    })
  
   })
  
})




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
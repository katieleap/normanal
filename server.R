library(shiny)
library(clusterPower)
library(ggplot2)
library(DT)
library(reshape2)
library(plyr)

shinyServer(function(input, output, session){
  ## allows us to hide and show model text by clicking ##
  observeEvent(input$model, {
    toggle("modeltext")
  })

  ## allows us to hide and show disclaimer text by clicking ##
  observeEvent(input$disclaimer, {
    toggle("disclaimertext")
  })  
  
### reactive values that we want to use whenever ### 
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
  ### end reactive values ###
  
  
### text output ###
  output$ap <- renderText(AP())
  output$dp <- renderText(DP())
  
  
### table output ###
  # difference n.clusters  n.per.cluster sigma.b sigma icc cv approx.power analytic.power
  
  # calctable runs when we ask it to calculate
  # it also needs to append to anything that has already been saved, but without becoming permanent
  t <- reactiveValues(data=NULL)
  s <- reactiveValues(data=NULL)
  power <- function(rho){
    df <- 2*(as.numeric(input$M)-1) # degrees of freedom
    deff <- 1+((CV()^2+1)*as.numeric(input$N)-1)*rho # correction factor with CV
    lambda <- (as.numeric(input$d)/as.numeric(input$sigma)) / sqrt(2*deff/(as.numeric(input$M)*as.numeric(input$N)))
    nullq <- qt(as.numeric(input$alpha)/2, df, ncp=0)  # quantile of the alpha error, under the null, ncp=0
    ap <- pt(nullq, df, lambda, lower.tail = TRUE) + pt(-nullq, df, lambda, lower.tail = FALSE)
    return(ap)
  }
  apr.power <- function(rho){
    power.t.test(n = as.numeric(input$N)*as.numeric(input$M)/(1+(as.numeric(input$N)-1)*rho),
                 delta=as.numeric(input$d), sig.level=as.numeric(input$alpha))$power
  }
  sbpower <- function(rho){
    df <- 2*(as.numeric(input$M)-1) # degrees of freedom
    ICC <- rho/(as.numeric(input$sigma)+rho)
    deff <- 1+((CV()^2+1)*as.numeric(input$N)-1)*ICC # correction factor with CV
    lambda <- (as.numeric(input$d)/as.numeric(input$sigma)) / sqrt(2*deff/(as.numeric(input$M)*as.numeric(input$N)))
    nullq <- qt(as.numeric(input$alpha)/2, df, ncp=0)  # quantile of the alpha error, under the null, ncp=0
    ap <- pt(nullq, df, lambda, lower.tail = TRUE) + pt(-nullq, df, lambda, lower.tail = FALSE)
    return(ap)
  }
  sbapr.power <- function(rho){
    ICC <- rho/(as.numeric(input$sigma)+rho)
    power.t.test(n = as.numeric(input$N)*as.numeric(input$M)/(1+(as.numeric(input$N)-1)*ICC),
                 delta=as.numeric(input$d), sig.level=as.numeric(input$alpha))$power
  }
  icccalc <- function(rho){rho/(as.numeric(input$sigma)+rho)}
  sbcalc <- function(rho){(rho*as.numeric(input$sigma))/(rho+1)}
  
  observeEvent(input$calc, {
    if((input$options == 'useCV' & length(input$options) == 1) || length(input$options) == 0 ) {
      t$data <- data.frame(delta=as.numeric(input$d), M=as.numeric(input$M),N=as.numeric(input$N), 
                           SB=round(SB(),4), sigma=as.numeric(input$sigma), ICC=round(ICC(),4), CV=CV(),
                           DP=DP(), AP=AP())
    } else {
      if(input$rhosigmab == 'ICC'){
        rho.values <- list(seq(as.numeric(input$rho1),as.numeric(input$rho2),length.out=as.numeric(input$numval)))
        rap <- laply(rho.values,power)
        rdp <- laply(rho.values,apr.power)
      t$data <- data.frame(delta=rep_len(as.numeric(input$d),length.out=as.numeric(input$numval)), 
                           M=rep_len(as.numeric(input$M),length.out=as.numeric(input$numval)),
                           N=rep_len(as.numeric(input$N),length.out=as.numeric(input$numval)), 
                           SB= laply(rho.values,sbcalc), 
                           sigma=rep_len(as.numeric(input$sigma),length.out=as.numeric(input$numval)), 
                           ICC=rho.values, CV=rep_len(CV(),length.out=as.numeric(input$numval)),
                           DP=rdp, AP=rap)
      colnames(t$data)[6] <- "ICC"
      } else {
        rho.values <- seq(as.numeric(input$rho1),as.numeric(input$rho2),length.out=as.numeric(input$numval))
        rap <- laply(rho.values,sbpower)
        rdp <- laply(rho.values,sbapr.power)
        t$data <- data.frame(delta=rep_len(as.numeric(input$d),length.out=as.numeric(input$numval)), 
                             M=rep_len(as.numeric(input$M),length.out=as.numeric(input$numval)),
                             N=rep_len(as.numeric(input$N),length.out=as.numeric(input$numval)), 
                             SB=rho.values, 
                             sigma=rep_len(as.numeric(input$sigma),length.out=as.numeric(input$numval)), 
                             ICC= laply(rho.values,icccalc), CV=rep_len(CV(),length.out=as.numeric(input$numval)),
                             DP=rdp, AP=rap)
        colnames(t$data)[4] <- "SB"
      }
    }
    })

  # savetable runs when we ask it to save and should append the new calculation to the old info
  # need a way to have calctable append to savetable if we run a new calculaton
  observeEvent(input$save, {
    s$data <- rbind(s$data,t$data)
    t$data <- NULL
  })
  
  # clearall makes an empty table and forgets everything we've done
  observeEvent(input$clearall, {
    t$data <- NULL
    s$data <- NULL
  })
  
  output$tablefiller <- renderText({
    validate(
      need(input$calc, 'Press the calculate button to create a table!'))
  })
    
  output$table <- DT::renderDataTable(rbind(s$data,t$data),options=list(paging=FALSE,searching=FALSE,
                                                           ordering=0, processing=0, info=0),
                                      class='compact hover row-border nowrap',
                                      colnames = c("Difference", "Number of Clusters", "Number per Cluster", "Sigma-B",
                                                              "Sigma", "ICC", "CV",
                                                              "Approximate Power", "Analytic Power"))
# we probably want a save as csv option, yeah?
  output$download <- downloadHandler(
    filename=function(){
      paste('poweranalysis', Sys.Date(),'.csv',sep='')
    },
    content = function(file) {
      write.csv(rbind(s$data,t$data),file)
    }
  )
  
  
### plot output ###
  output$plot <- renderPlot({
    validate(
      need(input$options == 'multipleICC' || length(input$options) > 1, "Input multiple values to plot!")
    )
    sigbrange <- function(){
    rho.values <- seq(as.numeric(input$rho1),as.numeric(input$rho2),length.out=as.numeric(input$numval))
    v.power <- Vectorize(sbpower)
    v.apr.power <- Vectorize(sbapr.power)
    power.df <- data.frame(rho.values,v.power(rho.values),v.apr.power(rho.values))
    colnames(power.df) <- c("Values", "Analytic", "Approximate")
    melt.power <- melt(power.df, id="Values")
    ggplot(data=melt.power,aes(x=Values, y = value, color = variable)) + geom_point() +
      labs(list(y="Power")) + scale_color_discrete(name="Power")
    }
    iccrange <- function(){
      rho.values <- seq(as.numeric(input$rho1),as.numeric(input$rho2),length.out=as.numeric(input$numval))
      v.power <- Vectorize(power)
      v.apr.power <- Vectorize(apr.power)
      power.df <- data.frame(rho.values,v.power(rho.values),v.apr.power(rho.values))
      colnames(power.df) <- c("Values", "Analytic", "Approximate")
      melt.power <- melt(power.df, id="Values")
      ggplot(data=melt.power,aes(x=Values, y = value, color = variable)) + geom_point() +
        labs(list(y="Power")) + scale_color_discrete(name="Power")
    }
    
    if (input$rhosigmab == 'ICC') iccrange() else sigbrange()
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
    conf <- binom.test(p.try$power*nsims, nsims)$conf.int[1:2]
    paste("Confidence Interval:", paste(round(conf,3),collapse="-"),sep=" ")
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
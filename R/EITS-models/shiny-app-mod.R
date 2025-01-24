#shinySIR: Interactive plotting for infectious disease models
#Sinead E. Morris (2020). shinySIR: Interactive Plotting for Mathematical Models of Infectious Disease Spread. R package version 0.1.2.
#https://CRAN.R-project.org/package=shinySIR #documentation
#https://github.com/SineadMorris/shinySIR #source code

#Last modified: 27 Nov 2024 by Madeline Lewis 
#Modified to allow for specification of compartments to not be plotted. The environment E can be large, 
#which precludes visualization of compartments representing individuals' states.

library(shinySIR)

#library(
#  dplyr,
#  ggplot2,
#  deSolve,
#  tidyr
#)

plot_model <- function(output, linesize, textsize, xlabel, ylabel, legend_title, levels, values, ...){
  
  output$variable <- factor(output$variable, levels = levels)
  
  ggplot(output, aes(x = time, y = value, colour = as.factor(variable))) +
    geom_line(size = linesize) +
    scale_colour_manual(legend_title, values = values, ...) +
    ylab(ylabel) + xlab(xlabel) +
    theme_bw() + theme(axis.text = element_text(size = textsize),
                       axis.title= element_text(size = textsize + 2),
                       legend.text = element_text(size = textsize),
                       legend.title = element_text(size = textsize + 2) )
  
}

solve_eqns <- function(eqns, ics, times, parms, remove_cols=c()){
  #print(ics)
  #print(parms)
  
  trySolve <- tryCatch(deSolve::lsoda(y = ics,
                                      times = times,
                                      func = eqns,
                                      parms = parms),
                       error = function(e) e,
                       warning = function(w) w)
  
  if (inherits(trySolve, "condition")) {
    print(paste("deSolve error:", trySolve$message))
    stop("ODE solutions are unreliable. Check model attributes e.g. equations, parameterization, and initial conditions.")
  } else {
    soln <- deSolve::lsoda(y = ics,
                           times = times,
                           func = eqns,
                           parms = parms)
  }
  data <- data.frame(soln)
  data[, remove_cols] <- list(NULL)
  head(data)
  output <- data %>% #tbl_df() %>%
    tidyr::gather(variable, value, 2:ncol(.))
  return(output)
}

run_shiny <- function(model = "EITS Model", neweqns = NULL,
                      ics = NULL,
                      tstart = 0, timestep = 1, tmax = 365,
                      parm0 = NULL,
                      parm_names = NULL,
                      parm_min = NULL,
                      parm_max = NULL,
                      sigfigs = 4,
                      showtable = TRUE,
                      linesize = 1.2, textsize = 14,
                      xlabel = "Time", ylabel = "Number of individuals",
                      legend_title = "Compartment",
                      slider_steps = NULL,
                      values = NULL, 
                      remove_cols = c(), ...
){
  
  # Get eqns & display name
  name <- model
  
  if (exists(model) & is.null(neweqns)) {
    print('using existing model')
    eqns <- get(model, mode = "function")
    name <- get_name(model)
  } else if (exists(model) & !is.null(neweqns)) {
    print("Your model name matches one of the built-in models. You can rename it using the 'model' argument.")
    eqns <- neweqns
  } else if (!exists(model) & is.null(neweqns)) {
    stop("Model name not recognized. Try one of the built-in models or specify your own with the 'neweqns' argument.")
  } else {
    print('using custom models')
    eqns <- neweqns
  }
  
  # Get parameters for built-in models
  if (is.null(parm0) & is.null(neweqns)) {
    params <- get_params(model)
    
    parm0 <- params$parm0
    parm_names <- params$parm_names
    
    parm_min <- params$parm_min
    parm_max <- params$parm_max
  }
  
  # Check parameter vectors when user-specified model is defined
  if (!is.null(neweqns) & ( is.null(parm0))) {
    stop("Missing parameter vector 'parm0'")
  }
  
  if (!is.null(neweqns) & ( is.null(parm_min))) {
    warning("Missing parameter vector 'parm_min': using 0.5 * parm0 as default.")
    parm_min <- parm0 * 0.5
  }
  
  if (!is.null(neweqns) & ( is.null(parm_max))) {
    warning("Missing parameter vector 'parm_max': using 1.5 * parm0 as default.")
    parm_max <- parm0 * 1.5
  }
  
  if (is.null(names(parm_min)) | is.null(names(parm_max)) | is.null(names(parm0))) {
    stop("parm0, parm_min, and parm_max must be named vectors.")
  }
  
  if (!is.null(neweqns) & ( is.null(parm_names))) {
    parm_names <- names(parm0)
    warning("Could not find names of parameters for interactive menu ('parm_names'). Using names of 'parm0' instead.")
  }
  
  # Check parameters appear in the same order in all vectors
  if ( !( all(sapply(list(names(parm0), names(parm_min), names(parm_max)), function(x) x == names(parm0)))) ){
    stop("The parameters in parm0, parm_min, and parm_max must have the same names, and appear in the same order.")
  }
  
  # Check parameter values
  if (any(parm_min >= parm_max)) {
    stop("All entries in parm_min must be less than their corresponding entries in parm_max.")
  }
  
  if (any(parm0 > parm_max) | any(parm0 < parm_min)) {
    warning("All entries in parm0 should be within the bounds of their corresponding entries in parm_min, parm_max.")
  }
  
  parm0 <- signif(parm0, sigfigs)
  parm_min <- signif(parm_min, sigfigs)
  parm_max <- signif(parm_max, sigfigs)
  
  
  # Get initial conditions
  if (is.null(ics) & is.null(neweqns)) {
    ics <- get_ics(model)
  } else if (is.null(ics) & !is.null(neweqns)) {
    stop("You must specify initial conditions for your own model using the 'ics' argument.")
  }
  if (is.null(names(ics))) {
    stop("ics must be a named vector.")
  }
  
  # get default ggplot colours
  if (is.null("values")) {
    gghues <- seq(15, 375, length = length(ics) + 1)
    values <- hcl(h = gghues, l = 65, c = 100)[1:length(ics)]
  }
  
  if (length(values) != (length(ics))) {
    warning("The length of the manual colour scale vector ('values') must equal the number of model variables. Using default ggplot colours instead.")
    
    gghues <- seq(15, 375, length = length(ics) + 1)
    values <- hcl(h = gghues, l = 65, c = 100)[1:length(ics)]
  }
  
  # User Interface (UI)
  ui <- pageWithSidebar(
    headerPanel(name),
    sidebarPanel(
      lapply(seq_along(parm0),
             function(x) numericInput(
               inputId = names(parm0)[x], 
               label = parm_names[x],
               value = parm0[x], 
               min = parm_min[x], 
               max = parm_max[x]
               #step = slider_steps[x]
               )
      )
    ),
    mainPanel(
      plotOutput("plot1"),
      br(), br(), br(),
      tableOutput("table1"),
    )
  )
  
  
  # Behind the scenes code (Server)
  server <- function(input, output){
    
    # Get initial population size (doesn't change with user input)
    START.N <- as.numeric(sum(ics))
    
    output$plot1 <- renderPlot({
      
      # Get parameters from user input
      parms_vector <- unlist(reactiveValuesToList(input))
      parms_vector <- c(parms_vector, N = START.N)
      
      # Time vector (total length input from shiny interface)
      times_vector <- seq(from = tstart, to = tmax, by = timestep)

      # Run ODE solver
      ODEoutput <- solve_eqns(eqns, ics, times = times_vector, parms = parms_vector, remove_cols = remove_cols)
      
      # Plot output
      plot_model(ODEoutput, linesize, textsize, xlabel, ylabel, legend_title, levels = names(ics), values, ...)
      
    })
    
    output$table1 <- renderTable({
      
      parms_vector <- unlist(reactiveValuesToList(input))
      parms_vector <- c(parms_vector, N = START.N)
      
      if ((model %in% c("SIR", "SIS")) & showtable == TRUE) {
        data.frame(
          Parameter = c("gamma", "beta"),
          Value = c(1/parms_vector["Ip"],
                    parms_vector["R0"] * (1/parms_vector["Ip"]) / parms_vector["N"])
        )
      } else if ((model %in% c("SIRS")) & showtable == TRUE) {
        data.frame(
          Parameter = c("gamma", "delta", "beta"),
          Value = c(1/parms_vector["Ip"],
                    1/parms_vector["Rp"],
                    parms_vector["R0"] * (1/parms_vector["Ip"]) / parms_vector["N"])
        )
      } else if ((model %in% c("SIRbirths", "SISbirths", "SIRvaccination")) & showtable == TRUE) {
        data.frame(
          Parameter = c("gamma", "beta"),
          Value = c(1/parms_vector["Ip"],
                    parms_vector["R0"] * (1/parms_vector["Ip"] + parms_vector["mu"]) / parms_vector["N"])
        )
      } else if ((model %in% c("SIRSbirths", "SIRSvaccination")) & showtable == TRUE) {
        data.frame(
          Parameter = c("gamma", "delta", "beta"),
          Value = c(1/parms_vector["Ip"],
                    1/parms_vector["Rp"],
                    parms_vector["R0"] * (1/parms_vector["Ip"] + parms_vector["mu"]) / parms_vector["N"])
        )
      }
      
    }, digits = -1, bordered = TRUE)
    
  }
  
  shinyApp(ui, server)
}

#END

library(shiny) # build the shiny app
library(shinyjs)
library(visNetwork) # multistate structure diagram
library(mstate) # data prep
library(DT) # show the table with function datatable(data, ...)
library(eha) # fit PH model with piecewise constant function for baseline hazard function
library(flexsurv) # fit PH model with natural cubic spline for baseline log cumulative hazard function
library(ggplot2) # for plots
library(scales)  # for plots
library(hesim)   # use params_surv() to store the fitted models
library(dplyr)
library(doParallel)

# load all the helper functions
r_files <- list.files(
  path = "helper_functions",
  pattern = "\\.R$",
  full.names = TRUE
)
invisible(lapply(r_files, source, local = FALSE))

# ui
ui <- fluidPage(
  
  theme = "bootstrap.css",
  shinyjs::useShinyjs(),
  withMathJax(),
  titlePanel("Toolkit for Multistate Disease Progression Simulation and Treatment Decision-Making Aid"),
  
  # Custom box style 
  tags$head(
    tags$style(HTML("
      .box {
        background: #ffffff;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 16px;
        margin-bottom: 16px;
        box-shadow: 0 2px 4px rgba(0,0,0,.05);
      }
      .box h4 {
        margin-top: 0;
        font-weight: 600;
      }
    "))
  ),
  
  tags$head(
    tags$style(HTML("
    .input-box { background: #f8f9fa; border-left: 4px solid #2b7bba; }
    .output-box { background: #ffffff; border-left: 4px solid #28a745; }
  "))
  ),
  
  tabsetPanel(
    id = "tabs",
    
    # ---- Tab 1: Multistate Structure ----
    tabPanel(
      "Multistate Structure",
      h5(HTML("In the section of Multistate Structure, you need to <br/>
      (1) define the number of states and the state names in the multistate model,<br/>
      (2) specify all possible transitions, and<br/>
      (3) mark which transitions are affected by the treatment.<br/>
              The diagram and transition matrix of the multistate model as the output will update automatically based on your selections.<br/>
              The treatment considered for treatment stratgy comparsion needs to be included as one of the states.")),
      navlistPanel(
        
        id = "MSMstruc",
        widths = c(2, 10),       
        
        tabPanel("States and transitions",
                 fluidRow(
                   column(
                     width = 3,
                     div(class = "box input-box",
                         h4("Define states"),
                         p(class = "text-muted",
                           "Define the number of states in the multistate model and provide custom names if desired."),
                         numericInput("numStates", "Number of States:", min = 3, max = 20, value = 3),
                         uiOutput("stateInputs")
                     )
                   ),
                   column(
                     width = 5,
                     div(class = "box input-box",
                         h4("Select all the possible transitions"),
                         uiOutput("transitionInputs")
                     )
                   ),
                   column(
                     width = 4,
                     div(class = "box input-box",
                         h4("Select transitions which are affected by the treatment"),
                         uiOutput("interventionTransitionsUI")
                     )
                   )
                 ),
                 br(),
                 div(style = "display: flex; justify-content: flex-end; margin-top: 20px;",
                     actionButton("goToMSMView", "Next: View the Transition Diagram and Matrix", class = "btn btn-primary")
                 )
        ),
        tabPanel("Transition diagram and matrix",
                 
                 fluidRow(
                   column(
                     width = 7,
                     div(class = "box output-box",
                         h4("Diagram of the multistate model"),
                         visNetworkOutput("transitionDiagram", height = "560px", width = "100%"),
                         div(class = "muted",
                             tags$em("Note: Blue edges represent normal transitions; black edges indicate treatment transitions.")
                        )
                     )
                   ),
                   column(
                     width = 5,
                     div(class = "box output-box",
                         h4("Transition Matrix"),
                         tableOutput("transitionMatrix")
                     )
                   )
                 ),
                 br(),
                 div(style = "display: flex; justify-content: flex-end; margin-top: 20px;",
                     actionButton("goToFitting", "Next: Multistate Model Fitting", class = "btn btn-primary")
                 )
        )
      )
    ),
    
    # ---- Tab 2: Multistate Fitting ----
    tabPanel(
      "Multistate Modeling",
      id="MSMFitting",
      tagList(
        h4("Multistate Modeling"),
        tags$p("In this section, we build the multistate model using flexible 
          parametric proportional hazards survival models 
          (Royston–Parmar models) for transitions not affected by the treatment."),
        tags$p("You can upload your dataset in the required format. If you do not have a real dataset, you can input the estimation 
          results of the multistate model manually. Plots of baseline hazard functions and hazard ratios are also presented.")
      ),
      
      navlistPanel(
        
        id="multi_fit",
        widths = c(2, 10),       
        
        tabPanel("Dataset", 
                 
                 h4("Would you like to upload a dataset for multistate modeling?"),
                 radioButtons(
                   "has_data", 
                   label = NULL,
                   choices = c("Yes" = "yes", "No" = "no"),
                   selected = character(0),
                   inline = TRUE) ,
                 conditionalPanel(
                   condition = "input.has_data == 'yes'",
                   h5("Please upload a dataset for multistate modeling in the CSV format."),
                   h5("The instructions for the data format are as follows:"),
                   tags$ul(
                     tags$li("Each row = one patient"),
                     tags$li(
                       "Required columns (case-insensitive matching is supported):",
                       tags$ul(
                         tags$li("Unique patient identifier (e.g., id / ID / Id)."),
                         tags$li("Time column per state: Column names must match the state names you defined 
               in the multistate structure (e.g., State1, State2, State3). Each column is the time of entering that state."),
                         tags$li("Status column per state: <state>.s (e.g., State1.s, State2.s, State3.s).
                                     Value should be 1 if the patient entered the corresponding state, 0 otherwise."),
                         tags$ul(tags$li("For a patient's status column(s) equal to 0, the time column(s) for the corresponding states should be the last observed time of the patient.")),
                         tags$li("Covariates: any additional columns for the covariates measured at baseline (time 0). Categorical covariates should be dummy-coded with a reference level.")
                       )
                     )
                   ),
                   uiOutput("uploadData"),
                   br(),
                   uiOutput("msdataTitle"),
                   DTOutput("msdataDT", height = "400px"),
                   br(),
                   fluidRow(
                     column(
                       width = 6,
                       uiOutput("transitionCountTitle"),
                       DTOutput("transitionCountDT")
                     )
                   ),
                   br(),
                   # Add one more dataset 
                   uiOutput("addSectionUI"),
                   br(),
                   actionButton("goToCovars", "Next: Covariate Assignment",
                                class = "btn btn-primary")
                 ),
                 # Show guidance if NO
                 conditionalPanel(
                   condition = "input.has_data == 'no'",
                   h5("Please proceed to the next subsection to input values for the model."),
                   tags$p("You will be asked to define covariate effects and baseline hazard functions manually."),
                   br(),
                   actionButton("goToCovars", "Next: Covariate Assignment", class = "btn btn-primary")
                 )
        ),
        
        tabPanel("Covariate Assignment",
                 conditionalPanel(
                   condition = "input.has_data == 'yes' || input.has_data == 'no'",
                   conditionalPanel(condition="input.has_data=='yes'",
                                    h4("Given the uploaded dataset, define covariates for transitions not affected by the treatment")),
                   conditionalPanel(condition="input.has_data=='no'",
                                    h4("Define covariates for transitions not affected by the treatment")),
                   uiOutput("covariateBlocks")
                 )
        ),
        
        tabPanel("Model Specification",
                 id="model_spec",
                 conditionalPanel(
                   condition = "input.has_data == 'yes' || input.has_data == 'no'",
                   h4("The flexible parametric proportional hazards models are fitted for all transitions not affected by the treatment."),
                   
                   radioButtons(
                     inputId = "baseline_hazard",
                     label = div(
                       style = "white-space: nowrap;",
                       "What model would you like to use to fit the baseline hazard function?"
                     ),
                     choices = c(
                       "Piecewise constant function" = "piecewise",
                       "Natural cubic spline (Royston–Parmar)" = "spline"
                     ),
                     selected = "piecewise",
                     inline = FALSE
                   ),
                   
                   # Conditional options for piecewise constant basleine hazard function
                   conditionalPanel(
                     condition = "input.baseline_hazard == 'piecewise'&input.has_data == 'yes'",
                     helpText("Set the maximum number of time intervals and the number of events in each interval for the piecewise constant function. 
                              The number of time intervals for each transition-specific model is obtained from the total number of events divided by the assigned number of events. 
                              If the obtained number of time intervals is greater than the assigned maximum number of intervals, then the maximum number of intervals will be used and the number
                              of events in each interval automatically increases accordingly."),
                     numericInput("max_intervals",
                                  "Maximum number of intervals:",
                                  value = 20, min = 1, max = 20, step = 1),
                     numericInput("num_events",
                                  "Number of events in each interval:",
                                  value = 20, min = 20, step = 1)
                   ),
                   
                   
                   conditionalPanel(
                     condition = "input.baseline_hazard == 'piecewise'",
                     br(),
                     conditionalPanel(
                       condition = "input.has_data == 'yes'",
                       actionButton("goToHR", "Next: Fit Transition-Specific Models", class = "btn-primary")
                     ),
                     conditionalPanel(
                       condition = "input.has_data == 'no'",
                       actionButton("goToHR", "Next: Specify Model Parameter Values", class = "btn-primary")
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.baseline_hazard == 'spline'",
                     br(),
                     conditionalPanel(
                       condition = "input.has_data == 'yes'",
                       actionButton("goToHR", "Next: Fit Transition-Specific Models", class = "btn-primary")
                     ),
                     conditionalPanel(
                       condition = "input.has_data == 'no'",
                       actionButton("goToHR", "Next: Specify Model Parameter Values", class = "btn-primary")
                     )
                   )
                   
                 )
        ),
        
        tabPanel("Parameter Specification",
                 
                 tabsetPanel(
                   id="param_spec",
                   tabPanel("Covariates",
                            conditionalPanel(
                              condition = "input.has_data == 'yes' || input.has_data == 'no'",
                              h4("Hazard Ratios of Covariates"),
                              conditionalPanel(
                                condition = "input.has_data == 'yes' ",
                                h5("The hazard ratios (HRs) of the covariates under transitions are estimated from the flexible parametric proportional hazards models. 
                                             They can be adjusted as needed.")),
                              conditionalPanel(
                                condition = "input.has_data == 'no'",
                                h5("Define the hazard ratios (HRs) of the covariates under transitions. The defaults are 1")
                              ),
                              h5(tags$i("Hazard ratios of covariates are the exponential of the corresponding parameter coefficients. The hazard ratios are greater than 0.", br(),
                                        "For example, a hazard ratio > 1 indicates that the greater the value of the covariate, the larger the hazards it has.",br(),
                                        "A hazard ratio < 1 indicates that the greater the value of the covariate, the smaller the hazards it has.")),
                              uiOutput("hr_inputs"),
                              br(),
                              div(style = "display: flex; justify-content: flex-end; margin-top: 20px;",
                                  actionButton("goToBaseHaz", "Next: Parameters in the baseline hazard functions", class = "btn btn-primary")
                              )
                            )
                   ),
                   tabPanel("Baseline Hazard Functions",
                            conditionalPanel(condition="input.baseline_hazard == 'spline'",
                                             conditionalPanel(
                                               condition = "input.has_data == 'yes' || input.has_data == 'no'",
                                               h4("Parameters of Baseline Cumulative Hazard Functions"),
                                               h5("Define the parameter coefficients of the baseline cumulative hazard function for each transition."),
                                               h5(tags$i("In the Royston–Parmar proportional hazards model, the logrithm of the baseline cumulative hazard function was modeled by a natural cubic spline.
                                                         While adjusting the numeric inputs of gammas and knots, the baseline hazard functions do not have negative values")),
                                               
                                             ),
                                             conditionalPanel(
                                               condition = "input.has_data == 'yes' ",
                                               h5(tags$i("The total number of knots are two boundary knots plus the number of interior knots. 
                                                   The number of interior knots is automatically determined during model fitting with the smallest the Bayesian information criteria. "))
                                             )
                            ),
                            uiOutput("baseline_hazards_inputs")
                   )
                 ),
        )
      )
    ),

    # ---- Tab 3: Treatment Strategies ----
    tabPanel(
      "Treatment Strategies",
      h4("Predict the restricted mean survival times (RMSTs) under different treatment strategies for a patient with given covariate values using the microsimulation method."),
      tagList(
        fluidRow(
          column(
            width = 12,
            numericInput(
              inputId = "rmst_time",
              label = strong("Time Horizon for Restricted Mean Survival Time"),
              value = 15,
              min = 0,
              step = 5
            ),
            numericInput(
              inputId = "num_sample_microsim",
              label = strong("Number of simulated patients for the microsimulation"),
              value = 1000,
              min = 0,
              step = 100
            ),
          )
       ),
       tags$hr(),
     ),
     uiOutput("interventionStrategiesUI"),
     tagList(
       tags$style(HTML("
    .btn-blue-1 { background-color:#0B2E4B; border-color:#0B2E4B; color:#fff; }
    .btn-blue-2 { background-color:#0D47A1; border-color:#0D47A1; color:#fff; }
    .btn-blue-3 { background-color:#1976D2; border-color:#1976D2; color:#fff; }
    .btn-blue-4 { background-color:#64B5F6; border-color:#64B5F6; color:#0B2E4B; }
    .btn-blue-5 { background-color:#BBDEFB; border-color:#BBDEFB; color:#0B2E4B; }

    .btn-blue-1:hover, .btn-blue-2:hover, .btn-blue-3:hover, .btn-blue-4:hover, .btn-blue-5:hover {
      filter: brightness(0.95);
    }
  ")),
       
       fluidRow(
         column(
           width = 12,
           h4("Start a new analysis with one of the following options:"),
           
           actionButton("reset_all",       "1) Re-run from the beginning",         class = "btn btn-blue-1"),
           actionButton("reset_data",      "2) Update the data",                   class = "btn btn-blue-2"),
           actionButton("reset_covariate", "3) Change covariate sets",             class = "btn btn-blue-3"),
           actionButton("reset_model",     "4) Change the model",                  class = "btn btn-blue-4"),
           actionButton("reset_micro",     "5) Reset the microsimulation setting", class = "btn btn-blue-5")
         )
       )
     )
     ),
   
   # ---- References ----
   tabPanel("References",
             h4("Royston, P. and Parmar, M. (2002). Flexible parametric proportional-hazards and proportionalodds models for censored survival data, with application to prognostic modelling and estimation of treatment effects. Statistics in Medicine 21(1):2175-2197."))
  ),
  
  br(),
  hr(
    HTML("Version 1.0 ---- &#169; The University of Texas MD Anderson Cancer Center"),
    p(
      "App developed by Can Xie, Xuelin Huang, Ruosha Li, and Nicholas Short",
      br(),
      a(href = "https://github.com/cxie19/TxMicroSim", "View source code at GitHub")
    )
  )
)

server <- function(input, output, session){
  # ---- Tab 1: Multistate Structure ----
  # reactiveValues 
  rv1 <- reactiveValues()
  
  rv1$stateNames <- reactive({
    n <- input$numStates %||% 3
    nm <- vapply(seq_len(n), function(i) input[[paste0("n", i)]] %||% paste("State", i), character(1))
    nm <- trimws(nm)
    nm[nm == ""] <- paste0("State ", which(nm == ""))
    make.unique(nm, sep = " ")
  })
  
  rv1$transitionEdges <- reactiveVal(data.frame(from=character(), to=character(), label=character(), stringsAsFactors=FALSE))
  rv1$interventionSel <- reactiveVal(character(0))
  rv1$transitionMatrix <- reactiveVal(NULL)
  
  # state name inputs
  output$stateInputs <- renderUI({
    n <- input$numStates %||% 3
    lapply(seq_len(n), function(i) {
      textInput(paste0("n", i), paste("State", i, "name"), value = paste0("State", i))
    })
  })
  
  # transition checkboxes: exclude self; default = only states after i
  output$transitionInputs <- renderUI({
    states <- rv1$stateNames()
    n <- length(states)
    if (n <= 1) return(NULL)
    
    lapply(seq_len(n), function(i) {
      to_choices <- setdiff(states, states[i])
      default_after <- if (i < n) states[(i + 1L):n] else character(0)
      
      checkboxGroupInput(
        inputId  = paste0("trans_from_", i),
        label    = paste0("From: ", states[i]),
        choices  = setNames(to_choices, paste0("To: ", to_choices)),
        selected = default_after,   # always reset to "later states"
        inline   = TRUE
      )
    })
  })
  
  
  # build transition ID matrix + edges from inputs; render matrix & intervention UI
  observe({
    req(rv1$stateNames())
    states <- rv1$stateNames()
    n <- length(states)
    
    id_mat <- matrix(NA_integer_, nrow = n, ncol = n, dimnames = list(states, states))
    edges  <- data.frame(from=character(), to=character(), label=character(), stringsAsFactors=FALSE)
    counter <- 1L
    
    for (i in seq_len(n)) {
      from_state <- states[i]
      to_states  <- input[[paste0("trans_from_", i)]] %||% character(0)
      to_states  <- setdiff(to_states, from_state)   # no self loops
      to_states  <- intersect(to_states, states)     # drop stale labels
      if (length(to_states)) {
        for (to_state in to_states) {
          id_mat[from_state, to_state] <- counter
          edges <- rbind(edges, data.frame(from = from_state, to = to_state,
                                           label = as.character(counter), stringsAsFactors = FALSE))
          counter <- counter + 1L
        }
      }
    }
    
    rv1$transitionEdges(edges)
    rv1$transitionMatrix(id_mat)
    output$transitionMatrix <- renderTable(id_mat, rownames = TRUE)
    
    if (nrow(edges)) {
      prev <- rv1$interventionSel()
      choice_names <- paste0(edges$label, ": ", edges$from, " \u2192 ", edges$to)
      keep_prev    <- intersect(prev, edges$label)
      output$interventionTransitionsUI <- renderUI({
        checkboxGroupInput(
          "intervention_transitions",
          label   = NULL,
          choices = setNames(edges$label, choice_names),
          selected = keep_prev
        )
      })
    } else {
      output$interventionTransitionsUI <- renderUI(NULL)
      rv1$interventionSel(character(0))
    }
  })
  
  observeEvent(input$intervention_transitions, {
    rv1$interventionSel(input$intervention_transitions %||% character(0))
  })
  
  observeEvent(input$goToMSMView, {
    updateNavlistPanel(
      session,
      inputId = "MSMstruc",   # the id of your tabsetPanel/navlistPanel
      selected = "Transition diagram and matrix"
    )
  })
  
  # network diagram 
  output$transitionDiagram <- renderVisNetwork({
    req(rv1$stateNames())
    nodes <- circle_layout_custom(
      labels      = rv1$stateNames(),
      radius      = 420,
      start_angle = 3*pi/4,   # top-left
      clockwise   = TRUE      # move clockwise
    )
    
    edges <- rv1$transitionEdges()
    
    if (nrow(edges)) {
      edges$color  <-  "#1f77b4"  # blue
      edges$width  <- 3
      edges$arrows <- "to"
      edges$smooth <- FALSE
      sel <- rv1$interventionSel()
      if (length(sel)) {
        idx <- which(edges$label %in% sel)
        edges$color[idx] <- "#000000"; edges$width[idx] <- 3
      }
    }
    
    visNetwork(nodes, edges, height = "560px") %>%
      visNodes(
        shape       = "circle",
        size        = 56,
        scaling     = list(min = 56, max = 56),
        borderWidth = 3,
        fixed       = list(x = TRUE, y = TRUE),
        font        = list(size = 22, face = "arial", vadjust = 0),
        color       = list(background = "#e8f1fb", border = "#2b7bba")
      ) %>%
      visEdges(arrows = "to", font = list(size = 26)) %>%
      visOptions(highlightNearest = TRUE,
                 nodesIdSelection = list(enabled = TRUE, main = "Select by state")) %>%
      visPhysics(enabled = FALSE)
  })
  
  
  observeEvent(input$goToFitting, {
    # rv_gate$data_ready <- TRUE
    updateTabsetPanel(
      session,
      inputId = "tabs",           # Outer tabsetPanel ID
      selected = "Multistate Modeling"  # The title of the outer tabPanel to switch to
    )
    
    updateTabsetPanel(
      session,
      inputId = "multi_fit",           
      selected = "Dataset"  
    )
  })
  
  # for later use
  states <- reactive(rv1$stateNames()) # state names, e.g., CP, AP, BP, HSCT, and Death
  tm     <- reactive(rv1$transitionMatrix()) # transition matrix
  edges <- reactive(rv1$transitionEdges()) # a dataframe with columns from, to, and label (1,2,3,..)
  tx_trans <- reactive(rv1$interventionSel()) # the sequence number for the transition affected by the treatment
  non_tx_trans <- reactive(setdiff(edges()$label,tx_trans())) # the sequence number for the transition not affected by the treatment
  tr_labels <- reactive(setNames(   # transition labels
    paste0(edges()$label, ": ", edges()$from, " → ", edges()$to),
    edges()$label
  ))
 
  # ----  Tab 2: Multistate Modeling ---- 
  rv_gate <- reactiveValues(#data_ready = FALSE,
    model_ready = FALSE,
    micro_ready = FALSE)
  rv2 <- reactiveValues()
  # reactiveValues
  rv2_1 <- reactiveValues()
  rv2_1$firstTransformed <- NULL  # the frozen, first transformed dataset
  rv2_1$extraList <- list()  # list of added transformed datasets
  rv2_1$combined <- reactiveVal(NULL)    # view that shows first + extras
  
  # covariate assignment
  rv2_2 <- reactiveValues()
  rv2_2$blockCountRV    <- reactiveVal(1)    # number of blocks
  rv2_2$blockCovNamesRV <- reactiveVal(list())  # covariate names (typed or selected)
  rv2_2$blockTargetsRV  <- reactiveVal(list())  # selected transitions
  
  # parameter specification
  rv2_4 <- reactiveValues()
  
  gamma_counts <- reactiveValues()
  knot_counts  <- reactiveValues()
  observer_flags <- reactiveValues()
  
  # ----  Tab 2-1: dataset ---- 
  output$uploadData <- renderUI({
    tagList(
      fileInput("dataFile", "Upload your dataset (CSV)", accept = ".csv"),
      DTOutput("uploadedDataDT", height = "400px")
    )
  })
  
  # If data is uploaded, read in the data 
  observeEvent(input$dataFile, {
    req(input$dataFile)
    # read data
    df <- read.csv(input$dataFile$datapath, check.names = FALSE)
    # load data
    output$uploadedDataDT <- renderDT({
      datatable(df, options = list(scrollX = TRUE, scrollY = 300))
    })
    # Check state name consistency 
    df_names <- names(df)
    
    if (is.null(tm())) {
      showModal(modalDialog(
        title = "Error: Transition Matrix Not Found",
        "Please define the states and transitions in the Multistate Structure tab before uploading data.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
      return()
    }
    
    required_time   <- states()
    required_status <- paste0(states(), ".s")
    
    # Lowercase matching
    df_names_lower <- tolower(df_names)
    required_time_lower <- tolower(required_time)
    required_status_lower <- tolower(required_status)
    
    missing_time <- required_time[!required_time_lower %in% df_names_lower]
    missing_status <- required_status[!required_status_lower %in% df_names_lower]
    
    if (length(missing_time) > 0 || length(missing_status) > 0) {
      msg <- paste(
        "The uploaded dataset does not have matching state names or required columns.",
        "\nMissing time columns:", paste(missing_time, collapse = ", "),
        "\nMissing status columns:", paste(missing_status, collapse = ", "),
        "\nPlease ensure the dataset follows the required format."
      )
      showModal(modalDialog(
        title = "Dataset Format Error",
        pre(msg), easyClose = TRUE, footer = modalButton("OK")
      ))
      return()
    }
    
    #  Proceed with msprep if valid 
    id_col <- resolve_single_ci(df_names, "id")
    time_cols   <- resolve_cols_ci(df_names, states())
    status_cols <- resolve_cols_ci(df_names, paste0(states(), ".s"))
    
    msdata <- tryCatch(
      mstate::msprep(
        time   = time_cols,
        status = status_cols,
        data   = df,
        trans  = tm(),
        id     = id_col,
        keep   = setdiff(df_names, c(id_col, time_cols, status_cols))
      ),
      error = function(e) {
        showNotification(paste("msprep failed:", e$message), type = "error")
        NULL
      }
    )
    
    if (is.null(msdata)) return()
    msdata <- msdata[msdata$time != 0, ]
    
    # Save the first dataset's transformed data (frozen)
    rv2_1$firstTransformed <- msdata
    # Initialize the combined view as the first dataset
    rv2_1$combined(msdata)
    
    # Display transformed data
    output$msdataTitle <- renderUI(h4("Transformed dataset for multistate modeling (first dataset)"))
    output$msdataDT <- renderDT({
      req(rv2_1$firstTransformed)
      datatable(rv2_1$firstTransformed, options = list(scrollX = TRUE, scrollY = 300))
    })
    
    # Display counts
    trans_counts <- updateTransitionCounts(rv2_1$firstTransformed, edges(), tm())
    output$transitionCountTitle <- renderUI(h4("Number of Patients per Transition"))
    output$transitionCountDT <- renderDT({
      datatable(trans_counts[, c("label", "n_patients")],
                colnames = c("Transition", "Number of Patients"),
                options = list(scrollX = TRUE, paging = FALSE, dom = "t"))
    })
    
    output$addSectionUI <- renderUI({
      tagList(
        br(),
        actionButton("addDataBtn", "Add one more data for specific transition(s)", class = "btn-primary"),
        conditionalPanel(
          condition = "input.addDataBtn > 0",
          uiOutput("addDataUI"),
          uiOutput("combinedDataTitle"),
          DTOutput("combinedDataDT", height = "400px"),
          fluidRow(
            column(
              width = 6,
              uiOutput("combinedCountTitle"),
              DTOutput("combinedCountDT")
            )
          ),
          br(),br(),
          h4("If there is more dataset to be uploaded, please scroll up for the button to add more dataset.")
        )
      )
    })
  })
  
  # Add one more data 
  observeEvent(input$addDataBtn, {
    
    output$addDataUI <- renderUI({
      fluidRow(
        column(
          width = 5,
          checkboxGroupInput(
            "transitionSelect",
            "Assign to transition(s):",
            choices = paste0(edges()$from, " → ", edges()$to), 
            inline = FALSE
          )
        ),
        column(
          width = 4,
          fileInput(
            "extraDataFile",
            "Upload an additional dataset (CSV) with the same column names as the previous one",
            accept = ".csv"
          )
        ),
        column(
          width = 3,
          div(
            style = "margin-top: 25px;",
            actionButton("combineBtn", "Combine this dataset", class = "btn-success")
          )
        )
      )
    })
  })
  
  #  Combine handler 
  observeEvent(input$combineBtn, {
    # 
    req(rv2_1$firstTransformed, input$extraDataFile, input$transitionSelect)
    
    df_extra <- read.csv(input$extraDataFile$datapath, check.names = FALSE)
    
    id_col      <- resolve_single_ci(names(df_extra), "id")
    time_cols   <- resolve_cols_ci(names(df_extra), states())
    status_cols <- resolve_cols_ci(names(df_extra), paste0(states(), ".s"))
    
    msdata_extra <- tryCatch(
      mstate::msprep(
        time   = time_cols,
        status = status_cols,
        data   = df_extra,
        trans  = tm(),
        id     = id_col,
        keep   = setdiff(names(df_extra), c(id_col, time_cols, status_cols))
      ),
      error = function(e) {
        showNotification(paste("msprep failed:", e$message), type = "error")
        NULL
      }
    )
    if (is.null(msdata_extra)) return()
    
    # Filter by selected transitions
    label_df <- data.frame(
      trans = seq(nrow(edges())),
      label = paste0(edges()$from, " → ",edges()$to))
    target_trans <- label_df$trans[label_df$label %in% input$transitionSelect]
    if (length(target_trans) == 0) {
      showNotification("Please select at least one transition.", type = "warning")
      return()
    }
    msdata_extra <- msdata_extra[msdata_extra$trans%in%target_trans,]
    
    # Append and rebuild combined from the frozen first dataset
    rv2_1$extraList[[length(rv2_1$extraList) + 1]] <- msdata_extra
    combined_data <- do.call(rbind, c(list(rv2_1$firstTransformed), rv2_1$extraList))
    rv2_1$combined(combined_data)
    
    # Update combined displays
    output$combinedDataTitle <- renderUI(h4("Combined dataset (first + added)"))
    output$combinedDataDT <- renderDT({
      datatable(rv2_1$combined(), options = list(scrollX = TRUE, scrollY = 300))
    })
    
    # counts for transitions
    trans_counts_combined <- updateTransitionCounts(rv2_1$combined(), edges(), tm())
    output$combinedCountTitle <- renderUI(h4("Number of Patients per Transition (combined data)"))
    output$combinedCountDT <- renderDT({
      datatable(trans_counts_combined[, c("label", "n_patients")],
                colnames = c("Transition", "Number of Patients"),
                options = list(scrollX = TRUE, paging = FALSE, dom = "t"))
    })
    
    showNotification(paste("✅ Added data for:", paste(input$transitionSelect, collapse = ", ")),
                     type = "message")
    
  })
  
  observeEvent(input$goToCovars, {
    updateNavlistPanel(
      session,
      inputId = "multi_fit",   
      selected = "Covariate Assignment"
    )
  })
  
  # for later use
  dat <- reactive(rv2_1$combined())
  
  # ---- Tab 2-2: Covariates ----
  nBlocks <- reactive(rv2_2$blockCountRV())
  savedCovs <- reactive(rv2_2$blockCovNamesRV())
  savedTrns <- reactive(rv2_2$blockTargetsRV())
  
  observe({
    if (nBlocks() <= 0) return()
    elig_edges <- edges()[edges()$label %in% non_tx_trans(), , drop = FALSE]
    if (!nrow(elig_edges)) return(NULL)
    
    choices_all <- setNames(
      elig_edges$label,
      paste0(elig_edges$label, ": ", elig_edges$from, " → ", elig_edges$to)
    )
    
    for (i in seq_len(nBlocks())) {
      local({
        idx <- i
        selectAllId  <- paste0("select_all_", idx)
        transInputId <- paste0("transitions_block_", idx)
        
        observeEvent(input[[selectAllId]], {
          if (isTRUE(input[[selectAllId]])) {
            updateCheckboxGroupInput(session, transInputId, selected = choices_all)
          } else {
            updateCheckboxGroupInput(session, transInputId, selected = character(0))
          }
        }, ignoreInit = TRUE)
      })
    }
  })
  
  
  observe({
    cov_list <- list()
    trn_list <- list()
    for (i in seq_len(nBlocks())) {
      val <- input[[paste0("covars_block_", i)]]
      if (!is.null(val) && length(val)==1 && grepl(",", val)) {
        val <- trimws(unlist(strsplit(val, ",")))
      }
      cov_list[[as.character(i)]] <- val %||% character(0)
      trn_list[[as.character(i)]] <- input[[paste0("transitions_block_", i)]] %||% character(0)
    }
    rv2_2$blockCovNamesRV(cov_list)
    rv2_2$blockTargetsRV(trn_list)
  })
  
  observeEvent(input$add_block, {rv2_2$blockCountRV(rv2_2$blockCountRV()+1)})
  
  output$covariateBlocks <- renderUI({
    
    if (nBlocks() <= 0) return(NULL)
    
    elig_edges <- edges()[edges()$label %in% non_tx_trans(), , drop = FALSE]
    if (!nrow(elig_edges)) return(NULL)
    
    choices_all <- setNames(
      elig_edges$label,
      paste0(elig_edges$label, ": ", elig_edges$from, " → ", elig_edges$to)
    )
    
    #already_assigned <- unique(unlist(savedTrns()))
    usedTransitions <- character(0)
    
    ui_blocks <- lapply(seq_len(nBlocks()), function(i) {
      sel_cov <- savedCovs()[[as.character(i)]] %||% character(0)
      sel_trn <- savedTrns()[[as.character(i)]] %||% character(0)
      # 
      # if (isTRUE(rv_gate$covariate_ready)) {
      #   sel_cov <- character(0)
      #   sel_trn <- character(0)
      # }
      # 
      availableChoices <- choices_all[!(choices_all %in% usedTransitions)]
      sel_trn <- intersect(sel_trn, availableChoices)
      usedTransitions <<- c(usedTransitions, sel_trn)
      
      # Covariate input depends on has_data
      if (input$has_data == "yes" && !is.null(dat())) {
        all_vars <- names(dat())
        excluded <- c("id","from","to","trans","Tstop","status","time",
                      states(), paste0(states(), ".s"))
        filter_cov <- setdiff(all_vars, excluded)
        sel_cov <- if (length(sel_cov) == 0) filter_cov else sel_cov
        covUI <- selectizeInput(
          paste0("covars_block_", i), "Select covariates",
          choices = filter_cov, selected = sel_cov, multiple = TRUE,
          options = list(closeAfterSelect = FALSE, openOnFocus = TRUE)
        )
      } else {
        covUI <- textAreaInput(
          paste0("covars_block_", i),
          "Enter covariate names separated by commas. For categorical covariates, write their dummy variable names.",
          placeholder = "e.g., sex, age_45_64, age_ge65",
          value = paste(sel_cov, collapse = ", ")
        )
      }
      
      # Select-all logic 
      transInputId <- paste0("transitions_block_", i)
      selectAllId  <- paste0("select_all_", i)
      selectAllValue <- length(sel_trn) == length(availableChoices)
      
      wellPanel(
        h5(paste("Covariate set", i)),
        covUI,
        checkboxGroupInput(
          transInputId,
          "Select transitions",
          choices = availableChoices,
          selected = sel_trn
        ),
        checkboxInput(selectAllId, "All transitions above", value = selectAllValue)
      )
    })
    
    # Dynamic disable button logic ---
    all_assigned <- all(elig_edges$label %in% usedTransitions)
    
    tagList(
      ui_blocks,
      # Build the "control" section (add button / disabled message)
      control_ui <- if (!all_assigned ) {
        div(style="margin-top:6px;",
            actionButton("add_block", "Add another covariate set",
                         class = "btn-outline-secondary"))
      } else {
        # The button "add one more covariate set" is hidden if all available transitions are selected
        tagList(div(style="margin-top:10px; color:#888;",
                    tags$em("All transitions are already selected in one covariate set — new blocks are disabled.")))
        # Extract all covariate names
        covariates <- unique(unlist(savedCovs()))
        desc_inputs <- lapply(covariates, function(cov) {
          textInput(
            inputId = paste0("desc_", cov),
            label = paste("Description for covariate:", cov),
            placeholder = paste("e.g.,", cov),
            value = cov  # default to the covariate name
          )
        })
        
        tagList( 
          div(style="margin-top:10px; font-style: italic; color:#555;",
              "All transitions not affected by the treatment have been assigned with covariate sets."
          ),
          br(),
          h5("Please provide a description for each covariate."),
          br(),
          div(desc_inputs),
          br(),
          actionButton("save_cov_desc", "Next: Model Specification", class = "btn btn-primary")
        )
      } 
    )
  })
  
  observeEvent(input$save_cov_desc, {
    updateNavlistPanel(
      session,
      inputId = "multi_fit",   
      selected = "Model Specification"
    )
  })
  
  # ---- Tab 2-3: Model Specification ----
  observeEvent(c(input$goToHR), {
    rv_gate$model_ready <- TRUE
    updateNavlistPanel(
      session,
      inputId = "multi_fit",   
      selected = "Parameter Specification"
    )
    
    updateTabsetPanel(
      session,
      inputId = "param_spec",   
      selected = "Covariates"   
    )
  })
  
  # ---- Tab 2-4 Parameter specification ---- 
  fittedModelsRV <- eventReactive(input$goToHR, {
    req(rv_gate$model_ready)
    req(dat())
    fits <- list()
    
    withProgress(message = "Fitting flexible parametric PH models...", value = 0, {
      n <- length(non_tx_trans())
      for (j in seq_along(non_tx_trans())) {
        tr <- non_tx_trans()[j]
        incProgress(1/n, detail = paste("Fitting model for transition", tr))
        
        df_tr <- subset(dat(), trans == as.integer(tr))
        if (nrow(df_tr) == 0) next
        
        assign_covs <- unlist(savedCovs()[sapply(savedTrns(), function(x) tr %in% x)])
        formula_str <- paste("Surv(time, status) ~", paste(assign_covs, collapse = " + "))
        f <- as.formula(formula_str)
        
        # Branch by model type
        if (input$baseline_hazard == "piecewise") {
          cuts <- make_event_intervals(df_tr, input$num_events, input$max_intervals,ceiling(max(dat()$time)))
          fit <- tryCatch(
            eha::pchreg(f, data = df_tr, cuts = cuts),
            error = function(e) { cat("Piecewise fit error:", e$message, "\n"); NULL }
          )
        } else if (input$baseline_hazard == "spline") {
          fit <- tryCatch(
            fit_flexsurvspline(df_tr, covariates = assign_covs),
            error = function(e) { cat("Spline fit error:", e$message, "\n"); NULL }
          )
        } else {
          fit <- NULL
        }
        
        if (!is.null(fit)) fits[[tr]] <- fit
      }
    })
    fits
  })
  
  observeEvent(input$has_data, {
    req(rv_gate$model_ready)
    if (input$has_data == "no") {
      fittedModelsRV <- reactiveVal(NULL)
      cat("Cleared fitted models due to has_data='no'\n")
    }
  })
  
  # ---- Tab 2-4: Hazard ratios ----
  # Build hazard ratio inputs & plots side-by-side
  output$hr_inputs <- renderUI({
    req(rv_gate$model_ready)
    req(input$has_data)
    
    # covariate → transitions mapping
    cov_tr_map <- list()
    for (i in seq_along(savedCovs())) {
      for (cov in savedCovs()[[i]]) {
        cov_tr_map[[cov]] <- unique(c(cov_tr_map[[cov]], savedTrns()[[i]]))
      }
    }
    
    # Step 1: Fit models once per transition
    coef_map <- list()
    if (input$has_data == "yes" && !is.null(dat())) {
      fits <- fittedModelsRV()
      coef_map <- lapply(fits, function(m) exp(coef(m)))
    }else if(input$has_data == "no" ){
      fits <- NULL
    }
    
    # Step 2: Determine global HR range (so vertical line aligns)
    all_HRs <- c()
    for (cov in names(cov_tr_map)) {
      trns <- cov_tr_map[[cov]]
      for (tr in trns) {
        if (!is.null(coef_map[[tr]]) && cov %in% names(coef_map[[tr]])) {
          all_HRs <- c(all_HRs, coef_map[[tr]][cov])
        } else {
          all_HRs <- c(all_HRs, 1)
        }
      }
    }
    hr_min <- min(all_HRs, na.rm = TRUE)
    hr_max <- max(all_HRs, na.rm = TRUE)
    
    # Step 3: Build UI per covariate
    ui_list <- list()
    
    for (cov in names(cov_tr_map)) {
      trns <- cov_tr_map[[cov]]
      if (length(trns) == 0) next
      
      # defaults
      default_HRs <- setNames(rep(1, length(trns)), trns)
      for (tr in trns) {
        if (!is.null(coef_map[[tr]]) && cov %in% names(coef_map[[tr]])) {
          default_HRs[tr] <- coef_map[[tr]][cov]
        }
      }
      
      trns_sorted <- trns[order(as.numeric(trns))]
      
      # user-defined description if available
      cov_label <- if (!is.null(input[[paste0("desc_", cov)]])) input[[paste0("desc_", cov)]] else cov
      
      # UI block for each covariate 
      cov_box <- wellPanel(
        style = "background-color: transparent; border: 1px solid #dee2e6; padding: 15px; border-radius: 8px;",
        h5(strong(paste("Hazard Ratios of", cov_label))),
        fluidRow(
          column(
            width = 3,
            tagList(lapply(trns_sorted, function(tr) {
              numericInput(
                inputId = paste0("HR_", cov, "_", tr),
                label   = tr_labels()[tr],
                value   = round(default_HRs[tr], 3),
                min     = 0,
                step    = 0.1, width = "120px"
              )
            }))
          ),
          column(
            width = 9,
            plotOutput(
              paste0("hr_plot_", cov),
              height = paste0(80 * length(trns_sorted), "px")
            )
          )
        )
      )
      
      ui_list[[length(ui_list)+1]] <- cov_box
      
      # render one combined forest plot per covariate 
      local({
        cov_local <- cov
        trns_local <- trns_sorted
        
        output[[paste0("hr_plot_", cov_local)]] <- renderPlot({
          # Pull user-entered HRs
          vals <- sapply(trns_local, function(tr) {
            input[[paste0("HR_", cov_local, "_", tr)]] %||% NA
          })
          
          # Initialize data for plotting
          df <- data.frame(
            transition = tr_labels()[trns_local],
            HR = vals,
            lower = NA_real_,
            upper = NA_real_,
            stringsAsFactors = FALSE
          )
          
          # Get model CIs if available
          if (input$has_data == "yes" && !is.null(fittedModelsRV())) {
            fits <- fittedModelsRV()
            
            for (tr in trns_local) {
              model <- fits[[tr]]
              
              if (!is.null(model) && cov_local %in% names(coef(model))) {
                coef_index <- which(names(coef(model)) == cov_local)
                est <- exp(coef(model)[coef_index])
                se <- if (input$baseline_hazard=="spline"){sqrt(vcov(model)[coef_index, coef_index])} else{sqrt(model$var[coef_index, coef_index])}
                ci_lower <- est * exp(-1.96 * se)
                ci_upper <- est * exp(1.96 * se)
                
                label <- tr_labels()[tr]
                user_val <- input[[paste0("HR_", cov_local, "_", tr)]]
                
                # Only show CI if user hasn't changed HR (rounded match)
                if (!is.null(user_val) && abs(user_val - round(est, 3)) < 1e-3) {
                  df[df$transition == label, "lower"] <- ci_lower
                  df[df$transition == label, "upper"] <- ci_upper
                }
              }
            }
          }
          
          # Ensure transition is a factor in correct order
          df$transition <- factor(df$transition, levels = rev(tr_labels()[trns_local]))
          
          # Plot
          ggplot(df, aes(x = HR, y = transition)) +
            geom_point(size = 4, color = "steelblue", na.rm = TRUE) +
            geom_errorbarh(
              aes(xmin = lower, xmax = upper),
              height = 0.2,
              na.rm = TRUE,
              color = "gray40"
            ) +
            geom_vline(xintercept = 1, linetype = "dashed") +
            scale_x_log10(
              labels = label_number(accuracy = 0.01), 
              #limits = c(hr_min, ),
              breaks = extended_breaks(n = 8)
            ) +
            theme_minimal(base_size = 12) +
            theme(
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x  = element_text(size = 14),  
              axis.text.y = element_text(size = 14),   
              axis.ticks.y = element_blank()
            )+ theme(
              plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 5.5)  # increase r
            )
          
        })
      })
    }
    
    do.call(tagList, ui_list)
  })
  
  # Collect edited HR values into data frame
  all_HRs <- reactive({
    df_list <- list()
    
    for (i in seq_along(savedCovs())) {
      covs <- savedCovs()[[i]]
      trns <- savedTrns()[[i]]
      for (cov in covs) {
        for (tr in trns) {
          input_id <- paste0("HR_", cov, "_", tr)
          val <- input[[input_id]]
          if (!is.null(val)) {
            df_list[[length(df_list)+1]] <- data.frame(
              transition = tr,
              covariate = cov,
              HR = val
            )
          }
        }
      }
    }
    if (length(df_list)) do.call(rbind, df_list) else NULL
  })
  
  observeEvent(input$goToBaseHaz, {
    updateNavlistPanel(
      session,
      inputId = "param_spec",   # the id of your tabsetPanel/navlistPanel
      selected = "Baseline Hazard Functions"
    )
  })
  
  # ---- Tab 2-4: Baseline hazard functions ----
  output$baseline_hazards_inputs <- renderUI({
    req(rv_gate$model_ready)
    req(input$has_data)
    
    # Determine default max time 
    max_time_default <- if (input$has_data == "yes" &&
                            !is.null(dat())) {
      ceiling(max(dat()$time, na.rm = TRUE))
    } else 20
    
    # --- Handle user adjustment vs reset on has_data change ---
    # If has_data changed, reset to default; otherwise, preserve previous input
    prev_xmax <- isolate(input$baseline_xmax)
    prev_has_data <- isolate(rv2$prev_has_data)  # store last known has_data
    
    if (is.null(prev_xmax) || is.null(prev_has_data) || input$has_data != prev_has_data) {
      xmax_value <- max_time_default
    } else {
      xmax_value <- prev_xmax
    }
    
    # Store current has_data state for next render
    rv2$prev_has_data <- input$has_data
    
    #  UI for max x- and y-axises
    if(input$has_data=="no"){
      max_x_set=NA
    }else if (input$has_data=="yes"){
      max_x_set= ceiling(max(dat()$time))
    }
    
    ui_list <<- list()
    ui_list[1] <<- list(
      fluidRow(
        column(
          width = 4,
          numericInput(
            "baseline_xmax",
            "Max x-axis value for time:",
            value = xmax_value,
            min = 1, step = 1, max = max_x_set
          )
        ),
        column(
          width = 4,
          numericInput(
            "baseline_ymax_h",
            "Max y-axis value for baseline hazard function:",
            value = 1, min = 1, step = 1, max = 50
          )
        ),
        column(
          width = 4,
          numericInput(
            "baseline_ymax_H",
            "Max y-axis value for baseline cumulative hazard function:",
            value = 5, min = 1, step = 1, max = 100
          )
        )
      )
    )
    
    ymax_h <- suppressWarnings(as.numeric(input$baseline_ymax_h))
    ymax_H <- suppressWarnings(as.numeric(input$baseline_ymax_H))
    
    xmax_value <- suppressWarnings(as.numeric(input$baseline_xmax))
    if (is.na(xmax_value) || length(xmax_value) == 0) {
      xmax_value <- if (!is.null(dat())) max(dat()$time, na.rm = TRUE) else 20
    }
    
    # get fits if available
    fits <- if (input$has_data == "yes") fittedModelsRV() else NULL
    
    # --- piecewise constant function ---
    if (input$baseline_hazard == "piecewise") {
      
      for (tr in non_tx_trans()) {
        local({
          tr_local <- tr
          fit_local <- if (!is.null(fits)) fits[[tr_local]] else NULL
          
          # Default cuts & hazards
          default_cuts <- if (!is.null(fit_local)) fit_local$cuts else seq(0, xmax_value, length.out = 5)
          default_haz  <- if (!is.null(fit_local)) fit_local$hazards else rep(0.01, length(default_cuts) - 1)
          
          # UI container
          ui_list[[length(ui_list) + 1]] <<- fluidRow(
            column(
              width = 4,
              h5(strong(paste("Transition", tr_labels()[tr_local]))),
              if (input$has_data == "no") {
                tagList(
                  numericInput(paste0("n_intervals_", tr_local),
                               "Number of intervals:",
                               value = length(default_haz),
                               min = 1, step = 1))
              } else {
                tagList(
                  paste0("The number of intervals is ", length(default_haz), "."),
                  br())
              },
              uiOutput(paste0("interval_inputs_", tr_local))
            ),
            column(width = 8, plotOutput(paste0("piecewise_plot_", tr_local), height = "300px")),
            tags$div(style = "height:30px;")
          )
          
          # Dynamic inputs
          output[[paste0("interval_inputs_", tr_local)]] <- renderUI({
            req(input$has_data)
            cuts <- default_cuts
            haz  <- default_haz
            # dynamically respond to number of intervals when has_data == "no"
            if (input$has_data == "no") {
              n_raw <- input[[paste0("n_intervals_", tr_local)]]
              
              validate(
                need(!is.null(n_raw), "Please choose number of intervals."),
                need(length(n_raw) == 1, "n_intervals must be a single number."),
                need(is.finite(n_raw) && n_raw >= 1, "n_intervals must be >= 1.")
              )
              
              n_intervals <- as.integer(n_raw)
              cuts <- seq(0, xmax_value, length.out = n_intervals + 1)
              haz <- rep(0.01, n_intervals)
            } else {
              n_intervals <- length(default_haz)
              cuts <- if (!is.null(fit_local)) fit_local$cuts else default_cuts
              haz <- if (!is.null(fit_local)) fit_local$hazards else default_haz
            }
            
            tagList(
              if (input$has_data == "no") {
                tagList(
                  h6(strong("Interval cut points:")),
                  lapply(seq_len(n_intervals + 1), function(i) {
                    numericInput(paste0("cut_", tr_local, "_", i),
                                 label = paste("Cut", i, ":"),
                                 value = round(cuts[i], 2), min = 0, step = 1)}))
              },
              h6(strong("Hazard rates (λ) for intervals:")),
              lapply(seq_len(n_intervals), function(i) {
                if (input$has_data == "yes"){
                  numericInput(paste0("haz_", tr_local, "_", i),
                               label = paste0("λ for", " time interval ",i, " [", round(cuts[i],2), ", ", round(cuts[i+1],2), ")"),
                               value = round(haz[i],4), min = 0, step = 0.01)
                }else{
                  numericInput(paste0("haz_", tr_local, "_", i),
                               label = paste0("λ for", " time interval ",i," [Cut ", i, ", Cut", i+1, ")"),
                               value = round(haz[i],4), min = 0, step = 0.01)
                }
              })
            )
          })
          
          pwc_inputs <- reactive({
            # determine interval count
            n_intervals <- if (input$has_data == "no") {
              input[[paste0("n_intervals_", tr_local)]]
            } else {
              length(default_haz)
            }
            req(n_intervals)
            
            # collect hazard inputs
            haz_list <- sapply(seq_len(n_intervals), function(i) {
              val <- input[[paste0("haz_", tr_local, "_", i)]]
              if (is.null(val)) NA_real_ else as.numeric(val)
            })
            haz_vals <- as.numeric(unlist(haz_list))
            
            # fallback for missing hazards
            if (length(haz_vals) == 0 || all(is.na(haz_vals))) {
              haz_vals <- if (!is.null(fit_local)) fit_local$hazards else default_haz
            }
            
            # collect cuts
            cuts_list <- lapply(seq_len(n_intervals + 1), function(i) {
              val <- input[[paste0("cut_", tr_local, "_", i)]]
              if (is.null(val)) NA_real_ else as.numeric(val)
            })
            cuts_vals <- as.numeric(unlist(cuts_list))
            
            # fallback for missing cuts
            if (length(cuts_vals) == 0 || all(is.na(cuts_vals))) {
              cuts_vals <- if (!is.null(fit_local)) fit_local$cuts else default_cuts
            }
            
            list(cuts = cuts_vals, haz = haz_vals)
          })
          
          # Plot
          output[[paste0("piecewise_plot_", tr_local)]] <- renderPlot({
            
            req(input$baseline_hazard, input$has_data) #rv2$ui_ready, 
            vals <- pwc_inputs()
            req(vals$cuts, vals$haz)
            
            # handle partially filled inputs gracefully
            if (anyNA(vals$cuts) || anyNA(vals$haz)) return(NULL)
            H_vals <- numeric(length(vals$cuts))
            for (i in seq_len(length(vals$cuts) - 1)) {
              H_vals[i + 1] <- H_vals[i] + (vals$cuts[i + 1] - vals$cuts[i]) * vals$haz[i]
            }
            # must satisfy length relation
            if (length(vals$cuts) != length(vals$haz) + 1) return(NULL)
            plot_piecewise_hazard(vals$cuts, vals$haz, tr_local, xmax_value,
                                  ylim_haz = ymax_h,
                                  ylim_cumhaz = ymax_H)
          })
        })
      }
    }
    
    
    # --- Natural cubic spline ---
    if(input$baseline_hazard == "spline"){
      
      for (tr in non_tx_trans()) {
        
        log_min <- log(0.1)
        log_max <- log(xmax_value)
        
        if (input$has_data == "no") {
          default_gammas <- c(-1, 1.5, 0.01)
          default_knots <- c(log_min, (log_min + log_max) / 2, log_max)
        }else if (input$has_data == "yes") {
          best_fit <- fits[[tr]]
          if (is.null(best_fit)) next
          default_gammas <- coef(best_fit)[grep("gamma", names(coef(best_fit)))]
          default_knots  <- best_fit$knots
        }
        
        if (is.null(gamma_counts[[tr]])) gamma_counts[[tr]] <- length(default_gammas)
        if (is.null(knot_counts[[tr]]))  knot_counts[[tr]]  <- length(default_knots)
        
        #  Current counts
        n_gammas <- gamma_counts[[tr]]
        n_knots  <- knot_counts[[tr]]
        
        # Gamma inputs 
        gamma_inputs <- tagList(
          lapply(seq_len(n_gammas), function(i) {
            id <- paste0("gamma_", tr, "_", i - 1)
            val <- if (i <= length(default_gammas)) round(default_gammas[i], 3) else 0.01
            label_text <- if (i == 1) {
              withMathJax("\\(\\gamma_{0}\\) (intercept)")
            } else if (i == 2) {
              withMathJax("\\(\\gamma_{1}\\) (for log(t))")
            } else {
              withMathJax(sprintf("\\(\\gamma_{%d}\\) (for knots)", i - 1))
            }
            numericInput(id, label = label_text, value = val, step = 0.1, width = "100%")
          })
        )
        
        K_internal <- max(0, n_knots - 2)
        internal <- if (K_internal > 0) {
          log_min + (1:K_internal)/(K_internal + 1) * (log_max - log_min)
        } else numeric(0)
        
        knots_current <- if (input$has_data=="no") c(log_min, internal, log_max) else default_knots
        
        # knot inputs 
        knot_inputs <- tagList(
          lapply(seq_along(knots_current), function(i) {
            lbl <- if (i == 1 || i == length(knots_current)) {
              paste("Knot", i, "(boundary)")
            } else paste("Knot", i, "(internal)")
            numericInput(
              paste0("knot_", tr, "_", i),
              label = lbl,
              value = round(knots_current[i], 3),
              step = 0.1,
              width = "100%"
            )
          })
        )
        
        # The add button 
        add_btn <- if (input$has_data=="no")
          actionButton(paste0("add_gamma_knot_", tr),
                       "Add Gamma & Knot",
                       class = "btn btn-sm btn-outline-primary")
        else NULL
        
        # Layout
        ui_list[[length(ui_list) + 1]] <<- fluidRow(
          column(
            width = 4,
            h5(strong(paste("Transition", tr_labels()[tr]))),
            fluidRow(
              column(width = 6, h5("Gammas"), gamma_inputs),
              column(width = 6, h5("Knots"),  knot_inputs)
            ),
            if (!is.null(add_btn)) div(style = "margin-top:10px;", add_btn)
          ),
          column(width = 8, plotOutput(paste0("spline_plot_", tr), height = "300px")),
          tags$div(style = "height:30px;")
        )
        
        # Reactive inputs for plotting 
        gamma_vals_re <- reactive({
          sapply(0:(n_gammas - 1), function(i)
            as.numeric(input[[paste0("gamma_", tr, "_", i)]]) %||% NA_real_)
        })
        knot_vals_re <- reactive({
          sapply(1:n_knots, function(i)
            as.numeric(input[[paste0("knot_", tr, "_", i)]]) %||% NA_real_)
        })
        
        # Plot 
        output[[paste0("spline_plot_", tr)]] <- renderPlot({
          
          req(input$baseline_xmax)
          times <- seq(0.1, input$baseline_xmax, length.out = 300)
          
          g <- gamma_vals_re()
          k <- knot_vals_re()
          
          if (is.null(g) || is.null(k) || any(is.na(g)) || any(is.na(k))) return(NULL)
          if (length(g) < 2 || length(k) < 2) return(NULL)
          
          logH <- function(gammas, time, knots) {
            # Defensive checks
            if (is.null(knots) || length(knots) < 2) return(NA_real_)
            knots <- suppressWarnings(as.numeric(knots))
            if (any(is.na(knots))) return(NA_real_)
            
            lowest.knot  <- knots[1]
            highest.knot <- knots[length(knots)]
            if (is.na(lowest.knot) || is.na(highest.knot)) return(NA_real_)
            if (!is.finite(lowest.knot) || !is.finite(highest.knot)) return(NA_real_)
            if (lowest.knot == highest.knot) return(NA_real_)
            
            # Nested cubic-spline basis function, fully guarded
            bfun <- function(time, knot) {
              tlog <- suppressWarnings(as.numeric(log(time)))
              knot  <- suppressWarnings(as.numeric(knot))
              if (any(is.na(c(tlog, knot, lowest.knot, highest.knot)))) return(NA_real_)
              if (!is.finite(tlog) || !is.finite(knot)) return(NA_real_)
              
              hk <- highest.knot
              lk <- lowest.knot
              denom <- hk - lk
              if (denom == 0 || is.na(denom)) return(NA_real_)
              
              t1 <- pmax(0, tlog - knot)^3
              t2 <- (hk - knot) * pmax(0, tlog - lk)^3 / denom
              t3 <- (knot - lk) * pmax(0, tlog - hk)^3 / denom
              
              as.numeric(t1 - t2 - t3)
            }
            
            # Evaluate basis
            bvalue <- if (length(knots) > 2) {
              sapply(knots[2:(length(knots) - 1)], function(x) bfun(time, x))
            } else numeric(0)
            
            # Safe combination
            if (any(is.na(bvalue))) return(NA_real_)
            g <- suppressWarnings(as.numeric(gammas))
            if (any(is.na(g))) return(NA_real_)
            if (length(g) < 2) return(NA_real_)
            
            as.numeric(g %*% c(1, log(time), bvalue))
          }
          
          logH_values <- suppressWarnings(sapply(times, function(x) logH(g, x, k)))
          if (any(is.na(logH_values))) return(NULL)
          H <- exp(logH_values)
          h <- c(diff(H) / diff(times), NA)
          
          par(mfrow = c(1, 2))
          plot(times, h, type = "l", lwd = 2, col = "darkblue",
               xlab = "Time", ylab = "Baseline hazard function",ylim=c(0,ymax_h),
               main = paste0("Transition ", tr, ": Baseline hazard"))
          plot(times, H, type = "l", lwd = 2, col = "darkred",
               xlab = "Time", ylab = "Baseline cumulative hazard",ylim=c(0,ymax_H),
               main = paste0("Transition ", tr, ": Baseline cumulative hazard"))
          par(mfrow = c(1, 1))
          
        })
      }
    }
    
    ui_list[[length(ui_list) + 1]] <- div(
      style = "text-align:right; margin-top:30px;",
      actionButton(
        "goToMicrosimulation",
        "Next: Treatment Strategies",
        class = "btn btn-primary"
      )
    )
    
    # Combine everything into one tagList
    do.call(tagList, ui_list)
  })
  
  observe({
    if (!length(non_tx_trans())) return()
    
    for (tr in non_tx_trans()) {
      # prevent duplicate observers
      if (isTRUE(observer_flags[[tr]])) next
      observer_flags[[tr]] <- TRUE
      
      # capture tr by value for this iteration
      local({
        tr_local <- tr
        btn_id   <- paste0("add_gamma_knot_", tr_local)
        
        observeEvent(input[[btn_id]], {
          isolate({
            gamma_counts[[tr_local]] <- as.integer((gamma_counts[[tr_local]] %||% 3)) + 1L
            knot_counts[[tr_local]]  <- as.integer((knot_counts[[tr_local]]  %||% 3)) + 1L
          })
        }, ignoreInit = TRUE)
      })
    }
  })
  
  observeEvent(input$has_data, {
    isolate({
      # Case 1: manual mode
      if (input$has_data == "no") {
        if (!length(non_tx_trans())) return()
        for (tr in non_tx_trans()) {
          gamma_counts[[tr]] <- 3
          knot_counts[[tr]]  <- 3
        }
        return()
      }
      
      # Case 2: model-based mode
      if (input$has_data == "yes") {
        fits <- tryCatch(fittedModelsRV(), error = function(e) NULL)
        # Skip if fits is not yet ready or invalid
        if (is.null(fits) || !is.list(fits)) return()
        if (length(fits) == 0) return()
        
        # Clear outdated transitions
        old_trs <- names(reactiveValuesToList(gamma_counts))
        for (tr in setdiff(old_trs, names(fits))) {
          gamma_counts[[tr]] <- NULL
          knot_counts[[tr]]  <- NULL
        }
        
        # Reset counts to match fitted models
        for (tr in names(fits)) {
          fit <- fits[[tr]]
          if (is.null(fit) || !is.list(fit)) next
          if (!"knots" %in% names(fit)) next
          coefs <- tryCatch(coef(fit), error = function(e) NULL)
          if (is.null(coefs) || !length(coefs)) next
          
          n_gamma <- length(grep("gamma", names(coefs)))
          n_knot  <- length(fit$knots)
          
          gamma_counts[[tr]] <- n_gamma
          knot_counts[[tr]]  <- n_knot
        }
      }
    })
  })
  
  observeEvent(input$goToMicrosimulation, {
    updateNavlistPanel(
      session,
      inputId = "tabs",   # the id of your tabsetPanel/navlistPanel
      selected = "Treatment Strategies"
    )
  })
  
  
  # ---- Tab 3: Treatment Strategies tab ----
  rv3 <- reactiveValues()
  rv3$covLabelsRV <- reactiveVal(NULL)
  sharedCovValues <- reactiveVal(list())  # numeric values of shared covariates
  interventionCovValuesRV <- reactiveVal(list()) # Store covariate values (shared + transition-specific) for each treatment transition
  
  output$interventionStrategiesUI <- renderUI({
    req(input$goToMicrosimulation)
    
    #  Basic checks
    if (length(tx_trans()) == 0)
      return(h5("No treatment transitions defined yet."))
    if (length(savedCovs()) == 0)
      return(h5("No covariates defined yet. Please fit the model first."))
    
    covariates <- unique(unlist(savedCovs()))
    covariates <- unique(c("Tstart", "StartState", covariates))
    
    # Build user-defined labels 
    cov_labels <- sapply(covariates, function(cov) {
      if (!is.null(input[[paste0("desc_", cov)]])) {
        input[[paste0("desc_", cov)]]
      } else if (cov == "StartState") {
        "Initial state of the microsimulation"
      } else if (cov == "Tstart") {
        "Treatment start time since prior state"
      } else {
        cov
      }
    }, USE.NAMES = TRUE)
    
    rv3$covLabelsRV(cov_labels)
    
    choices_named <- setNames(as.list(covariates), cov_labels)
    
    tagList(
      withMathJax(),  
      h4("Input the values for the covariates for each treatment transition. Plot the bar plots for treatment (yes vs no) at the transitions."),
      h4("Step 1: Choose Shared Covariates"),
      helpText("Select the covariates that have the same values across all treatment transitions."),
      checkboxGroupInput(
        inputId = "shared_covs",
        label = NULL,
        choices = choices_named,
      ),
      uiOutput("covariateValueInputs"),
      tags$hr(),
      h4("Step 2: Define Covariate Values for Each Treatment Transition"),
      helpText("For non-shared covariates, specify their values for each treatment transition."),
      uiOutput("transitionCovariateInputs"),
      if (input$has_data=="no"){
        helpText("One bar
shows an RMST; an error bar indicates \\(\\pm 1\\) standard deviation truncated at 0 and the prespecified time horizon.")
      },
      br(),
      div(
        style = "text-align:center;",
        actionButton("run_all_transitions", "Run microsimulation for all transitions", class = "btn btn-success")
      ),
      br()
    )
  })
  
  observe({
    req(input$shared_covs)
    
    for (cov in input$shared_covs) {
      local({
        cov_local <- cov
        id_local  <- paste0("covval_shared_", cov_local)
        
        observeEvent(input[[id_local]], {
          tmp <- sharedCovValues() %||% list()
          tmp[[cov_local]] <- input[[id_local]]
          sharedCovValues(tmp)
        }, ignoreInit = TRUE)
      })
    }
  })
  
  output$covariateValueInputs <- renderUI({
    req(input$goToMicrosimulation)
    
    cov_labels <- rv3$covLabelsRV() %||% list()
    covariates <- unique(unlist(savedCovs()))
    covariates <- unique(c("Tstart", "StartState", covariates))
    
    shared_covs <- input$shared_covs %||% character(0)
    if (length(shared_covs) == 0) return(NULL)
    
    start_states <- unique(edges()[edges()$label %in% tx_trans(),"from" ]) 
    
    wellPanel(
      h5("Shared Covariate Values"),
      fluidRow(
        lapply(shared_covs, function(cov) {
          label <- cov_labels[[cov]] %||% cov
          
          inputTag <- if (cov == "StartState") {
            selectInput(
              inputId = paste0("covval_shared_", cov),
              label = label,
              choices = start_states,
              selected = start_states[1] 
            )
          } else {
            numericInput(
              inputId = paste0("covval_shared_", cov),
              label = label,
              value = 0,  
              step = 0.1
            )
          }
          
          column(width = 6, inputTag)
        })
      )
    )
  })
  
  # Render per-transition covariate inputs (non-shared only)
  output$transitionCovariateInputs <- renderUI({
    req(input$goToMicrosimulation)
    
    start_states <-  unique(edges()[edges()$label %in% tx_trans(),"from" ])
    
    covariates <- unique(unlist(savedCovs()))
    
    shared_covs <- input$shared_covs %||% character(0)
    non_shared_covs <- setdiff(unique(c("Tstart", "StartState", covariates)), shared_covs)
    
    if (length(non_shared_covs) == 0) {
      h5("All covariates are shared — no transition-specific inputs needed.")
    }
    
    # Retrieve covariate labels
    cov_labels <- rv3$covLabelsRV() %||% setNames(covariates, covariates)
    
    # Build transition panels 
    n_trs <- length(tx_trans())
    cols_per_row <- 3
    n_rows <- ceiling(n_trs / cols_per_row)
    
    rows <- lapply(seq_len(n_rows), function(r) {
      start_idx <- (r - 1) * cols_per_row + 1
      end_idx <- min(r * cols_per_row, n_trs)
      
      cols <- lapply(start_idx:end_idx, function(i) {
        tr <- tx_trans()[i]
        tr_label <- paste("Transition", tr_labels()[tr])
        
        # --- Each transition shown as vertical list of numeric inputs ---
        inputs <- lapply(non_shared_covs, function(cov) {
          label <- cov_labels[[cov]] %||% cov
          if (cov == "StartState") {
            selectInput(
              inputId = paste0("covval_", cov, "_", tr),
              label = label,
              choices = if (length(start_states) > 0) start_states else "None",
              selected = if (length(start_states) > 0) start_states[1] else "None"
            )
          } else {
            numericInput(
              inputId = paste0("covval_", cov, "_", tr),
              label = label,
              value = 0,
              step = 0.1
            )
          }
        })
        
        column(
          width = 4,
          wellPanel(
            h5(strong(tr_label)),
            tagList(inputs),
            div(style = "text-align:center; margin-top:10px;",
                actionButton(paste0("run_tr_", tr), "Run microsimulation", class = "btn btn-primary btn-sm")),
            # Each transition gets its own plot
            plotOutput(paste0("microsimPlot_", tr), height = "250px")
          )
        )
      })
      fluidRow(cols)
    })
    
  })
  
  build_intervention_cov_values <- function() {
    if (length(tx_trans()) == 0) return()
    
    # Get list of all state names in order
    all_states <- unique(c(edges()$from, edges()$to))
    state_map <- setNames(seq_along(all_states), all_states)
    
    cov_labels <- rv3$covLabelsRV() %||% list()
    
    all_covs <- names(cov_labels) %||% character(0)
    shared_covs <- input$shared_covs %||% character(0)
    non_shared_covs <- setdiff(all_covs, shared_covs)
    
    vals_existing <- list()
    
    for (tr in tx_trans()) {
      # Collect shared values
      shared_vals <- lapply(shared_covs, function(cov) sharedCovValues()[[cov]] %||% 0)
      names(shared_vals) <- shared_covs
      
      # Collect non-shared values 
      non_shared_vals <- lapply(non_shared_covs, function(cov) {
        input[[paste0("covval_", cov, "_", tr)]] %||% NA_real_
      })
      names(non_shared_vals) <- non_shared_covs
      
      cov_values <- c(shared_vals, non_shared_vals)
      
      # Convert StartState to numeric 
      if ("StartState" %in% names(cov_values)) {
        start_char <- as.character(cov_values[["StartState"]])
        if (start_char %in% names(state_map)) {
          cov_values[["StartState"]] <- state_map[[start_char]]
        } else if (is.numeric(cov_values[["StartState"]])) {
          # already numeric, leave as is
        } else {
          cov_values[["StartState"]] <- NA_real_
        }
      }
      
      vals_existing[[tr]] <- cov_values
    }
    
    interventionCovValuesRV(vals_existing)
  }
  
  # Create a dynamic observer for each treatment transition button
  observe({
    for (tr in tx_trans()) {
      local({
        tr_local <- tr
        observeEvent(input[[paste0("run_tr_", tr_local)]], {
          isolate({
            req(input$goToMicrosimulation)
            build_intervention_cov_values()
          })
        }, ignoreInit = TRUE)
      })
    }
  })
  
  # create the list for the parameters for the transitions
  paramsSurvListByInterventionRV <- reactive({
    
    req(input$goToMicrosimulation)
    req(rv_gate$micro_ready)
    
    if (length(tx_trans()) == 0) return(NULL)
    
    all_HRs_df <- all_HRs() %||% data.frame()
    cov_values_all <- interventionCovValuesRV()
    surv_lists_all <- list()
    
    cat("\n==== Building params_surv_list for each treatment transition ====\n")
    
    # Loop over each treatment transition as "active"
    for (tr_active in tx_trans()) {
      cat("\n--- Active treatment transition:", tr_active, "---\n")
      surv_list <- list()
      
      for (tr in edges()$label){
        hr_sub <- subset(all_HRs_df, transition == tr)
        # ---------- Non-treatment transitions ----------
        if (tr %in% non_tx_trans()) {
          cat("Processing non-treatment:", tr, "\n")
          if(input$has_data=="yes"){
            req(fittedModelsRV())
            fits <- fittedModelsRV()
            beta_fit <- tryCatch(coef(fits[[tr]]), error = function(e) NULL)
            beta_fit <- beta_fit[hr_sub$covariate]
            if (all(!is.na(beta_fit))) {
              hr_fit <- exp(as.numeric(beta_fit))        # original (unrounded) HRs
              # replace where 3-decimal rounding matches hr_sub$HR
              same3 <- round(hr_fit, 3) == round(hr_sub$HR, 3)
              hr_sub$HR[same3] <- hr_fit[same3]
            }
          }
          if (input$baseline_hazard == "piecewise") {
            if (input$has_data == "no") {
              cut_ids <- grep(paste0("^cut_", tr, "_"), names(input), value = TRUE)
              cuts <- as.numeric(sapply(cut_ids, function(id) input[[id]]))
              cuts <- cuts[!is.na(cuts)]
              cuts <- sort(unique(cuts))
              haz <- sapply(1:(length(cuts) - 1),
                            function(i) input[[paste0("haz_", tr, "_", i)]] %||% NA_real_)
            } else {
              req(fittedModelsRV())
              fits <- fittedModelsRV()
              cuts <- fits[[tr]]$cuts
              haz <- sapply(1:(length(cuts) - 1),
                            function(i) input[[paste0("haz_", tr, "_", i)]] %||% NA_real_)
              
              haz_fit <- tryCatch(fits[[tr]]$hazards, error = function(e) NULL)
              if(!is.null(haz_fit)){
                same3 <- !is.na(haz) & !is.na(haz_fit) & (round(haz_fit, 3) == round(haz, 3))
                haz[same3] <- haz_fit[same3]
              }
            }
            
            cuts <- as.numeric(cuts)
            haz <- as.numeric(haz)
            cuts <- cuts[!is.na(cuts)]
            haz <- haz[!is.na(haz)]
            
            if (length(haz) == 0 || length(cuts) == 0 || length(haz) != (length(cuts) - 1)) {
              next
            }
            
            coefs <- lapply(haz, function(x) as.matrix(x))
            names(coefs) <- paste0("h0", seq_along(coefs))
            for (i in 1:length(coefs)){
              colnames(coefs[[i]]) <- "intercept"
            }
            coefs[[1]] <- matrix(c(coefs[[1]][1, 1], log(hr_sub$HR)), nrow = 1)
            colnames(coefs[[1]]) <- c("intercept", hr_sub$covariate)
            
            surv_list[[tr]] <- hesim::params_surv(
              coefs = coefs,
              dist = "pwexp",
              aux = list(time = cuts[-1])
            )
          }
          
          if (input$baseline_hazard == "spline") {
            # build exact expected input IDs based on counts
            ng <- gamma_counts[[tr]] %||% 0
            nk <- knot_counts[[tr]]  %||% 0
            gamma_ids <- paste0("gamma_", tr, "_", 0:(ng - 1))
            knot_ids  <- paste0("knot_",  tr, "_", 1:nk)
            gammas <- suppressWarnings(as.numeric(unlist(lapply(gamma_ids, function(id) input[[id]]))))
            knots  <- suppressWarnings(as.numeric(unlist(lapply(knot_ids,  function(id) input[[id]]))))
            gammas <- gammas[!is.na(gammas)]; knots <- knots[!is.na(knots)]
            if (length(gammas) == 0 || length(knots) == 0) next
            
            if (input$has_data=="yes"){
              req(fittedModelsRV())
              fits <- fittedModelsRV()
              gammas_fit <- coef(fits[[tr]])[grep("gamma", names(coef(fits[[tr]])))]
              knots_fit  <- fits[[tr]]$knots
              same3_gammas <- round(gammas_fit, 3) == round(gammas, 3)
              gammas[same3_gammas] <- gammas_fit[same3_gammas]
              same3_knots <- round(knots_fit, 3) == round(knots, 3)
              knots[same3_knots] <- knots_fit[same3_knots]
            }
            
            coefs <- lapply(gammas, function(x) as.matrix(x))
            names(coefs) <- paste0("gamma", seq_along(coefs))
            for (i in 1:length(coefs)){
              colnames(coefs[[i]]) <- "intercept"
            }
            coefs[[1]] <- matrix(c(coefs[[1]][1, 1], log(hr_sub$HR)), nrow = 1)
            colnames(coefs[[1]]) <- c("intercept", hr_sub$covariate)
            
            surv_list[[tr]] <- hesim::params_surv(
              coefs = coefs,
              dist = "survspline",
              aux = list(knots = knots, scale = "log_cumhazard", timescale = "log")
            )
          }
          
        }else if (tr %in% tx_trans()) { #  treatment transitions 
          cov_vals <- cov_values_all[[tr]] %||% list()
          
          # Assign Tstart logic
          if (tr == tr_active) {
            tstart_val <- cov_vals[["Tstart"]] %||% 0
          } else {
            tstart_val <- 1e10
          }
          
          coefhsct <- as.matrix(tstart_val)
          colnames(coefhsct) <- "Tstart"
          
          surv_list[[tr]] <- hesim::params_surv(
            coefs = list(est = coefhsct),
            dist = "fixed"
          )
          
          cat("  Transition", tr, "Tstart =", tstart_val, "\n")
        }
      }
      
      surv_lists_all[[tr_active]] <- hesim::params_surv_list(surv_list)
    }
    
    # Add one more scenario: all Tstart = 1e10 
    cat("\n--- Building 'all-inactive' treatment scenario ---\n")
    surv_list_all_inactive <- list()
    for (tr in edges()$label) {
      if (tr %in% non_tx_trans()) {
        # reuse first non-treatment structure if exists
        surv_list_all_inactive[[tr]] <- surv_lists_all[[1]][[tr]]
      } else if (tr %in% tx_trans()) {
        coefhsct <- as.matrix(1e10)
        colnames(coefhsct) <- "Tstart"
        surv_list_all_inactive[[tr]] <- hesim::params_surv(
          coefs = list(est = coefhsct),
          dist = "fixed"
        )
        cat("  Transition", tr, "Tstart = 1e10\n")
      }
    }
    surv_lists_all[["all_inactive"]] <- hesim::params_surv_list(surv_list_all_inactive)
    
    cat("\n✅ Built", length(surv_lists_all), "params_surv_list objects (including all_inactive)\n")
    surv_lists_all
  })
  
  observe({
    req(paramsSurvListByInterventionRV())
    cat("\n=== paramsSurvListByInterventionRV() updated ===\n")
    print(paramsSurvListByInterventionRV())
  })
  
  run_single_microsim <- function(tr_active) {
    req(paramsSurvListByInterventionRV())
    req(is.numeric(input$rmst_time), input$rmst_time > 0)
    req(is.numeric(input$num_sample_microsim), input$num_sample_microsim>0)
    
    all_states <- unique(c(edges()$from, edges()$to))
    state_map <- setNames(seq_along(all_states), all_states)
    
    cat("interventionCovValuesRV()")
    print(interventionCovValuesRV())
    
    cat("\n--- Running microsimulation for:", tr_active, "---\n")
    
    transmod_lists <- paramsSurvListByInterventionRV()
    cov_values_all <- interventionCovValuesRV()
    
    transmod_params_yes <- isolate(transmod_lists[[tr_active]])
    transmod_params_no  <- isolate(transmod_lists[["all_inactive"]])
    cov_vals <- cov_values_all[[tr_active]] %||% list()
    
    start_state_num <- suppressWarnings(as.numeric(cov_vals[["StartState"]]))
    if (is.null(start_state_num) || length(start_state_num) == 0L || is.na(start_state_num)) start_state_num <- 1
    
    cov_vals_for_data <- cov_vals
    cov_vals_for_data[["StartState"]] <- NULL
    cov_vals_for_data[["Tstart"]] <- NULL
    
    sim_data <- data.frame(intercept = 1)
    for (nm in names(cov_vals_for_data)) {
      v <- cov_vals_for_data[[nm]]
      if (length(v) == 1L && (is.numeric(v) || suppressWarnings(!is.na(as.numeric(v))))) {
        sim_data[[nm]] <- suppressWarnings(as.numeric(v))
      }
    }
    
    tstart_yes <- suppressWarnings(as.numeric(cov_vals[["Tstart"]]))
    strategies_yes <- data.frame(Tstart = tstart_yes, strategy_id = 1)
    strategies_no  <- data.frame(Tstart = 1e10, strategy_id = 2)
    
    one_ok <- TRUE
    sim_yes <- sim_no <- NULL
    prior_intervention_state <- edges()[tr_active,"from"]
    intervention_state <- edges()[tr_active,"to"]
    tryCatch({
      set.seed(1)
      sim_yes <- gen_trans_time(paramsurv_list=transmod_params_yes, trans_mat= tm(), sim_data = sim_data, 
                                start_state=start_state_num, 
                                prior_intervention_state=state_map[[prior_intervention_state]],intervention_state=state_map[[intervention_state]],
                                thor=input$rmst_time, n_sample=input$num_sample_microsim, model=input$baseline_hazard)
      set.seed(2)
      sim_no <- gen_trans_time(paramsurv_list=transmod_params_no, trans_mat= tm(), sim_data = sim_data, 
                               start_state=start_state_num, 
                               prior_intervention_state=state_map[[prior_intervention_state]],intervention_state=state_map[[intervention_state]],
                               thor=input$rmst_time, n_sample=input$num_sample_microsim, model=input$baseline_hazard)
      if (nrow(sim_yes) == 0L || nrow(sim_no) == 0L) one_ok <- FALSE
    }, error = function(e) {
      one_ok <<- FALSE
      message("⚠️ Skipped ", tr_active, " due to error: ", e$message)
    })
    
    if (!one_ok) return(NULL)
    
    sim_yes$strategy_id <- 1
    sim_no$strategy_id  <- 2
    sim_all <- rbind(sim_yes, sim_no)
    
    rmst <- sim_all %>%
      group_by(strategy_id, id) %>%
      slice_max(total_time, with_ties = FALSE) 
    
    rmst_summary <- rmst %>%
      group_by(strategy_id) %>%
      summarise(
        mean_rmst = mean(total_time, na.rm = TRUE),
        sd_rmst   = sd(total_time,   na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      mutate(
        transition = tr_active,
        strategy   = ifelse(strategy_id == 1, "Yes", "No")
      )
    
    if (input$has_data=="no"){
      rmst_summary <- rmst_summary %>%
        group_by(strategy_id) %>%
        mutate( # truncated standard deviation bar
          lb = pmax(0, mean_rmst - sd_rmst),
          ub = pmin(input$rmst_time, mean_rmst + sd_rmst))
    }
    
    rmst_summary$strategy <- factor(rmst_summary$strategy, levels = c("Yes", "No"))
    
    output[[paste0("microsimPlot_", tr_active)]] <- renderPlot({
      req(rv_gate$micro_ready)
      
      ggplot(rmst_summary, aes(x = strategy, y = mean_rmst, fill = strategy)) +
        geom_col(width = 0.6) +
        geom_errorbar(
          aes(ymin = lb, ymax = ub),
          width = 0.15, linewidth = 0.8
        ) +
        geom_text(
          aes(label = number(mean_rmst, accuracy = 0.01)),
          nudge_x = 0.10,                         
          nudge_y = 0.02 * input$rmst_time,       
          hjust = 0, vjust = 0,
          size = 5
        ) +
        coord_cartesian(clip = "off") +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("Treatment at state", edges()[tr_active, "from"]),
          x = NULL, y = "Restricted Mean Survival Time", fill = ""
        ) +
        theme(
          legend.position = "none",
          axis.text.x  = element_text(size = 16),
          axis.text.y  = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.margin = margin(5.5, 30, 5.5, 5.5)  # extra right margin for labels
        ) +
        scale_y_continuous(
          limits = c(0, input$rmst_time),
          expand = expansion(mult = c(0, 0.08))
        )
    })
    
    cat("✅ Completed single microsimulation for", tr_active, "\n")
  }
  
  observeEvent(input$run_all_transitions, {
    
    rv_gate$micro_ready <- TRUE
    req(input$goToMicrosimulation)
    
    build_intervention_cov_values()
    req(paramsSurvListByInterventionRV())
    cat("\n=== Running all microsimulations ===\n")
    for (tr_active in tx_trans()) {
      run_single_microsim(tr_active)
    }
    cat("\n✅ All microsimulations completed.\n")
  }, ignoreInit = TRUE)
  
  
  observe({
    req(paramsSurvListByInterventionRV())
    # For each treatment transition, create a separate observer
    lapply(tx_trans(), function(tr_active) {
      observeEvent(input[[paste0("run_tr_", tr_active)]], {
        rv_gate$micro_ready <- TRUE
        run_single_microsim(tr_active)
      }, ignoreInit = TRUE)
    })
  })
  
  observeEvent(input$reset_all, {
    updateTabsetPanel(
      session,
      inputId = "tabs",           
      selected = "Multistate Structure"  
    )
    
    updateTabsetPanel(
      session,
      inputId = "MSMstruc",           
      selected = "States and transitions"  
    )
  })
  
  observeEvent(input$reset_data, {
    updateTabsetPanel(
      session,
      inputId = "tabs",           
      selected = "Multistate Modeling"  
    )
    
    updateTabsetPanel(
      session,
      inputId = "multi_fit",           
      selected = "Dataset"  
    )
  })
  
  observeEvent(input$reset_covariate, {
    updateTabsetPanel(
      session,
      inputId = "tabs",           
      selected = "Multistate Modeling"  
    )
    
    updateTabsetPanel(
      session,
      inputId = "multi_fit",           
      selected = "Covariate Assignment"  
    )
  })
  
  observeEvent(input$reset_model, {
    updateTabsetPanel(
      session,
      inputId = "tabs",           
      selected = "Multistate Modeling"  
    )
    
    updateTabsetPanel(
      session,
      inputId = "multi_fit",           
      selected = "Model Specification"  
    )
  })
  
  observeEvent(input$reset_micro, {
    updateTabsetPanel(
      session,
      inputId = "tabs",           
      selected = "Treatment Strategies" 
    )
  })
  
  # ---- RESET ---- #
  # tab 1: multistate structure
  reset_struc_func <- function(){
    rv1$stateNames <- reactive({
      n <- input$numStates %||% 3
      nm <- vapply(seq_len(n), function(i) input[[paste0("n", i)]] %||% paste("State", i), character(1))
      nm <- trimws(nm)
      nm[nm == ""] <- paste0("State ", which(nm == ""))
      make.unique(nm, sep = " ")
    })
    rv1$transitionEdges <- reactiveVal(data.frame(from=character(), to=character(), label=character(), stringsAsFactors=FALSE))
    rv1$interventionSel <- reactiveVal(character(0))
    rv1$transitionMatrix <- reactiveVal(NULL)
    updateNumericInput(session, "numStates", value = 3)
  }
  
  # tab2_1:data
  reset_data_func <- function(){
    rv2_1$firstTransformed <- NULL
    rv2_1$extraList <- list()
    rv2_1$combined(NULL)
    
    shinyjs::reset("dataFile")
    shinyjs::reset("extraDataFile")
    shinyjs::reset("transitionSelect")
    
    output$uploadedDataDT        <- renderDT(NULL)
    output$msdataTitle           <- renderUI(NULL)
    output$msdataDT              <- renderDT(NULL)
    output$transitionCountTitle  <- renderUI(NULL)
    output$transitionCountDT     <- renderDT(NULL)
    output$addDataUI             <- renderUI(NULL)
    output$combinedDataTitle     <- renderUI(NULL)
    output$combinedDataDT        <- renderDT(NULL)
    output$combinedCountTitle    <- renderUI(NULL)
    output$combinedCountDT       <- renderDT(NULL)
    output$addSectionUI          <- renderUI(NULL)
  }
  
  # tab 2-2: covariate
  reset_covariate_func <- function(){
    rv2_2$blockCountRV(1)
    rv2_2$blockCovNamesRV(list(NULL))
    rv2_2$blockTargetsRV(list(NULL))
  }
  
  # tab 2-3: model
  reset_model_func <- function(){
    updateNumericInput(session, "baseline_hazard", value = "piecewise")
    updateNumericInput(session, "max_intervals", value = 20)
    updateNumericInput(session, "num_events", value = 20)
    rv_gate$model_ready <- FALSE
  }
  
  # tab 3: microsimulation
  reset_micro_func <- function(){
    updateNumericInput(session, "rmst_time", value = 15)
    updateNumericInput(session, "num_sample_microsim", value = 1000)
    updateNumericInput(session, "shared_covs", value = character(0))
    rv_gate$micro_ready <- FALSE
  }
  
  # reset_all (<- tab 1)
  observeEvent(input$reset_all, {
    reset_struc_func() # tab 1
    reset_data_func() # tab 2-1
    reset_covariate_func() # tab 2-2
    reset_model_func() # tab 2-3
    reset_micro_func()# tab 3
  }, ignoreInit = TRUE)
  
  # reset_data (<- tab 2-1)
  observeEvent(input$reset_data, {
    reset_data_func() # tab 2-1
    reset_covariate_func() # tab 2-2
    reset_model_func() # tab 2-3
    reset_micro_func()# tab 3
  }, ignoreInit = TRUE)
  
  # reset_covariate (<- tab 2-2)
  observeEvent(input$reset_covariate, {
    reset_covariate_func() # tab 2-2
    reset_model_func() # tab 2-3
    reset_micro_func()# tab 3
  }, ignoreInit = TRUE)
  
  # reset_model (<- tab 2-3)
  observeEvent(input$reset_model, {
    reset_model_func() # tab 2-3
    reset_micro_func()# tab 3
  }, ignoreInit = TRUE)
  
  # reset_model (<- tab 3)
  observeEvent(input$reset_micro, {
    reset_micro_func()# tab 3
  }, ignoreInit = TRUE)
  
}

shinyApp(ui = ui, server = server)

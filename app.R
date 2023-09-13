# app.R
# Author : Brooks Benard, bbenard@stanford.edu
# Date: 06/22/2023
# updated 09/13/2023
# Description: This script creates the ui and server code required to run the DepMap Shiny app that performs mutation and cancer enrichment analyses based on gene effect distributions

# Define ui ----
ui <- fluidPage(tags$head(tags$style(
  HTML(
    "
        /* CSS styles here */
        .main-content {
          max-height: 100vh;
          overflow-y: auto;
        }
        .tab-content {
          width: 100%; /* Adjust the width as needed */
        }
        .plot-container {
          width: 100%;
          height: calc(80vh - 70px); /* Adjust the offset value if needed */
          margin-top: 0px !important; /* Adjust the top margin */
          margin-bottom: 0px !important; /* Adjust the bottom margin */
        }
        /* Increase the size of plot elements */
        .plot-container svg {
          font-size: 20px;
        }
        /* Increase the size of axis labels */
        .plot-container .y.axis text,
        .plot-container .x.axis text {
          font-size: 20px;
        }
        /* Increase the size of axis tick labels */
        .plot-container .y.axis .tick text,
        .plot-container .x.axis .tick text {
          font-size: 20px;
        }
        "
  )
)),
navbarPage(
  inverse = TRUE,
  title = NULL,
  tabPanel(
    "DepMap Shiny",
    titlePanel("DepMap analysis interface"),
    sidebarLayout(
      sidebarPanel(
        textInput(
          "cancer_type",
          label = "Select Cancer Of Interest:",
          value = "",
          placeholder = "Pan-Cancer, AML, etc."
        ),
        textInput(
          "genes_of_interest",
          label = "Enter Gene Of Interest:",
          value = "",
          placeholder = "HUGO ID (all caps)"
        ),
        actionButton("analysis_start", "Run Analysis")
      ),
      mainPanel(tabsetPanel(
        type = "tabs",
        tabPanel(
          "Mutation Enrichment",
          id = "majorTab1",
          tabsetPanel(
            tabPanel(
              "Effect distribution",
              div(class = "plot-container", plotOutput("effect", height = "90%")),
              downloadButton("download_effect", "Download Plot"),
              downloadButton("download_effect_data", "Download Plot Data")
            ),
            tabPanel(
              "Distribution by cancer type",
              div(class = "plot-container", plotOutput("effect_group", height = "90%")),
              downloadButton("download_effect_group", "Download Plot"),
              downloadButton("download_effect_group_data", "Download Plot Data")
            ),
            tabPanel(
              "Mutation enrichment",
              div(class = "plot-container", plotOutput("enrichment", height = "90%")),
              downloadButton("download_enrichment", "Download Plot"),
              downloadButton("download_enrichment_data", "Download Plot Data")
            ),
            tabPanel(
              "Top differentially mutated genes",
              div(class = "plot-container", plotOutput("top_genes", height = "90%")),
              downloadButton("download_top_genes", "Download Plot"),
              downloadButton("download_top_genes_data", "Download Plot Data")
            ),
            tabPanel(
              "Mutation NES plots",
              div(class = "plot-container", plotOutput("nes_plots", height = "90%")),
              downloadButton("download_nes_plots", "Download Plot"),
              downloadButton("download_nes_plots_data", "Download Plot Data")
            ),
            tabPanel(
              "Mutation distributions",
              div(class = "plot-container", plotOutput("mut_dist", height = "90%")),
              downloadButton("download_mut_dist", "Download Plot"),
              downloadButton("download_mut_dist_data", "Download Plot Data")
            )
          )
        ),
        tabPanel(
          "Cancer Enrichment",
          id = "majorTab2",
          tabsetPanel(
            tabPanel(
              "Effect distribution",
              div(class = "plot-container", plotOutput("effect_cancer", height = "90%")),
              downloadButton("download_effect_cancer", "Download Plot"),
              downloadButton("download_effect_cancer_data", "Download Plot Data")
            ),
            tabPanel(
              "Distribution by cancer type",
              div(class = "plot-container", plotOutput("effect_group_cancer", height = "90%")),
              downloadButton("download_effect_cancer_group", "Download Plot"),
              downloadButton("download_effect_cancer_group_data", "Download Plot Data")
            ),
            tabPanel(
              "Cancer enrichment",
              div(class = "plot-container", plotOutput("cancer_enrichment", height = "90%")),
              downloadButton("download_cancer_enrichment", "Download Plot"),
              downloadButton("download_cancer_enrichment_data", "Download Plot Data")
            ),
            tabPanel(
              "Top differentially enriched cancers",
              div(class = "plot-container", plotOutput("top_cancer_plots", height = "90%")),
              downloadButton("download_top_cancer_plots", "Download Plot"),
              downloadButton("download_top_cancer_plots_data", "Download Plot Data")
            ),
            tabPanel(
              "Cancer type distributions",
              div(class = "plot-container",plotOutput("cancer_distribution_plots", height = "90%")),
              downloadButton("download_cancer_distribution_plots", "Download Plot"),
              downloadButton(
                "download_cancer_distribution_plots_data",
                "Download Plot Data"
              )
            ),
            tabPanel(
              "Cancer NES plots",
              div(class = "plot-container", plotOutput("cancer_nes_plots", height = "90%")),
              downloadButton("download_cancer_nes_plots", "Download Plot"),
              downloadButton("download_cancer_nes_plots_data", "Download Plot Data")
            )
          )
        )
      ))
    )
  )
))

# Define server logic ----
server <- function(input, output, session) {
  # Store the result of run_enrichment() in a reactive value
  enriched_plots <- eventReactive(input$analysis_start, {
    req(input$genes_of_interest, input$cancer_type)
    
    genes_of_interest <-
      toupper(input$genes_of_interest)  # Convert to uppercase
    cancer_type <-
      toupper(input$cancer_type)  # Convert to uppercase
    
    crispr_dependency_enrichment(genes_of_interest = genes_of_interest,
                                 cancer_type = cancer_type)
  })
  
  # Define reactive values for the input variables
  cancer_type <- reactiveVal()
  genes_of_interest <- reactiveVal()
  
  observe({
    cancer_type(input$cancer_type)
    genes_of_interest(input$genes_of_interest)
  })
  
  # Function to extract the plots from the reactive expression
  get_plots <- function() {
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      return(plots)
    } else {
      return(list(NULL))
    }
  }
  
  # Function to extract the data from the reactive expression
  get_data <- function() {
    data <- enriched_plots()$data
    if (!is.null(data)) {
      return(data)
    } else {
      return(list(NULL))
    }
  }
  
  # Render the plots based on the selected plot index
  output$effect <- renderPlot({
    plots <- get_plots()
    if (!is.null(plots)) {
      plots[[1]]
    }
  })
  
  output$effect_cancer <- renderPlot({
    plots <- get_plots()
    if (!is.null(plots)) {
      plots[[1]]
    }
  })
  
  output$effect_group <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[2]]
    }
  })
  
  output$effect_group_cancer <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[2]]
    }
  })
  
  output$enrichment <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[3]]
    }
  })
  
  output$top_genes <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[4]]
    }
  })
  
  output$mut_dist <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[5]]
    }
  })
  
  output$nes_plots <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[6]]
    }
  })
  
  output$cancer_enrichment <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[7]]
    }
  })
  
  output$top_cancer_plots <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[8]]
    }
  })
  
  output$cancer_distribution_plots <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[9]]
    }
  })
  
  output$cancer_nes_plots <- renderPlot({
    plots <- enriched_plots()$plots
    if (!is.null(plots)) {
      plots[[10]]
    }
  })
  
  # Download effect distribution info ----
  # Download plot for "Effect distribution" tab
  # mutations
  output$download_effect <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_distribution.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[1]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[1]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Effect distribution" tab
  # mutations
  output$download_effect_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_distribution_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[1]])) {
        write.csv(data[[1]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # cancers
  output$download_effect_cancer <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_distribution_cancer.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[1]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[1]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Effect distribution" tab
  # cancers
  output$download_effect_cancer_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_distribution_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[1]])) {
        write.csv(data[[1]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download cancer effect distribution info ----
  # Download plot for "Distribution by cancer type" tab
  # mutations
  output$download_effect_group <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_group.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[2]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[2]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Distribution by cancer type" tab
  # mutations
  output$download_effect_group_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_group_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[1]])) {
        write.csv(data[[1]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download plot for "Distribution by cancer type" tab
  # cancers
  output$download_effect_cancer_group <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_group_cancer.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[2]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[2]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Distribution by cancer type" tab
  # cancers
  output$download_effect_cancer_group_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_effect_group_cancer_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[1]])) {
        write.csv(data[[1]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download mutation enrichment info ----
  # Download plot for "Mutation enrichment" tab
  output$download_enrichment <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_mutation_enrichment.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[3]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[3]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Mutation enrichment" tab
  output$download_enrichment_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_mutation_enrichment_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[2]])) {
        write.csv(data[[2]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download top enriched mutations info ----
  # Download plot for "Top differentially mutated genes" tab
  output$download_top_genes <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_top_genes.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[4]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[4]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Top differentially mutated genes" tab
  output$download_top_genes_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_top_genes_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[3]])) {
        write.csv(data[[3]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download mutation NES info ----
  # Download plot for "Mutation NES plots" tab
  output$download_nes_plots <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_mutation_nes_plots.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[6]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[6]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Mutation NES plots" tab
  output$download_nes_plots_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_mutation_nes_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[5]])) {
        write.csv(data[[5]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download mutation distribution info ----
  # Download plot for "Mutation distributions" tab
  output$download_mut_dist <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_mut_dist.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[5]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[5]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Mutation distributions" tab
  output$download_mut_dist_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_mut_dist_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[4]])) {
        write.csv(data[[4]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download cancer enrichment info ----
  # Download plot for "Cancer enrichment" tab
  output$download_cancer_enrichment <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_cancer_enrichment.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[7]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[7]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Cancer enrichment" tab
  output$download_cancer_enrichment_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_cancer_enrichment_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[6]])) {
        write.csv(data[[6]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download top cancers enriched info ----
  # Download plot for "Top differentially enriched cancers" tab
  output$download_top_cancer_plots <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_top_cancers.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[8]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[8]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Top differentially enriched cancers" tab
  output$download_top_cancer_plots_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_top_cancers_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[7]])) {
        write.csv(data[[7]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
  
  # Download cancer distribution info ----
  # Download plot for "Cancer type distributions" tab
  output$download_cancer_distribution_plots <-
    downloadHandler(
      filename = function() {
        paste0(cancer_type(),
               "_",
               genes_of_interest(),
               "_cancer_distribution.png")
      },
      content = function(file) {
        plots <- get_plots()
        if (!is.null(plots[[9]])) {
          png(file,
              width = 800,
              height = 600,
              res = 300)
          print(plots[[9]])
          dev.off()
        }
      },
      contentType = "image/png"
    )
  
  # Download plot data for "Cancer type distribution" tab
  output$download_cancer_distribution_plots_data <-
    downloadHandler(
      filename = function() {
        paste0(cancer_type(),
               "_",
               genes_of_interest(),
               "_cancer_distribution_data.csv")
      },
      content = function(file) {
        data <- get_data()
        if (!is.null(data[[8]])) {
          write.csv(data[[8]], file, row.names = FALSE)
        }
      },
      contentType = "text/csv"
    )
  
  # Download cancer NES info ----
  # Download plot for "Cancer NES plots" tab
  output$download_cancer_nes_plots <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_cancer_nes.png")
    },
    content = function(file) {
      plots <- get_plots()
      if (!is.null(plots[[10]])) {
        png(file,
            width = 800,
            height = 600,
            res = 300)
        print(plots[[10]])
        dev.off()
      }
    },
    contentType = "image/png"
  )
  
  # Download plot data for "Cancer NES plots" tab
  output$download_cancer_nes_plots_data <- downloadHandler(
    filename = function() {
      paste0(cancer_type(),
             "_",
             genes_of_interest(),
             "_cancer_nes_data.csv")
    },
    content = function(file) {
      data <- get_data()
      if (!is.null(data[[9]])) {
        write.csv(data[[9]], file, row.names = FALSE)
      }
    },
    contentType = "text/csv"
  )
}
shinyApp(ui = ui, server = server)
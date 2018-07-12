library(shiny)

# Required libraries
packagelist <- ("insect")

# Check which packages are missing
missing <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]

# Install missing packages with dependencies
if(length(missing)) install.packages(missing, dependencies = TRUE)

library("insect")

ui <- fluidPage(
   titlePanel(h1("symulator")),
   titlePanel(h4("clade-level identification for Symbiodinium ITS2 sequences")),
   fluidRow(column(12, br(), textInput("sequence", "Paste your sequence here:", 
                                       width = "200%"))),
   fluidRow(column(12, br(),
                   sliderInput(
                     "threshold","Akaike weight threshold", min = 0.5, max = 1, value = 0.9))),
   #fluidRow(column(12, br(), fileInput("fasta", "or Load ITS2 .fasta file", multiple = FALSE, accept = NULL))),
   fluidRow(column(6, actionButton("do", "Go", icon = icon("barcode", lib = "glyphicon")), br(), br(), br())),
   fluidRow(
     column(12,
            mainPanel(
              tabsetPanel(
                tabPanel("Results", br(), h2(textOutput("text1"))),
                tabPanel("Report", br(), plotOutput("plot1"))
                #tabPanel("More info", br(), plotOutput("plot1"))
              ),
              width = 12
            )
     )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  tree <- readRDS("./data/classification_tree2.rds")
  observeEvent(input$do,{
    myseq <- toupper(input$sequence)
    myseq <- gsub("[ -]", "", myseq)
    myseq <- gsub("U", "T", myseq)
    thresh <- input$threshold
    if(grepl("[^ACGT]", myseq)){
      output$text1 <- renderText({
        "ERROR: SEQUENCE CONTAINS NON STANDARD CHARACTERS"
      })
    }else{
      myseq <- char2dna(myseq)
      names(myseq) <- "S1"
      myseq[[1]] <- shave(myseq[[1]], motif = attr(tree, "model"))
      classif <- classify(myseq, tree, threshold = thresh)
      if(is.na(classif$score[1])) classif$score[1] <- 0.999
      if(classif$score[1] == 1) classif$score[1] <- 0.999
      if(classif$taxon[1] == "root") classif$taxon[1] <- "unknown"
      output$text1 <- renderText({
        paste0("Clade ", classif$taxon[1], " with probability ", classif$score[1], "\n")
        #out  <- data.frame(taxon = paste0("Clade ", classif$taxon[1]), probability = classif$score[1]) 
        #out
      })
    }

    # output$plot1 <- renderPlot({
    #   plot(0:1, 0:1, type = "n")
    # })
  })

}

# Run the application 
shinyApp(ui = ui, server = server)


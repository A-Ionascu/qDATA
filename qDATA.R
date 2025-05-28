

# List packages to check installation
packages <- c("ggplot2", "dplyr","tidyr","shiny","shinythemes","DT",
              "gridExtra","ggthemes","nortest","tseries","reporter","qqconf","qqplotr",
              "shinycssloaders","shinybusy","shinycustomloader","car", "rstatix",
              "fontawesome","RColorBrewer","ggtext","coin","exactRankTests","bslib","rlist",
              "htmltools")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))





if(interactive()){
  
  
  # just make ui nice
  
  ui <- shinyUI( navbarPage( 
    #theme = shinytheme("cerulean"),
    theme = bs_theme(version = 5, bootswatch = "cosmo",
                     base_font = "Times New Roman", heading_font = "Times New Roman"),
    title="qDATA",
    main_page <-  tabPanel(title=HTML(paste0("Livak method" 
                                             #"2", tags$sup("-\u0394\u0394Ct") ," 
                                             #method"
                                             )), value="livak",
                           titlePanel(HTML(paste0("Livak 2", tags$sup("-\u0394\u0394Ct")," method"))),
                           sidebarLayout(
                             sidebarPanel(#style = "position:fixed;width:20%;",
                                          title="File input",
                                          fileInput("file", "Choose .csv file",
                                                    accept = c('.csv','.xls', ".xlsx"), buttonLabel=list(icon("upload"), paste("Browse"))),
                                          width = 3,
                                          br(),
                                          checkboxGroupInput("col_n","Data to analyse:",choices = c() ),
                                          br(),
                                          selectizeInput("typeofdata","Choose type of data:", 
                                                         choices = list("Linear form of \u0394Ct values" = "2^_transformed", 
                                                                        "\u0394Ct values" = "untransformed")),
                                          br(),
                                          numericInput("pvalue",label=HTML(paste0("Statistical significance threshold (\u03B1):")),value = c(0.05)),
                                          actionButton("update2","Update results", icon=icon("redo")),
                                          br(),br(),br(),
                                          numericInput("decNum","Number of decimals (not applicable to p values):",value = c(3)),
                                          actionButton("update","Update results",  icon=icon("redo")),
                                          br(),br(),br(),
                                          selectInput(inputId = "plot_orientation", label="Choose orientation for X-axis labels:", 
                                                      choices = list("Horizontal" = "horizontal", "30 degrees angle" = "angle30", "45 degree angle" = "angle45", 
                                                                     "60 degree angle" = "angle60", "Vertical" = "vertical"),
                                                      selected = "horizontal"),
                                          br(),
                             ),
                             
                             mainPanel(
                               tabsetPanel(
                                 tabPanel(
                                   title="Original file", value="input",
                                          br(),
                                          h4("Upload data using the `Browse` button."),
                                          br(),
                                          h4("Prepare your data table as in the the downloadable input model:"),
                                          downloadButton("input_model_download","Download model of input table"),
                                          br(),br(),
                                          DT::dataTableOutput('original'),
                                          br(),br(),br(),br()),
                                 
                                 # Summary statistics
                                 tabPanel(title="Summary statistics", value="summary",
                                          br(),
                                          h3("Summary statistics"),
                                          h4(HTML(paste0("Table and plots of descriptive statistics parameters corresponding to either 2", 
                                                         tags$sup("-\u0394Ct"), " values (linear form of \u0394Ct values) or \u0394Ct values", 
                                                         " calculated from each experimental group of interest."))),
                                          br(),
                                          downloadButton("summary_download","Export summary table"),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput('summary'), type="html", loader="dnaspin"),
                                          br(),br(),
                                          h3("Boxplot"),
                                          h4("Boxplot representation of ", HTML(paste0("2", tags$sup("-\u0394Ct")))," values or \u0394Ct values."),
                                          withLoader(plotOutput(outputId = "boxplot"), type="html", loader="dnaspin"),
                                          downloadButton("download_boxplots", "Export boxplots"),
                                          br(),br(),
                                          h3("Violin plot"),
                                          h4("Violin plots of ", HTML(paste0("2", tags$sup("-\u0394Ct")))," values or \u0394Ct values."),
                                          withLoader(plotOutput(outputId = "violin"), type="html", loader="dnaspin"),
                                          downloadButton("download_violin_plots", "Export violin plots"),
                                          br(),br(),
                                          h3("Histograms"),
                                          h4("Histograms of ", HTML(paste0("2", tags$sup("-\u0394Ct"))), " values or \u0394Ct values corresponding to each experimental group of interest.
                                           Mean value is shown as a continuous line and median values is shown as a dashed line."),
                                          withLoader(plotOutput(outputId = "histogram", height = "800px"), type="html", loader="dnaspin"),
                                          downloadButton("download_histograms", "Export histograms"),
                                          br(),br(),br(),br()),
                                 
                                 # Normality analysis
                                 tabPanel(title="Normality statistics", value="normality",
                                          br(),
                                          h3("Normality statistics"),
                                          h4(HTML(paste0("Table of normality tests (including Kolmogorov-Smirnov, Shapiro-Wilk, Anderson-Darling and Jarque-Bera) corresponding to either 2", 
                                                         tags$sup("-\u0394Ct"), " values (linear form of \u0394 values) or \u0394Ct values calculated from each experimental group of interest."))),
                                          br(),
                                          downloadButton("normality_download", "Export normality table"),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput('normality'), type="html", loader="dnaspin"),
                                          br(),br(),
                                          h3("Q-Q plots"),
                                          h4("Quantile-quantile (Q-Q) plots of ", HTML(paste0("2", tags$sup("-\u0394Ct")))," values or \u0394Ct values."),
                                          br(),
                                          withLoader(plotOutput(outputId = "qqplot"), type="html", loader="dnaspin"),
                                          downloadButton("downloadQQ", "Export Q-Q plots"),
                                          br(),br(),br(),br()),
                                 
                                 # Two samples testing
                                 tabPanel(title="Two samples testing", value="pairs",
                                          br(),
                                          h3(HTML(paste0("Two sample groups statistical testing"))),
                                          h4(HTML(paste0("Statistical testing between two sample groups of interest using either 2", 
                                                         tags$sup("-\u0394Ct"), " values (linear form of \u0394 values) or \u0394Ct values."))),
                                          br(),
                                          h4("There are three different types of statistical tests available. 
                                           Select the most suitable test based on the corresponding assumptions such as data normality and equal/unequal variances between groups.
                                           The assumptions of each test are presented below:"),
                                          br(),
                                          dataTableOutput("types_of_paired_tests"),
                                          br(),br(),
                                          h4("Select the type of statistical test:"),
                                          selectInput("select_paired", label=NULL, 
                                                      choices = list("Student's two-sample t-test"="ttest",
                                                                     "Student's t-test Welch corrected"="welchttest",
                                                                     "Mann-Whitney U test (Wilcoxon)"="wilcoxon")),
                                          br(),br(),
                                          h3("Two sample groups statistical testing"),
                                          downloadButton("paired_download","Export calculated table"),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput("paired"), type="html", loader="dnaspin"),
                                          br(),br(),br(),br()),
                                 
                                 # Multiple samples testing
                                 tabPanel(title="Multiple samples testing", value="multiple",
                                          br(),
                                          h3("Multiple sample groups statistical testing"),
                                          h4(HTML(paste0("Statistical testing between multiple sample groups of interest using either 2", 
                                                         tags$sup("-\u0394Ct"), " values (linear form of \u0394 values) or \u0394Ct values."))),
                                          br(),
                                          h4("There are three different types of statistical tests available. 
                                           Select the most suitable test based on the corresponding assumptions such as data normality and equal/unequal variances between groups.
                                           The assumptions of each test are presented below:"),
                                          br(),
                                          dataTableOutput("types_of_multiple_tests"),
                                          br(),br(),
                                          h4("Select the type of statistical test:"),
                                          selectInput("select_multiple",label=NULL, 
                                                      choices = list("One-way ANOVA"="anova",
                                                                     "Welch's ANOVA"="welchanova",
                                                                     "Kruskal-Wallis ANOVA"="kruskalwallis")),
                                          br(),
                                          h3("Multiple sample groups testing"),
                                          downloadButton("multiple_download","Export calculated table"),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput("multiple_table"), type="html", loader="dnaspin"),
                                          br(),br(),br(),
                                          h3("Post-hoc analysis"),
                                          downloadButton("posthoc_download","Export post-hoc table"),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput("posthoc_table"), type="html", loader="dnaspin"),
                                          br(),br(),br(),br()),
                                 
                                 # Fold Change
                                 tabPanel(title="Fold Change", value="FC",h3("Fold change"),
                                          br(),
                                          h4("Barplot of FC values calculated between each experimental group of interest using the Livak ",
                                             HTML(paste0("2", tags$sup("-\u0394\u0394Ct"))), " formula. Whiskers are shown as standard error (SE) for each comparison."),
                                          br(),
                                          downloadButton("downloadFC", "Export FC plot"),
                                          withLoader(plotOutput("plotFC", height = 700), type="html", loader="dnaspin"),
                                          br(),br(),br(),
                                          h4("Summary table of fold-change (FC = ", HTML(paste0("average 2",tags$sup("-\u0394\u0394Ct"),")")), 
                                             " values between each pair of experimental group of interest using the Livak formula."),
                                          br(),
                                          downloadButton("fc_download","Export FC table"),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput('fc'), type="html", loader="dnaspin"),
                                          br(),br(),br(),br()),
                                 
                                 # Log FC
                                 tabPanel(title=HTML(paste0("Log",tags$sub("2")," Fold Change")), value="log2FC",
                                          br(),
                                          h3(HTML(paste0("Log",tags$sub("2")," Fold Change"))),
                                          br(),
                                          h4("Barplot of ", HTML(paste0("log",tags$sub("2"),"FC")), " values calculated between each experimental group of interest using the Livak ",
                                             HTML(paste0("2", tags$sup("-\u0394\u0394Ct"))), " formula. Whiskers are shown as standard error (SE) for each comparison."),
                                          br(),
                                          downloadButton("downloadlogFC", HTML(paste0("Export log",tags$sub("2"),"FC plot"))),
                                          withLoader(plotOutput("plotlogFC", height = 700), type="html", loader="dnaspin"),
                                          br(),br(),br(),
                                          h4("Summary table of", HTML(paste0("log",tags$sub("2"),"FC")), 
                                             HTML(paste0("( = average log",tags$sub("2"),"2",tags$sup("-\u0394\u0394Ct"),")")), 
                                             " values between each pair of experimental group of interest using the Livak formula."),
                                          br(),
                                          downloadButton("log2fc_download",HTML(paste0("Export Log", tags$sub("2"), "FC table"))),
                                          br(),br(),
                                          withLoader(DT::dataTableOutput('log2fc'), type="html", loader="dnaspin"),
                                          br(),br(),br(),br()),
                                 
                                 # Export tab
                                 tabPanel(title="Export files", value="download",
                                          br(),
                                          h3("Export calculated values in a simple tabelar format to be used as input in other bioinformatics and statistics tools."),
                                          br(),
                                          h4(HTML(paste0("Export the \u0394Ct data:"))),
                                          downloadButton("outputDctvalues", HTML(paste0("Export \u0394Ct table"))),
                                          br(),br(),
                                          h4(HTML(paste0("Export the 2", tags$sup("-\u0394Ct"), " data:"))),
                                          downloadButton("output2Dctvalues", HTML(paste0("Export 2", tags$sup("-\u0394Ct"), " table"))),
                                          br(),br(),
                                          h4(HTML(paste0("Export the 2", tags$sup("-\u0394\u0394Ct"), " FC data:"))),
                                          downloadButton("outputDDctvalues", HTML(paste0("Export 2", tags$sup("-\u0394\u0394Ct"), " table"))),
                                          br(),br(),br(),br(),
                                          h4(HTML(paste0("Export all two sample tests"))),
                                          downloadButton("all_two_tests", HTML(paste0("Export two sample tests table"))),
                                          br(),br(),
                                          h4(HTML(paste0("Export all ANOVA tests"))),
                                          downloadButton("all_anova_tests", HTML(paste0("Export ANOVA tests table"))),
                                          br(),br(),
                                          h4(HTML(paste0("Export all post-hoc sample tests"))),
                                          downloadButton("all_posthoc_tests", HTML(paste0("Export post-hoc tests table"))),
                                          br(),br(),br(),br())
                               )
                             )
                           )),
    
    about_page <- tabPanel(title="About",
                           titlePanel("About"),
                           h3("qRTdatA (quick Real-Time PCR Data Analyser) is a tool created with R Shiny in the Drosophila Laboratory, Department of Genetics, Faculty of Biology, University of Bucharest."),
                           h5("March-July 2023"),
                           br(),
                           h3("For citing our tool, please use: Ionascu, A., Ecovoiu, A.A., Chifiriuc, M.C., Ratiu, A.C. (2023). qDATA - an R application implementing a practical framework for analyzing quantitative Real-Time PCR data. bioRxiv, 2023-11."),
                           h3("Access qDATA article from", tags$a(href="https://doi.org/10.1101/2023.11.29.569183", "BioRxiv")),
                           h3("For additional information, suggestions or inquiries, please contact us via email: a.ionascu20@s.bio.unibuc.ro"),
                           br(),br(),br(),
                           h3("This tool uses several abbreviations which are explained here:"),
                           br(),
                           h4("Input table and `Original file` tab: "),
                           h5("BR = biological replicate"),
                           h5("TR = technical replicate"),
                           h5("Ct = cycle threshold from the qRT-PCR machine"),
                           br(),
                           h4("`Summary statistics` tab:"),
                           h5("SD = standard deviation"),
                           h5("SE = standard error"),
                           h5(HTML(paste0("CI-lower/upper = confidence interval values calculated as mean &#177 SE"))),
                           h5("IQR = interquartile range"),
                           br(),
                           h4("`Two samples testing` tab:"),
                           h5("df = degrees of freedom"),
                           br(),
                           h4("`Multiple sample testing` tab:"),
                           h5("DFg = degrees of freeedom for compared groups"),
                           h5("DFv = degrees of freedom for values coresponding to compared groups"),
                           h5("p Bonferroni adjusted = Bonferroni correction applied to calculated p values"),
                           h5("p FDR adjusted = false discovery rate (FDR) correction applied to calculated p values")
    )
    
  ))
  
  
  
  server <- function(input, output, session){ 
    
    #####
    # Read input table and prepate it for further processing
    
    # prepare input_model
    input_model <- reactive({
      input_table <- read.csv("input_model.csv", header=TRUE, sep=",")
      input_table
    })
    output$input_model_download <- downloadHandler(
      filename = function(){paste("input_model_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(input_model(), fname) })
    
    # read table input from user
    myData <- reactive({
      inFile <- input$file
      if(is.null(inFile)) return(NULL)
      data_user <- read.csv(inFile$datapath, header = TRUE)
      colnames(data_user) <- c("group", "case", "type", "gene", "RB", "RT", "Ct")
      data_user
    })
    
    # modify columns names of input from user
    myData_original_file <- reactive({
      inFile <- input$file
      if(is.null(inFile)) return(NULL)
      data_user <- read.csv(inFile$datapath, header = TRUE)
      colnames(data_user) <- c("Sample Group", "Case", "Type", "Gene", "BR", "TR", "Ct")
      data_user
    })
    
    # output the original user inpput table
    output$original <- DT::renderDataTable( 
      DT::datatable(myData_original_file(),
                    filter="top", options = list(pageLength = -1,
                                                 lengthMenu = list(c(10,50,100,-1),c('10','50','100','All')),
                                                 columnDefs = list(list(className = 'dt-center', targets = 0:7)),
                                                 dom = "ft")))
    
    
    #############################################################################
    #####
    # Create custom functions
    
    # declare SE function
    SE <- function(x){
      se <- sd(x)/sqrt(length(x))
      return(se)}
    
    # declare CI function
    # CI lower
    CI_lower <- function(x){
      ci_lower <- mean(x) - SE(x)
      return(ci_lower)}
    # CI upper
    CI_upper <- function(x){
      ci_upper <- mean(x) + SE(x)
      return(ci_upper)} 
    
    # declare decimals function
    decimals <- function(x, k) trimws(format(round(x, k), nsmall=k))
    
    
    #############################################################################
    #####
    # Calculate Dct values using a list of tables 
    
    processedTable_Dct <- reactive({
      inFile <- input$file
      #extension <- tools::file_ext(inFile$name)
      #filepath <- inFile$datapath
      #data_raw <- switch(extension,
      #                   csv = readr::read_csv(filepath),
      #                   xls = readxl::read_xls(filepath, sheet=1),
      #                   xlsx = readxl::read_xlsx(filepath, sheet=1))
      
      data_raw <- read.csv(inFile$datapath, header = TRUE)
      data <- as.data.frame(data_raw)
      colnames(data) <- c("group", "case", "type", "gene", "RB", "RT", "Ct")
      
      # declare Dct function
      Dct <- function(x){
        genelist <<- list()
        s <- 1
        GOI <<- data %>%
          filter(group == x, type == "GOI")
        HK <<- data %>%
          filter(group == x, type == "HK")
        for(g in unique(GOI$gene)){
          geneGOI <<- GOI %>%
            filter(gene == g)
          if(max(geneGOI$RB) == max(HK$RB)){
            iteration <<- data.frame(matrix(NA, nrow=max(geneGOI$RT)*max(HK$RT), ncol=max(geneGOI$RB)+3))
            colnames(iteration) <<- c("group","case","gene", rep(c("Dct.RB."), times=max(geneGOI$RB)))}
          if(max(geneGOI$RB) > max(HK$RB)){
            iteration <<- data.frame(matrix(NA, nrow=max(gene$RT)*max(HK$RT), ncol=max(HK$RB)+3))
            colnames(iteration) <<- c("group","case","gene", rep(c("Dct.RB."), times=max(HK$RB)))}
          if(max(geneGOI$RB) < max(HK$RB)){
            iteration <<- data.frame(matrix(NA, nrow=max(geneGOI$RT)*max(HK$RT), ncol=max(geneGOI$RB)+3))
            colnames(iteration) <<- c("group","case","gene", rep(c("Dct.RB."),times=max(geneGOI$RB)))}
          for(b in seq(1,max(geneGOI$RB),1)){
            if( any(HK$RB == b)){
              RB_GOI <<- geneGOI %>%
                filter(geneGOI$RB == b)
              RB_HK <<- HK %>%
                filter(HK$RB == b)
              r <- 1
              for(t_goi in seq(1,max(RB_GOI$RT),1)){
                for(t_hk in seq(1,max(RB_HK$RT),1)){
                  iteration[r,3+b] <<- RB_GOI$Ct[t_goi] - RB_HK$Ct[t_hk]
                  iteration$group[r] <<- RB_GOI$group[1]
                  iteration$case[r] <<- RB_GOI$case[1]
                  iteration$gene[r] <<- RB_GOI$gene[1]
                  r <- r+1 }}}
            else{break}
            colnames(iteration)[3+b] <<- paste(colnames(iteration)[3+b],b,sep="")}
          genelist[[s]] <<- iteration
          s <- s+1 }}
      
      # function for calculating the number of rows required for comparisons
      r <- c(0)
      for(i in unique(data$group)){
        iteration <- data %>%
          filter(group == i, type == "GOI")
        r <- (r + length(unique(iteration$gene)))}
      
      ### 
      # Loop Dct function nicely
      listDct <- list()
      listDct_ordered <- list()
      s <- c(1)
      m <- c(1)
      for(i in unique(data$group)){
        Dct(i)
        listDct[[m]] <- genelist 
        for(t in seq(1,length(listDct[[m]]),1)){
          Dct_table <- data.frame(listDct[[m]][[t]]$group, listDct[[m]][[t]]$case, listDct[[m]][[t]]$gene)
          colnames(Dct_table) <- c("group","case","gene")
          for(p in seq(4,ncol(listDct[[m]][[t]]),1)){
            Dct_table <- cbind(Dct_table, listDct[[m]][[t]][p])
            colnames(Dct_table)[p] <- colnames(listDct[[m]][[t]][p])}
          listDct_ordered[[s]] <- Dct_table
          s <- s+1
          Dct2 <- c() 
        }
        m <- m+1 
      }
      #listDct_ordered
      
      listDct_raw <- listDct_ordered
      listDct_ordered <- list()
      for(l in seq(1,length(listDct_raw),1)){
        if( isTRUE(grep("experimental", unique(listDct_raw[[l]]$case)) == 1 )) {
          listDct_ordered <- list.append(listDct_ordered , listDct_raw[[l]]) }}
      for(l in seq(1,length(listDct_raw),1)){
        if( isTRUE(grep("control", unique(listDct_raw[[l]]$case)) == 1 )) {
          listDct_ordered <- list.append(listDct_ordered , listDct_raw[[l]]) }}
      
      listDct_ordered
    })
    
    
    #############################################################################
    #####
    # Create 2^(-Dct) values using a list of tables
    
    calculationTable_2Dct <- reactive({
      table_2Dct <- processedTable_Dct()
      
      for(i in seq(1,length(table_2Dct),1)){
        for(j in seq(4,ncol(table_2Dct[[i]]),1)){
          table_2Dct[[i]][,j] <- 2^(-table_2Dct[[i]][,j])
        }
      }
      
      table_2Dct
    })
    
    
    #############################################################################
    ######
    # Custom functions for various future tables
    
    # Calculate the number of columns (comparisons) for future background tables
    coloane <- reactive({
      coloane <- c(0)
      data <- as.data.frame(myData())
      
      for(i in unique(data$group)){
        iteration <- data %>%
          filter(group == i, type == "GOI")
        coloane <- (coloane + length(unique(iteration$gene)))} 
      coloane
    })
    
    # Calculate the number of biological replicates
    RB <- reactive({
      RB <- c(0)
      data <- as.data.frame(myData())
      RB <- max(data$RB)
      RB
    })
    
    # Calculate the number of technical replicates
    RT <- reactive({
      RT <- c(0)
      data <- as.data.frame(myData())
      RT <- max(data$RT)
      RT
    })
    
    # Find the number of decimals to work with
    decVal <- reactive({
      if(input$update == 0){
        decVal <- c(3)
        return(decVal) }
      
      else{
        a <- eventReactive(input$update,{
          decVal2 <- input$decNum
          decVal2 }) }
      
      if(is.null(a())){
        decVal <- c(3) }
      else{ decVal <- a() }
      
      decVal
    })
    
    # Ask for input the significance threshold
    pval <- reactive({
      if(input$update2 == 0){
        pval <- c(0.05)
        return(pval) }
      
      else{
        a <- eventReactive(input$update2, {
          pval2 <- input$pvalue
          pval2 }) }
      
      if(is.null(a())){
        pval <- c(0.05) }
      else{pval <- a()}
      
      pval
    })
    
    
    #############################################################################
    #####
    # Create output table of Dct values
    outputDct <- reactive({
      listDct_ordered <- processedTable_Dct()
      Dct_vector <- c()
      if(RB() == 1){
        outputDct_table <- data.frame(matrix(NA,nrow=(RT())^2,ncol=coloane())) }
      else{
        outputDct_table <- data.frame(matrix(NA,nrow=(RB()^RT())^2,ncol=coloane())) }
      for(i in seq(1, length(listDct_ordered),1)){
        colnames(outputDct_table)[i] <- paste(listDct_ordered[[i]]$group[1]," ", listDct_ordered[[i]]$gene[1],sep="")
        for(j in seq(4, ncol(listDct_ordered[[i]]),1)){
          Dct_vector <- append(Dct_vector,listDct_ordered[[i]][,j]) }
        
        outputDct_table[,i] <- Dct_vector[seq_len(nrow(outputDct_table))]
        Dct_vector <- c()
      }  
      
      outputDct_table <- outputDct_table[rowSums(is.na(outputDct_table)) != ncol(outputDct_table), ]
      outputDct_table
    })
    
    
    #############################################################################
    #####
    # Create names for checkboxes
    observeEvent(input$file, {
      req(outputDct())
      updateCheckboxGroupInput(
        session,
        "col_n",
        choices = colnames(outputDct()),
        selected = colnames(outputDct()))
    })
    
    
    #############################################################################
    #####
    # Prepare output Dct table (from `Download files` tab)
    outputDct_decimals <- reactive({
      req(input$col_n)
      outputDct_decimals <- outputDct() %>% select_at(input$col_n)
      for(i in seq(1,ncol(outputDct_decimals),1)){
        outputDct_decimals[,i] <- decimals(outputDct_decimals[,i],decVal()) }
      outputDct_decimals
    })
    
    output$outputDctvalues <- downloadHandler(
      filename = function(){paste("Dct_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(outputDct_decimals(), fname) } )
    
    
    #############################################################################
    #####
    # Prepare output table of 2^(-Dct) values
    output2Dct <- reactive({
      req(input$col_n)
      outputDct_table <- outputDct() %>% select_at(input$col_n)
      output2Dct_table <- 2^(-outputDct_table)
      output2Dct_table
    })
    
    output2Dct_decimals <- reactive({
      output2Dct_decimals <- output2Dct()
      for(i in seq(1,ncol(output2Dct_decimals),1)){
        output2Dct_decimals[,i] <- decimals(output2Dct_decimals[,i],decVal()) }
      output2Dct_decimals
    })
    
    output$output2Dctvalues <- downloadHandler(
      filename = function(){paste("2^-Dct_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(output2Dct_decimals(), fname) })
    
    
    #############################################################################
    #####
    # Fix colour based on the number of comparisons
    custom_colours <- reactive({
      output2Dct_table <- output2Dct()
      nume_coloane <- colnames(output2Dct_table)
      culori <- brewer.pal(ncol(output2Dct_table), "Set1")
      names(culori) <- levels(nume_coloane)
      custom_colours <- scale_colour_manual(values=culori)
      custom_colours
    })
    
    
    #############################################################################
    #####
    # Summary statistics
    
    # Create summary table according to selected input
    summary2Dct <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{
        output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      summary2Dct <- data.frame(matrix(NA, nrow=ncol(output2Dct_table), ncol=11))
      colnames(summary2Dct) <- c("Sample group","Mean","Median","SD","Variance","SE",
                                 "CI-lower","CI-upper","IQR","Minimum","Maximum")
      
      for(i in seq(1, ncol(output2Dct_table), 1)){
        summary2Dct$`Sample group`[i] <- colnames(output2Dct_table[i])
        summary2Dct$Mean[i] <- decimals(mean(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$Median[i] <- decimals(median(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$SD[i] <- decimals(sd(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$Variance[i] <- decimals(var(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$SE[i] <- decimals(SE(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$`CI-lower`[i] <- decimals(CI_lower(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$`CI-upper`[i] <- decimals(CI_upper(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$IQR[i] <- decimals(IQR(c(output2Dct_table[,i]), na.rm = TRUE),decVal())
        summary2Dct$Minimum[i] <- decimals(min(c(na.omit(output2Dct_table[,i]))),decVal())
        summary2Dct$Maximum[i] <- decimals(max(c(na.omit(output2Dct_table[,i]))),decVal())
      }
      
      outputDct <- outputDct()
      my_colours <- brewer.pal(ncol(outputDct), "Set1")
      names(my_colours) <- levels(factor(c(colnames(outputDct))))
      my_scale <- scale_fill_manual(values=my_colours)
      
      # Declare boxplot function 
      boxplot <- reactive({
        boxplot_2Dct <- na.omit(stack(output2Dct_table)) %>%
          ggplot(aes(x=ind, y=values, fill=ind)) +
          geom_boxplot(width=0.5, color="black") +
          geom_jitter(width=0.2, color="black", size=2) +
          theme_bw() +
          theme(
            legend.position = "none",
            plot.title = element_text(size=15, face="bold"),
            axis.title.y = element_text(size=20, face="bold", margin=margin(r=20)),
            axis.text.y = element_text(size=15, face="bold") ) +
          #scale_x_discrete(guide = guide_axis(angle = 0)) +
          my_scale
        
        boxplot_2Dct
      })
      
      # Declare violin plot function
      violin <- reactive({
        violin_2Dct <- na.omit(stack(output2Dct_table)) %>%
          ggplot( aes(x=ind, y=values, fill=ind)) +
          geom_violin(trim=FALSE, width=0.75) +
          geom_boxplot(width=0.1, color="black") +
          geom_jitter(width=0.25, color="black", size=2, alpha=0.9) +
          theme_bw() +
          theme(
            legend.position = "none",
            plot.title = element_text(size=15, face="bold"),
            axis.title.y = element_text(size=20, face="bold", margin=margin(r=20)),
            axis.text.y = element_text(size=15, face="bold") ) +
          my_scale
        
        violin_2Dct
      })
      
      # Declare histogram function
      histogram <- reactive({
        histogram_2Dct <- na.omit(stack(output2Dct_table))
        histogram_list <- vector('list', length = ncol(histogram_2Dct))
        j<-1
        for(i in unique(histogram_2Dct$ind)){
          iteration <- histogram_2Dct %>%
            filter(histogram_2Dct$ind == i)
          histo <- ggplot(iteration, aes(x=values)) + 
            geom_histogram(aes(y=after_stat(density)), position = "identity", binwidth = abs(max(iteration$values))/10, color="black", fill="skyblue")+
            geom_density(color="red", fill="coral1", alpha=0.5) +
            geom_vline(xintercept = mean(iteration$values), color="black", linewidth=1.25) +
            geom_vline(xintercept = median(iteration$values), color="black", linewidth=1.25, linetype="dashed")+
            theme_classic() +
            theme(
              plot.title = element_text(size=20, face="bold"),
              axis.title.x = element_text(size=20, face="bold", margin=margin(t=10)),
              axis.text.x = element_text(size=20, face="bold"),
              axis.title.y = element_text(size=20, face="bold", margin=margin(r=20)),
              axis.text.y = element_text(size=15, face="bold") ) +
            ggtitle(paste(i))
          
          if(input$typeofdata == "2^_transformed"){
            histo <- histo + labs( x = expression(bold("2"^"-\u0394Ct"~"values")), y = "Density" ) }
          else{
            histo <- histo + labs( x = expression(bold("\u0394Ct"~"values")), y = "Density" ) }
          
          histogram_list[[j]] <- histo
          j<-j+1 }
        
        do.call("grid.arrange", c(histogram_list, ncol=2))
      })
      
      # Adapt to user preference of plot orientation
      req(input$plot_orientation)
      if(input$plot_orientation == "horizontal"){
        boxplot2 <- reactive({
          plot <- boxplot() + theme(axis.text.x=element_text(size=20, face="bold", colour="black", hjust = 0.5, margin = margin(t=10), vjust=0.5))
          plot })
        violin2 <- reactive({
          plot <- violin() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5, margin = margin(t=10), vjust=0.5))
          plot })
      }
      if(input$plot_orientation == "vertical"){
        boxplot2 <- reactive({
          plot <- boxplot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 90, vjust = 0.5))
          plot })
        violin2 <- reactive({
          plot <- violin() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 90, vjust = 0.5))
          plot })
      }
      if(input$plot_orientation == "angle30"){
        boxplot2 <- reactive({
          plot <- boxplot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 30, vjust = 0.5))
          plot })
        violin2 <- reactive({
          plot <- violin() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 30, vjust = 0.5))
          plot })
      }
      if(input$plot_orientation == "angle45"){
        boxplot2 <- reactive({
          plot <- boxplot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 45, vjust = 0.5))
          plot })
        violin2 <- reactive({
          plot <- violin() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 45, vjust = 0.5))
          plot })
      }
      if(input$plot_orientation == "angle60"){
        boxplot2 <- reactive({
          plot <- boxplot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 60, vjust = 0.5))
          plot })
        violin2 <- reactive({
          plot <- violin() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 60, vjust = 0.5))
          plot })
      }
      
      
      
      # Output boxplots, violin plots and histograms according to selected input  
      
      if(input$typeofdata == "2^_transformed"){
        
        # boxplots
        boxplot_updated <- reactive({
          boxplot3 <- boxplot2() + labs(y = expression(bold("2"^"-\u0394Ct"~"values")), x = "")
          boxplot3 })
        output$boxplot <- renderPlot({ boxplot_updated() })
        output$download_boxplots <- downloadHandler(
          filename = function(){paste("boxplots_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, boxplot_updated(), dpi = 600, width = 16, height=10 ) })
        
        # violin plots
        violin_updated <- reactive({
          violin3 <- violin2() + labs(y = expression(bold("2"^"-\u0394Ct"~"values")), x = "")
          violin3 })
        output$violin <- renderPlot({ violin_updated() })
        output$download_violin_plots <- downloadHandler(
          filename = function(){paste("violin_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, violin_updated(), dpi = 600, width = 16, height=10 ) })
        
        # histograms
        output$histogram <- renderPlot({ histogram() })
        output$download_histograms <- downloadHandler(
          filename = function(){paste("histograms_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, histogram(), dpi = 600, width = 16, height=10 ) })
      }
      else{
        
        # boxplots
        boxplot_updated <- reactive({
          boxplot3 <- boxplot2() + labs(y = expression(bold("\u0394Ct"~"values")), x = "")
          boxplot3 })
        output$boxplot <- renderPlot({ boxplot_updated() })
        output$download_boxplots <- downloadHandler(
          filename = function(){paste("boxplots_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, boxplot_updated(), dpi = 600, width = 16, height=10 ) })
        
        # violin plots
        violin_updated <- reactive({
          violin3 <- violin2() + labs(y = expression(bold("\u0394Ct"~"values")), x = "")
          violin3 })
        output$violin <- renderPlot({ violin_updated() })
        output$download_violin_plots <- downloadHandler(
          filename = function(){paste("violin_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, violin_updated(), dpi = 600, width = 16, height=10 ) })
        
        # histograms
        output$histogram <- renderPlot({ histogram() })
        output$download_histograms <- downloadHandler(
          filename = function(){paste("histograms_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, histogram(), dpi = 600, width = 16, height=10 ) })
      }
      
      summary2Dct
    })
    
    # Output summary table 
    output$summary <- DT::renderDataTable(
      DT::datatable(summary2Dct(), rownames = FALSE, filter="top", options = list( pageLength = -1,
                                                                                   lengthMenu = list(c(10,50,100,-1),c('10','50','100','All')),
                                                                                   columnDefs = list(list(className = 'dt-center', targets = 1:10)),
                                                                                   dom="ft")))
    
    # Download summary table according to selected input
    observe(
      if(input$typeofdata == "2^_transformed"){
        output$summary_download <- downloadHandler(
          filename = function(){paste("summary_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
          content = function(fname){ write.csv(summary2Dct(), fname) })
      }
      else{
        output$summary_download <- downloadHandler(
          filename = function(){paste("summary_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
          content = function(fname){ write.csv(summary2Dct(), fname) })
      }
    )
    
    
    
    
    #############################################################################
    #####
    # Normality statistics
    
    # Prepare normality table and Q-Q plots
    normality <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{
        output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      normality_table <- data.frame(matrix(NA, nrow=ncol(output2Dct_table)*3, ncol=5))
      colnames(normality_table) <- c("Sample group","Normality test","Statistic","p","Type")
      rand <- c(seq(1,ncol(output2Dct_table)*4,4))
      
      for(i in seq(1, ncol(output2Dct_table),1)){
        normality_table[rand[i],1] <- paste(colnames(output2Dct_table[i]))
        #normality_table[rand[i]+1,1] <- paste(colnames(output2Dct_table[i]))
        #normality_table[rand[i]+2,1] <- paste(colnames(output2Dct_table[i]))
        #normality_table[rand[i]+3,1] <- paste(colnames(output2Dct_table[i]))
        
        iteration <- c(output2Dct_table[,i])
        iteration <- iteration[!is.na(iteration)]
        
        #a <- rnorm(100,mean=mean(iteration),sd=sd(iteration))
        ks <- ks.test(unique(floor(iteration)), "pnorm", mean=mean(iteration), sd=sd(iteration))
        sw <- shapiro.test(iteration)
        ad <- ad.test(iteration)
        #jb <- jarque.bera.test(iteration)
        
        # Kolmogorov-Smirnov
        normality_table[rand[i],2] <- paste(ks$method)
        normality_table[rand[i],3] <- paste("D = ", decimals(ks$statistic,decVal()))
        normality_table[rand[i],4] <- paste(decimals(ks$p.value,5))
        if(ks$p.value < pval()){ normality_table[rand[i],5] <- paste("Parametric") }
        else{ normality_table[rand[i],5] <- paste("Nonparametric") }
        
        # Shapiro-Wilk
        normality_table[rand[i]+1,2] <- paste(sw$method)
        normality_table[rand[i]+1,3] <- paste("W = ", decimals(sw$statistic,decVal()))
        normality_table[rand[i]+1,4] <- paste(decimals(sw$p.value,5))
        if(sw$p.value > pval()){ normality_table[rand[i]+1,5] <- paste("Parametric") }
        else{ normality_table[rand[i]+1,5] <- paste("Nonparametric") }
        
        # Anderson-Darling
        normality_table[rand[i]+2,2] <- paste(ad$method)
        normality_table[rand[i]+2,3] <- paste("A" %p% supsc("2")," = ", decimals(ad$statistic,decVal()))
        normality_table[rand[i]+2,4] <- paste(decimals(ad$p.value,5))
        if(sw$p.value > pval()){  normality_table[rand[i]+2,5] <- paste("Parametric") }
        else{ normality_table[rand[i]+2,5] <- paste("Nonparametric") }
        
        # Jarque-Bera
        #normality_table[rand[i]+3,2] <- paste(jb$method)
        #normality_table[rand[i]+3,3] <- paste("X" %p% supsc("2")," = ", decimals(jb$statistic,decVal()))
        #normality_table[rand[i]+3,4] <- paste(decimals(jb$p.value,5))
        #if(sw$p.value > 0.05){  normality_table[rand[i]+3,5] <- paste("Parametric") }
        #else{normality_table[rand[i]+3,5] <- paste("Nonparametric")}
      }
      
      # Create Q-Q plots
      qqplot <- reactive({
        QQplot_list <- vector('list', length = ncol(output2Dct_table))
        qqplot_table <- na.omit(stack(output2Dct_table))
        j<-1
        for(i in unique(qqplot_table$ind)){
          iteration <- qqplot_table %>%
            filter(qqplot_table$ind == i)
          qqindividual <- ggplot(data = iteration, mapping = aes(sample = values)) +
            stat_qq_band(color="blue") +
            stat_qq_line() +
            stat_qq_point() +
            theme_classic() +
            labs( title = paste(i), x = "Theoretical Quantiles", y = "Sample Quantiles") +
            theme(
              plot.title = element_text(size=15, face="bold"),
              axis.title.x = element_text(size=15, face="bold", margin=margin(t=20)),
              axis.text.x = element_text(size=20, face="bold"),
              axis.title.y = element_text(size=15, face="bold", margin=margin(r=20)),
              axis.text.y = element_text(size=15, face="bold") ) 
          
          QQplot_list[[j]] <- qqindividual
          j <- j + 1
        }
        do.call("grid.arrange", c(QQplot_list, ncol=2))
      })
      
      # Ouput and doanload Q-Q plots
      output$qqplot <- renderPlot({ qqplot() })
      
      if(input$typeofdata == "2^_transformed"){
        output$downloadQQ <- downloadHandler(
          filename = function(){paste("Q-Q_plots_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){
            ggsave(file, qqplot(), dpi = 600, width = 16, height=10 ) }) }
      else{
        output$downloadQQ <- downloadHandler(
          filename = function(){paste("Q-Q_plots_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
          content = function(file){ ggsave(file, qqplot(), dpi = 600, width = 16, height=10 ) }) }
      
      
      
      normality_table
    })
    
    # Output normality table 
    output$normality <- DT::renderDataTable(
      DT::datatable(normality(), rownames = FALSE, filter="top", options = list(pageLength = -1,
                                                                                lengthMenu = list(c(10,50,100,-1),c('10','50','100','All')),
                                                                                columnDefs = list(list(className = 'dt-center', targets = 2:4)),
                                                                                dom="ft")))
    
    
    # Download normality statistics table
    observe(
      if(input$typeofdata == "2^_transformed"){
        output$normality_download <- downloadHandler(
          filename = function(){paste("normality_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
          content = function(fname){ write.csv(normality(), fname) }) }
      else{
        output$normality_download <- downloadHandler(
          filename = function(){paste("normality_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
          content = function(fname){  write.csv(normality(), fname) }) }
    )
    
    
    
    #############################################################################
    #####
    # Two samples tab 
    
    # Prepare the always-on assumptions table
    teste <- reactive({
      teste <- data.frame(matrix(NA, nrow=3,ncol=3))
      colnames(teste) <- c("Test type", "Data normality", "Equal variances")
      teste[1,1] <- c("Student's two-sample t-test")
      teste[2,1] <- c("Student's two-sample t-test Welch corrected")
      teste[3,1] <- c("Mann-Whitney U test (Wilcoxon rank sum)")
      
      teste[1,2] <- paste0(fa(name="circle-check", fill="green"))
      teste[1,3] <- paste0(fa(name="circle-check", fill="green"))
      
      teste[2,2] <- paste0(fa(name="circle-check", fill="green"))
      teste[2,3] <- paste0(fa(name="circle-xmark", fill="red"))
      
      teste[3,2] <- paste0(fa(name="circle-xmark", fill="red"))
      teste[3,3] <- paste0(fa(name="circle-check", fill="green"))
      
      teste
    })
    
    output$types_of_paired_tests <- renderDataTable(
      teste(), escape = FALSE, options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                              columnDefs = list(list(className = 'dt-center', targets = 2:3)),
                                              dom="t"))
    
    ### Student's two-way t-test 
    ttests <- reactive({
      TTEST <- function(x,y){
        t.test(x,y, alternative=c("two.sided"), paired=FALSE, var.equal=TRUE, conf.level=0.95)}
      
      # Create t-test table according to selected input
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      Ttest <- data.frame(matrix(NA, nrow=c(ncol(output2Dct_table)^2 - ncol(output2Dct_table))/2, ncol=8))
      colnames(Ttest) <- c("Comparison","Significance","t","p","df","95%CI-lower","95%CI-upper","Method")
      
      s <- 1
      sqsum <- c()
      for(i in seq(1,ncol(output2Dct_table),1)){
        first <- c(output2Dct_table[,i])
        for(j in seq(1, ncol(output2Dct_table),1)){
          second <- c(output2Dct_table[,j])
          if(i != j){  
            iteration <- i^2 + j^2
            if(any(sqsum == iteration)){ next }
            else{
              sqsum <- append(sqsum, iteration)
              STT <- TTEST(x=first, y=second)
              Ttest$Comparison[s] <- paste0(colnames(output2Dct_table[i])," vs ",colnames(output2Dct_table[j]))
              
              if(STT$p.value <= pval()){
                Ttest$Significance[s] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
              else{
                Ttest$Significance[s] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
              
              Ttest$t[s] <- decimals(STT$statistic,decVal())
              Ttest$p[s] <- STT$p.value
              Ttest$df[s] <- floor(STT$parameter)
              Ttest$`95%CI-lower`[s] <- decimals(STT$conf.int[1],decVal())
              Ttest$`95%CI-upper`[s] <- decimals(STT$conf.int[2],decVal())
              Ttest$Method[s] <- STT$method
              
              s <- s+1
            }}}}
      
      Ttest 
    })
    
    
    ### Welch t-test
    welchttests <- reactive({
      WELCH <- function(x,y){
        t.test(x,y, alternative=c("two.sided"), paired=FALSE,var.equal=FALSE, conf.level=0.95)}
      
      # Create Welch t-test table according to selected input
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      welchTtest <- data.frame(matrix(NA, nrow=c(ncol(output2Dct_table)^2 - ncol(output2Dct_table))/2, ncol=8))
      colnames(welchTtest) <- c("Comparison","Significance","t","p","df","95%CI-lower","95%CI-upper","Method")
      
      s <- 1
      sqsum <- c()
      for(i in seq(1,ncol(output2Dct_table),1)){
        first <- c(output2Dct_table[,i])
        for(j in seq(1, ncol(output2Dct_table),1)){
          second <- c(output2Dct_table[,j])
          if(i != j){ 
            iteration <- i^2 + j^2
            if(any(sqsum == iteration)){ next }
            else{
              sqsum <- append(sqsum, iteration)
              STT <- WELCH(x=first, y=second)
              welchTtest$Comparison[s] <- paste0(colnames(output2Dct_table[i])," vs ",colnames(output2Dct_table[j]))
              
              if(STT$p.value <= pval()){
                welchTtest$Significance[s] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
              else{
                welchTtest$Significance[s] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
              
              welchTtest$t[s] <- decimals(STT$statistic,decVal())
              welchTtest$p[s] <- STT$p.value
              welchTtest$df[s] <- floor(STT$parameter)
              welchTtest$`95%CI-lower`[s] <- decimals(STT$conf.int[1],decVal())
              welchTtest$`95%CI-upper`[s] <- decimals(STT$conf.int[2],decVal())
              welchTtest$Method[s] <- STT$method
              
              s <- s+1
            }}}}
      
      welchTtest 
    })
    
    
    ### Wilcoxon test 
    wilcoxontests <- reactive({
      WILCOXON <- function(x,y){
        wilcox.exact(x,y, alternative=c("two.sided"), paired=FALSE,conf.int=TRUE, conf.level=0.95, exact = TRUE)}
      
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      wilcoxontest <- data.frame(matrix(NA, nrow=c(ncol(output2Dct_table)^2 - ncol(output2Dct_table))/2, ncol=7))
      colnames(wilcoxontest) <- c("Comparison","Significance","W","p","95%CI-lower","95%CI-upper","Method")
      
      s <- 1
      sqsum <- c()
      for(i in seq(1,ncol(output2Dct_table),1)){
        first <- c(output2Dct_table[,i])
        for(j in seq(1, ncol(output2Dct_table),1)){
          second <- c(output2Dct_table[,j])
          if(i != j){ 
            iteration <- i^2 + j^2
            if(any(sqsum == iteration)){ next }
            else{
              sqsum <- append(sqsum, iteration)
              STT <- WILCOXON(x=first, y=second)
              wilcoxontest$Comparison[s] <- paste0(colnames(output2Dct_table[i])," vs ",colnames(output2Dct_table[j]))
              
              if(STT$p.value <= pval()){
                wilcoxontest$Significance[s] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
              else{
                wilcoxontest$Significance[s] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
              
              wilcoxontest$W[s] <- decimals(STT$statistic,decVal())
              wilcoxontest$p[s] <- STT$p.value
              wilcoxontest$`95%CI-lower`[s] <- decimals(STT$conf.int[1],decVal())
              wilcoxontest$`95%CI-upper`[s] <- decimals(STT$conf.int[2],decVal())
              wilcoxontest$Method[s] <- STT$method
              
              s <- s+1
            }}}}
      
      wilcoxontest 
    })
    
    
    # Output two samples tests
    observe(
      if(input$typeofdata == "2^_transformed"){
        if(input$select_paired == "ttest"){
          paired_tests <- reactive({
            paired <- ttests()
            paired })
          output$paired <- DT::renderDataTable(
            DT::datatable(paired_tests(),filter = "top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)), dom="ft"),escape=FALSE ))
          download_paired_tests <- reactive({
            download_paired_tests <- paired_tests()
            for(i in seq(1, nrow(download_paired_tests),1)){
              if(download_paired_tests$p[i] <= pval()){
                download_paired_tests$Significance[i] <- c("Yes") }
              else{ download_paired_tests$Significance[i] <- c("No")}}
            download_paired_tests })
          output$paired_download <- downloadHandler(
            filename = function(){paste("ttest_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_paired_tests(), fname) }) }
        if(input$select_paired == "welchttest"){
          paired_tests <- reactive({
            paired <- welchttests()
            paired })
          output$paired <- DT::renderDataTable(
            DT::datatable(paired_tests(),filter = "top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)), dom="ft"),escape=FALSE ))
          download_paired_tests <- reactive({
            download_paired_tests <- paired_tests()
            for(i in seq(1, nrow(download_paired_tests),1)){
              if(download_paired_tests$p[i] <= pval()){
                download_paired_tests$Significance[i] <- c("Yes") }
              else{ download_paired_tests$Significance[i] <- c("No")}}
            download_paired_tests })
          output$paired_download <- downloadHandler(
            filename = function(){paste("welch_ttest_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_paired_tests(), fname) }) }
        if(input$select_paired == "wilcoxon"){
          paired_tests <- reactive({
            paired <- wilcoxontests()
            paired })
          output$paired <- DT::renderDataTable(
            DT::datatable(paired_tests(),filter = "top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:7)), dom="ft"),escape=FALSE ))
          download_paired_tests <- reactive({
            download_paired_tests <- paired_tests()
            for(i in seq(1, nrow(download_paired_tests),1)){
              if(download_paired_tests$p[i] <= pval()){
                download_paired_tests$Significance[i] <- c("Yes") }
              else{ download_paired_tests$Significance[i] <- c("No")}}
            download_paired_tests })
          output$paired_download <- downloadHandler(
            filename = function(){paste("wilcoxon_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_paired_tests(), fname) }) }
      }
      else{
        if(input$select_paired == "ttest"){
          paired_tests <- reactive({
            paired <- ttests()
            paired })
          output$paired <- DT::renderDataTable(
            DT::datatable(paired_tests(),filter = "top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)), dom="ft"),escape=FALSE ))
          download_paired_tests <- reactive({
            download_paired_tests <- paired_tests()
            for(i in seq(1, nrow(download_paired_tests),1)){
              if(download_paired_tests$p[i] <= pval()){
                download_paired_tests$Significance[i] <- c("Yes") }
              else{ download_paired_tests$Significance[i] <- c("No")}}
            download_paired_tests })
          output$paired_download <- downloadHandler(
            filename = function(){paste("ttest_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_paired_tests(), fname) }) }
        if(input$select_paired == "welchttest"){
          paired_tests <- reactive({
            paired <- welchttests()
            paired })
          output$paired <- DT::renderDataTable(
            DT::datatable(paired_tests(),filter = "top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)),dom="ft"),escape=FALSE ))
          download_paired_tests <- reactive({
            download_paired_tests <- paired_tests()
            for(i in seq(1, nrow(download_paired_tests),1)){
              if(download_paired_tests$p[i] <= pval()){
                download_paired_tests$Significance[i] <- c("Yes") }
              else{ download_paired_tests$Significance[i] <- c("No")}}
            download_paired_tests })
          output$paired_download <- downloadHandler(
            filename = function(){paste("welch_ttest_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_paired_tests(), fname) }) }
        if(input$select_paired == "wilcoxon"){
          paired_tests <- reactive({
            paired <- wilcoxontests()
            paired })
          output$paired <- DT::renderDataTable(
            DT::datatable(paired_tests(),filter = "top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:7)),dom="ft"),escape=FALSE ))
          download_paired_tests <- reactive({
            download_paired_tests <- paired_tests()
            for(i in seq(1, nrow(download_paired_tests),1)){
              if(download_paired_tests$p[i] <= pval()){
                download_paired_tests$Significance[i] <- c("Yes") }
              else{ download_paired_tests$Significance[i] <- c("No")}}
            download_paired_tests })
          output$paired_download <- downloadHandler(
            filename = function(){paste("wilcoxon_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){ write.csv(download_paired_tests(), fname) }) }
      }
    )
    
    
    ##############################################################################
    #####
    # Multiple samples testing
    
    # Prepare the always-on assumptions table
    teste2 <- reactive({
      teste <- data.frame(matrix(NA, nrow=3,ncol=4))
      colnames(teste) <- c("Test type", "Post-hoc test", "Data normality", "Equal variances")
      teste[1,1] <- c("One-way ANOVA")
      teste[2,1] <- c("Welch's ANOVA")
      teste[3,1] <- c("Kruskal-Wallis ANOVA")
      
      teste[1,2] <- c("Tukey HSD test")
      teste[2,2] <- c("Games-Howell test")
      teste[3,2] <- c("Dunn test")
      
      teste[1,3] <- paste0(fa(name="circle-check", fill="green"))
      teste[1,4] <- paste0(fa(name="circle-check", fill="green"))
      
      teste[2,3] <- paste0(fa(name="circle-check", fill="green"))
      teste[2,4] <- paste0(fa(name="circle-xmark", fill="red"))
      
      teste[3,3] <- paste0(fa(name="circle-xmark", fill="red"))
      teste[3,4] <- paste0(fa(name="circle-check", fill="green"))
      
      teste
    })
    
    output$types_of_multiple_tests <- renderDataTable(
      teste2(), escape = FALSE, options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                               columnDefs = list(list(className = 'dt-center', targets = 3:4)),
                                               dom="t"))
    
    
    ### one-way ANOVA
    ANOVA <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      pre_anova <- na.omit(stack(output2Dct_table))
      anova_table  <- data.frame(matrix(NA, nrow = 1, ncol=6))
      colnames(anova_table) <- c("Test","Significance","F","p", "DFg", "DFv")
      # where DFg = degrees of freedom for groups
      # where DFv = degrees of freedom for values
      
      res.aov <- aov(values~ind, data=pre_anova)
      s_anova <- anova_summary(res.aov)
      
      anova_table$Test[1] <- paste("ANOVA between chosen groups")
      anova_table$DFg[1] <- s_anova$DFn
      anova_table$DFv[1] <- s_anova$DFd
      anova_table$F[1] <- s_anova$F
      anova_table$p[1] <- s_anova$p
      
      if(s_anova$p <= pval()){
        anova_table$Significance[1] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
      else{
        anova_table$Significance[1] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
      
      anova_table
    })
    
    ### Tukey HSD test (post-hoc)
    TUKEY <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      pre_anova <- na.omit(stack(output2Dct_table))
      res.aov <- aov(values~ind, data=pre_anova)
      
      TK <- TukeyHSD(res.aov)
      pre_tukey <- data.frame(TK$ind)
      Tukey_table <- data.frame(matrix(NA,nrow=c(ncol(output2Dct_table)^2 - ncol(output2Dct_table))/2, ncol=7))
      colnames(Tukey_table) <- c("Comparison","Significance","p","Difference","Lower","Upper","Method")
      
      comparisons <- c(rownames(pre_tukey))
      Tukey_table$Comparison <- comparisons
      
      for(i in seq(4,ncol(Tukey_table),1)){
        Tukey_table[,i] <- decimals(pre_tukey[,i-3],decVal()) }
      Tukey_table$p <- pre_tukey$p.adj
      for(s in seq(1,nrow(Tukey_table),1)){
        if(Tukey_table$p[s] <= pval()){
          Tukey_table$Significance[s] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
        else{
          Tukey_table$Significance[s] <- paste0(c("No"), "   ", fa(name="circle", fill="red", prefer_type = "solid")) }
      }
      Tukey_table$Method <- c("Tukey multiple comparisons of means")
      
      Tukey_table
    })
    
    
    ### Welch ANOVA
    WELCHANOVA <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      pre_anova <- na.omit(stack(output2Dct_table))
      anova_table  <- data.frame(matrix(NA, nrow = 1, ncol=6))
      colnames(anova_table) <- c("Test","Significance","F","p", "DFg", "DFv")
      # where DFg = degrees of freedom for groups
      # where DFv = degrees of freedom for values
      
      res.aov <- oneway.test(values~ind, data=pre_anova, var.equal = FALSE)
      s_anova <- res.aov
      
      anova_table$Test[1] <- paste("Welch's ANOVA between chosen groups")
      anova_table$DFg[1] <- decimals(s_anova$parameter[1],decVal())
      anova_table$DFv[1] <- decimals(s_anova$parameter[2],decVal())
      anova_table$F[1] <- decimals(s_anova$statistic,decVal())
      anova_table$p[1] <- s_anova$p.value
      
      if(s_anova$p.value <= pval()){
        anova_table$Significance[1] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
      else{
        anova_table$Significance[1] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
      
      anova_table
    })
    
    ### Games-Howell test (post-hoc)
    GAMESHOWELL <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      pre_anova <- na.omit(stack(output2Dct_table))
      GM <- games_howell_test(values~ind, data=pre_anova, conf.level=0.95, detailed = TRUE)
      
      GM_table <- data.frame(matrix(NA,nrow=c(ncol(output2Dct_table)^2 - ncol(output2Dct_table))/2, ncol=8))
      colnames(GM_table) <- c("Comparison","Significance","t","p","Mean difference","95% CI lower","95% CI upper","Method")
      
      for(i in seq(1, nrow(GM_table),1)){
        GM_table[i,1] <- paste0(GM$group1[i], "vs", GM$group2[i]) }
      
      GM_table$t <- decimals(GM$statistic,decVal())
      GM_table$p <- GM$p.adj
      GM_table$`Mean difference` <- decimals(GM$estimate,decVal())
      GM_table$`95% CI lower` <- decimals(GM$conf.low, decVal())
      GM_table$`95% CI upper` <- decimals(GM$conf.high, decVal())
      GM_table$Method <- GM$method
      
      for(s in seq(1,nrow(GM_table),1)){
        if(GM_table$p[s] <= pval()){
          GM_table$Significance[s] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
        else{
          GM_table$Significance[s] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
      }
      
      GM_table
    })
    
    
    ### Kruskal-Wallis ANOVA
    KRUSKALWALLIS <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      pre_KW <- na.omit(stack(output2Dct_table))
      KW_table  <- data.frame(matrix(NA, nrow = 1, ncol=5))
      colnames(KW_table) <- c("Test","Significance","Statistic","p", "df")
      
      KW <- kruskal.test(values~ind, data=pre_KW)
      
      KW_table$Test[1] <- KW$method
      KW_table$Statistic[1] <- decimals(KW$statistic,decVal())
      KW_table$p[1] <- KW$p.value
      KW_table$df[1] <- KW$parameter
      
      if(KW_table$p <= pval()){
        KW_table$Significance[1] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
      else{
        KW_table$Significance[1] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
      
      KW_table
    })
    
    ### Dunn test (post-hoc)
    DUNN <- reactive({
      req(input$col_n)
      req(input$typeofdata)
      if(input$typeofdata == "2^_transformed"){
        output2Dct_table <- output2Dct() %>% select_at(input$col_n) }
      else{ output2Dct_table <- outputDct() %>% select_at(input$col_n) }
      
      pre_kruskal <- na.omit(stack(output2Dct_table))
      dunn <- dunn_test(data=pre_kruskal, values~ind, detailed = TRUE)
      dunn_bonferroni <- dunn_test(data=pre_kruskal, values~ind, detailed = TRUE, p.adjust.method = "bonferroni")
      dunn_fdr <- dunn_test(data=pre_kruskal, values~ind, detailed = TRUE, p.adjust.method = "fdr")
      
      dunn_table <- data.frame(matrix(NA,nrow=c(ncol(output2Dct_table)^2 - ncol(output2Dct_table))/2, ncol=8))
      colnames(dunn_table) <- c("Comparison","Significance","Z","p","p Bonferroni adjusted","p FDR adjusted", "Mean rank difference","Method")
      
      for(i in seq(1, nrow(dunn_table),1)){
        dunn_table[i,1] <- paste0(dunn$group1[i], "vs", dunn$group2[i]) }
      
      dunn_table$Z <- decimals(dunn$statistic,decVal())
      dunn_table$p <- dunn$p
      dunn_table$`p Bonferroni adjusted` <- dunn_bonferroni$p.adj
      dunn_table$`p FDR adjusted` <- dunn_fdr$p.adj
      dunn_table$`Mean rank difference` <- decimals(dunn$estimate,decVal())
      dunn_table$Method <- dunn$method
      
      for(s in seq(1,nrow(dunn_table),1)){
        if(dunn_table$p[s] <= pval()){
          dunn_table$Significance[s] <- paste0(c("Yes"),"   ",fa(name="circle", fill="green", prefer_type = "solid")) }
        else{
          dunn_table$Significance[s] <- paste0(c("No"), "   ",fa(name="circle", fill="red", prefer_type = "solid")) }
      }
      
      dunn_table
    })
    
    
    # Output multiple samples tests
    observe(
      if(input$typeofdata == "2^_transformed"){
        if(input$select_multiple == "anova"){
          multiple_test <- reactive({
            multiple <- ANOVA()
            multiple })
          posthoc_tests <- reactive({
            posthoc <- TUKEY()
            posthoc })
          output$multiple_table <- DT::renderDataTable(
            DT::datatable(multiple_test(),
                          options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                         columnDefs = list(list(className = 'dt-center', targets = 2:6)),dom="ft"),escape=FALSE ))
          download_multiple_test <- reactive({
            download_multiple_tests <- multiple_test()
            for(i in seq(1, nrow(download_multiple_tests),1)){
              if(download_multiple_tests$p[i] <= pval()){
                download_multiple_tests$Significance[i] <- c("Yes") }
              else{ download_multiple_tests$Significance[i] <- c("No")}}
            download_multiple_tests })
          output$multiple_download <- downloadHandler(
            filename = function(){paste("anova_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){ write.csv(download_multiple_test(), fname) })
          output$posthoc_choice <- renderText({ paste0(h3("Tukey HSD post-hoc test")) })
          output$posthoc_table <- DT::renderDataTable(
            DT::datatable(posthoc_tests(),filter="top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:7)),dom="ft"),escape=FALSE))
          download_posthoc_tests <- reactive({
            download_posthoc_tests <- posthoc_tests()
            for(i in seq(1, nrow(download_posthoc_tests),1)){
              if(download_posthoc_tests$p[i] <= pval()){
                download_posthoc_tests$Significance[i] <- c("Yes") }
              else{ download_posthoc_tests$Significance[i] <- c("No")}}
            download_posthoc_tests })
          output$posthoc_download <- downloadHandler(
            filename = function(){paste("tukey_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_posthoc_tests(), fname) })
        }
        if(input$select_multiple == "welchanova"){
          multiple_test <- reactive({
            multiple <- WELCHANOVA()
            multiple })
          posthoc_tests <- reactive({
            posthoc <- GAMESHOWELL()
            posthoc })
          output$multiple_table <- DT::renderDataTable(
            DT::datatable(multiple_test(),
                          options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                         columnDefs = list(list(className = 'dt-center', targets = 2:6)),dom="ft"),escape=FALSE))
          download_multiple_test <- reactive({
            download_multiple_tests <- multiple_test()
            for(i in seq(1, nrow(download_multiple_tests),1)){
              if(download_multiple_tests$p[i] <= pval()){
                download_multiple_tests$Significance[i] <- c("Yes") }
              else{ download_multiple_tests$Significance[i] <- c("No")}}
            download_multiple_tests })
          output$multiple_download <- downloadHandler(
            filename = function(){paste("welch_anova_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){ write.csv(download_multiple_test(), fname) })
          output$posthoc_choice <- renderText({ paste0(h3("Games-Howell post-hoc test")) })
          output$posthoc_table <- DT::renderDataTable(
            DT::datatable(posthoc_tests(),filter="top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)),dom="ft"),escape=FALSE))
          download_posthoc_tests <- reactive({
            download_posthoc_tests <- posthoc_tests()
            for(i in seq(1, nrow(download_posthoc_tests),1)){
              if(download_posthoc_tests$p[i] <= pval()){
                download_posthoc_tests$Significance[i] <- c("Yes") }
              else{ download_posthoc_tests$Significance[i] <- c("No")}}
            download_posthoc_tests })
          output$posthoc_download <- downloadHandler(
            filename = function(){paste("games_howell_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){ write.csv(download_posthoc_tests(), fname) })
        }
        if(input$select_multiple == "kruskalwallis"){
          multiple_test <- reactive({
            multiple <- KRUSKALWALLIS()
            multiple })
          posthoc_tests <- reactive({
            posthoc <- DUNN()
            posthoc })
          output$multiple_table <- DT::renderDataTable(
            DT::datatable(multiple_test(),
                          options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                         columnDefs = list(list(className = 'dt-center', targets = 2:5)),dom="ft"),escape=FALSE ))
          download_multiple_test <- reactive({
            download_multiple_tests <- multiple_test()
            for(i in seq(1, nrow(download_multiple_tests),1)){
              if(download_multiple_tests$p[i] <= pval()){
                download_multiple_tests$Significance[i] <- c("Yes") }
              else{ download_multiple_tests$Significance[i] <- c("No")}}
            download_multiple_tests })
          output$multiple_download <- downloadHandler(
            filename = function(){paste("kw_anova_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_multiple_test(), fname) })
          output$posthoc_choice <- renderText({ paste0(h3("Dunn post-hoc test")) })
          output$posthoc_table <- DT::renderDataTable(
            DT::datatable(posthoc_tests(),filter="top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)),dom="ft"),escape=FALSE))
          download_posthoc_tests <- reactive({
            download_posthoc_tests <- posthoc_tests()
            for(i in seq(1, nrow(download_posthoc_tests),1)){
              if(download_posthoc_tests$p[i] <= pval()){
                download_posthoc_tests$Significance[i] <- c("Yes") }
              else{ download_posthoc_tests$Significance[i] <- c("No")}}
            download_posthoc_tests })
          output$posthoc_download <- downloadHandler(
            filename = function(){paste("dunn_table_linear_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_posthoc_tests(), fname) })
        }
      }
      else{
        if(input$select_multiple == "anova"){
          multiple_test <- reactive({
            multiple <- ANOVA()
            multiple })
          posthoc_tests <- reactive({
            posthoc <- TUKEY()
            posthoc })
          output$multiple_table <- DT::renderDataTable(
            DT::datatable(multiple_test(),
                          options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                         columnDefs = list(list(className = 'dt-center', targets = 2:6)),dom="ft"),escape=FALSE ))
          download_multiple_test <- reactive({
            download_multiple_tests <- multiple_test()
            for(i in seq(1, nrow(download_multiple_tests),1)){
              if(download_multiple_tests$p[i] <= pval()){
                download_multiple_tests$Significance[i] <- c("Yes") }
              else{ download_multiple_tests$Significance[i] <- c("No")}}
            download_multiple_tests })
          output$multiple_download <- downloadHandler(
            filename = function(){paste("anova_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_multiple_test(), fname) })
          output$posthoc_choice <- renderText({ paste0(h3("Tukey HSD post-hoc test")) })
          output$posthoc_table <- DT::renderDataTable(
            DT::datatable(posthoc_tests(),filter="top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:7)),dom="ft"),escape=FALSE))
          download_posthoc_tests <- reactive({
            download_posthoc_tests <- posthoc_tests()
            for(i in seq(1, nrow(download_posthoc_tests),1)){
              if(download_posthoc_tests$p[i] <= pval()){
                download_posthoc_tests$Significance[i] <- c("Yes") }
              else{ download_posthoc_tests$Significance[i] <- c("No")}}
            download_posthoc_tests })
          output$posthoc_download <- downloadHandler(
            filename = function(){paste("tukey_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_posthoc_tests(), fname) })
        }
        if(input$select_multiple == "welchanova"){
          multiple_test <- reactive({
            multiple <- WELCHANOVA()
            multiple })
          posthoc_tests <- reactive({
            posthoc <- GAMESHOWELL()
            posthoc })
          output$multiple_table <- DT::renderDataTable(
            DT::datatable(multiple_test(),
                          options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                         columnDefs = list(list(className = 'dt-center', targets = 2:6)),dom="ft"),escape=FALSE ))
          download_multiple_test <- reactive({
            download_multiple_tests <- multiple_test()
            for(i in seq(1, nrow(download_multiple_tests),1)){
              if(download_multiple_tests$p[i] <= pval()){
                download_multiple_tests$Significance[i] <- c("Yes") }
              else{ download_multiple_tests$Significance[i] <- c("No")}}
            download_multiple_tests })
          output$multiple_download <- downloadHandler(
            filename = function(){paste("welch_anova_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_multiple_test(), fname) })
          output$posthoc_choice <- renderText({ paste0(h3("Games-Howell post-hoc test")) })
          output$posthoc_table <- DT::renderDataTable(
            DT::datatable(posthoc_tests(),filter="top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)),dom="ft"),escape=FALSE))
          download_posthoc_tests <- reactive({
            download_posthoc_tests <- posthoc_tests()
            for(i in seq(1, nrow(download_posthoc_tests),1)){
              if(download_posthoc_tests$p[i] <= pval()){
                download_posthoc_tests$Significance[i] <- c("Yes") }
              else{ download_posthoc_tests$Significance[i] <- c("No")}}
            download_posthoc_tests })
          output$posthoc_download <- downloadHandler(
            filename = function(){paste("games_howell_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){ write.csv(download_posthoc_tests(), fname) })
          
        }
        if(input$select_multiple == "kruskalwallis"){
          multiple_test <- reactive({
            multiple <- KRUSKALWALLIS()
            multiple })
          posthoc_tests <- reactive({
            posthoc <- DUNN()
            posthoc })
          output$multiple_table <- DT::renderDataTable(
            DT::datatable(multiple_test(),
                          options = list(lengthChange = FALSE, searching=FALSE, paging=FALSE,
                                         columnDefs = list(list(className = 'dt-center', targets = 2:5)),dom="ft"),escape=FALSE ))
          download_multiple_test <- reactive({
            download_multiple_tests <- multiple_test()
            for(i in seq(1, nrow(download_multiple_tests),1)){
              if(download_multiple_tests$p[i] <= pval()){
                download_multiple_tests$Significance[i] <- c("Yes") }
              else{ download_multiple_tests$Significance[i] <- c("No")}}
            download_multiple_tests })
          output$multiple_download <- downloadHandler(
            filename = function(){paste("kw_anova_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_multiple_test(), fname) })
          output$posthoc_choice <- renderText({ paste0(h3("Dunn post-hoc test")) })
          output$posthoc_table <- DT::renderDataTable(
            DT::datatable(posthoc_tests(),filter="top",
                          options = list(pageLength=-1,lengthMenu=list(c(10,50,100,-1),c('10','50','100','All')),
                                         columnDefs = list(list(className = 'dt-center', targets = 2:8)),dom="ft"),escape=FALSE))
          download_posthoc_tests <- reactive({
            download_posthoc_tests <- posthoc_tests()
            for(i in seq(1, nrow(download_posthoc_tests),1)){
              if(download_posthoc_tests$p[i] <= pval()){
                download_posthoc_tests$Significance[i] <- c("Yes") }
              else{ download_posthoc_tests$Significance[i] <- c("No")}}
            download_posthoc_tests })
          output$posthoc_download <- downloadHandler(
            filename = function(){paste("dunn_table_untransformed_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
            content = function(fname){  write.csv(download_posthoc_tests(), fname) })
        }
      }
    )
    
    
    #############################################################################
    #####
    # Prepare Livak function and DDct table
    DDct_table <- reactive({
      req(input$col_n)
      outputDct_table <- outputDct() %>% select_at(input$col_n)
      DDct <- c()
      Livak <- function(x,y){
        for(a in seq(1, nrow(outputDct_table[x]),1)){
          for(b in seq(1, nrow(outputDct_table[y]),1)){
            iteration <- c(outputDct_table[a,x] - outputDct_table[b,y])
            DDct <<- append(DDct,iteration)
            DDct <<- na.omit(DDct)
          }}}
      
      # Create raw DDct table 
      if(RB() == 1){
        DDct_table <- data.frame(matrix(NA, nrow=c(RT()^4), ncol=c(ncol(outputDct_table)^2-ncol(outputDct_table))/2)) }
      else{
        DDct_table <- data.frame(matrix(NA, nrow=c((RB()^RT())^2), ncol=c(ncol(outputDct_table)^2-ncol(outputDct_table))/2)) }
      
      
      relevant_comparisons <- data.frame(matrix(NA, nrow=c((ncol(outputDct_table)^2-ncol(outputDct_table))/2), ncol=4))
      colnames(relevant_comparisons) <- c("comp1", "comp2", "n1", "n2")
      p <- 1
      for(i in seq(1,ncol(outputDct_table),1)){
        for(j in seq(1,ncol(outputDct_table),1)){
          if(i != j){
            if(isTRUE(any(relevant_comparisons$n2 == i))){
              if(isTRUE(any(relevant_comparisons$n1 == j))){ next }
              }
            relevant_comparisons$comp1[p] <- colnames(outputDct_table)[i]
            relevant_comparisons$comp2[p] <- colnames(outputDct_table)[j]
            relevant_comparisons$n1[p] <- i
            relevant_comparisons$n2[p] <- j
            p <- p + 1
          }}}
      
      
      
      DDct <- c()
      for(r in seq(1, nrow(relevant_comparisons),1)){
        i <- relevant_comparisons$n1[r]
        j <- relevant_comparisons$n2[r]
        
        Livak(x=i, y=j)
        colnames(DDct_table)[r] <- paste(colnames(outputDct_table)[i], " vs ", colnames(outputDct_table)[j])
        livak <- 2^(-DDct)
        if (length(livak) == length(DDct_table[,r])) {
          DDct_table[,r] <- c(livak) }
        else {
          DDct_table[,r] <- c(livak, rep(NA, length(DDct_table[,r]) - length(livak))) }
        DDct <- c()
        }
      
      DDct_table
      
    })
    
    # Update 2^-DDct table to number of selected decimals
    DDct_table_decimals <- reactive({
      DDct_table_decimals <- DDct_table()
      for(i in seq(1,ncol(DDct_table_decimals),1)){
        DDct_table_decimals[,i] <- decimals(DDct_table_decimals[,i],decVal()) }
      DDct_table_decimals
    })
    
    # Download 2^-DDct table
    output$outputDDctvalues <- downloadHandler(
      filename = function(){paste("2^-DDct_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){  write.csv(DDct_table_decimals(), fname) })
    
    
    #############################################################################
    ##### 
    #Fold Change tab 
    
    # FC table based on 2^-DDct values
    foldchange <- reactive({
      DDct_table <- DDct_table()
      DDct_table <- na.omit(DDct_table)
      
      # Plotting summary DDct
      FCtable <- data.frame(matrix(NA, nrow=ncol(DDct_table), ncol=12))
      colnames(FCtable) <- c("Comparison","Type","FC","Median","SD","Variance","SE","CI-lower","CI-upper","IQR","Minimum","Maximum")
      
      for(i in seq(1, ncol(DDct_table),1)){
        FCtable$Comparison[i] <- paste(colnames(DDct_table)[i])
        FCtable$FC[i] <- mean(na.omit(DDct_table[,i]))
        FCtable$Median[i] <- median(na.omit(DDct_table[,i]))
        FCtable$SD[i] <- sd(na.omit(DDct_table[,i]))
        FCtable$Variance[i] <- var(na.omit(DDct_table[,i]))
        FCtable$SE[i] <- SE(na.omit(DDct_table[,i]))
        FCtable$`CI-lower`[i] <- as.numeric(FCtable$FC[i] - FCtable$SE[i])
        FCtable$`CI-upper`[i] <- as.numeric(FCtable$FC[i] + FCtable$SE[i])
        FCtable$IQR[i] <- IQR(DDct_table[,i], na.rm = TRUE)
        FCtable$Minimum[i] <- min(na.omit(DDct_table[,i]))
        FCtable$Maximum[i] <- max(na.omit(DDct_table[,i]))
        
        if(FCtable$FC[i] > 1){
          FCtable$Type[i] <- paste0(c("Up-regulated"),"   ",fa(name="arrow-up", fill="blue", prefer_type = "solid")) }
        else{
          FCtable$Type[i] <- paste0(c("Down-regulated"),"   ",fa(name="arrow-down", fill="red", prefer_type = "solid")) }
      }
      
      for(i in seq(3,ncol(FCtable),1)){
        FCtable[,i] <- decimals(FCtable[,i],decVal()) }
      
      # Prepare FC plot
      FCtable_plot <- FCtable
      FCtable_plot$Comparison <- factor(FCtable_plot$Comparison, levels=FCtable_plot$Comparison)
      FC_plot <- reactive({
        req(input$col_n)
        ggplot(FCtable_plot, aes(x=Comparison, y=as.numeric(FC)), stat="identity") +
          geom_bar(stat="identity", aes(fill = FC < 1 )) +
          scale_fill_manual(guide=FALSE, breaks = c(TRUE,FALSE), values=c("tomato1", "dodgerblue2")) +
          geom_errorbar(aes(x=Comparison, ymin=c(as.numeric(FC)-as.numeric(SE)), ymax=c(as.numeric(FC)+as.numeric(SE))),
                        width=0.8, colour="black", alpha=0.9, linewidth=1) +
          geom_hline(yintercept = 1, color="black", linewidth=1, linetype = "dashed") +
          scale_x_discrete(labels = gsub('\\s vs','\n vs', gsub('vs \\s', 'vs \n' , FCtable_plot$Comparison))) +
          theme_bw() +
          theme(
            legend.position = "none", 
            plot.title = element_text(size=20, face="bold"),
            axis.title.x = element_text(size=20, face="bold"),
            axis.title.y = element_text(size=20, face="bold"),
            axis.text.y = element_text(size=15, face="bold") ) +
          labs(y = expression(bold("FC",'\n')), x = "\n\n Comparison" ) 
      })
      
      # Adapt to user preference of plot orientation
      req(input$plot_orientation)
      if(input$plot_orientation == "horizontal"){
        FC_plot_updated <- reactive({
          plot <- FC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5, margin = margin(t=10)))
          plot }) }
      if(input$plot_orientation == "vertical"){
        FC_plot_updated <- reactive({
          plot <- FC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 90, vjust = 0.5))
          plot }) }
      if(input$plot_orientation == "angle30"){
        FC_plot_updated <- reactive({
          plot <- FC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 30, vjust = 0.5))
          plot }) }
      if(input$plot_orientation == "angle45"){
        FC_plot_updated <- reactive({
          plot <- FC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 45, vjust = 0.5))
          plot }) }
      if(input$plot_orientation == "angle60"){
        FC_plot_updated <- reactive({
          plot <- FC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 60, vjust = 0.5))
          plot }) }
      
      # Output and download FC plot
      output$plotFC <- renderPlot({ FC_plot_updated() })
      output$downloadFC <- downloadHandler(
        filename = function(){paste("FC_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
        content = function(file){ ggsave(file, plot=FC_plot_updated(), dpi = 600, width = 16, height=10) } )
      
      FCtable
    })
    
    # output FC table
    output$fc <- DT::renderDataTable(
      DT::datatable(foldchange(), filter = "top",
                    options = list( pageLength = -1,
                                    lengthMenu = list(c(10,50,100,-1),c('10','50','100','All')),
                                    columnDefs = list(list(className = 'dt-center', targets = 2:12)),dom="ft"), escape = FALSE))
    
    # download FC table
    download_foldchange <- reactive({
      download_FC_table <- foldchange()
      for(i in seq(1,nrow(download_FC_table),1)){
        if(download_FC_table$FC[i] > 1){
          download_FC_table$Type[i] <- c("Up-regulated") }
        else{
          download_FC_table$Type[i] <- c("Down-regulated") }
      }
      download_FC_table
    })
    output$fc_download <- downloadHandler(
      filename = function(){paste("FC_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(download_foldchange(), fname) } )
    
    
    
    
    #############################################################################  
    #####
    #Log2FC tab
    
    # log2FC table based on 2^-DDct values
    log2foldchange <- reactive({
      DDct_table <- DDct_table()
      DDct_table <- na.omit(DDct_table)
      
      # Plotting summary log2DDct
      log2FC_table <- data.frame(matrix(NA, nrow=ncol(DDct_table), ncol=5))
      colnames(log2FC_table) <- c("Comparison","Type","log2FC","SE-lower","SE-upper")
      
      for(i in seq(1, ncol(DDct_table),1)){
        log2FC_table$Comparison[i] <- paste(colnames(DDct_table)[i])
        log2FC_table$log2FC[i] <- log(mean(na.omit(DDct_table[,i])),2)
        log2FC_table$`SE-lower`[i] <- log(mean(na.omit(DDct_table[,i]))-SE(na.omit(DDct_table[,i])),2)
        log2FC_table$`SE-upper`[i] <- log(mean(na.omit(DDct_table[,i]))+SE(na.omit(DDct_table[,i])),2)
        
        if(log2FC_table$log2FC[i] > 0){
          log2FC_table$Type[i] <- paste0(c("Up-regulated"),"   ",fa(name="arrow-up", fill="blue", prefer_type = "solid")) }
        else{
          log2FC_table$Type[i] <- paste0(c("Down-regulated"),"   ",fa(name="arrow-down", fill="red", prefer_type = "solid")) }
      }
      
      log2FCtable_plot <- log2FC_table
      
      for(i in seq(3,ncol(log2FC_table),1)){
        log2FC_table[,i] <- decimals(log2FC_table[,i],decVal())
      }
      
      # Prepare log2FC plot
      log2FCtable_plot$Comparison <- factor(log2FCtable_plot$Comparison, levels=log2FCtable_plot$Comparison)
      logFC_plot <- reactive({
        req(input$col_n)
        ggplot(log2FCtable_plot, aes(x=Comparison, y=log2FC), stat="identity") +
          geom_bar(stat="identity", aes(fill = log2FC < 0 )) +
          scale_fill_manual(guide=FALSE, breaks = c(TRUE,FALSE), values=c("tomato1", "dodgerblue2")) +
          geom_errorbar(aes(x=Comparison, ymin=`SE-lower`, ymax=`SE-upper`),
                        width=0.8, colour="black", alpha=0.9, linewidth=1) +
          scale_x_discrete(labels = gsub('\\s vs','\n vs', gsub('vs \\s', 'vs \n' , log2FCtable_plot$Comparison))) +
          theme_bw() +
          theme(
            legend.position = "none", 
            plot.title = element_text(size=20, face="bold"),
            axis.title.x = element_text(size=20, face="bold"),
            axis.title.y = element_text(size=20, face="bold"),
            axis.text.y = element_text(size=15, face="bold") ) +
          labs(y = expression(bold("log"["2"]~"FC"~'\n')),x = "\n\nComparison" ) 
      })
      
      # Adapt to user preference of plot orientation
      req(input$plot_orientation)
      if(input$plot_orientation == "horizontal"){
        logFC_plot_updated <- reactive({
          plot <- logFC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5, margin = margin(t=10), vjust=0.5))
          plot }) }
      if(input$plot_orientation == "vertical"){
        logFC_plot_updated <- reactive({
          plot <- logFC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 90, vjust = 0.5))
          plot }) }
      if(input$plot_orientation == "angle30"){
        logFC_plot_updated <- reactive({
          plot <- logFC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 30, vjust = 0.5))
          plot }) }
      if(input$plot_orientation == "angle45"){
        logFC_plot_updated <- reactive({
          plot <- logFC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 45, vjust = 0.5))
          plot }) }
      if(input$plot_orientation == "angle60"){
        logFC_plot_updated <- reactive({
          plot <- logFC_plot() + theme(axis.text.x=element_text(size=15, face="bold", colour="black", hjust = 0.5,angle = 60, vjust = 0.5))
          plot }) }
      
      # Output and download logFC plot
      output$plotlogFC <- renderPlot({ logFC_plot_updated() })
      output$downloadlogFC <- downloadHandler(
        filename = function(){paste("log2FC_plot_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".png", sep="")},
        content = function(file){  ggsave(file, plot=logFC_plot_updated(), dpi = 600, width = 16, height=10) } )
      
      log2FC_table
    })
    
    # Output log2FC table
    output$log2fc <- DT::renderDataTable(
      DT::datatable(log2foldchange(), filter = "top",
                    options = list(pageLength = -1,
                                   lengthMenu = list(c(10,50,100,-1),c('10','50','100','All')),
                                   columnDefs = list(list(className = 'dt-center', targets = 2:5)),dom="ft"), escape=FALSE))
    
    # Download log2FC table
    download_log2foldchange <- reactive({
      download_log2FC_table <- log2foldchange()
      for(i in seq(1,nrow(download_log2FC_table),1)){
        if(download_log2FC_table$log2FC[i] > 0){
          download_log2FC_table$Type[i] <- c("Up-regulated") }
        else{
          download_log2FC_table$Type[i] <- c("Down-regulated") }
      }
      download_log2FC_table
    })
    output$log2fc_download <- downloadHandler(
      filename = function(){paste("log2FC_table_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(download_log2foldchange(), fname) } )
    
    
    
  
    
    
    export_all_two_sample_tests <- reactive({
      
      initial_ttest <- ttests()
      initial_welch <- welchttests()
      initial_wilcoxon <- wilcoxontests()
      
      export <- data.frame(matrix(NA, nrow=nrow(initial_ttest)+1, ncol=7))
      colnames(export) <- c("Comparison", "Two Sample t-test", "", 
                            "Welch Two Sample t-test", "",
                            "Exact Wilcoxon rank sum test", "")
      export[1,2] <- c("t") ; export[1,3] <- c("p")
      export[1,4] <- c("t") ; export[1,5] <- c("p")
      export[1,6] <- c("W") ; export[1,7] <- c("p")
      
      export$Comparison[2:nrow(export)] <- initial_ttest$Comparison
      
      export[2:nrow(export), 2] <- initial_ttest$t
      export[2:nrow(export), 3] <- initial_ttest$p
      
      export[2:nrow(export), 4] <- initial_welch$t
      export[2:nrow(export), 5] <- initial_welch$p
      
      export[2:nrow(export), 6] <- initial_wilcoxon$W
      export[2:nrow(export), 7] <- initial_wilcoxon$p
      
      export
    })
    
    output$all_two_tests <- downloadHandler(
      filename = function(){paste("two_sample_tests_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(export_all_two_sample_tests(), fname) } )
    
    
    
    
    export_all_anova_tests <- reactive({
      
      initial_anova <- ANOVA()
      initial_welch_anova <- WELCHANOVA()
      initial_kruskal <- KRUSKALWALLIS()
      
      export <- data.frame(matrix(NA, nrow=3, ncol=5))
      colnames(export) <- c("Test", "Statistic", "p", "DFg", "DFv")
      
      export[1,1] <- initial_anova[1,1] ; export[1,2:5] <- initial_anova[1,3:6]
      export[2,1] <- initial_welch_anova[1,1] ; export[2,2:5] <- initial_welch_anova[1,3:6]
      export[3,1] <- initial_kruskal[1,1] ; export[3,2:4] <- initial_kruskal[1,3:5] ; export[3,5] <- c("")
      
      export
    })
    
    output$all_anova_tests <- downloadHandler(
      filename = function(){paste("anova_tests_",gsub(":","-",format(Sys.time(),'%d-%m-%Y_%H-%M-%S')), ".csv", sep="")},
      content = function(fname){ write.csv(export_all_anova_tests(), fname) } )
    
    
    
    
    export_all_posthoc_tests <- reactive({
      
      initial_tukey <- TUKEY()
      initial_gameshowell <- GAMESHOWELL()
      initial_dunn <- DUNN()
      
      export <- data.frame(matrix(NA, ))
      
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    }
}



# shiny app
shinyApp(ui=ui, server=server)

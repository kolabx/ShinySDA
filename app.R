library(shiny)
library(shinyWidgets)
library(rclipboard)
library(shinydashboard)
library(ggplot2)
library(data.table)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(grid)
library(gridExtra) 
library(dplyr)
library(Seurat)
library(shinyFiles)

source("fxs.R")


ui <- dashboardPage(skin="red",
                    dashboardHeader(title = "ShinySDA"),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Main Tab", tabName = "MainDash", icon = icon("dashboard"),
                                 badgeLabel = "underconst.", badgeColor = "yellow"),
                        menuItem("QC plots", tabName = "QCplots", icon = icon("affiliatetheme"),
                                 badgeLabel = "soon", badgeColor = "red"),
                        menuItem("Batch removal", tabName = "BatchRemove", icon = icon("allergies"),
                                 badgeLabel = "soon", badgeColor = "red"),
                        menuItem("Tab 4. (dev)", tabName = "FourthDash", icon = icon("arrows-alt"),
                                 badgeLabel = "soon", badgeColor = "red"),
                        menuItem("Enrichment Analysis", tabName = "Enrichment", icon = icon("dashboard"),
                                 badgeLabel = "underconst.", badgeColor = "yellow"),
                        menuItem("Source code", icon = icon("file-code-o"), 
                                 href = "https://")
                      )
                    ),
                    
                    dashboardBody(
                      tags$head(
                        tags$style(HTML("
                                        .content-wrapper {
                                        background-color: black !important;
                                        }
                                        .main-sidebar {
                                        background-color: black !important;
                                        }
                                        .multicol .shiny-options-group{
                            -webkit-column-count: 5; /* Chrome, Safari, Opera */
                                        -moz-column-count: 5;    /* Firefox */
                                        column-count: 5;
                                        -moz-column-fill: balanced;
                                        -column-fill: balanced;
                                        }
                                        .checkbox{
                                        margin-top: 0px !important;
                                        -webkit-margin-after: 0px !important; 
                                        }
                                        "))),
                      tabItems(
                        # First tab content
                        tabItem(tabName = "MainDash",
                                h2("Main Dashboard of ShinySDA"),
                                fluidRow(
                                  valueBoxOutput("InfoBox", width = 6),
                                  
                                  box(textInput("SDAroot", "Path to SDA folders. Expects dimnames in one dir up.", 
                                                value ="../../../Expts/TetCombo2/data/sda_results"), #/Tet_SDADGE_DropSim_mi10000_nc40_N21509_rep1
                                    uiOutput("select.folder"),
                                      actionButton("loadSDA", "Load SDA"),
                                      actionButton("getGeneAnn", "Get Gene Annotations"),
                                      actionButton("runtSNE", "Run tSNE (cs-all)"),
                                      actionButton("runtSNEQCfilt", "Run tSNE (cs-qc)"),
                                     # textInput("loadSDAmsg", "File Status", "not loaded"),
                                      width = 10
                                      #fileInput("SDAin", "Browse")
                                      ),
                                  box(
                                    title = "QC_MaxScore_filt", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqcMaxScorefilt"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE CS All", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("tSNE_CS_all"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE_CS_QC", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("tSNE_CS_qc"), 
                                    width = 5, background = "black"
                                  )
                                  
                                  
                                  
                                ) 
                        ),
                        
                        # QC plots 
                        tabItem(tabName = "QCplots",
                                h2("Full QC plots content"),
                                fluidRow(
                                  box(
                                    title = "QC1", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc1"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC2", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc2"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC3", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc3"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC4", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc4"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC5", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc5"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC6", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc6"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Convergence", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("convergence"), 
                                    width = 5, background = "black"
                                  ) ,
                                  box(
                                    title = "Gene Loading Hist.", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("loadhist"), 
                                    width = 5, background = "black"
                                  ) ,
                                  box(
                                    title = "Score Dist.", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("scoredist"), 
                                    width = 5, background = "black"
                                  ) ,
                                  #The PIP is the mean of the posterior. You can think of it as a measure of how likely it is that this variable is in the true model.
                                  box(
                                    title = "Post. inclus. prob. dist", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("pipdist"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Slack-Slab Prior", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("slackslabprior"), 
                                    width = 10, background = "black"
                                  )
                                  
                                  
                                  
                                  ,
                                  box(
                                    title = "main", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("main"), 
                                    width = 5, background = "black"
                                  )
                                  
                                )
                        ),
                        
                        # ThirdDash content
                        tabItem(tabName = "BatchRemove",
                                h2("Batch Removal"),
                                fluidRow(
                                  box(
                                    title = "SDA projected on tSNE", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    actionButton("prevSDA_br", "Prev comp"),
                                    actionButton("nextSDA_br", "Next comp"),
                                    plotOutput("SDAtsne_br1"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Batch Removal Selection", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, 
                                    column(10,
                                           uiOutput("CompBatchCheckBoxSelect")
                                    ),
                                    actionButton("save_batch_selection", "Save Selection"),
                                    actionButton("load_batch_selection", "Load last Selection"),
                                    actionButton("reset_batch_selection", "Reset last Selection"),
                                    width=5, background = "black"),
                                  box(
                                    title = "Scores order by Meta", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresAcross"), 
                                    width = 10, background = "black"
                                  ),
                                  
                                box(
                                  title = "Run tSNE Batch Removed", status = "primary", solidHeader = TRUE,
                                  collapsible = TRUE, 
        
                                  selectInput("tSNEiter", "tSNE n-iter:",
                                              c("Fast-500" = "500",
                                                "Med-2000" = "2000",
                                                "Med-5000" = "5000",
                                                "Robust-10000" = "10000"), selected = "500"),
                                  actionButton("run_tSNE_CS_batch", "Run tSNE (batch-removed)"),
                                  width = 10, background = "black",
                                ),
                                box(
                                  title = "tSNE Batch removed", status = "primary", solidHeader = TRUE,
                                  collapsible = TRUE,
                                  plotOutput("tSNE_CS_batch"), 
                                  width = 10, background = "black"
                                )
                                )
                        ),
                        
                        # FourthDash tab content
                        tabItem(tabName = "FourthDash",
                                h2("FourthDash tab content")
                        ),
                        
                        # Enrichment tab content
                        tabItem(tabName = "Enrichment",
                                fluidRow(
                                  
                                )
                        )
                        
                      )
                        ))


server <- function(input, output, session) {
  
  shinyFileChoose(input,'file', session=session,roots=c(wd='.'))
  
  envv=reactiveValues(y=NULL)
  envv$InfoBox_sub = "Load in SDA"
  

  output$select.folder <-
    renderUI(expr = selectInput(inputId = 'folder.name',
                                label = 'Folder Name',
                                choices = list.dirs(path = input$SDAroot,
                                                    full.names = FALSE,
                                                    recursive = FALSE)))
  
  
  # observeEvent(input$folder.name, {
  #   
  #   
  #   # print(head(input$folder.name))
  #   # print(head(input$SDAroot))
  #   
  #   envv$path2SDA_dyn <- paste0(input$SDAroot, "/", input$folder.name)
  #   
  #   # print(head( envv$path2SDA_dyn))
  #   
  # })


output$InfoBox <- renderValueBox({
    valueBox(
      value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
      subtitle = envv$InfoBox_sub,
      icon = icon("area-chart"),
      color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
    )
  })
  

## ObserveEvents--------------------------------------
observeEvent(input$loadSDA, {
    
    # print(head(paste0(input$SDAroot, "/", input$folder.name)))
  
  envv$path2SDA_dyn <- paste0(input$SDAroot, "/", input$folder.name)
  
    if(file.exists(envv$path2SDA_dyn)) {

      # loadloc = "../../../Expts/TetCombo2/data/sda_results/Tet_SDADGE_DropSim_mi10000_nc40_N21509_rep1"
      # input <- list(); envv <- list()
      # envv$path2SDA_dyn = loadloc
      
      # updateTextInput(session, "loadSDAmsg", value = "Loading SDA")
      envv$InfoBox_sub = "Loading SDA"
      
      SDAres <- SDAtools::load_results(
        results_folder = envv$path2SDA_dyn,
        data_path =  stringr::str_split(envv$path2SDA_dyn, "sda_results")[[1]][1])


      #update the names
      colnames(SDAres$scores) <- paste("SDA", 1:ncol(SDAres$scores), sep="")
      rownames(SDAres$loadings[[1]]) <- paste("SDA", 1:ncol(SDAres$scores), sep="")

      #adds some stats
      SDAres      <- AddCompStats(SDAres)

      envv$SDAres <- SDAres
      
      # updateTextInput(session, "loadSDAmsg", value = "SDA Loaded")
      
      
      envv$MaxScore.thr <- quantile(SDAres$component_statistics$max_score, c(.95))
      envv$QC_components  <- SDAres$component_statistics[max_score<envv$MaxScore.thr][order(Component)]$Component
      envv$QC_compIter <- min(envv$QC_components)
      
      if(file.exists(paste0(stringr::str_split(envv$path2SDA_dyn, "sda_results")[[1]][1], "ComboSeuratObj_AgSpT_MetaDF.rds"))){
        MetaDF <- readRDS(paste0(stringr::str_split(envv$path2SDA_dyn, "sda_results")[[1]][1], "ComboSeuratObj_AgSpT_MetaDF.rds"))
        envv$MetaDF <- MetaDF
        envv$InfoBox_sub = "SDA & Meta Loaded"
      } else {
        envv$InfoBox_sub = "SDA Loaded - Meta not found"
      }
      

    } else { 
      # updateTextInput(session, "loadSDAmsg", value = "File not found")
    
      
      }
  })


## get biomart
observeEvent(input$getGeneAnn, {
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    library(biomaRt)
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    SDAres <- envv$SDAres
    
    
    
    
    if(!file.exists(paste0(base.path, head.path, "_biomaRt_gene_loc_human.rds"))){
      
      gene_locations <- get.location(gene.symbols=colnames(SDAres$loadings[[1]]),
                                     data_set = "hsapiens_gene_ensembl",
                                     gene_name = "external_gene_name")
      
      saveRDS(gene_locations, paste0(base.path, head.path, "_biomaRt_gene_loc_human.rds"))
      
    } else {
      gene_locations <- readRDS(paste0(base.path, head.path, "_biomaRt_gene_loc_human.rds"))
    }
    
    if(!file.exists(paste0(base.path, head.path, "_SDAtools_GeneLoc_human.rds"))){
      
      GeneLoc         <- SDAtools::load_gene_locations(path = base.path,
                                                       genes = colnames(SDAres$loadings[[1]]),
                                                       organism = "hsapiens_gene_ensembl",
                                                       name="human")
      chromosome.lengths <- SDAtools::load_chromosome_lengths(organism = "hsapiens_gene_ensembl")
      
      saveRDS(GeneLoc, paste0(base.path, head.path, "_SDAtools_GeneLoc_human.rds"))
      saveRDS(chromosome.lengths, paste0(base.path, head.path, "_SDAtools_chromosome.lengths.rds"))
      
    } else {
      GeneLoc                <- readRDS(paste0(base.path, head.path, "_SDAtools_GeneLoc_human.rds"))
      chromosome.lengths     <- readRDS(paste0(base.path, head.path, "_SDAtools_chromosome.lengths.rds"))
    }
    
    envv$chromosome.lengths <- chromosome.lengths
    envv$GeneLoc <- GeneLoc
    envv$gene_locations <- gene_locations
    
  }


  
  
})


## run raw tSNE
observeEvent(input$runtSNE, {
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    SDAres <- envv$SDAres
   
    
    if(!file.exists(paste0(base.path, head.path, "_tSNE_CellScore_AllComps.rds"))){
      
      envv$InfoBox_sub = "Starting tSNE with cell scores - all comps.. wait"
      tsne_CS_all <- Rtsne::Rtsne(SDAres$scores[,], verbose=TRUE, pca=FALSE, perplexity = 40, 
                                  max_iter=1000, num_threads = 8)
    
      
      envv$InfoBox_sub = "Saving tSNE with cell scores - all comps.. wait"
      saveRDS(tsne_CS_all, file=paste0(base.path, head.path, "_tSNE_CellScore_AllComps.rds"))
      
      
      
    } else {
      tsne_CS_all <- readRDS(paste0(base.path, head.path, "_tSNE_CellScore_AllComps.rds"))
    }
    envv$InfoBox_sub = "raw tSNE with cell scores complete"
    envv$tsne_CS_all <- tsne_CS_all
    
  }
  
 
 
})

##Qc tsne
observeEvent(input$runtSNEQCfilt, {
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    SDAres <- envv$SDAres
    
    suffix <- paste(setdiff(1:as.numeric(SDAres$command_arguments$num_comps), envv$QC_components), collapse = "_")
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    
    tSNE_n.iter = 1000
    
    if(!file.exists(paste0(base.path, head.path, "_tSNE_CellScore_QCfil_",tSNE_n.iter,"_", suffix, ".rds"))){
      
      
      
      
      
      
      envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
      tsne_CS_qc <- Rtsne::Rtsne(SDAres$scores[,envv$QC_components], verbose=TRUE, pca=FALSE, perplexity = 40, 
                                 max_iter=tSNE_n.iter, num_threads = 8)
      
      envv$InfoBox_sub = "Saving tSNE with cell scores - qc comps.. wait"
      saveRDS(tsne_CS_qc, file=paste0(base.path, head.path, "_tSNE_CellScore_QCfil_",tSNE_n.iter,"_", suffix, ".rds"))
      
      
      
    } else {
      tsne_CS_qc <- readRDS(paste0(base.path, head.path, "_tSNE_CellScore_QCfil_",tSNE_n.iter,"_", suffix, ".rds"))
    }
    envv$InfoBox_sub = "qc tSNE with cell scores complete"
    envv$tsne_CS_qc <- tsne_CS_qc
    
  }
  
  
  
  
})


## Main tab--------------------------------------

#### QC threshold plot

output$SDAqcMaxScorefilt <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    MaxScore.thr <- quantile(SDAres$component_statistics$max_score, c(.95))
    MaxScore.thr9 <- quantile(SDAres$component_statistics$max_score, c(.9))
    MaxScore.thr75 <- quantile(SDAres$component_statistics$max_score, c(.75))
    
    sum(SDAres$component_statistics$max_score > MaxScore.thr)
    
    ggplot(SDAres$component_statistics, aes(y=max_score)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                          outlier.size=2) + 
      geom_abline(slope=0, intercept = MaxScore.thr, colour="red") + 
      geom_abline(slope=0, intercept = MaxScore.thr9, colour="dodgerblue")+ 
      geom_abline(slope=0, intercept = MaxScore.thr75, colour="navy")+ 
      theme_bw() + ggtitle(paste0("Max_Score quantile-thresholds remove:\n 95th perc: ", sum(SDAres$component_statistics$max_score > MaxScore.thr),
                                  "\n 90th perc: ", sum(SDAres$component_statistics$max_score > MaxScore.thr9),
                                  "\n 75th perc: ", sum(SDAres$component_statistics$max_score > MaxScore.thr75)))
    
    
    }
  
  
})

output$tSNE_CS_all <- renderPlot({
  if(is.null(envv$tsne_CS_all)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_all$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_all", "tSNE2_all")
    tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
    tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
    
    
    ggplot(tsneDF, aes(tSNE1_all, tSNE2_all, color=(SumScore))) +
      geom_point(size=0.5) + theme_bw() + 
      scale_color_distiller(palette = "Spectral") +
      theme(legend.position = "bottom") +
      ggtitle("tSNE SDA qc Components\n Sum absolute-cell-scores normalized by its mean \n ")
    
  }
})

output$tSNE_CS_qc <- renderPlot({
  if(is.null(envv$tsne_CS_qc)){
    plot(x=0, y=0, main="Load an SDA")

  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_qc$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_qc", "tSNE2_qc")
    
    
    
    
    if(is.null(envv$MetaDF)){
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
      
      ggplot(tsneDF, aes(tSNE1_qc, tSNE2_qc, color=(SumScore))) +
        geom_point(size=0.5) + theme_bw() + 
        scale_color_distiller(palette = "Spectral") +
        theme(legend.position = "bottom") +
        ggtitle("tSNE SDA qc Components\n  Sum absolute-cell-scores normalized by its mean \n ")
      
    } else {
      MetaDF <- envv$MetaDF
      tsneDF$library_size <- MetaDF[rownames(tsneDF), ]$library_size
      
      ggplot(tsneDF, aes(tSNE1_qc, tSNE2_qc, color=log10(library_size))) +
        geom_point(size=0.5) + theme_bw() + 
        scale_color_distiller(palette = "Spectral") +
        theme(legend.position = "bottom") +
        ggtitle("tSNE SDA qc Components\n log10 library size \n ")
    }
    

    
  }
})
  
## QC tabs--------------------------------------

### QC plots 

output$SDAqc1 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
  SDAres <- envv$SDAres
  ggplot(SDAres$component_statistics, aes(max_score, max_loading,
                                          label = Component_name_plot)) +
    geom_point() + geom_label_repel() + theme_bw()+ ggtitle("")}
  
  
})

output$SDAqc2 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
  SDAres <- envv$SDAres
  ggplot(SDAres$component_statistics, aes(max_score, mean_score,
                                          label = Component_name_plot)) +
    geom_point() + geom_label_repel() + theme_bw()+ ggtitle("")}


})

output$SDAqc3 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
  SDAres <- envv$SDAres
  ggplot(SDAres$component_statistics, aes(mean_score, mean_loading,
                                          label = Component_name_plot)) +
    geom_point() + geom_label_repel() + theme_bw() + ggtitle("")}


})

output$SDAqc4 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
  SDAres <- envv$SDAres
  ggplot(SDAres$component_statistics, aes(sd_score, sd_loading,
                                          label = Component_name_plot)) +
    geom_point() + geom_label_repel() + theme_bw()}


})

output$SDAqc5 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
  SDAres <- envv$SDAres
  ggplot(SDAres$component_statistics, aes(ssqrd_score, ssqrd_loading,
                                          label = Component_name_plot)) +
    geom_point() + geom_label_repel() + theme_bw() + ggtitle("")
  }


})

output$SDAqc6 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    ggplot(SDAres$component_statistics, aes(max_loading, mean_loading,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw() + ggtitle("")
  }
  
  
  
})

output$convergence <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    hist(as.numeric(SDAres$loadings[[1]]), 
         main = "Overall distribution of gene loadings\nrep1", 
         xlab = "Gene Loading", breaks = 300, xlim = range(-.15,.15), col="skyblue")
  }
  
  
  
})

output$loadhist <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    SDAtools::check_convergence(SDAres)
  }
  
  
  
})

output$scoredist <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    qplot(asinh(as.numeric(SDAres$scores)), 
          binwidth = 0.01, main = "Overall distribution of individual scores", 
          xlab = "asinh(Score)")+ scale_y_log10()
  }
  
  
  
})

output$pipdist <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    qplot(as.numeric(SDAres$pips[[1]]), geom = "histogram", 
          binwidth = 0.005) + xlab("PIP") + ylab("Count") + scale_y_log10()
  }
  
  
  
})

output$slackslabprior <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    SDAres <- envv$SDAres
    
    
    library(ggforce)
    
    density_data <- rbind(data.table(density=as.data.table(density(SDAres$loadings[[1]][SDAres$pips[[1]]>0.5], bw=1e-4, n=5000)[c("x","y")]), type="Slab (PIP>0.5)"),
                          data.table(density=as.data.table(density(SDAres$loadings[[1]][SDAres$pips[[1]]<0.5],bw=1e-4, n=5000)[c("x","y")]), type="Spike (PIP<0.5)"))
    
    setnames(density_data, c("gene_loading","density","type"))
    
    
    sparsity_plot <- ggplot(density_data, aes(gene_loading, density, colour=type)) +
      geom_line() +
      facet_zoom(xy = density<40 & abs(gene_loading)<0.1) +
      scale_color_brewer(palette = "Set1") + 
      theme_bw() +
      theme(legend.title=element_blank()) +
      labs(x="Gene Loading",y="Density") +
      scale_x_continuous(labels = function(x) as.character(x)); 
    sparsity_plot

  }


})

## Batch removal tab

observeEvent(input$nextSDA_br, {
  SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
  if(which(SDAorder==envv$QC_compIter) < length(SDAorder)){
    envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) + 1]
  }
  
   
  
 
})

observeEvent(input$prevSDA_br, {
  
  SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
  if(which(SDAorder==envv$QC_compIter) < length(SDAorder)){
    envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) - 1]
  }
  
  
})


output$SDAtsne_br1 <- renderPlot({
  
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    zN = envv$QC_compIter
    SDAres <- envv$SDAres
    tempDFX <- as.data.frame(envv$tsne_CS_qc$Y)
    colnames(tempDFX) <- c("tSNE1_qc", "tSNE2_qc")
    
    tempDFX$SDAComp <- SDAres$scores[,paste0("SDA", zN, sep="")]
    
    ggplot(tempDFX, aes(tSNE1_qc, tSNE2_qc,  color=cut(asinh(SDAComp), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
      geom_point(size=0.1) + theme_bw() +
      scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
      theme(legend.position = "bottom", aspect.ratio=1) + 
      simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) + ggtitle(paste0("SDA", zN, sep=""))
    
  }
})

observeEvent(input$CompBatchCheckBoxSelect, {
  
  envv$Remove_comps <- input$CompBatchCheckBoxSelect
  
  
})

observeEvent(input$run_tSNE_CS_batch, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    SDAres <- envv$SDAres
    suffix <- paste(envv$Remove_comps, collapse = "_")
    tSNE_n.iter <- as.numeric(input$tSNEiter)
    
    choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    selected <- setdiff(choice, envv$Remove_comps)
    
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    if(!file.exists(paste0(base.path, head.path, "_tSNE_CellScore_QCfil_",tSNE_n.iter,"_", suffix, ".rds"))){
      
      envv$InfoBox_sub = "Starting tSNE with cell scores - batch-removal.. wait"
      tsne_CS_batch <- Rtsne::Rtsne(SDAres$scores[,selected], verbose=TRUE, pca=FALSE, perplexity = 40, 
                                 max_iter=tSNE_n.iter, num_threads = 8)
      
      envv$InfoBox_sub = "Saving tSNE with cell scores - batch-removal.. wait"
      saveRDS(tsne_CS_batch, file=paste0(base.path, head.path, "_tSNE_CellScore_QCfil_",tSNE_n.iter,"_", suffix, ".rds"))
      
      
      
    } else {
      tsne_CS_batch <- readRDS(paste0(base.path, head.path, "_tSNE_CellScore_QCfil_",tSNE_n.iter,"_", suffix, ".rds"))
    }
    envv$InfoBox_sub = "batch-removal tSNE with cell scores complete"
    envv$tsne_CS_batch <- tsne_CS_batch
    
  }
  
  # 
  
  
    

  
  
  
  
})


observeEvent(input$save_batch_selection, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    SDAres <- envv$SDAres
    
  
    # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    selected <- envv$Remove_comps
    # print(head(selected))
    saveRDS(selected, file=paste0(base.path, head.path, "_BatchSelectedComps", ".rds"))
    
    
  }
  })

observeEvent(input$load_batch_selection, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    if(file.exists(paste0(base.path, head.path, "_BatchSelectedComps", ".rds"))){
      selected <- readRDS(paste0(base.path, head.path, "_BatchSelectedComps", ".rds"))
    } else {
      choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) 
      selected <- setdiff(choice, envv$QC_components)
     
    }
    
    
    # print(head(selected))
    envv$Remove_comps <- selected
  }
})

observeEvent(input$reset_batch_selection, {
  choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) 
  selected <- setdiff(choice, envv$QC_components)
  
    envv$Remove_comps <- selected

})



output$CompBatchCheckBoxSelect <- renderUI({
  
  choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
  if(is.null(envv$Remove_comps)){
    selected <- setdiff(choice, envv$QC_components) #setdiff(choice, paste0("SDA",envv$QC_components))
    
  } else {
    selected <- envv$Remove_comps
    
  }
  
  tags$div(align = 'left',
           class = 'multicol',
           checkboxGroupInput("CompBatchCheckBoxSelect", "components",
                              choices=choice, selected = selected))

})


output$SDAScoresAcross <- renderPlot({
  
  # ColFac_DONR.ID <- CDID()
  
  if(is.null(envv$MetaDF)){
    print("No Comp")
  } else {
    SDAScores <- envv$SDAres$scores
    ComponentN <- as.numeric(envv$QC_compIter)
    MetaDF <- envv$MetaDF
    MetaDF <- MetaDF[rownames(SDAScores),]

    
    # colnames(MetaDF)
    
    
    
    
    
    
    ggplot(data.table(cell_index = 1:nrow(SDAScores), 
                             score = asinh((SDAScores[, paste0("SDA", ComponentN)])^3), 
                             experiment = MetaDF$EXP.ID, 
                             ColFac = MetaDF$BarcodePrefix), 
                  aes(cell_index, score, colour = ColFac)) + 
      geom_point(size = 0.5, stroke = 0) + 
      xlab("Cell Index") + ylab("asinh(Score^3)") + 
      #scale_color_brewer(palette = "Paired") + 
      theme_bw() + 
      theme(legend.position = "bottom") + 
      guides(colour = guide_legend(ncol = 4, override.aes = list(size=2, alpha=1))) +
      # scale_colour_manual(values =(col_vector),
      #                     guide = guide_legend(nrow=2)) +
      # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
      ggtitle(paste0("SDAV", ComponentN))
    
 
    
  }
  
})





output$tSNE_CS_batch <- renderPlot({
  
  if(is.null(envv$tsne_CS_batch)){
    plot(x=0, y=0, main="tsne CS Batch not found")

  } else {

    tsneDF <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_batch", "tSNE2_batch")




    if(is.null(envv$MetaDF)){
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)

      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=(SumScore))) +
        geom_point(size=0.5) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n ")+
        theme(legend.position = "bottom", aspect.ratio=1)

    } else {
      MetaDF <- envv$MetaDF
      tsneDF$library_size <- MetaDF[rownames(tsneDF), ]$library_size

      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=log10(library_size))) +
        geom_point(size=0.5) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n log10 library size \n ")+
        theme(legend.position = "bottom", aspect.ratio=1)
    }
  }
 
})


  
}

shinyApp(ui, server)


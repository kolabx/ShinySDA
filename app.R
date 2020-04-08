

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

library("BiocParallel")
register(MulticoreParam(4))


source("fxs.R")


ui <- dashboardPage(skin="red",
                    dashboardHeader(title = "ShinySDA"),
                    #https://rstudio.github.io/shinydashboard/appearance.html#icons
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Main Tab", tabName = "MainDash", icon = icon("dashboard")),
                        menuItem("QC plots", tabName = "QCplots", icon = icon("wrench")),
                        menuItem("Batch removal", tabName = "BatchRemove", icon = icon("toolbox")),
                        menuItem("DGE Batch-Removed", tabName = "DGEsda", icon = icon("autoprefixer")),
                        menuItem("Gene Explorer", tabName = "GeneExplorer", icon = icon("dna")),
                        menuItem("Save Out", tabName = "SaveOut", icon = icon("save"))
                        # menuItem("Source code", icon = icon("file-code-o"), 
                        #          href = "https://")
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
                                                value ="../../../Expts/TetCombo4/data/sda_results"), #/Tet_SDADGE_DropSim_mi10000_nc40_N21509_rep1
                                      uiOutput("select.folder"),
                                      actionButton("loadSDA", "Load SDA"),
                                      actionButton("getGeneAnn", "Get Gene Annotations"),
                                      actionButton("getSDAGo", "Get SDA GO Enrichments"),
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
                        
                        # Batch removal content
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
                                    title = "Run tSNE Batch Removed", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, 
                                    
                                    selectInput("tSNEiter", "tSNE n-iter:",
                                                c("Fast-500" = "500",
                                                  "Fast2-1000" = "1000",
                                                  "Med1-2000" = "2000",
                                                  "Med2-5000" = "5000",
                                                  "Robust-10000" = "10000",
                                                  "OverKill-20000"="20000"), selected = "500"),
                                    selectInput("tSNEpp", "tSNE prplx:",
                                                c("1" = "1",
                                                  "10" = "10",
                                                  "50-default" = "50",
                                                  "100-med" = "100",
                                                  "200-max" = "200",
                                                  "300-insane!"  = "300"), selected = "50"),
                                    actionButton("run_tSNE_CS_batch", "Run tSNE (batch-removed)"),
                                    actionButton("SDAScoresChi_clus", "Show/Hide Pairwise Clustering"),
                                    width = 5, background = "black",
                                  ),
                                  box(
                                    title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresChiPos", height = 300), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "ChiSqrRes Scores Neg cellscores", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresChiNeg", height = 300), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Scores order by Meta", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresAcross"), 
                                    width = 10, background = "black"
                                  ),
                                  box(
                                    title = "Scores order by magnitude", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresAcross2", height = 800), 
                                    width = 10, background = "black"
                                  ),
                                  
                                  
                                  box(
                                    title = "tSNE Batch removed", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("DGE_SDA_tSNE"), 
                                    width = 10, background = "black"
                                  )
                                )
                        ),
                        
                        # DGE-SDA-BatchRemove tab content
                        tabItem(tabName = "DGEsda",
                                h2("DGE_SDA Batched Removed DGE"),
                                fluidRow(
                                  box(
                                    title = "tSNE Meta Exploration", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    
                                    plotOutput("tSNE_CS_batch1"),
                                    selectInput("Metaselect1", "Meta select:",
                                                c("SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "ExpID" = "ExpID",
                                                  "EXP.ID" = "EXP.ID",
                                                  "SingleR_Labels" = "SingleR_Labels",
                                                  "SingleR_Labels_Fine" = "SingleR_Labels_Fine",
                                                  "Phase" = "Phase"), selected = "SampleDate"),
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE Meta Exploration", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    
                                    plotOutput("tSNE_CS_batch2"),
                                    selectInput("Metaselect2", "Meta select:",
                                                c("SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "ExpID" = "ExpID",
                                                  "EXP.ID" = "EXP.ID",
                                                  "SingleR_Labels" = "SingleR_Labels",
                                                  "SingleR_Labels_Fine" = "SingleR_Labels_Fine",
                                                  "Phase" = "Phase"), selected = "SampleDate"),
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE Meta Exploration", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    
                                    plotOutput("tSNE_CS_batch3"),
                                    selectInput("Metaselect3", "Meta select:",
                                                c("SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "ExpID" = "ExpID",
                                                  "EXP.ID" = "EXP.ID",
                                                  "SingleR_Labels" = "SingleR_Labels",
                                                  "SingleR_Labels_Fine" = "SingleR_Labels_Fine",
                                                  "Phase" = "Phase"), selected = "SampleDate"),
                                    width = 5, background = "black"
                                  ),
                                  box(title = "tSNE SDA projection", status = "primary", solidHeader = TRUE,
                                      collapsible = TRUE,
                                      actionButton("prevSDA_br2", "Prev comp"),
                                      actionButton("nextSDA_br2", "Next comp"),
                                      textInput("SDAVn", "SDA comp. No:"),
                                      plotOutput("SDAtsne_br2"),
                                      width=5, background = "black"
                                  ),
                                  box(title = "SDA score tabulation", status = "primary", solidHeader = TRUE,
                                      collapsible = TRUE,
                                      plotOutput("SDAtsne_br2Tab"),
                                      width=10, background = "black"
                                  )
                                ),
                                box(
                                  title = "Pos. Loadings GO", status = "primary", solidHeader = TRUE,
                                  collapsible = TRUE,
                                  plotOutput("GOpos"), #plotlyOutput
                                  width = 5, background = "black"
                                ),
                                box(
                                  title = "Neg. Loadings GO", status = "primary", solidHeader = TRUE,
                                  collapsible = TRUE,
                                  plotOutput("GOneg"),
                                  width = 5, background = "black"
                                ),
                                
                                box(title = "Pos. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                                    tableOutput("packageTablePos")
                                ),
                                box(title = "Neg. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                                    tableOutput("packageTableNeg")
                                )
                        ),
                        
                        # Batch-Removed Exploration
                        tabItem(tabName = "GeneExplorer",
                                h2("DGE_SDA Batched Removed DGE Explorer"),
                                fluidRow(
                                  box(
                                    title = "Inputs", status = "warning", solidHeader = TRUE,
                                    "Multiple formatting of gene sets accepted", 
                                    br(), "List can be seperated by comma e.g. from ", 
                                    br(), "   or spaces e.g. from Excel", 
                                    br(), "Also, single or double quotes or not",
                                    #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                                    textInput("GeneSet", "A set of genes", "'CD19', 'CD20', 'MS4A1', 'IGHM', 'IGHA2'"),
                                    width = 5
                                  ),
                                  
                                  box(
                                    title = "BatchRemoved-DGE Expr", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExprSDAtSNE"),
                                    width = 5
                                  ),
                                  
                                  
                                  box(
                                    title = "Positive Loadings", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GenesEnrichSDAPos"),
                                    width = 10
                                  ),
                                  
                                  box(
                                    title = "Negative Loadings", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GenesEnrichSDANeg"),
                                    width = 10
                                  )
                                  
                                )
                        ),
                        
                        # Batch-Removed Exploration
                        tabItem(tabName = "SaveOut",
                                h2("Save the results for downstream analysis"),
                                fluidRow(
                                  box(
                                    title = "Save as Seurat Object", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    actionButton("SaveAsSerObj", "Save as Seurat Obj"),
                                    width = 5, background = "black"
                                  )
                                  
                                )
                        )
                        
                      )
                        ))


server <- function(input, output, session) {
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  
  shinyFileChoose(input,'file', session=session,roots=c(wd='.'))
  
  envv=reactiveValues(y=NULL)
  envv$InfoBox_sub = "Load in SDA"
  
  
  output$select.folder <-
    renderUI(expr = selectInput(inputId = 'folder.name',
                                label = 'Folder Name',
                                choices = list.dirs(path = input$SDAroot,
                                                    full.names = FALSE,
                                                    recursive = FALSE)))
  
  
  
  
  ## Main tab--------------------------------------
  
  
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
      # loadloc = "../../../../Conrad/R/Utah/sda_results/Testis_Run2/TestisV2_150Kgenes_SDA_mi10000_nc50_N26093_rep1"
      # input <- list(); envv <- list()
      # envv$path2SDA_dyn = loadloc
      # 
      # updateTextInput(session, "loadSDAmsg", value = "Loading SDA")
      envv$InfoBox_sub = "Loading SDA"
      
      
      if(sum(grepl("_dimnames", list.files(dirname(envv$path2SDA_dyn), recursive = T)))==0){
        DimNamesPath <- paste0(dirname(dirname(envv$path2SDA_dyn)), "/")
      } else {
        DimNamesPath <- paste0(dirname(envv$path2SDA_dyn), "/")
        
      }
      # print(DimNamesPath)
      
      
      
      
      
      SDAres <- SDAtools::load_results(
        results_folder = envv$path2SDA_dyn,
        data_path =  DimNamesPath)
      
      
      #update the names
      colnames(SDAres$scores) <- paste("SDA", 1:ncol(SDAres$scores), sep="")
      rownames(SDAres$loadings[[1]]) <- paste("SDA", 1:ncol(SDAres$scores), sep="")
      
      #adds some stats
      SDAres      <- AddCompStats(SDAres)
      
      envv$SDAres <- SDAres
      
      # updateTextInput(session, "loadSDAmsg", value = "SDA Loaded")
      
      QuantThr = 0.95
      
      envv$MaxScore.thr <- quantile(SDAres$component_statistics$max_score, c(QuantThr))
      envv$QC_components  <- SDAres$component_statistics[max_score<envv$MaxScore.thr][order(Component)]$Component
      envv$QC_compIter <- min(envv$QC_components)
      
      NComps <- SDAres$command_arguments$num_comps
      
      TopN = 150
      envv$TopN <- TopN
      
      SDA_TopNpos <- (as.data.frame(lapply(1:NComps, function(xN){
        as.data.frame(print_gene_list(results = SDAres, i=xN, PosOnly = T, NegOnly = F, TopN = TopN))[1:TopN,1]
      })))
      colnames(SDA_TopNpos) <- paste0("SDAV", 1:NComps)
      
      SDA_TopNneg <- (as.data.frame(lapply(1:NComps, function(xN){
        as.data.frame(print_gene_list(results = SDAres, i=xN, PosOnly = F, NegOnly = T, TopN = TopN))[1:TopN,1]
      })))
      colnames(SDA_TopNneg) <- paste0("SDAV", 1:NComps)
      
      envv$SDA_TopNpos <- SDA_TopNpos
      envv$SDA_TopNneg <- SDA_TopNneg
      
      MetaPath.base <- stringr::str_split(envv$path2SDA_dyn, "sda_results")[[1]][1]
      # MetaPath.base = "../../../Expts/TetCombo3/data/"
      MetaPath <- list.files(MetaPath.base, full.names = T)[grepl("_MetaDF", list.files(MetaPath.base))][1]
      
      
      if(file.exists(MetaPath)){
        # Maybe this can be a dropdown too..
        #for now, 
        
        MetaDF <- readRDS(MetaPath)
        
        if(!is.null(MetaDF$SubjectId)) {
          MetaDF$SubjectId <- paste0("Rh", MetaDF$SubjectId)
        }
        
        envv$MetaDF <- MetaDF
        
        sum(rownames(envv$SDAres$scores) %in% rownames(envv$MetaDF))
        
        envv$SDAres$scores <- envv$SDAres$scores[!grepl("remove", rownames(envv$SDAres$scores)),]
        
        
        envv$MetaDF <- envv$MetaDF[rownames(envv$SDAres$scores), ]
        
        envv$SDAres$n$individuals <- nrow(envv$SDAres$scores)
        envv$SDAres$command_arguments$N <- as.character(nrow(envv$SDAres$scores))
        
        
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
      
      head.path <- gsub("/", "", head.path)
      
      print(head.path)
      print(envv$path2SDA_dyn)
      
      SDAres <- envv$SDAres
      
      # stringr::str_split("../../../../Conrad/R/Utah/sda_results/Testis_Run2", "sda_results/")
      
      # basename("../../../../Conrad/R/Utah/sda_results/Testis_Run2")
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_human.rds"))){
        
        gene_locations <- get.location(gene.symbols=colnames(SDAres$loadings[[1]]),
                                       data_set = "hsapiens_gene_ensembl",
                                       gene_name = "external_gene_name")
        
        saveRDS(gene_locations, paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_human.rds"))
        
      } else {
        gene_locations <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_biomaRt_gene_loc_human.rds"))
      }
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_human.rds"))){
        
        GeneLoc         <- SDAtools::load_gene_locations(path = base.path,
                                                         genes = colnames(SDAres$loadings[[1]]),
                                                         organism = "hsapiens_gene_ensembl",
                                                         name="human")
        chromosome.lengths <- SDAtools::load_chromosome_lengths(organism = "hsapiens_gene_ensembl")
        
        saveRDS(GeneLoc, paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_human.rds"))
        saveRDS(chromosome.lengths, paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_chromosome.lengths.rds"))
        
      } else {
        GeneLoc                <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_human.rds"))
        chromosome.lengths     <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_chromosome.lengths.rds"))
      }
      
      envv$chromosome.lengths <- chromosome.lengths
      envv$GeneLoc <- GeneLoc
      envv$gene_locations <- gene_locations
      
    }
    
    
    
    
  })
  
  ## get GO
  observeEvent(input$getSDAGo, {
    
    if(is.null(envv$SDAres)){ 
      envv$InfoBox_sub = "Load SDA"
    } else {
      
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      NComps <- as.numeric(envv$SDAres$command_arguments$num_comps)
      
      head.path <- gsub("/", "", head.path)
      
      
      SDAres <- envv$SDAres
      
      envv$InfoBox_sub <- paste0(NComps, " comps, ~30 sec per comp")
      
      library(AnnotationHub) # source("https://bioconductor.org/biocLite.R"); biocLite("AnnotationHub")
      library(clusterProfiler) # source("https://bioconductor.org/biocLite.R"); biocLite("clusterProfiler")
      
      hub <- AnnotationHub()
      
      
      query(hub, "org.Mmu.eg.db")
      HuMAN.names <- query(hub, "org.Mmu.eg.db")#org.Hs.eg.db  org.MM.eg
      
      
      HuMAN <- hub[["AH75746"]] #AH75742 #query(hub, "org.MM.eg"), AH52234, AH57184.... new human AH73986? or AH70572? old? AH66156 query(hub, "org.Hs.eg.db")
      
      
      
      envv$GOAnn <- list(HuMAN = HuMAN, HuMAN.names = HuMAN.names)
      
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps.rds"))){
        
        GO_data <- list()
        
        for (i in 1:NComps){
          print(i)
          envv$InfoBox_sub <- paste0("Getting GO: %", round(i/NComps,3)*100)
          
          print("...negatives")
          GO_data[[paste0("V",i,"N")]] <- GO_enrichment(results =SDAres, i, side="N", geneNumber = 100, threshold=0.05, OrgDb =HuMAN)
          print("...positives")
          GO_data[[paste0("V",i,"P")]] <- GO_enrichment(results =SDAres, i, side="P", geneNumber = 100, threshold=0.05, OrgDb =HuMAN)
        }
        
        saveRDS(GO_data, paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps.rds"))
        
      } else{
        print("Loaded previously saved GO annotations.")
        
        GO_data <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps.rds"))
      }
      
      envv$InfoBox_sub <- "GO data loaded"
      envv$GO_data  <- GO_data
      
    }
    
    
    
    
  })
  
  
  ## run raw tSNE
  observeEvent(input$runtSNE, {
    
    if(is.null(envv$SDAres)){ 
      envv$InfoBox_sub = "Load SDA"
    } else {
      
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
      SDAres <- envv$SDAres
      
      # length(which(duplicated(SDAres$scores[,])))
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"))){
        
        # print(nrow(SDAres$scores[,]))
        tsnepp <- round((nrow(SDAres$scores[,]) - 1)/3.2)
        tsnepp <- ifelse(tsnepp < 50, tsnepp, 50)
        # print(tsnepp)
        
       
        
        envv$InfoBox_sub = "Starting tSNE with cell scores - all comps.. wait"
        tsne_CS_all <- Rtsne::Rtsne(SDAres$scores[,], verbose=TRUE, pca=FALSE, 
                                    perplexity = tsnepp, 
                                    max_iter=1000, 
                                    num_threads = 8, 
                                    check_duplicates = F)
        
        
        envv$InfoBox_sub = "Saving tSNE with cell scores - all comps.. wait"
        saveRDS(tsne_CS_all, file=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"))
        
        envv$InfoBox_sub = paste0("tSNE complete; pp=", tsnepp)
        if(length(which(duplicated(SDAres$scores[,])))>0) envv$InfoBox_sub = paste0("tSNE complete : dups=", length(which(duplicated(SDAres$scores[,]))))
        
        
      } else {
        tsne_CS_all <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"))
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
      
      head.path <- gsub("/", "", head.path)
      
      tsnepp <- round((nrow(SDAres$scores[,]) - 1)/3.2)
      tsnepp <- ifelse(tsnepp < 50, tsnepp, 50)
      
      tSNE_n.iter = 1000
      tsnepp
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"))){
        
        
        
        
        
        
        envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
        tsne_CS_qc <- Rtsne::Rtsne(SDAres$scores[,envv$QC_components], verbose=TRUE, pca=FALSE, 
                                   perplexity = tsnepp, 
                                   max_iter=tSNE_n.iter, num_threads = 8, check_duplicates = F)
        
        envv$InfoBox_sub = "Saving tSNE with cell scores - qc comps.. wait"
        saveRDS(tsne_CS_qc, file=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"))
        
        
        
      } else {
        tsne_CS_qc <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"))
      }
      envv$InfoBox_sub = "qc tSNE with cell scores complete"
      envv$tsne_CS_qc <- tsne_CS_qc
      
    }
    
    
    
    
  })
  
  
  
  
  
  
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
  
  ## Batch removal tab--------------------------------------
  
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
  
  observeEvent(input$SDAScoresChi_clus, {
    
   if(is.null(envv$SDAScoresChi_clusBTN)) {
     envv$SDAScoresChi_clusBTN = "ON"
   } else if( envv$SDAScoresChi_clusBTN == "OFF"){
     envv$SDAScoresChi_clusBTN = "ON"
   } else if(envv$SDAScoresChi_clusBTN == "ON"){
     envv$SDAScoresChi_clusBTN = "OFF"
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
      
      ggplot(tempDFX, aes(tSNE1_qc, tSNE2_qc,  color=cut(asinh(SDAComp^3), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
        geom_point(size=0.1) + theme_bw() +
        scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        theme(legend.position = "bottom", aspect.ratio=1) + 
        simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) + ggtitle(paste0("SDA", zN, sep=""))+ylab("asinh(SDAscore^3)")
      
    }
  })
  
  observeEvent(input$CompBatchCheckBoxSelect, {
    
    envv$Remove_comps <- input$CompBatchCheckBoxSelect
    
    
  })
  
  observeEvent(input$run_tSNE_CS_batch, {
    
    
    if(is.null(envv$SDAres)){ 
      envv$InfoBox_sub = "Load SDA"
    } else {
      # envv$Remove_comps = c(13, 40)
      SDAres <- envv$SDAres
      if(length(envv$Remove_comps) > 50) {
        # suffix <- paste(envv$Remove_comps, collapse = "")
        
        suffix <- paste0("LargeSetOfComps_", length(envv$Remove_comps))
        
      } else if(length(envv$Remove_comps) > 30 & length(envv$Remove_comps) <= 50 ) {
        
        suffix <- paste0(findIntRuns(as.numeric(unlist(strsplit(envv$Remove_comps, ",")))), collapse="")
        
      } else if (length(envv$Remove_comps) <= 30) {
        suffix <- paste(envv$Remove_comps, collapse = "_")
      }
      
      tSNE_n.iter <- as.numeric(input$tSNEiter) # tSNE_n.iter = 1000
      tSNE_pp <- as.numeric(input$tSNEpp) # tSNE_pp = 50
      
      choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
      selected <- setdiff(choice, envv$Remove_comps)
      
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tSNE_pp,"", suffix, ".rds"))){
        
        envv$InfoBox_sub = "Starting tSNE with cell scores - batch-removal.. wait"
        tsne_CS_batch <- Rtsne::Rtsne(SDAres$scores[,selected], verbose=TRUE, pca=FALSE, 
                                      perplexity = tSNE_pp, 
                                      max_iter=tSNE_n.iter, num_threads = 8, check_duplicates = F)
        
        envv$InfoBox_sub = "Saving tSNE with cell scores - batch-removal.. wait"
        saveRDS(tsne_CS_batch, file=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tSNE_pp,"", suffix, ".rds"))
        
        
        
      } else {
        tsne_CS_batch <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tSNE_pp,"", suffix, ".rds"))
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
      
      head.path <- gsub("/", "", head.path)
      
      
      SDAres <- envv$SDAres
      
      
      # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
      selected <- envv$Remove_comps
      # print(head(selected))
      saveRDS(selected, file=paste0(envv$path2SDA_dyn, "/", head.path,"_BatchSelectedComps", ".rds"))
      
      
    }
  })
  
  observeEvent(input$load_batch_selection, {
    
    
    if(is.null(envv$SDAres)){ 
      envv$InfoBox_sub = "Load SDA"
    } else {
      
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
      
      if(file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_BatchSelectedComps", ".rds"))){
        selected <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_BatchSelectedComps", ".rds"))
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
  
  
  
  output$SDAScoresChiPos <- renderPlot({
    
    # ColFac_DONR.ID <- CDID()
    
    if(is.null(envv$MetaDF)){
      print("No Comp")
    } else {
      SDAScores <- envv$SDAres$scores
      ComponentN <- as.numeric(envv$QC_compIter)
      MetaDF <- envv$MetaDF
      MetaDF <- MetaDF[rownames(SDAScores),]
      
      
      
      
      PosCompsDF <- as.data.frame(lapply(levels(factor(MetaDF$EXP.ID)), function(CondX){
        
        apply(SDAScores[rownames(MetaDF)[which(MetaDF$EXP.ID == CondX)], ], 2, 
              function(x){
                round(sum(x>0)/nrow(SDAScores)*100, 2)
              })
      }))
      
      colnames(PosCompsDF) <- levels(factor(MetaDF$EXP.ID))
      
      # print(head(PosCompsDF))
      # print(min(PosCompsDF))
      

      PosCompsDF <- PosCompsDF[gtools::mixedsort(rownames(PosCompsDF)),]
      PosCompsDF$SDA <- factor(rownames(PosCompsDF), levels=rownames(PosCompsDF))
      
      

      ChiT <- chisq.test(PosCompsDF[,1:(ncol(PosCompsDF)-1)])
      
      ChiTres <- ChiT$residuals
      ChiTres[which(is.na(ChiTres))] = 0
      
      ChiResSD <- round(apply(ChiTres, 1, sd),2)
      ChiResSD[which(is.na(ChiResSD))] <- 0
      ChiResSD[ChiResSD < 0.2] <- ""
      
      if(is.null(envv$SDAScoresChi_clusBTN)) {
        clustStat = F
      } else {
        clustStat <- ifelse(envv$SDAScoresChi_clusBTN=="ON", T, F)
      }
      
      
      pheatmap::pheatmap((t(ChiT$residuals)),
                         cluster_cols = clustStat, cluster_rows = clustStat,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10), 
                         labels_col = paste0(rownames(PosCompsDF), " sd_", ChiResSD)
                         )
      
      
 
      
      
      
      
      
      
    }
    
  })
  
  
  output$SDAScoresChiNeg <- renderPlot({
    
    # ColFac_DONR.ID <- CDID()
    
    if(is.null(envv$MetaDF)){
      print("No Comp")
    } else {
      SDAScores <- envv$SDAres$scores
      ComponentN <- as.numeric(envv$QC_compIter)
      MetaDF <- envv$MetaDF
      MetaDF <- MetaDF[rownames(SDAScores),]
      
      
      
      NegCompsDF <- as.data.frame(lapply(levels(factor(MetaDF$EXP.ID)), function(CondX){
        
        apply(SDAScores[rownames(MetaDF)[which(MetaDF$EXP.ID == CondX)], ], 2, 
              function(x){
                round(sum(x<0)/nrow(SDAScores)*100, 1)
              })
      }))
      
      colnames(NegCompsDF) <- levels(factor(MetaDF$EXP.ID))
      
      # print(head(NegCompsDF))
      
      
  
      
      NegCompsDF <- NegCompsDF[gtools::mixedsort(rownames(NegCompsDF)),]
      NegCompsDF$SDA <- factor(rownames(NegCompsDF), levels=rownames(NegCompsDF))
      
      ChiT <- chisq.test(NegCompsDF[,1:(ncol(NegCompsDF)-1)])
      
      ChiTres <- ChiT$residuals
      ChiTres[which(is.na(ChiTres))] = 0
      
      ChiResSD <- round(apply(ChiTres, 1, sd),2)
      ChiResSD[which(is.na(ChiResSD))] <- 0
      ChiResSD[ChiResSD < 0.2] <- ""
      # print(ChiResSD)
      
      # print(head(ChiTres))
            
      if(is.null(envv$SDAScoresChi_clusBTN)) {
        clustStat = F
      } else {
        clustStat <- ifelse(envv$SDAScoresChi_clusBTN=="ON", T, F)
      }
      
      pheatmap::pheatmap((t(ChiTres)),
                         cluster_cols = clustStat, cluster_rows = clustStat,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10), 
                         labels_col = paste0(rownames(NegCompsDF), " sd_", ChiResSD)
                         )
      
      
   
      
      
    }
    
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
        scale_colour_manual(values =(col_vector),
                            guide = guide_legend(nrow=2)) +
        # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
        ggtitle(paste0("SDAV", ComponentN))
      
      
      
    }
    
  })
  
  
  output$SDAScoresAcross2 <- renderPlot({
    
    # ColFac_DONR.ID <- CDID()
    
    if(is.null(envv$MetaDF)){
      print("No Comp")
    } else {
      SDAScores <- envv$SDAres$scores
      ComponentN <- as.numeric(envv$QC_compIter)
      MetaDF <- envv$MetaDF
      MetaDF <- MetaDF[rownames(SDAScores),]
      
      
      # colnames(MetaDF)
      
      
      
      
      tempDF <- data.frame(cell_index = 1:nrow(SDAScores), 
                           score = asinh((SDAScores[, paste0("SDA", ComponentN)])^3), 
                           experiment = MetaDF$EXP.ID, 
                           ColFac = MetaDF$BarcodePrefix)
      
      tempDF <- tempDF[order(tempDF$score),]
      tempDF$cell_index <- 1:nrow(tempDF)
      
      # tempDF$cell_index_cut <- cut(tempDF$cell_index, quantile(tempDF$cell_index))
      
      # print(levels(tempDF$cell_index_cut))
      
      ggplot(tempDF, 
             aes(cell_index, score, colour = ColFac)) + 
        geom_jitter(size=1, width=0, height = 3, alpha = .6) +
        # geom_boxplot(aes(x= factor(cut(cell_index, quantile(tempDF$cell_index))), y=score, 
                         # colour = ColFac), outlier.colour = "red", outlier.shape = 8) +
        # geom_point(size = 0.5, stroke = 0) + 
        xlab("Cell Index") + ylab("asinh(Score^3)") + 
        #scale_color_brewer(palette = "Paired") + 
        theme_bw() + 
        theme(legend.position = "none") + 
        guides(colour = guide_legend(ncol = 4, override.aes = list(size=2, alpha=1))) +
        scale_colour_manual(values =(col_vector),
                            guide = guide_legend(nrow=2)) +
        # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
        ggtitle(paste0("SDAV", ComponentN)) + 
        geom_smooth(method = "lm", formula = y ~ x, size = 2, colour="red")+ 
        geom_smooth(method = "loess", formula = y ~ x, size = 2, colour="dodgerblue")  +
        facet_wrap(~experiment, ncol = 2)
      
      
      
    }
    
  })
  
  
  output$DGE_SDA_tSNE <- renderPlot({
    
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
  
  ## Batch Removed DGE --------------------------------------
  
  
  output$tSNE_CS_batch1 <- renderPlot({
    
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
          ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
          theme(legend.position = "bottom", aspect.ratio=1)
        
      } else {
        MetaDF <- envv$MetaDF
        
        tsneDF$Meta <- MetaDF[rownames(tsneDF), input$Metaselect1]
        
        ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=Meta)) +
          geom_point(size=0.1, alpha=.4)+ theme_bw() +
          theme(legend.position = "bottom", aspect.ratio=1) +
          ggtitle(paste0("tSNE - batch removed cell scores\n", input$Metaselect1)) +
          scale_color_manual(values = col_vector
                             #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
          ) + 
          guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol=5))
        
        
      }
    }
    
  })
  
  output$tSNE_CS_batch2 <- renderPlot({
    
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
          ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
          theme(legend.position = "bottom", aspect.ratio=1)
        
      } else {
        MetaDF <- envv$MetaDF
        
        tsneDF$Meta <- MetaDF[rownames(tsneDF), input$Metaselect2]
        
        ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=Meta)) +
          geom_point(size=0.1, alpha=.4)+ theme_bw() +
          theme(legend.position = "bottom", aspect.ratio=1) +
          ggtitle(paste0("tSNE - batch removed cell scores\n", input$Metaselect2)) +
          scale_color_manual(values = col_vector
                             #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
          ) + 
          guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol=5))
        
        
      }
    }
    
  })
  
  output$tSNE_CS_batch3 <- renderPlot({
    
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
          ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
          theme(legend.position = "bottom", aspect.ratio=1)
        
      } else {
        MetaDF <- envv$MetaDF
        
        tsneDF$Meta <- MetaDF[rownames(tsneDF), input$Metaselect3]
        
        ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=Meta)) +
          geom_point(size=0.1, alpha=.4)+ theme_bw() +
          theme(legend.position = "bottom", aspect.ratio=1) +
          ggtitle(paste0("tSNE - batch removed cell scores\n", input$Metaselect3)) +
          scale_color_manual(values = col_vector
                             #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
          ) + 
          guides(colour = guide_legend(override.aes = list(size=2, alpha=1), ncol=5))
        
        
      }
    }
    
  })
  
  
  observeEvent(input$nextSDA_br2, {
    
    
    SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
    
    
    
    if(envv$QC_compIter + 1 %in% SDAorder){
      # envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) + 1]
      # updateTextInput(session, "SDAVn", value = envv$QC_compIter)
      updateTextInput(session, "SDAVn", value = SDAorder[which(SDAorder==envv$QC_compIter) + 1])
      
    }
    
    
    
    
  })
  
  observeEvent(input$prevSDA_br2, {
    
    # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    # SDAorder <- setdiff(choice, envv$Remove_comps)
    
    SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
    
    
    if(envv$QC_compIter - 1 %in% SDAorder){
      #envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) - 1]
      # updateTextInput(session, "SDAVn", value = envv$QC_compIter)
      updateTextInput(session, "SDAVn", value = SDAorder[which(SDAorder==envv$QC_compIter) - 1])
    }
    
    
    
    
  })
  
  
  observeEvent(input$SDAVn, {
    
    
    
    if(!is.null(envv$QC_compIter)){
      
      SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
      
      if(is.null(input$SDAVn) | (!(as.numeric(input$SDAVn) %in% SDAorder))) {
        
        
        print("Its null")
        updateTextInput(session, "SDAVn", value = envv$QC_compIter)
      } else {
        envv$QC_compIter = as.numeric(input$SDAVn)
      }
    }
    
    
    
    # # SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
    # if(!is.null(envv$QC_compIter)){
    #   
    #   if(!is.null(input$SDAVn)) {
    #     
    #     
    #   } else {
    #     # input$SDAVn <- envv$QC_compIter
    #   }
    # } else {
    #     
    #   }
    
    # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    # SDAorder <- setdiff(choice, envv$Remove_comps)
    # 
    # input$SDAVn <- as.character(envv$QC_compIter)
    # 
    # if(as.numeric(input$SDAVn) %in% SDAorder) {
    # 
    #   envv$QC_compIter = as.numeric(input$SDAVn)
    # 
    # } else {
    
    
    
    
    
  })
  
  
  output$SDAtsne_br2 <- renderPlot({
    
    if(is.null(envv$SDAres)){
      plot(x=0, y=0, main="Load an SDA")
    } else {
      
      zN = envv$QC_compIter
      
      SDAres <- envv$SDAres
      tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
      rownames(tempDFX)  <- rownames(envv$SDAres$scores)
      colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")
      
      if(zN %in% envv$Remove_comps) RemoveTag = "removed" else RemoveTag = "kept"
      
      
      tempDFX$SDAComp <- cut(asinh(SDAres$scores[,paste0("SDA", zN, sep="")]^3), 
                             breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf))
      
      tempDFX$SDAComp <- factor(tempDFX$SDAComp, 
                                levels = c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) )
      # print(tempDFX$SDAComp)
      # print(factor(tempDFX$SDAComp))
      # print(factor(tempDFX$SDAComp, 
      #              levels = c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) ))
      
      ggplot(tempDFX, aes(tSNE1_batch, tSNE2_batch,  color=tempDFX$SDAComp)) +
        geom_point(size=0.1) + theme_bw() +
        scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        theme(legend.position = "bottom", aspect.ratio=1) + 
        simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) +
        ggtitle(paste0("SDA", zN, " :: ", RemoveTag))+
        ylab("asinh(SDAscore^3)")
      
      
      
      
      
    }
  })
  
  output$SDAtsne_br2Tab <- renderPlot({
    
    if(is.null(envv$SDAres)){
      plot(x=0, y=0, main="Load an SDA")
    } else {
      
      zN = envv$QC_compIter
      
      SDAres <- envv$SDAres
      tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
      rownames(tempDFX)  <- rownames(envv$SDAres$scores)
      colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")
      
      if(zN %in% envv$Remove_comps) RemoveTag = "removed" else RemoveTag = "kept"
      
      
      tempDFX$SDAComp <- SDAres$scores[,paste0("SDA", zN, sep="")]
      
      
      
      if(!is.null(envv$MetaDF)){
        MetaDF <- envv$MetaDF
        tempDFX$Meta <- MetaDF[rownames(tempDFX), input$Metaselect3]
        tempDFX <- table(cut(asinh(tempDFX$SDAComp^3), 
                             c(-Inf, -1, -.5, 0, .5, 1, Inf)), tempDFX$Meta)
        # print(tempDFX)
        # print(rownames(tempDFX))
        
        
        tempDFX <-  tempDFX[rowSums(tempDFX)!=0, ]
        ppg2 <- ggplot(reshape2::melt(tempDFX)) +
          geom_bar(aes(x=Var2, y=value, fill=factor(Var1, levels=c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) )), 
                   stat="identity", width = 0.7, position="fill") +
          theme_bw()  + scale_fill_manual(values=rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue"))) +
          theme(legend.position="bottom",
                legend.direction="horizontal",
                legend.title = element_blank(),
                axis.text.x = element_text(angle = 90)) +
          scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                             labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
          ggtitle(paste0("Relative Contribution\n","SDA", zN, " :: ", RemoveTag)) + ylab("Relative % cells")
        
        print(ppg2)
        
        
      } else {
        plot(x=0, y=0, main="No Meta")
      }
      
      
      
    }
  })
  
  output$GOpos <- renderPlot({
    
    
    
    if(is.null(envv$GO_data)){
      plot(x=0, y=0, main="Load an SDA the GO data")
      
    } else {
      GO_data  <- envv$GO_data
      zN = envv$QC_compIter
      go_volcano_plot(x=GO_data, 
                      component = paste("V", zN, "P", sep=""))+ 
        theme_bw()+ theme(aspect.ratio = 1)
      
    }
    
  })
  
  output$GOneg <- renderPlot({
    
    if(is.null(envv$GO_data)){
      plot(x=0, y=0, main="Load an SDA the GO data")
      
    } else {
      GO_data  <- envv$GO_data
      zN = envv$QC_compIter
      go_volcano_plot(x=GO_data, component = paste("V", zN, "N", sep=""))+ 
        theme_bw()+ theme(aspect.ratio = 1)
      
    }
    
  })
  
  output$packageTablePos <- renderTable({
    
    print_gene_list(results=envv$SDAres, as.numeric(envv$QC_compIter), PosOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(50))
  }, digits = 1)
  
  output$packageTableNeg <- renderTable({
    print_gene_list(results=envv$SDAres, as.numeric(envv$QC_compIter), NegOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(50))
  }, digits = 1)
  
  
  ## Enrich N Explore-----------
  
  
  output$GenesEnrichSDAPos <- renderPlot({
    
    SDAres <- envv$SDAres
    SDA_TopNpos <- envv$SDA_TopNpos
    
    envv$TopN
    
    # N = total number of genes (usually not entire genome, since many have unk func)
    N=length(colnames(SDAres$loadings[[1]]))
    # k = number of genes submitted, top N
    k = envv$TopN
    
    GeneSet <- input$GeneSet
    
    
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }
      
      
    } else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      
    }
    
    GeneSetNot <- GeneSet[!GeneSet %in% colnames(SDAres$loadings[[1]][,])]
    
    print("length of your genes:")
    print(length(GeneSet))
    GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]][,])]
    print("length of your genes in this dataset:")
    print(length(GeneSet))
    
    
    
    
    # print("length of your genes in this dataset:")
    # print(length(GeneSet))
    
    plotEnrich(GeneSetsDF=SDA_TopNpos, 
               GeneVec = GeneSet, 
               plotTitle= paste0("Gene-set enrichment\n SDA top ", k, " pos loadings\nGene universe size: ", N, "\n Hypergeometric test: * adj.p < 0.01 \n Genes not found: ",
                                 paste0(GeneSetNot, collapse = ", ")),
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  
  output$GenesEnrichSDANeg <- renderPlot({
    
    SDAres <- envv$SDAres
    SDA_TopNneg <- envv$SDA_TopNneg
    
    # envv$TopN
    
    # N = total number of genes (usually not entire genome, since many have unk func)
    N=length(colnames(SDAres$loadings[[1]]))
    # k = number of genes submitted, top N
    k = envv$TopN
    
    GeneSet <- input$GeneSet
    
    
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }
      
      
    } else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      
    }
    
    
    GeneSetNot <- GeneSet[!GeneSet %in% colnames(SDAres$loadings[[1]][,])]
    
    print("length of your genes:")
    print(length(GeneSet))
    GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]][,])]
    print("length of your genes in this dataset:")
    print(length(GeneSet))
    
    
    
    
    # print("length of your genes in this dataset:")
    # print(length(GeneSet))
    
    plotEnrich(GeneSetsDF=SDA_TopNneg, 
               GeneVec = GeneSet, 
               plotTitle= paste0("Gene-set enrichment\n SDA top ", k, " neg loadings\nGene universe size: ", N, "\n Hypergeometric test: * adj.p < 0.01 \n Genes not found: ",
                                 paste0(GeneSetNot, collapse = ", ")),
               xLab = "SDA Comps",
               N=N,
               k=k)
    
    
  })
  
  
  
  
  output$GeneExprSDAtSNE <- renderPlot({
    
    GeneSet <- input$GeneSet
    
    
    if(length(grep(",", GeneSet)) == 0){
      
      if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
        GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
      } else {
        GeneSet <- unlist(strsplit(GeneSet, " "))
      }
      
      
    } else {
      GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
      
    }
    
    
    
    if(is.null(envv$SDAres)){
      plot(x=0, y=0, main="Load an SDA")
    } else {
      
      tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
      rownames(tempDFX)  <- rownames(envv$SDAres$scores)
      colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")
      
      tempDFX$GeneExpr <- rep(0, nrow(tempDFX))
      
      
      SDAres <- envv$SDAres
      
      GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]])]
      
      if(length(GeneSet)>1){
        
        GeneExpr <- SDAres$scores %*% SDAres$loadings[[1]][,as.character(GeneSet)]
        GeneExpr <- as.data.frame(rowSums(GeneExpr))
        
        TitleX = paste0("Sum-Expr of :", paste(GeneSet, collapse = "_") )
        
        
      } else if(length(GeneSet)==1){
        GeneExpr <- SDAres$scores %*% SDAres$loadings[[1]][,as.character(GeneSet)]
        TitleX = paste0("Expr of :", GeneSet )
        
        
      } else if(!length(GeneSet)>=1)  {
        GeneExpr <- SDAres$scores %*% rep(0, nrow(SDAres$loadings[[1]]))
        TitleX = "No genes in input"
      }
      
      
      
      tempDFX[rownames(GeneExpr), ]$GeneExpr <- GeneExpr[,1]
      
      
      
      
      
      
      # tempDFX <- (envv$tSNEGEx_br)
      print(head(tempDFX))
      # TitleX <- envv$tSNEGEx_tit
      
      
      ggplot(tempDFX, aes(tSNE1_batch, tSNE2_batch,  color=cut(asinh(GeneExpr^3),
                                                               breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
        geom_point(size=.5) + theme_bw() +
        scale_color_manual("Expr", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) +
        guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
        theme(legend.position = "bottom", aspect.ratio=1) +
        simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) +
        ggtitle(paste0("SDA-Batch-removed DGE\n", TitleX))+
        ylab("asinh(GeneExpr^3)")
      
      
      
      
    }
    
    
    
    
  })
  
  ## SAve out ------------
  
  observeEvent(input$SaveAsSerObj, {
    
    
    if(is.null(envv$SDAres)){ 
      envv$InfoBox_sub = "Load SDA"
    } else {
      
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
      
      SDAres <- envv$SDAres
      
      
      
      # print(names(envv))
      library(Seurat)
      
      ## create an empty Seurat object
      Mat1 <- abs(Matrix::rsparsematrix(20, nrow(SDAres$scores), density = .8))
      
      
      colnames(Mat1) <- rownames(SDAres$scores)
      rownames(Mat1) <- paste0("Gene", 1:nrow(Mat1))
      
      
      
      
      Ser1 <- CreateSeuratObject(Mat1)
      Ser1 <- NormalizeData(Ser1) #needed for FindVariableFeatures
      Ser1 <- FindVariableFeatures(Ser1) #needed for ScaleData
      Ser1 <- ScaleData(Ser1, features = rownames(x = Ser1)) #needed for RunPCA
      Ser1 <- RunPCA(Ser1, npcs = 5, verbose = T) # creates @ reduction object
      
      
      
      reduction.data <- CreateDimReducObject(
        embeddings = (SDAres$scores),
        loadings = t(SDAres$loadings[[1]]),
        assay = "RNA",
        stdev = apply(SDAres$scores, 2, sd),
        key = "SDA_",
        misc = list(command.args = SDAres$command_arguments,
                    n.cells = SDAres$n,
                    pip.mat = SDAres$pips[[1]],
                    pip.frac = SDAres$pip_fraction,
                    iterations = seq_len(ncol(SDAres$free_energy)) * as.numeric(SDAres$command$free_freq),
                    free.energu = SDAres$free_energy[1, ])
      )
      
      # Ser0@reductions[["SDA"]] <- reduction.data
      
      Ser1@reductions[["SDA"]] <- reduction.data
      
      rownames(envv$tsne_CS_batch$Y) <- rownames((SDAres$scores))
      
      tsne.reduction_CS_BR <- CreateDimReducObject(
        embeddings = (envv$tsne_CS_batch$Y),
        key = "tSNECSBR_",
        assay = "RNA"
      )
      
      Ser1@reductions[["tSNECSBR"]] <- tsne.reduction_CS_BR
    }
    
    Ser1@misc$SDA_processing_results <- list(GOAnn = envv$GOAnn,
                                             GO_data = envv$GO_data,
                                             MaxScore.thr = envv$MaxScore.thr,
                                             path2SDA_dyn = envv$path2SDA_dyn,
                                             SDA_TopNneg = envv$SDA_TopNneg,
                                             TopN = envv$TopN,
                                             SDA_TopNpos = envv$SDA_TopNpos,
                                             Remove_comps = envv$Remove_comps,
                                             QC_components = envv$QC_components,
                                             tsne_CS_raw = envv$tsne_CS_all,
                                             MetaDF = envv$MetaDF,
                                             chromosome.lengths = envv$chromosome.lengths,
                                             gene_locations = envv$gene_locations)
    
    
    print("Saving...")
    saveRDS(Ser1, file=paste0(envv$path2SDA_dyn, "/", head.path,"_FinalSerObj", ".rds"))
    print("Save complete...")
    
    
    #   1] "MaxScore.thr"  "y"             "GOAnn"         "path2SDA_dyn"  "SDAres"        "SDA_TopNneg"   "TopN"         
    # [8] "tsne_CS_all"   "tsne_CS_qc"    "Remove_comps"  "tsne_CS_batch" "QC_components" "QC_compIter"   "SDA_TopNpos"  
    # [15] "InfoBox_sub"   "MetaDF"        "GO_data"  
  })
  
  
  
}

shinyApp(ui, server)


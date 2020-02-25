#this is a work-in progress area. These fx will eventually be generalized and put in OOSAP. 



simplify <- theme(legend.position = "none",
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank())
simplify2 <- theme(legend.position = "right",
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank())

AddCompStats <- function(SDAres, sdThrsh = 0.04, maxscoreThrsh=20, maxloadThrsh = 1, redoCalc=T){
  if (redoCalc)   SDAres$component_statistics <- NULL
  if (is.null(SDAres$component_statistics) ) {
    
    SDAres$component_statistics <- as.data.frame(data.table(
      Component = 1:SDAres$n$components, 
      Component_name = dimnames(SDAres$scores)[[2]], 
      max_score = apply(abs(SDAres$scores),  2, max),
      max_loading = apply(abs(SDAres$loadings[[1]]), 1, max),
      mean_score = apply(abs(SDAres$scores),  2, mean),
      mean_loading = apply(abs(SDAres$loadings[[1]]), 1, mean),
      sd_score = apply(abs(SDAres$scores),  2, sd),
      sd_loading = apply(abs(SDAres$loadings[[1]]), 1, sd),
      ssqrd_score = apply(SDAres$scores,  2, function(x) sum(x^2)),
      ssqrd_loading = apply(SDAres$loadings[[1]], 1, function(x) sum(x^2))
    )[order(-Component)])
    
  }
  
  SDAres$component_statistics$Component_namev2    <- SDAres$component_statistics$Component_name
  SDAres$component_statistics$Component_name_plot <- SDAres$component_statistics$Component_name
  
  SDAres$component_statistics[which(SDAres$component_statistics$sd_loading<sdThrsh),]$Component_namev2         <- rep("", length(which(SDAres$component_statistics$sd_loading<sdThrsh)))
  
  SDAres$component_statistics[which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh),]$Component_name_plot <- rep("", length(which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh)))
  
  SDAres$component_statistics <- data.table(SDAres$component_statistics)
  return(SDAres)
  
}


convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  humanx <- data.frame(genesV2)
  #print(head(humanx))
  return(humanx)
}

print_loadings_scores <- function(SDAResult = NA, ComponentN=NA, ColFac = NA, Prefix="SDAV", GeneLoc=NA){
  library(ggthemes)
  library(scales)
  if(!is.factor(ColFac)) ColFac <- factor(ColFac)
  
  SDAScores    <- SDAResult$scores
  SDALoadings <- SDAResult$loadings[[1]]
  
  
  #cowplot::plot_grid(
  g1 <- ggplot(data.table(cell_index = 1:nrow(SDAScores), 
                          score = SDAScores[, paste0(Prefix, ComponentN)], experiment = gsub("_.*", 
                                                                                             "", gsub("[A-Z]+\\.", "", rownames(SDAScores))), ColFac = ColFac), 
               aes(cell_index, score, colour = ColFac)) + 
    geom_point(size = 0.5, stroke = 0) + 
    xlab("Cell Index") + ylab("Score") + 
    #scale_color_brewer(palette = "Paired") + 
    
    
    theme_bw() + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
    scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(length(levels(ColFac))),
                        guide = guide_legend(nrow=2)) +
    #guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) + 
    ggtitle(paste0(Prefix, ComponentN)) 
  #,
  
  g2 <- genome_loadings(SDALoadings[ComponentN,], 
                        label_both = T, 
                        max.items = 10, 
                        gene_locations =   GeneLoc)
  #, ncol = 1)
  print(g1)
  print(g2)
}


get.location <- function(gene.symbols, data_set, gene_name){
  
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = data_set)
  
  mapTab <- getBM(attributes = c(gene_name,'chromosome_name','start_position'),
                  filters = gene_name, values = gene.symbols, mart = ensembl, uniqueRows=TRUE)
  
  mapTab <- as.data.table(mapTab)
  
  setnames(mapTab, gene_name,"gene_symbol")
  
  # Remove duplicate genes!!
  # first which genes are duplicated
  duplicates <- mapTab[duplicated(mapTab, by="gene_symbol")]$gene_symbol
  # order by chr so patch versions to go bottom, then choose first unique by name
  unduplicated <- unique(mapTab[gene_symbol %in% duplicates][order(chromosome_name)], by="gene_symbol")
  
  # remove duplicates and replace with unique version
  mapTab <- mapTab[!gene_symbol %in% duplicates]
  mapTab <- rbind(mapTab, unduplicated)
  
  # change all patch chr names to NA
  mapTab[!chromosome_name %in% c(1:22,"X","Y","MT")]$chromosome_name <- NA # mice actually 19
  mapTab[is.na(chromosome_name)]$start_position <- NA
  return(mapTab)
}


GO_enrichment <- function(results = NULL, component, geneNumber = 100, threshold=0.01, side="N", OrgDb = NULL){
  # results = SDAres; component = 1; geneNumber = 100; threshold=0.01; side="N"; OrgDb = HuMAN
  
  if(side=="N"){
    top_genes <- data.table(as.matrix(results$loadings[[1]][component, ]), keep.rownames = TRUE)[order(V1)][1:geneNumber]$rn
  }else{
    top_genes <- data.table(as.matrix(results$loadings[[1]][component, ]), keep.rownames = TRUE)[order(-V1)][1:geneNumber]$rn
  }
  
  gene_universe <- data.table(as.matrix(results$loadings[[1]][component,]), keep.rownames = TRUE)$rn
  
  ego <- enrichGO(gene = top_genes,
                  universe = gene_universe,
                  OrgDb = OrgDb,
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio)/frac_to_numeric(ego@result$BgRatio)
  
  ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  
  return(ego@result)
  
}

frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))


go_volcano_plot <- function(x=GO_data, component="V5N", extraTitle=""){
  if(extraTitle=="") extraTitle = paste("Component : ", component, sep="")
  
  #print(
  ggplot(data.table(x[[component]]), aes(GeneOdds/BgOdds, -log(pvalue), size=Count)) +
    geom_point(aes(colour=p.adjust<0.05)) +
    scale_size_area() +
    geom_label_repel(data = data.table(x[[component]])[order(p.adjust)][1:30][p.adjust<0.7], aes(label = Description, size=0.25), size = 3, force=2) + 
    ggtitle(paste("",extraTitle, sep="\n") ) +
    xlab("Odds Ratio") +
    scale_x_log10(limits=c(1,NA), breaks=c(1,2,3,4,5,6,7,8))
  #)
}

print_gene_list <- function(results = results, i, PrintRes = F, PosOnly = F, NegOnly = F, AbsLoad = T, TopN = 150) {
  # results = SDAres; i=1; PosOnly = T; NegOnly = F; TopN = TopN
  
  if(AbsLoad)  tmp <- data.table(as.matrix(results$loadings[[1]][i,]), keep.rownames = TRUE)[order(-abs(V1))][1:TopN]
  
  if(PosOnly)  tmp <- data.table(as.matrix(results$loadings[[1]][i,]), keep.rownames = TRUE)[order(-(V1))][1:TopN]
  
  if(NegOnly)  tmp <- data.table(as.matrix(results$loadings[[1]][i,]), keep.rownames = TRUE)[order((V1))][1:TopN]
  
  
  setnames(tmp, c("Gene.Name","Loading"))
  setkey(tmp, Gene.Name)
  
  
  # Display Result
  if(PrintRes) print(tmp[order(-abs(Loading))]) else return(tmp[order(-abs(Loading))])
}



plotEnrich <- function(GeneSetsDF, GeneVec, plotTitle="", xLab="", N = NULL, k = NULL, ReturnPval=T, BiPlot = F){
  
  
  GeneVec_overlap <- apply(GeneSetsDF, 2, function(x){
    length(which(x %in% GeneVec))
  })
  table(GeneVec_overlap)
  
  
  
  m = length(GeneVec)
  n = N - m
  
  ## Random expectation
  marked.proportion <- m / N; marked.proportion
  exp.x <- k * marked.proportion; exp.x
  
  x = GeneVec_overlap
  
  ## Fold enrichment, as computed by David
  fold.enrichment <-  (x / k ) / (m / N)
  
  # barplot(fold.enrichment, las=2)
  
  
  
  # ggplot(data=data.frame(x=names(fold.enrichment),
  #            y=fold.enrichment), aes(x=x, y=y)) + theme_bw()  +
  #   geom_bar(stat="identity") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Fold-enrichment \n Sertoli cell-only syndrome vs control DE genes")
  
  p.value <-  phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
  # p.value <-  p.adjust(p.value, "BH")
  p.value <-  p.adjust(p.value, "fdr")
  
  fold.enrichment <- fold.enrichment[gtools::mixedsort(names(fold.enrichment))]
  
  p.value <-  p.value[names(fold.enrichment)]
  
  tDF <- data.frame(x=(fold.enrichment),
                    y=p.value)
  rownames(tDF) <- names(fold.enrichment)
  # tDF[tDF$y <0.01,]
  
  # print(head(p.value))
  # print(head(names(fold.enrichment)))
  # print(names(p.value)[p.value<0.01])
  # print(paste(names(p.value)[p.value<0.01]))
  # print(rownames(tDF[tDF$y <0.01,]))
  
  tempSigComps <- rownames(tDF[tDF$y <0.01,])
  
  tempLen <- length(tempSigComps)
  
  tempKeepLen <- tempLen - length(grep("Removed", tempSigComps))
  
  if(tempLen> 16){
    plotTitle <- paste0(plotTitle, "\nSig Comps ", tempKeepLen, "/", tempLen, ":", 
                        paste(rownames(tDF[tDF$y <0.01,])[1:7], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[8:15], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[16:tempLen], collapse = ", "))
    
  } else  if(tempLen> 8){
    plotTitle <- paste0(plotTitle, "\nSig Comps: ", tempKeepLen, "/", tempLen, ":",
                        paste(rownames(tDF[tDF$y <0.01,])[1:7], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[8:tempLen], collapse = ", "))
    
  } else {
    plotTitle <- paste0(plotTitle, "\nSig Comps: ", tempKeepLen, "/", tempLen, ":",
                        paste(rownames(tDF[tDF$y <0.01,]), collapse = ", "))
    
  }
  
  if(BiPlot) print(ggplot(data=tDF, aes(x=x, y=y, label = ifelse(y < 0.01, "*", ""))) + theme_bw()  +
                     geom_point() + geom_line() +
                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                     geom_text(vjust = 0) +
                     ggtitle(plotTitle) +
                     xlab("Fold-enrichment") +
                     ylab("1 - P-value = 1 - P(X>=x) = P(X<x)"))
  
  tDF <- data.frame(x=factor(names(fold.enrichment), levels = names(fold.enrichment)),
                    y=fold.enrichment,
                    p=p.value)
  
  rownames(tDF) <- names(fold.enrichment)
  
  
  
  print(ggplot(data=tDF, aes(x=x, y=y, label = ifelse(p < 0.01, "*", ""))) + theme_bw()  +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          geom_text(vjust = 0) + 
          ggtitle(plotTitle) + 
          xlab(xLab) + 
          ylab("Fold-enrichment"))
  
  
  if(ReturnPval) return(1-p.value)
  
}



findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) %in% 1:2) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

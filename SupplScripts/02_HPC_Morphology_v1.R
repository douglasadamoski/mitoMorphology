#!/usr/bin/env /mnt/nfs/home/douglas/anaconda3/envs/R413/bin/Rscript --vanilla
#
#  02_HPC_Morphology_v1.R
#
#
#

#
args = commandArgs(trailingOnly=TRUE)

# Available folders
#ActualRData <- "A99.28_endogena_parkin_controle_minMitoLabel_300.RData"
#ActualRData <- "A99.29_ectopica_parkin_controle_minMitoLabel_300.RData"
#ActualRData <- "A99.28_ectopica_LC3_controle_minMitoLabel_300.RData"


args[1] -> ActualRData

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) > 1) {
  # default output file
  stop("Only one argument must be supplied (input file).n", call.=FALSE)
}

print(ActualRData)

# Get rootFolder
rootFolder <- getwd()

NumberPositiveMito <- 1
minMaxDProteinLabel <- 0


library(pbapply)
library(dplyr)
library(ggplot2)
library(ggmosaic)
library(ggfortify)
library(ggExtra)

# If loop to run everything at once
# utils::memory.limit(size=32000)

# Create main data
dir.create(path = paste0(rootFolder, "/", "04-Morphology"), showWarnings = FALSE)

###### FUNCTIONS BLOCK
### START
# Create function
ColumnRetrieval <- function(tableNow, columnNow){
  # print(fileName)
  return(tableNow[,columnNow])
}

### END
###### FUNCTIONS BLOCK




##########
### Create folder
ActualRData_name <- gsub(".RData", "", ActualRData)

# Add mito inducer label
if(grepl("LC3", ActualRData_name)){
  mitoInducerNow <- "CHQ"
}
if(grepl("parkin", ActualRData_name)){
  mitoInducerNow <- "CCCP"
}


# Create dir to receive data
dir.create(path = paste0(rootFolder, "/", "04-Morphology", "/", ActualRData_name, "/"),
           showWarnings = FALSE)


####### LABEL POSITIVE MITOCHONDRIA
### Create function
GateProteinaMitocondria_function <- function(tableNow){
  
  if(dim(tableNow)[1] > 0){
    
    # Limpa a marcação atual
    tableNow$PositiveCell <- FALSE
    tableNow$PositiveMitochondria <- FALSE
    
    # ADDED TO REMOVE ABNORMAL MITOCHONDRIA
    tableNow <- tableNow[which(!is.na(tableNow$D)),]
    tableNow <- tableNow[which(!is.na(tableNow$Perim.)),]
    
    # Marca as mitocôndrias positivas
    tableNow$PositiveMitochondria[which(
      tableNow$Protein_Mean > minProteinLabel &
        (tableNow$Protein_Max * tableNow$D) > minMaxDProteinLabel &
        tableNow$Protein_Max > minProteinMaxLabel &
        tableNow$Protein_Median > minProteinMedianLabel &
        tableNow$Protein_StdDev > minProteinStdDevLabel
    )] <- TRUE
    
    
    # Marca as CÉLULAS positivas
    # Tentativa com DPLYR para agilizar
    TempPositiveTable_plyr_toLabel <- tableNow %>%
      group_by(FileOrigin_base) %>%
      summarize(positiveMito = sum(PositiveMitochondria, na.rm=TRUE),
                totalMito = n()
      )
    tableNow[(tableNow$FileOrigin_base %in% TempPositiveTable_plyr_toLabel$FileOrigin_base[TempPositiveTable_plyr_toLabel$positiveMito > NumberPositiveMito]), "PositiveCell"] <- TRUE
    
    # Add values
    TempPositiveTable_plyr_toLabel$FractionMitoPositive <- 100*(TempPositiveTable_plyr_toLabel$positiveMito/TempPositiveTable_plyr_toLabel$totalMito)
    TempPositiveTable_plyr_toLabel$negativeMito <- (TempPositiveTable_plyr_toLabel$totalMito - TempPositiveTable_plyr_toLabel$positiveMito)
    tableNow[,c("totalMito", "positiveMito", "negativeMito", "FractionMitoPositive")] <- NA
    tableNow[,c("totalMito", "positiveMito", "negativeMito", "FractionMitoPositive")] <-  TempPositiveTable_plyr_toLabel[match(tableNow$FileOrigin_base, TempPositiveTable_plyr_toLabel$FileOrigin_base), c("totalMito", "positiveMito", "negativeMito", "FractionMitoPositive") ]
    return(tableNow)
  }
  
}

#####
# Create function
PerimeterThreshold_function <- function(tableNow){
  # print(fileName)
  # Ordena a tabela
  tableNow <- tableNow[order(tableNow$FileOrigin_base), ]
  tableNow$MitoCat_Perim <- NA
  tableNow$CellCat_Perim <- NA
  
  tableNow$CellCat_Perim_negative <- NA
  
  
  # Add intermediate
  tableNow[which(
    !is.na(tableNow$Perim.)
  ),"MitoCat_Perim"] <- "2-Intermediate"
  # Add fragmented
  tableNow[which(
    !is.na(tableNow$Perim.) &
      tableNow$Perim. < minPerimNow
  ),"MitoCat_Perim"] <- "1-Fragmented"
  # Add tubular
  tableNow[which(
    !is.na(tableNow$Perim.) &
      tableNow$Perim. > maxPerimNow
  ),"MitoCat_Perim"] <- "3-Tubular"
  #
  tableNow <- tableNow[which(!is.na(tableNow$MitoCat_Perim)),]
  
  ##### Add CellCat for:
  # Cells labelled otherwise:
  # Only FinalResults$PositiveMitochondria TRUE
  # 
  TempTable_plyr_toLabel <- tableNow[which(
    #tableNow$Info01 == "01-Mock" | 
    tableNow$PositiveMitochondria & tableNow$PositiveCell
  ), ] %>%
    group_by(FileOrigin_base) %>%
    count(MitoCat_Perim)
  TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
  
  # Loop to maximize
  SummarizedData <- data.frame(row.names=unique(TempTable_plyr_toLabel$FileOrigin_base))
  
  for(i in unique(TempTable_plyr_toLabel$FileOrigin_base)  ){
    tempNow <- TempTable_plyr_toLabel[TempTable_plyr_toLabel$FileOrigin_base == i, ]
    #
    SummarizedData[i,"Category"] <- tempNow$MitoCat_Perim[ which.max(tempNow$n) ]
    remove(tempNow)
  }
  
  # Add to column the cell information
  if(is.null(SummarizedData[tableNow$FileOrigin_base, "Category"])){
    
  } else {
    tableNow$CellCat_Perim <- SummarizedData[tableNow$FileOrigin_base, "Category"]
  }
  remove(SummarizedData)
  
  
  
  ##### Add CellCat for:
  # Only FinalResults$PositiveMitochondria FALSE
  # 
  TempTable_plyr_toLabel <- tableNow[which(
    tableNow$Info01 == "01-Mock" & !tableNow$PositiveMitochondria | 
      !tableNow$PositiveMitochondria & tableNow$PositiveCell
  ), ] %>%
    group_by(FileOrigin_base) %>%
    count(MitoCat_Perim)
  TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
  
  # Loop to maximize
  SummarizedData <- data.frame(row.names=unique(TempTable_plyr_toLabel$FileOrigin_base))
  
  for(i in unique(TempTable_plyr_toLabel$FileOrigin_base)  ){
    tempNow <- TempTable_plyr_toLabel[TempTable_plyr_toLabel$FileOrigin_base == i, ]
    #
    SummarizedData[i,"Category"] <- tempNow$MitoCat_Perim[which.max(tempNow$n)]
    remove(tempNow)
  }
  # Add to column the cell information
  if(is.null(SummarizedData[tableNow$FileOrigin_base, "Category"])){
    
  } else {
    tableNow$CellCat_Perim_negative <- SummarizedData[tableNow$FileOrigin_base, "Category"]
  }
  remove(SummarizedData)
  
  # 
  #
  tableNow <- tableNow[which(!(is.na(tableNow$CellCat_Perim) & is.na(tableNow$CellCat_Perim_negative))),]
  
  #tableNow[,c("CellCat_Perim",
  #            "CellCat_Perim_negative")]
  
  return(tableNow)
}


###### D CLASS DEFINITION

#####
#####
# Create function
DThreshold_function <- function(tableNow){
  # print(fileName)
  # Ordena a tabela
  tableNow <- tableNow[order(tableNow$FileOrigin_base), ]
  tableNow$MitoCat_D <- NA
  tableNow$CellCat_D <- NA
  
  tableNow$CellCat_D_negative <- NA
  
  
  # Add intermediate
  tableNow[which(
    !is.na(tableNow$D)
  ),"MitoCat_D"] <- "2-Intermediate"
  # Add fragmented
  tableNow[which(
    !is.na(tableNow$D) &
      tableNow$D < minDNow
  ),"MitoCat_D"] <- "1-Fragmented"
  # Add tubular
  tableNow[which(
    !is.na(tableNow$D) &
      tableNow$D > maxDNow
  ),"MitoCat_D"] <- "3-Tubular"
  #
  tableNow <- tableNow[which(!is.na(tableNow$MitoCat_D)),]
  
  ##### Add CellCat for:
  # Cells labelled otherwise:
  # Only FinalResults$PositiveMitochondria TRUE
  # 
  TempTable_plyr_toLabel <- tableNow[which(
    #tableNow$Info01 == "01-Mock" | 
    tableNow$PositiveMitochondria & tableNow$PositiveCell
  ), ] %>%
    group_by(FileOrigin_base) %>%
    count(MitoCat_D)
  TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
  
  # Loop to maximize
  SummarizedData <- data.frame(row.names=unique(TempTable_plyr_toLabel$FileOrigin_base))
  
  for(i in unique(TempTable_plyr_toLabel$FileOrigin_base)  ){
    tempNow <- TempTable_plyr_toLabel[TempTable_plyr_toLabel$FileOrigin_base == i, ]
    #
    SummarizedData[i,"Category"] <- tempNow$MitoCat_D[which.max(tempNow$n)]
    remove(tempNow)
  }
  # Add to column the cell information
  if(is.null(SummarizedData[tableNow$FileOrigin_base, "Category"])){
    
  } else {
    tableNow$CellCat_D <- SummarizedData[tableNow$FileOrigin_base, "Category"]
  }
  remove(SummarizedData)
  
  
  
  ##### Add CellCat for:
  # Only FinalResults$PositiveMitochondria FALSE
  # 
  TempTable_plyr_toLabel <- tableNow[which(
    tableNow$Info01 == "01-Mock" & !tableNow$PositiveMitochondria | 
      !tableNow$PositiveMitochondria & tableNow$PositiveCell
  ), ] %>%
    group_by(FileOrigin_base) %>%
    count(MitoCat_D)
  TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
  
  # Loop to maximize
  SummarizedData <- data.frame(row.names=unique(TempTable_plyr_toLabel$FileOrigin_base))
  
  for(i in unique(TempTable_plyr_toLabel$FileOrigin_base)  ){
    tempNow <- TempTable_plyr_toLabel[TempTable_plyr_toLabel$FileOrigin_base == i, ]
    #
    SummarizedData[i,"Category"] <- tempNow$MitoCat_D[which.max(tempNow$n)]
    remove(tempNow)
  }
  # Add to column the cell information
  if(is.null(SummarizedData[tableNow$FileOrigin_base, "Category"])){
    
  } else {
    tableNow$CellCat_D_negative <- SummarizedData[tableNow$FileOrigin_base, "Category"]
  }
  remove(SummarizedData)
  
  # 
  #
  tableNow <- tableNow[which(!(is.na(tableNow$CellCat_D) & is.na(tableNow$CellCat_D_negative))),]
  
  tableNow[,c("CellCat_D",
              "CellCat_D_negative")]
  
  return(tableNow)
}










# fake param
minProteinLabel <- 300
minProteinMaxLabel <- 0
minProteinMedianLabel <- 0
minProteinStdDevLabel <- 0

# RawData_size_mitol_list 
for(minProteinStdDevLabel in c(0)){
  for(minProteinMedianLabel in c(0)){

      for(minProteinMaxLabel in c(0, 1000)){ #, 1000, 2000, 2500
        for(minProteinLabel in c(0, 700, 1000)){ #0, 300, 500, 
          #
          
          #Do no all zero
          if(sum(c(minProteinStdDevLabel,
                   minProteinMedianLabel,
                   minProteinMaxLabel,
                   minProteinLabel)) > 0 # & minProteinMaxLabel > minProteinLabel
				   ){
            
            
            print(paste0("minProteinLabel now: ", minProteinLabel))
            print(paste0("minProteinMaxLabel now: ", minProteinMaxLabel))
            print(paste0("minProteinMedianLabel now: ", minProteinMedianLabel))
            print(paste0("minProteinStdDevLabel now: ", minProteinStdDevLabel))
            
            ### COLLECT RDATA
            print(paste0("Coletando os dados novamente"))
            load(paste0(rootFolder, "/", "02-PlatesRData", "/", ActualRData))
            print(paste0("Dados coletados"))
            
            
            # Trimm table
            RawData_size_mitol_mitosig_list <- pblapply(X=RawData_size_mitol_list,
                                                        FUN=GateProteinaMitocondria_function)
            
            #
            remove(RawData_size_mitol_list)
            gc()
            
            
            # Create actual loop dir
            dir.create(path = paste0(rootFolder, "/",
                                     "04-Morphology", "/",
                                     ActualRData_name, "/",
                                     "MinL_", minProteinLabel,
                                     "minMaxD_",minMaxDProteinLabel,
                                     "MinMax_", minProteinMaxLabel,
                                     "MinMedian_", minProteinMedianLabel,
                                     "MinStdDev_", minProteinStdDevLabel,
                                     "/"), showWarnings = FALSE)
            
            
            # PERIMETER DEFINITION
            minPerimNow <- 4.3
            maxPerimNow <- 7.3
            
            # D DEFINITION
            minDNow <- 0.6
            maxDNow <- 0.9
            
            
            print(paste0("Iniciando loops de perimetro para esse corte intensidade"))
            
            ### ONLY PERIMETER
            for(minPerimNow in c(3.7, 4.3)   ){ #seq(4,6,by=0.1) #3.7, 4.0, 4.2, 4.3
              for(maxPerimNow in c(6.5, 7.3) ){ #6.5, 7.0, 7.1, 7.3

                
                #
                print(paste0("minPerimNow:", minPerimNow, " maxPerimNow:", maxPerimNow))   
                ####
                # PERIMETER
                RawData_size_mitol_mitosig_clas1_list <- pblapply(X=RawData_size_mitol_mitosig_list,
                                                                  FUN=PerimeterThreshold_function)
                # Remove old list
                #remove(RawData_size_mitol_mitosig_list)
                RawData_size_mitol_mitosig_clas1_list -> RawData_size_mitol_mitosig_classF_list
                remove(RawData_size_mitol_mitosig_clas1_list)
                gc()
                
                
                ################################ PLOTS PER CONDITION
                ########### Per condition
                #
                print(paste0("coletando as colunas")) 
                # Retrieve column with data
                OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                                 FUN=ColumnRetrieval,
                                                                 columnNow="FileOrigin_base")))
                
                # Retrieve column with data
                OverallUniqueTag <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                           FUN=ColumnRetrieval,
                                                           columnNow="UniqueTag")))
                
                # Retrieve column with data
                OverallUniqueConditionTag <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                                    FUN=ColumnRetrieval,
                                                                    columnNow="UniqueConditionTag")))
                
                # Retrieve column with data
                OverallMitoCat_Perim <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                               FUN=ColumnRetrieval,
                                                               columnNow="MitoCat_Perim")))
                
                # Retrieve column with data
                OverallCellCat_Perim <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                               FUN=ColumnRetrieval,
                                                               columnNow="CellCat_Perim")))
                
                
                
                # Retrieve column with data
                OverallCellCat_Perim_negative <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                                        FUN=ColumnRetrieval,
                                                                        columnNow="CellCat_Perim_negative")))
                
                
                
                OverallCellLine <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                          FUN=ColumnRetrieval,
                                                          columnNow="CellLine"))) 
                OverallInfo03 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="Info03"))) 
                OverallInfo01 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="Info01"))) 
                OverallInfo04 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="Info04"))) 
                OverallInfo02 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="Info02"))) 
                OverallLabeling <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                          FUN=ColumnRetrieval,
                                                          columnNow="Labeling"))) 
                
                OverallPerim <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                          FUN=ColumnRetrieval,
                                                          columnNow="Perim."))) 
                OverallD <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                       FUN=ColumnRetrieval,
                                                       columnNow="D"))) 
                OverallProtein_Kurt <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                   FUN=ColumnRetrieval,
                                                   columnNow="Protein_Kurt"))) 
                OverallPositiveMitochondria <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_classF_list,
                                                              FUN=ColumnRetrieval,
                                                              columnNow="PositiveMitochondria"))) 
                
                #
                TemporaryData <- (data.frame(OverallFileOrigin_base,
                                             OverallUniqueTag,
                                             OverallUniqueConditionTag,
                                             OverallMitoCat_Perim,
                                             OverallCellCat_Perim,
                                             
                                             
                                             OverallCellCat_Perim_negative,
                                             OverallPerim,
                                             OverallD,
                                             OverallProtein_Kurt,
                                             OverallCellLine,
                                             OverallInfo01,
                                             OverallInfo02,
                                             OverallInfo03,
                                             OverallInfo04,
                                             OverallPositiveMitochondria,
                                             OverallLabeling))
                
                remove(OverallFileOrigin_base)
                remove(OverallUniqueTag)
                remove(OverallUniqueConditionTag)
                remove(OverallMitoCat_Perim)
                remove(OverallCellCat_Perim)
                #remove(OverallMitoCat_D)
                #remove(OverallCellCat_D)
                remove(OverallCellCat_Perim_negative)
                #remove(OverallCellCat_D_negative)
                remove(OverallPositiveMitochondria)
                remove(OverallCellLine)
                remove(OverallInfo01)
                remove(OverallInfo02)
                remove(OverallInfo03)
                remove(OverallInfo04)
                remove(OverallLabeling)
                remove(OverallProtein_Kurt)
                remove(OverallD)
                remove(OverallPerim)
                
                print(paste0("Salvando tabela")) 
                
                write.csv2(TemporaryData,
                           file=paste0(rootFolder, "/",
                                       "04-Morphology", "/",
                                       ActualRData_name, "/",
                                       "MinL_", minProteinLabel,
                                       "minMaxD_",minMaxDProteinLabel,
                                       "MinMax_", minProteinMaxLabel,
                                       "MinMedian_", minProteinMedianLabel,
                                       "MinStdDev_", minProteinStdDevLabel,
                                       "/", "03_Table_MinP_", minPerimNow,
                                       "_MaxP_", maxPerimNow,
                                       ".csv"))
                
                TemporaryData <- TemporaryData[with(TemporaryData, order(
                  OverallInfo03,
                  OverallCellLine,
                  OverallInfo01,
                  OverallInfo04,
                  OverallInfo02
                )), ]
                TemporaryData$OverallLabeling <- factor(TemporaryData$OverallLabeling, levels=unique(TemporaryData$OverallLabeling))
                
                
                
                ########################
                
                parameterMitoNow <- "OverallD"
                parameterGACNow <- "OverallProtein_Kurt"
                condNow <- unique(TemporaryData$OverallUniqueConditionTag)[6]
                
                mitoStatus <- TRUE
                MitoCatNow <- unique(TemporaryData$OverallMitoCat_Perim)[2]
                
                parameterMitoNow <- "OverallD"
                parameterGACNow <- "OverallProtein_Kurt"
                
                
                    for(parameterMitoNow in c("OverallD",
                                              "OverallPerim")){
                      for(parameterGACNow in c("OverallProtein_Kurt")){
                        for(condNow in unique(TemporaryData$OverallUniqueConditionTag)){
                          for(mitoStatus in c(TRUE, FALSE)){
                            for(MitoCatNow in "3-Tubular"){ #unique(TemporaryData$OverallMitoCat_Perim)
                          
                          PlotData <- TemporaryData[which(TemporaryData$OverallUniqueConditionTag %in% condNow &
                                                            TemporaryData$OverallPositiveMitochondria %in% mitoStatus &
                                                            TemporaryData$OverallMitoCat_Perim %in% MitoCatNow),
                                                    c(parameterMitoNow, parameterGACNow, "OverallUniqueConditionTag")]
                          
                          PlotData <- PlotData[which(is.finite(PlotData[,parameterMitoNow]) &
                                                       is.finite(PlotData[,parameterGACNow]) ),]
                          
                          if(dim(PlotData)[1] > 5){
                            
                            PlotData$Sampled <- FALSE
                            PlotData$Sampled[sample(1:dim(PlotData)[1], min(c(5000, dim(PlotData)[1])))] <- TRUE
                            
                            
                            fit <- lm(formula=get(parameterGACNow)~get(parameterMitoNow), data=PlotData)
                            
                            
                            p <- ggplot(PlotData[PlotData$Sampled,], aes(x = get(parameterMitoNow),
                                                                         y = (get(parameterGACNow)),
                                                                         color = OverallUniqueConditionTag)) +
                              geom_point(alpha=0.4, shape=16) +
                              theme_classic() +
                              theme(legend.position = "bottom") +
                              coord_cartesian(ylim=c(quantile(PlotData[PlotData$Sampled, parameterGACNow], 0.01, na.rm=TRUE),
                                                     quantile(PlotData[PlotData$Sampled, parameterGACNow], 0.97, na.rm=TRUE)),
                                              xlim=c(quantile(PlotData[PlotData$Sampled, parameterMitoNow], 0.01, na.rm=TRUE),
                                                     quantile(PlotData[PlotData$Sampled, parameterMitoNow], 0.99, na.rm=TRUE))
                              ) +
                              guides(color=guide_legend(nrow=2,byrow=FALSE)) +
                              scale_color_manual(values=c(rgb(77,146,33,  220, max=255),
                                                          rgb(197,27,125, 200, max=255),
                                                          rgb(215,48,39, 200, max=255),
                                                          rgb(69,117,180, 200, max=255) )  ) +
                              # stat_summary(fun.data= mean_cl_normal) + 
                              geom_smooth(method="lm",
                                          se = TRUE) +
                              xlab(paste0(parameterMitoNow))+
                              ylab(paste0(parameterGACNow)) +
                              ggtitle(paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                             
                                             " Slope = ",signif(fit$coef[[2]], 5)))
                            #scale_x_continuous(limits = c(00, 1.6))
                            #formula = y ~ poly(x, 2)) 
                            
                            
                            
                            
                            # Marginal histograms by group
                            p <- ggMarginal(p,
                                            type = "density", #c("density", "histogram", "boxplot", "violin", "densigram")
                                            groupColour = TRUE,
                                            margins = "both",
                                            groupFill = TRUE)
                            
                            png(file=paste0(rootFolder, "/",
                                            "04-Morphology", "/",
                                            ActualRData_name, "/",
                                            "MinL_", minProteinLabel,
                                            "minMaxD_",minMaxDProteinLabel,
                                            "MinMax_", minProteinMaxLabel,
                                            "MinMedian_", minProteinMedianLabel,
                                            "MinStdDev_", minProteinStdDevLabel,
                                            "/", "06_",
                                            "P_", minPerimNow,
                                            "-", maxPerimNow,
                                            "Mito", mitoStatus,
                                            "Cat_",MitoCatNow,
                                            "-_", parameterMitoNow,
                                            "vs", parameterGACNow, "_", condNow, ".png"),
                                width=2000, height=2000, res=400)
                            print(p)
                            dev.off()
                          }
                          
                          
                          
                          
                        }
                      }
                    }
                    
                    
                    
                  }
                }

                
                
                ##########
                ########## END
                ##########
                ########## 06 - MORPHOLOGY PLOTS
                ##########
                
                
              }
              }
            
          }}}}}

  
#!/usr/bin/env /mnt/nfs/home/douglas/anaconda3/envs/R413/bin/Rscript --vanilla
#
#  03_HPC_Mitophagy_v1.R
#
#
#

#
args = commandArgs(trailingOnly=TRUE)

# Available folders
#ActualRData <- "A99.29_ectopica_parkin_controle_minMitoLabel_300.RData"

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


# If loop to run everything at once
# utils::memory.limit(size=32000)

# Create main data
dir.create(path = paste0(rootFolder, "/", "09-MitophagyTubular"), showWarnings = FALSE)

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
dir.create(path = paste0(rootFolder, "/", "09-MitophagyTubular", "/", ActualRData_name, "/"),
           showWarnings = FALSE)


  ####### LABEL POSITIVE MITOCHONDRIA
  ### Create function
  GateProteinaMitocondria_function <- function(tableNow){
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
  
  # Create function
  MitophagyMeanAdd_function <- function(tableNow){
    
    ##### Define mitophagy ratio
    tableNow$PositiveMitophagySensor <- FALSE
    tableNow$RatioMitophagy_Cellwise <- (tableNow$Mitophagy_Mean / tableNow$CellMeanMitophagy)
    # Calcula as razões na célula entre média
    
    
    # Loop to each cell
    # Dplyr para gerar os valores
    MinMitoLabel_plyr <- tableNow %>%
      group_by(FileOrigin_base) %>%
      dplyr::summarize(CellMitophagy_MinMean = mean(Mitophagy_Min, na.rm=TRUE),
                       TotalGACs = n()
      )
    data.frame(MinMitoLabel_plyr) -> MinMitoLabel_plyr
    MinMitoLabel_plyr$FileOrigin_base -> rownames(MinMitoLabel_plyr)
    
    ####
    tableNow[,c("CellMitophagy_MinMean")] <- MinMitoLabel_plyr[tableNow$FileOrigin_base, c("CellMitophagy_MinMean") ]
    
    tableNow$RatioMitophagy_Cellwise_MinMean <- (tableNow$Mitophagy_Mean / tableNow$CellMitophagy_MinMean)
    
    # Define positive cell
    # Marca as células positivas
    tableNow$PositiveMitophagySensor <- FALSE
    tableNow$PositiveMitophagySensor[tableNow$CellMeanMitophagy > minMitophagyLabel] <- TRUE
    
    return(tableNow)
  }
  
  
  
  
  
  # fake param
  minMitophagyLabel <- 100
  minProteinLabel <- 2000
  minProteinMaxLabel <- 0
  minProteinMedianLabel <- 0
  minProteinStdDevLabel <- 0
  
  # RawData_size_mitol_list 
  for(minProteinStdDevLabel in c(0)){
    for(minProteinMedianLabel in c(0)){
      for(minMitophagyLabel in c(100)){
        for(minProteinMaxLabel in c(0)){
          for(minProteinLabel in c(2000)){
            #
            
            #Do no all zero
            if(sum(c(minProteinStdDevLabel,
                     minProteinMedianLabel,
                     minProteinMaxLabel,
                     minProteinLabel)) > 0){
              
              print(paste0("minProteinLabel now: ", minProteinLabel))
              print(paste0("minProteinMaxLabel now: ", minProteinMaxLabel))
              print(paste0("minProteinMedianLabel now: ", minProteinMedianLabel))
              print(paste0("minProteinStdDevLabel now: ", minProteinStdDevLabel))
              
              ### COLLECT RDATA
              print(paste0("Coletando os dados novamente"))
              load(paste0(rootFolder, "/", "02-PlatesRData", "/", ActualRData))
              print(paste0("Dados coletados"))
              
              
              print(paste0("Fazendo gate de proteína"))
              # Trimm table
              RawData_size_mitol_list <- pblapply(X=RawData_size_mitol_list,
                                                  FUN=GateProteinaMitocondria_function)

              print(paste0("Fazendo gate de mitofagia"))
              # Apply dataset-wise
              RawData_size_mitol_list <- pblapply(X=RawData_size_mitol_list,
                                                  FUN=MitophagyMeanAdd_function)

              #
              gc()
              RawData_size_mitol_mitosig_list <- RawData_size_mitol_list
              remove(RawData_size_mitol_list)
              gc()
              #

              # Create actual loop dir
              dir.create(path = paste0(rootFolder, "/",
                                       "09-MitophagyTubular", "/",
                                       ActualRData_name, "/",
                                       "MinL_", minProteinLabel,
                                       "minMaxD_",minMaxDProteinLabel,
                                       "MinMax_", minProteinMaxLabel,
                                       "MinMedian_", minProteinMedianLabel,
                                       "MinStdDev_", minProteinStdDevLabel,
                                       "MinMitophagy_", minMitophagyLabel,
                                       "/"), showWarnings = FALSE)
              
              
              ##########
              ########## 01 - MITOCHONDRIA POSITIVITY PLOTS
              ##########
              ########## START
              ##########
              
              
              ########## PLOT FRACTION POSITIVE CELLS PER WELL
              # Retrieve column with data
              OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                               FUN=ColumnRetrieval,
                                                               columnNow="FileOrigin_base")))
              
              OverallUniqueConditionTag <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                                  FUN=ColumnRetrieval,
                                                                  columnNow="UniqueConditionTag")))
              
              OverallPositiveCell <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                            FUN=ColumnRetrieval,
                                                            columnNow="PositiveCell")))
              
              OverallFractionMitoPositive <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                                    FUN=ColumnRetrieval,
                                                                    columnNow="FractionMitoPositive"))) 
              
              OveralltotalMito <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                         FUN=ColumnRetrieval,
                                                         columnNow="totalMito"))) 
              
              OverallCellMeanProtein <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                               FUN=ColumnRetrieval,
                                                               columnNow="CellMeanProtein"))) 
              
              
              
              
              OverallUniqueTag <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                         FUN=ColumnRetrieval,
                                                         columnNow="UniqueTag"))) 
              
              
              OverallCellLine <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="CellLine"))) 
              OverallInfo03 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info03"))) 
              OverallInfo01 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info01"))) 
              OverallInfo04 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info04"))) 
              OverallInfo02 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info02"))) 
              OverallLabeling <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="Labeling"))) 
              
              
              
              TemporaryData <- (data.frame(OverallFileOrigin_base,
                                           OverallUniqueConditionTag,
                                           OverallPositiveCell,
                                           OveralltotalMito,
                                           OverallFractionMitoPositive,
                                           OverallUniqueTag,
                                           OverallCellMeanProtein,
                                           OverallCellLine,
                                           OverallInfo01,
                                           OverallInfo02,
                                           OverallInfo03,
                                           OverallInfo04,
                                           OverallLabeling))
              
              # Making unique
              remove(OverallFileOrigin_base)
              remove(OverallUniqueConditionTag)
              remove(OverallPositiveCell)
              remove(OverallFractionMitoPositive)
              remove(OverallUniqueTag)
              remove(OveralltotalMito)
              remove(OverallCellLine)
              remove(OverallInfo01)
              remove(OverallInfo02)
              remove(OverallInfo03)
              remove(OverallInfo04)
              remove(OverallLabeling)
              remove(OverallCellMeanProtein)
              
              gc()
              
              print("Making unique Temporary Data for boxplot")
              TemporaryData <- unique(TemporaryData)
              print("Finished unique Temporary Data for boxplot")
              
              
              print(paste0("Ordenando temporary data para quantidade mito"))
              TemporaryData <- TemporaryData[with(TemporaryData, order(
                OverallInfo03,
                OverallCellLine,
                OverallInfo01,
                OverallInfo04,
                OverallInfo02
              )), ]
              TemporaryData$OverallLabeling <- factor(TemporaryData$OverallLabeling, levels=unique(TemporaryData$OverallLabeling))
              
              #
              gc()
              
              
              print(paste0("Plotting Boxplot quantidade mito"))
              
              p <- ggplot(data=unique(TemporaryData[!is.na(TemporaryData$OverallFractionMitoPositive), c("OverallLabeling", "OverallUniqueConditionTag", "OverallFractionMitoPositive")]),
                          aes(x=OverallLabeling,
                              y=OverallFractionMitoPositive)) + 
                geom_boxplot(outlier.shape = NA) +
                theme_classic() + geom_jitter(shape=16,
                                              position=position_jitter(0.2),
                                              size=0.3,
                                              alpha=0.1) +
                theme(axis.text.x = element_text(angle = 0),
                      axis.title.x=element_text(size=5, color="black")) +
                coord_cartesian(clip = "off",
                                ylim = c(0,100)) +
                ylab("Mean mitochondrial positivity\nper cell (%)")+
                xlab(paste0(ActualRData_name, "\n",
                            "MinL_", minProteinLabel,
                            "minMaxD_",minMaxDProteinLabel,
                            "MinMax_", minProteinMaxLabel,
                            "\n",
                            "MinMedian_", minProteinMedianLabel,
                            "MinStdDev_", minProteinStdDevLabel,
                            "MinMitophagy_", minMitophagyLabel))
              
              
              png(paste0(rootFolder, "/",
                         "09-MitophagyTubular", "/",
                         ActualRData_name, "/",
                         "MinL_", minProteinLabel,
                         "minMaxD_",minMaxDProteinLabel,
                         "MinMax_", minProteinMaxLabel,
                         "MinMedian_", minProteinMedianLabel,
                         "MinStdDev_", minProteinStdDevLabel,
                         "MinMitophagy_", minMitophagyLabel,
                         "/", "01-PercentPositity_boxplot.png"),
                  width=4600, height=1800, res=600)
              print(p+geom_text(color="black",
                                x = 0.4,
                                y = -10,
                                label= MyLabelsNow,
                                hjust = 1,
                                vjust = 1,
                                size = 2.4)
              )
              dev.off()
              

              print(paste0("Plotting Boxplot média mitos mito")) 
              
              p <- ggplot(data=unique(TemporaryData[which(!is.na(TemporaryData$OverallCellMeanProtein)), c("OverallLabeling", "OverallUniqueConditionTag", "OverallCellMeanProtein")]),  # & (TemporaryData$OverallInfo01 == "01-Mock")
                          aes(x=OverallLabeling,
                              y=log2(OverallCellMeanProtein+1))) + 
                geom_boxplot(outlier.shape = NA) +
                theme_classic() + geom_jitter(shape=16,
                                              position=position_jitter(0.2),
                                              size=0.3,
                                              alpha=0.1) +
                theme(axis.text.x = element_text(angle = 0),
                      axis.title.x=element_text(size=5, color="black")) +
                coord_cartesian(clip = "off") +
                ylab("Mean mitochondrial\nmKO2 labelling") +
                xlab(paste0(ActualRData_name, "\n",
                            "MinL_", minProteinLabel,
                            "minMaxD_",minMaxDProteinLabel,
                            "MinMax_", minProteinMaxLabel,
                            "\n",
                            "MinMedian_", minProteinMedianLabel,
                            "MinStdDev_", minProteinStdDevLabel,
                            "MinMitophagy_", minMitophagyLabel))
              
              
              png(paste0(rootFolder, "/",
                         "09-MitophagyTubular", "/",
                         ActualRData_name, "/",
                         "MinL_", minProteinLabel,
                         "minMaxD_",minMaxDProteinLabel,
                         "MinMax_", minProteinMaxLabel,
                         "MinMedian_", minProteinMedianLabel,
                         "MinStdDev_", minProteinStdDevLabel,
                         "MinMitophagy_", minMitophagyLabel,
                         "/", "02-MeanIntensity_boxplot.png"),
                  width=4600, height=1800, res=600)
              print(p+geom_text(color="black",
                                x = 0.4,
                                y = -10,
                                label= MyLabelsNow,
                                hjust = 1,
                                vjust = 1,
                                size = 2.4)
              )
              dev.off()
              
              
              
              
              
              
              
              # SAVE TABLES
              
              write.csv2(TemporaryData,
                         file=paste0(rootFolder, "/",
                                     "09-MitophagyTubular", "/",
                                     ActualRData_name, "/",
                                     "MinL_", minProteinLabel,
                                     "minMaxD_",minMaxDProteinLabel,
                                     "MinMax_", minProteinMaxLabel,
                                     "MinMedian_", minProteinMedianLabel,
                                     "MinStdDev_", minProteinStdDevLabel,
                                     "MinMitophagy_", minMitophagyLabel,
                                     "/", "02-PercentPositity_boxplot.csv")
              )
              
              
              print(paste0("Preparando tabela com quantidade de celulas"))
              ####################### BAR PLOT
              # Dplyr para gerar os valores
              TempPositiveTable_plyr <- TemporaryData[TemporaryData$OverallPositiveCell ,] %>%
                group_by(OverallUniqueTag, OverallUniqueConditionTag) %>%
                dplyr::summarize(Frequency = mean(OverallFractionMitoPositive, na.rm=TRUE),
                                 TotalCell = n()
                )
              
              write.csv2(TempPositiveTable_plyr,
                         file=paste0(rootFolder, "/",
                                     "09-MitophagyTubular", "/",
                                     ActualRData_name, "/",
                                     "MinL_", minProteinLabel,
                                     "minMaxD_",minMaxDProteinLabel,
                                     "MinMax_", minProteinMaxLabel,
                                     "MinMedian_", minProteinMedianLabel,
                                     "MinStdDev_", minProteinStdDevLabel,
                                     "MinMitophagy_", minMitophagyLabel,
                                     "/", "02-PercentPositity_percell.csv")
              )
              
              #
              remove(p)
              remove(TemporaryData)
              remove(TempPositiveTable_plyr)
              
              ##########
              ########## END
              ##########
              ########## 05 - MITOCHONDRIA POSITIVITY PLOTS
              ##########
              
              RawData_size_mitophagy_list <- RawData_size_mitol_mitosig_list
              remove(RawData_size_mitol_mitosig_list)
              gc()
              #
              
              #
              # Retrieve column with data
              OverallRatioMitophagy_Cellwise <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                                       FUN=ColumnRetrieval,
                                                                       columnNow="RatioMitophagy_Cellwise")))
              
              OverallRatioMitophagy_Cellwise_MinMean <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                                               FUN=ColumnRetrieval,
                                                                               columnNow="RatioMitophagy_Cellwise_MinMean")))
              
              
              OverallUniqueConditionTag<- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                                 FUN=ColumnRetrieval,
                                                                 columnNow="UniqueConditionTag")))
              OverallUniqueTag<- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="UniqueTag")))
              OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                               FUN=ColumnRetrieval,
                                                               columnNow="FileOrigin_base")))
              
              OverallPositiveCell<- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                           FUN=ColumnRetrieval,
                                                           columnNow="PositiveCell")))
              
              OverallPositiveMitochondria<- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                                   FUN=ColumnRetrieval,
                                                                   columnNow="PositiveMitochondria")))
              
              OverallPositiveMitophagySensor <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                                       FUN=ColumnRetrieval,
                                                                       columnNow="PositiveMitophagySensor")))
              
              OverallCellLine <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                        FUN=ColumnRetrieval,
                                                        columnNow="CellLine")))
              
              OverallInfo01 <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info01")))
              
              OverallInfo02 <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info02")))
              
              OverallInfo03 <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info03")))
              
              
              OverallInfo04 <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                      FUN=ColumnRetrieval,
                                                      columnNow="Info04")))
              
              
              OverallLabeling<- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                       FUN=ColumnRetrieval,
                                                       columnNow="Labeling")))
              
              #
              OverallPerim <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                       FUN=ColumnRetrieval,
                                                       columnNow="Perim.")))
              
              
              
              #
              TemporaryData <- (data.frame(OverallFileOrigin_base,
                                           OverallUniqueTag,
                                           OverallUniqueConditionTag,
                                           OverallLabeling,
                                           OverallCellLine,
                                           OverallInfo01,
                                           OverallInfo02,
                                           OverallInfo03,
                                           OverallInfo04,
                                           OverallPerim,
                                           OverallRatioMitophagy_Cellwise,
                                           OverallPositiveMitophagySensor,
                                           OverallPositiveCell,
                                           OverallPositiveMitochondria,
                                           OverallRatioMitophagy_Cellwise_MinMean
              ))
              # get only cell with sensor
              TemporaryData <- TemporaryData[TemporaryData$OverallPositiveMitophagySensor,]
              
              # drop column
              TemporaryData$OverallPositiveMitophagySensor <- NULL
              
              remove(OverallFileOrigin_base)
              remove(OverallUniqueTag)
              remove(OverallUniqueConditionTag)
              remove(OverallRatioMitophagy_Cellwise)
              remove(OverallPositiveMitochondria)
              remove(OverallPositiveMitophagySensor)
              remove(OverallPositiveCell)
              remove(OverallRatioMitophagy_Cellwise_MinMean)
              remove(OverallCellLine)
              remove(OverallInfo01)
              remove(OverallInfo02)
              remove(OverallInfo03)
              remove(OverallInfo04)
              remove(OverallLabeling)
              remove(OverallPerim)
              
              # Block positive mitochondria in NEGATIVE cells
              # V0.3 addition
              TemporaryData[which(!TemporaryData$OverallPositiveCell), "OverallPositiveMitochondria"] <- FALSE
              
              
              # Remove a lista pra liberar memória
              remove(RawData_size_mitophagy_list)
              gc()
              gc()
              
              ####
              TemporaryData <- TemporaryData[with(TemporaryData, order(
                OverallCellLine,
                OverallInfo03,
                OverallInfo01,
                OverallInfo04,
                OverallInfo02
              )), ]
              TemporaryData$OverallLabeling <- factor(TemporaryData$OverallLabeling, levels=unique(TemporaryData$OverallLabeling))
              gc()
              
              
              #####OverallPerim
              
              # Create outside copy for perimeter thresholds
              TemporaryData -> TemporaryDataOutside
              remove(TemporaryData)
              
              TemporaryDataOutside$MitoCat_Perim <- NA
              
              ##
              RatioMitophagyCutoff <- 1.41
              
              ##
              possibleMitophagyRatios <- c(seq(from=1.41, to=1.41, by=0.1))
              
              ## 
              
              ##
              
              #### MITOPHAGY CUTS
              for(RatioMitophagyCutoff in possibleMitophagyRatios){
                
                ##
                print(paste0("Estou no corte ", RatioMitophagyCutoff))
                
                #######
                ####### OVER MEAN OF MINIMUMS CALCULATIONS
                #######
                ####### START
                #######
        
                
                minPerimNow <- 4.3
                maxPerimNow <- 7.3
                
                for(minPerimNow in c(4.3)){
                  for(maxPerimNow in c(7.3)){
                    
                    #
                    print(paste("Minimum perimeter: ", minPerimNow, "; maximum perimeter: ", maxPerimNow))
                    
                    #
                    #TemporaryData$OverallPerim
                    
                    #minPerim
                    #maxPerim
                    #TemporaryDataOutside[,]
                    #TemporaryData
                    TemporaryDataOutside$MitoCat_Perim <- NA
                    
                    # Add intermediate
                    TemporaryDataOutside[which(
                      !is.na(TemporaryDataOutside$OverallPerim)
                    ),"MitoCat_Perim"] <- "2-Intermediate"
                    # Add fragmented
                    TemporaryDataOutside[which(
                      !is.na(TemporaryDataOutside$OverallPerim) &
                        TemporaryDataOutside$OverallPerim < minPerimNow
                    ),"MitoCat_Perim"] <- "1-Fragmented"
                    # Add tubular
                    TemporaryDataOutside[which(
                      !is.na(TemporaryDataOutside$OverallPerim) &
                        TemporaryDataOutside$OverallPerim > maxPerimNow
                    ),"MitoCat_Perim"] <- "3-Tubular"
                    #
                    TemporaryDataOutside <- TemporaryDataOutside[which(!is.na(TemporaryDataOutside$MitoCat_Perim)),]

                    
                    list() -> StoredByMitoType
                    list() -> StoredByMitoType_toBar
                    # SELECT ONLY TUBULAR
                    # mitoTypeNow <- "3-Tubular"
                    for(mitoTypeNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
                      #
                      TemporaryData <- TemporaryDataOutside[which(TemporaryDataOutside$MitoCat_Perim == mitoTypeNow ), ]

                      #
                      
                      
                      #
                      print(paste0("Fazendo calculos com a razao sobre a média das mitocondrias"))
                      
                      TemporaryData$PositiveMitophagy <- FALSE
                      # Define positive mitochondrias using Cutoff
                      TemporaryData$PositiveMitophagy[(
                        TemporaryData$OverallRatioMitophagy_Cellwise_MinMean > RatioMitophagyCutoff
                      )] <- TRUE
                      
                      
                      #
                      print(paste0("Preparando dplyr"))
                      
                      
                      #####
                      ##### Add Number of Mitophagic dots for POSITIVE MITOCHONDRIA
                      #  
                      TempTable_plyr_toLabel_positive <- TemporaryData[which(TemporaryData$OverallPositiveMitochondria & TemporaryData$OverallPositiveCell
                      ), ] %>%
                        group_by(OverallFileOrigin_base, OverallUniqueTag, OverallUniqueConditionTag) %>%
                        count(PositiveMitophagy)
                      
                      TempTable_plyr_toLabel_positive <- data.frame(TempTable_plyr_toLabel_positive)
                      
                      #####
                      ##### Add Number of Mitophagic dots for NEGATIVE MITOCHONDRIA
                      #
                      TempTable_plyr_toLabel_negative <- TemporaryData[which(!TemporaryData$OverallPositiveMitochondria
                      ), ] %>%
                        group_by(OverallFileOrigin_base, OverallUniqueTag, OverallUniqueConditionTag) %>%
                        count(PositiveMitophagy)
                      TempTable_plyr_toLabel_negative <- data.frame(TempTable_plyr_toLabel_negative)
                      
                      # Join tables
                      colnames(TempTable_plyr_toLabel_positive)[5] <- "PositiveProteinMitophagy"
                      colnames(TempTable_plyr_toLabel_negative)[5] <- "NegativeProteinMitophagy"
                      TempTable_plyr_toLabel <- full_join(TempTable_plyr_toLabel_negative,
                                                          TempTable_plyr_toLabel_positive,
                                                          by = c("OverallFileOrigin_base", "OverallUniqueTag", "OverallUniqueConditionTag", "PositiveMitophagy"))
                      
                      remove(TempTable_plyr_toLabel_positive)
                      remove(TempTable_plyr_toLabel_negative)
                      
                      # Add zeroes
                      TempTable_plyr_toLabel[is.na(TempTable_plyr_toLabel$NegativeProteinMitophagy), "NegativeProteinMitophagy"] <- 0
                      TempTable_plyr_toLabel[is.na(TempTable_plyr_toLabel$PositiveProteinMitophagy), "PositiveProteinMitophagy"] <- 0
                      
                      # Loop to add total mitochondria
                      TempTable_plyr_toLabel$TotalConsideredMitochondria <- NA
                      TempTable_plyr_toLabel$TotalProtPosMitochondria <- NA
                      TempTable_plyr_toLabel$TotalProtNegMitochondria <- NA
                      
                      print(paste0("List splitting a dplyr"))
                      
                      # Perform LIST split
                      TempTable_plyr_toLabelSplitted <- split(TempTable_plyr_toLabel,
                                                              f = TempTable_plyr_toLabel[,"OverallFileOrigin_base"])
                      
                      #
                      
                      OverallFixingMitophagy_function <- function(tableNow){
                        if(dim(tableNow)[1] == 2){
                          # nada
                        } else {
                          #print(baseNow)
                          tempLine <- tableNow[1, ]
                          
                          if(tempLine$PositiveMitophagy){
                            # Se for uma linha TRUE, add a false
                            tempLine$PositiveMitophagy <- FALSE
                          } else {
                            # Se for uma linha FALSE, add a true
                            tempLine$PositiveMitophagy <- TRUE
                          }
                          tempLine[,c("NegativeProteinMitophagy", "PositiveProteinMitophagy")] <- 0
                          tableNow <- rbind(tableNow, tempLine)
                          remove(tempLine)
                        }
                        tableNow[, c("TotalConsideredMitochondria")] <-
                          sum(tableNow[,c("NegativeProteinMitophagy", "PositiveProteinMitophagy")])
                        
                        tableNow[, c("TotalProtPosMitochondria")] <-
                          sum(tableNow[,c("PositiveProteinMitophagy")])
                        
                        tableNow[, c("TotalProtNegMitochondria")] <-
                          sum(tableNow[,c("NegativeProteinMitophagy")])
                        return(tableNow)
                      }
                      
                      print(paste0("operando na lista dplyr"))
                      
                      # Trimm table
                      TempTable_plyr_toLabelPerformed <- pblapply(X=TempTable_plyr_toLabelSplitted,
                                                                  FUN=OverallFixingMitophagy_function)
                      
                      ###
                      TempTable_plyr_toLabel <- bind_rows(TempTable_plyr_toLabelPerformed) #, .id = "column_label"
                      print(paste0("Fim operando na lista dplyr"))
                      
                      # Remove old list
                      remove(TempTable_plyr_toLabelSplitted)
                      remove(TempTable_plyr_toLabelPerformed)
                      gc()
                      
                      # ordena
                      TempTable_plyr_toLabel <- TempTable_plyr_toLabel[order(TempTable_plyr_toLabel$OverallFileOrigin_base),]
                      
                      rownames(TempTable_plyr_toLabel) <- 1:dim(TempTable_plyr_toLabel)[1]
                      
                      
                      TempTable_plyr_toLabel$FractionMitophagyNegative <- (100*(TempTable_plyr_toLabel$NegativeProteinMitophagy/TempTable_plyr_toLabel$TotalProtNegMitochondria))
                      TempTable_plyr_toLabel$FractionMitophagyPositive <- (100*(TempTable_plyr_toLabel$PositiveProteinMitophagy/TempTable_plyr_toLabel$TotalProtPosMitochondria))
                      
                      TempTable_plyr_toLabel$FractionMitophagyPositive[is.nan(TempTable_plyr_toLabel$FractionMitophagyPositive)] <- NA
                      
                      
                      
                      
                      # Return +/- labelling
                      TempTable_plyr_toLabel[, c("OverallCellLine", "OverallInfo01", "OverallInfo02", "OverallInfo03", "OverallInfo04",
                                                 "OverallLabeling") ] <-
                        TemporaryData[match(TempTable_plyr_toLabel$OverallUniqueTag,
                                            TemporaryData$OverallUniqueTag), c("OverallCellLine", "OverallInfo01", "OverallInfo02", "OverallInfo03", "OverallInfo04",
                                                                               "OverallLabeling")]
                      
                      
                      
                      # Save actual data
                      write.csv2(TempTable_plyr_toLabel[TempTable_plyr_toLabel$PositiveMitophagy,],
                                 file=paste0(rootFolder, "/",
                                             "09-MitophagyTubular", "/",
                                             ActualRData_name, "/",
                                             "MinL_", minProteinLabel,
                                             "minMaxD_",minMaxDProteinLabel,
                                             "MinMax_", minProteinMaxLabel,
                                             "MinMedian_", minProteinMedianLabel,
                                             "MinStdDev_", minProteinStdDevLabel,
                                             "MinMitophagy_", minMitophagyLabel,
                                             "/", "03_CellWise_minMean_Pmin",
                                             minPerimNow*10,
                                             "Pmax_",maxPerimNow*10,
                                             "_", mitoTypeNow,
                                             "_CutRatio_",
                                             RatioMitophagyCutoff,".csv"))
                      
                      
                      
                      
                      
                      # Change data format
                      # Retain only positive mitophagy line
                      TempTable_plyr_toLabel <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$PositiveMitophagy),]
                      
                      # Sample some points
                      # Label used points
                      TempTable_plyr_toLabel$SampleNow <- FALSE
                      for(condNow in unique(TempTable_plyr_toLabel$OverallUniqueConditionTag) ){
                        
                        TempTable_plyr_toLabel$SampleNow[sample(x=which(TempTable_plyr_toLabel$OverallUniqueConditionTag == condNow),
                                                                size=min(500, min(table(TempTable_plyr_toLabel$OverallUniqueConditionTag))),
                                                                replace=FALSE)] <- TRUE
                        
                      }
                      
                      # Remove some columns
                      TempTable_plyr_toLabel[,c("PositiveMitophagy",
                                                "NegativeProteinMitophagy",
                                                "PositiveProteinMitophagy",
                                                "TotalProtPosMitochondria",
                                                "TotalProtNegMitochondria",
                                                "TotalConsideredMitochondria")] <- NULL
                      
                      
                      # make a copy
                      TempTable_plyr_toLabel2 <- TempTable_plyr_toLabel 
                      TempTable_plyr_toLabel$Classification <- "0-Negative"
                      TempTable_plyr_toLabel2$Classification <- "1-Positive"
                      
                      
                      TempTable_plyr_toLabel$FractionMitophagyNegative -> TempTable_plyr_toLabel$FractionMitophagy
                      TempTable_plyr_toLabel$FractionMitophagyNegative <- NULL
                      TempTable_plyr_toLabel$FractionMitophagyPositive <- NULL
                      
                      TempTable_plyr_toLabel2$FractionMitophagyPositive -> TempTable_plyr_toLabel2$FractionMitophagy
                      TempTable_plyr_toLabel2$FractionMitophagyNegative <- NULL
                      TempTable_plyr_toLabel2$FractionMitophagyPositive <- NULL
                      
                      TempTable_plyr_toLabel <- rbind(TempTable_plyr_toLabel, TempTable_plyr_toLabel2)
                      remove(TempTable_plyr_toLabel2)
                      
                      #
                      TempTable_plyr_toLabel$OverallUniqueConditionTag_Content <- sapply(1:dim(TempTable_plyr_toLabel)[1], function(w){
                        paste0(TempTable_plyr_toLabel$OverallUniqueConditionTag[w], "_", as.character(TempTable_plyr_toLabel$Classification[w]))
                      })
                      
                      #
                      
                      TempTable_plyr_toBar <- TempTable_plyr_toLabel[, ] %>%
                        group_by(OverallUniqueTag,
                                 OverallUniqueConditionTag,
                                 Classification,
                                 OverallUniqueConditionTag_Content,
                                 OverallLabeling,
                                 OverallCellLine,
                                 OverallInfo01,
                                 OverallInfo02,
                                 OverallInfo03,
                                 OverallInfo04) %>%
                        summarize(Fraction_mean = mean(FractionMitophagy, na.rm=TRUE),
                                  Fraction_median = median(FractionMitophagy, na.rm=TRUE),
                                  totalMito = n()
                        )
                      data.frame(TempTable_plyr_toBar) -> TempTable_plyr_toBar
                      
                      #TempTable_plyr_toBar$OverallLabeling
                      
                      #
                      ## CREATE LISTS FOR THE PLOTS
                      
                      GACList <- list()
                      #for(i in 1:length( unique(TempTable_plyr_toLabel$OverallInfo01)  ) ){
                      #  GACList[[i]] <- unique(TempTable_plyr_toLabel$OverallInfo01)[i]
                      #}
                      #GACList[[length(GACList)+1]] <- unique(TempTable_plyr_toLabel$OverallInfo01)
                      #GACList[[length(GACList)]] <- GACList[[length(GACList)]][GACList[[length(GACList)]] != "01-Mock"]
                      GACList[[1]] <- unique(TempTable_plyr_toLabel$OverallInfo01)
                      # GACList[[1]] <- GACList[[length(GACList)]][GACList[[length(GACList)]] != "01-Mock"]
                      
                      #
                      if("02-GAC.wt" %in% GACList[[1]]){
                        GACList[[1]] <- GACList[[1]][GACList[[1]]!="01-Mock"]
                      }
                      
                      CellList <- list()
                      for(i in 1:length( unique(TempTable_plyr_toLabel$OverallCellLine)  ) ){
                        CellList[[i]] <- unique(TempTable_plyr_toLabel$OverallCellLine)[i]
                      }
                      # CellList[[length(CellList)+1]] <- unique(TempTable_plyr_toLabel$OverallCellLine)
                      
                      
                      
                      #
                      MitophagyList <- list()
                      for(i in 1:length( unique(TempTable_plyr_toLabel$OverallInfo03)  ) ){
                        MitophagyList[[i]] <- unique(TempTable_plyr_toLabel$OverallInfo03)[i]
                      }
                      # MitophagyList[[length(MitophagyList)+1]] <- unique(TempTable_plyr_toLabel$OverallInfo03)
                      
                      
                      colorValues_noalpha_all <- c(rgb(255,255,255,255,max=255),
                                                   rgb(255,120,180,255,max=255),
                                                   
                                                   rgb(255,255,255,255,max=255),
                                                   rgb(255,0,0,255,max=255),
                                                   
                                                   rgb(255,255,255,255,max=255),
                                                   rgb(0,210,140,255,max=255),
                                                   
                                                   rgb(255,255,255,255,max=255),
                                                   rgb(0,160,255,255,max=255))
                      
                      
                      
                      colorValues_border_noalpha_all <- c(colorValues_noalpha_all[2],
                                                          colorValues_noalpha_all[2],
                                                          
                                                          colorValues_noalpha_all[4],
                                                          colorValues_noalpha_all[4],
                                                          
                                                          colorValues_noalpha_all[6],
                                                          colorValues_noalpha_all[6],
                                                          
                                                          colorValues_noalpha_all[8],
                                                          colorValues_noalpha_all[8])
                      
                      
                      
                      
                      
                      

                      
                      #### ADD LOOP HERE
                      #GACNow <- GACList[[5]]
                      #MitophagyNow <- MitophagyList[[3]]
                      for(CellNOW in CellList){
                        for(GACNow in GACList){
                          for(MitophagyNow in MitophagyList ){
                            
                            print(paste0("MinL_", minProteinLabel,
                                         "minMaxD_",minMaxDProteinLabel,
                                         "MinMax_", minProteinMaxLabel,
                                         "MinMedian_", minProteinMedianLabel,
                                         "MinStdDev_", minProteinStdDevLabel,
                                         "MinMitophagy_", minMitophagyLabel,
                                         "/",
                                         "Groups_", paste0(GACNow, collapse="_"),
                                         "Cell_", paste0(CellNOW, collapse="_"),
                                         "_Mitophagy_", paste0(MitophagyNow, collapse="_"),
                                         "_ratio_", RatioMitophagyCutoff,
                                         "Pmin", minPerimNow*10,
                                         "Pmax_",maxPerimNow*10,
                                         "_", mitoTypeNow
                            ))
                            
                            
                            # Color selection
                            if(length(GACNow) == 4){
                              colorValues_noalpha_all[c(1,3,5,7,2,4,6,8)] -> colorValues_noalpha
                              colorValues_border_noalpha_all[c(1,3,5,7,2, 4,6,8)] -> colorValues_border_noalpha
                            } else {
                              if(length(GACNow) == 3){
                                colorValues_noalpha_all[c(3,5,7,4,6,8)] -> colorValues_noalpha
                                colorValues_border_noalpha_all[c(3,5,7,4,6,8)] -> colorValues_border_noalpha
                              } else {
                                if(GACNow == "02-GAC.wt"){
                                  colorValues_noalpha_all[3:4] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[3:4] -> colorValues_border_noalpha
                                }
                                if(GACNow == "03-GAC.K320A"){
                                  colorValues_noalpha_all[5:6] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[5:6] -> colorValues_border_noalpha
                                }
                                if(GACNow == "04-GAC.R382D"){
                                  colorValues_noalpha_all[7:8] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[7:8] -> colorValues_border_noalpha
                                }
                                if(GACNow == "01-Mock"){
                                  colorValues_noalpha_all[1:2] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[1:2] -> colorValues_border_noalpha
                                }
                              }
                            }
                            
                            ####
                            TempTable_plyr_toLabel[which(
                              TempTable_plyr_toLabel$OverallInfo01 %in% GACNow &
                                TempTable_plyr_toLabel$OverallCellLine %in% CellNOW &
                                TempTable_plyr_toLabel$OverallInfo03 %in% MitophagyNow
                            ),] -> DataNow
                            
                            
                            
                            # DEFINE AXIS LIMITS
                            yAxisLimits <- c(0, 
                                             floor(max(
                                               
                                               sapply(unique(DataNow$OverallUniqueConditionTag), function(w){
                                                 quantile(DataNow$FractionMitophagy[which(DataNow$OverallUniqueConditionTag == w)],
                                                          0.75, na.rm=TRUE)
                                               }), na.rm=TRUE)+10))
                            ####
                            
                            #### 
                            p <- ggplot(data = DataNow,
                                        aes(x = OverallLabeling, #OverallUniqueConditionTag
                                            y = FractionMitophagy,
                                            colour = interaction(OverallInfo01, Classification, sep=':'),
                                            fill= interaction(OverallInfo01, Classification, sep=':'))) + #, fill = Classification
                              geom_violin(trim=TRUE,
                                          position=position_dodge(width=0.8),
                                          lwd=0.35)+
                              # scale_fill_manual(values=colorValues_noalpha) +
                              scale_fill_manual(values = colorValues_noalpha,
                                                aesthetics = c("fill")) +
                              scale_colour_manual(values = colorValues_border_noalpha,
                                                  aesthetics = c("colour")) +
                              theme_classic() +
                              coord_cartesian(clip = "off",
                                              ylim = yAxisLimits) +
                              theme(plot.margin = unit(c(1,1,1,3), "lines"),
                                    axis.text.x=element_text(size=10, face="bold", color="black"),
                                    axis.text.y=element_text(size=8, color="black"),
                                    axis.title.y=element_text(size=8, color="black"),
                                    axis.title.x=element_text(size=5, color="black"),
                                    legend.position="none" ) +
                              ylab("Mitochondria in mitophagy\nper cell (%)") + 
                              ggnewscale::new_scale_colour()+
                              ggnewscale::new_scale_fill() +
                              ggforce::geom_sina(data = DataNow[DataNow$SampleNow,],
                                                 aes(x = OverallLabeling, #OverallUniqueConditionTag
                                                     y = FractionMitophagy,
                                                     colour=FractionMitophagy
                                                 ),
                                                 scale = "count", # area"count" "width"
                                                 method = "density", # density # counts
                                                 #maxwidth,
                                                 position=position_dodge(width=0.8),
                                                 alpha = 0.5,
                                                 shape=16,
                                                 size=0.07) +
                              scale_colour_steps(low=rgb(100,100,100, 255,max=255),
                                                 high="black",
                                                 breaks=c(0,0.1) ) +
                              ggnewscale::new_scale_colour()+
                              ggnewscale::new_scale_fill()+
                              scale_colour_manual(values = c(rep("black", 100)),
                                                  aesthetics = c("colour")) +
                              geom_boxplot(position=position_dodge(width=0.8),
                                           width=0.5,
                                           fill="transparent",
                                           color=rgb(50,50,50, 200, max=255),
                                           outlier.shape = NA,
                                           coef=0) +
                              xlab(paste0(ActualRData_name, "\n",
                                          "MinL_", minProteinLabel,
                                          "minMaxD_",minMaxDProteinLabel,
                                          "MinMax_", minProteinMaxLabel,
                                          "MinMedian_", minProteinMedianLabel,
                                          "\n",
                                          "MinStdDev_", minProteinStdDevLabel,
                                          "MinMitophagy_", minMitophagyLabel,
                                          "_ratio_", RatioMitophagyCutoff,
                                          "\n",
                                          "minP_", minPerimNow, "maxP_", maxPerimNow, " - ", mitoTypeNow))
                            
                            # 
                            png(paste0(rootFolder, "/",
                                       "09-MitophagyTubular", "/",
                                       ActualRData_name, "/",
                                       "MinL_", minProteinLabel,
                                       "minMaxD_",minMaxDProteinLabel,
                                       "MinMax_", minProteinMaxLabel,
                                       "MinMedian_", minProteinMedianLabel,
                                       "MinStdDev_", minProteinStdDevLabel,
                                       "MinMitophagy_", minMitophagyLabel,
                                       "/",
                                       "03-", paste0(GACNow, collapse="_"),
                                       "_", paste0(CellNOW, collapse="_"),
                                       "_My_", paste0(MitophagyNow, collapse="_"),
                                       "_R_", RatioMitophagyCutoff,
                                       "P", minPerimNow*10,
                                       "-",maxPerimNow*10,
                                       "_", mitoTypeNow,
                                       ".png"),
                                width=(760+(187* length(unique(as.character(DataNow$OverallLabeling)))  )), height=1600, res=600)
                            print(
                              p+geom_text(color="black",
                                          x = 0.7,
                                          y = -5,
                                          label= MyLabelsNow,
                                          hjust = 1,
                                          vjust = 1,
                                          size = 2.7)
                            )
                            dev.off()
                            
                            remove(p)
                            
                            
                            
                            
                            
                            
                          }
                        }
                      }
                      
                      
                      ### Save table
                      ########
                      TempTable_plyr_toLabel$MitoCat <- mitoTypeNow
                      StoredByMitoType[[(length(StoredByMitoType)+1)]] <- TempTable_plyr_toLabel
                      
                      TempTable_plyr_toBar$MitoCat <- mitoTypeNow
                      StoredByMitoType_toBar[[(length(StoredByMitoType_toBar)+1)]] <- TempTable_plyr_toBar
                      
                      
                      
                      ### mito cat end
                    }
                    
                    
                    
                    
                    
                    
                    
                    
                    #### Create general plot
                    TempTable_plyr_toLabel <- bind_rows(StoredByMitoType)
                    StoredByMitoType_toBar <- bind_rows(StoredByMitoType_toBar)
                    
                    ## CREATE LISTS FOR THE PLOTS
                    
                    GACList <- list()
                    for(i in 1:length( unique(TempTable_plyr_toLabel$OverallInfo01)  ) ){
                      GACList[[i]] <- unique(TempTable_plyr_toLabel$OverallInfo01)[i]
                    }
                    #GACList[[length(GACList)+1]] <- unique(TempTable_plyr_toLabel$OverallInfo01)
                    #GACList[[length(GACList)]] <- GACList[[length(GACList)]][GACList[[length(GACList)]] != "01-Mock"]
                    # GACList[[1]] <- unique(TempTable_plyr_toLabel$OverallInfo01)
                    # GACList[[1]] <- GACList[[length(GACList)]][GACList[[length(GACList)]] != "01-Mock"]
                    

                    CellList <- list()
                    for(i in 1:length( unique(TempTable_plyr_toLabel$OverallCellLine)  ) ){
                      CellList[[i]] <- unique(TempTable_plyr_toLabel$OverallCellLine)[i]
                    }
                    # CellList[[length(CellList)+1]] <- unique(TempTable_plyr_toLabel$OverallCellLine)
                    
                    
                    
                    #
                    MitophagyList <- list()
                    for(i in 1:length( unique(TempTable_plyr_toLabel$OverallInfo03)  ) ){
                      MitophagyList[[i]] <- unique(TempTable_plyr_toLabel$OverallInfo03)[i]
                    }
                    # MitophagyList[[length(MitophagyList)+1]] <- unique(TempTable_plyr_toLabel$OverallInfo03)
                    
                    
                    colorValues_noalpha_all <- c(rgb(255,255,255,255,max=255),
                                                 rgb(255,120,180,255,max=255),
                                                 
                                                 rgb(255,255,255,255,max=255),
                                                 rgb(255,0,0,255,max=255),
                                                 
                                                 rgb(255,255,255,255,max=255),
                                                 rgb(0,210,140,255,max=255),
                                                 
                                                 rgb(255,255,255,255,max=255),
                                                 rgb(0,160,255,255,max=255))
                    
                    
                    
                    colorValues_border_noalpha_all <- c(colorValues_noalpha_all[2],
                                                        colorValues_noalpha_all[2],
                                                        
                                                        colorValues_noalpha_all[4],
                                                        colorValues_noalpha_all[4],
                                                        
                                                        colorValues_noalpha_all[6],
                                                        colorValues_noalpha_all[6],
                                                        
                                                        colorValues_noalpha_all[8],
                                                        colorValues_noalpha_all[8])
                    
                    
                    
                    
                    
                    
                    
                    
                    #### ADD LOOP HERE
                    #GACNow <- GACList[[5]]
                    #MitophagyNow <- MitophagyList[[3]]
                    for(CellNOW in CellList){
                      for(GACNow in GACList){
                        for(MitophagyNow in MitophagyList ){
                          
                          for(PositivityNow in unique(TempTable_plyr_toLabel$Classification) ){
                            
                            
                            print(paste0("JOINED MinL_", minProteinLabel,
                                         "minMaxD_",minMaxDProteinLabel,
                                         "MinMax_", minProteinMaxLabel,
                                         "MinMedian_", minProteinMedianLabel,
                                         "MinStdDev_", minProteinStdDevLabel,
                                         "MinMitophagy_", minMitophagyLabel,
                                         "/",
                                         "Groups_", paste0(GACNow, collapse="_"),
                                         "Cell_", paste0(CellNOW, collapse="_"),
                                         "_Mitophagy_", paste0(MitophagyNow, collapse="_"),
                                         "_ratio_", RatioMitophagyCutoff,
                                         "Pmin", minPerimNow*10,
                                         "Pmax_",maxPerimNow*10
                            ))
                            
                            
                            # Color selection
                            if(length(GACNow) == 4){
                              colorValues_noalpha_all[c(1,3,5,7,2,4,6,8)] -> colorValues_noalpha
                              colorValues_border_noalpha_all[c(1,3,5,7,2, 4,6,8)] -> colorValues_border_noalpha
                            } else {
                              if(length(GACNow) == 3){
                                colorValues_noalpha_all[c(3,5,7,4,6,8)] -> colorValues_noalpha
                                colorValues_border_noalpha_all[c(3,5,7,4,6,8)] -> colorValues_border_noalpha
                              } else {
                                if(GACNow == "02-GAC.wt"){
                                  colorValues_noalpha_all[3:4] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[3:4] -> colorValues_border_noalpha
                                }
                                if(GACNow == "03-GAC.K320A"){
                                  colorValues_noalpha_all[5:6] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[5:6] -> colorValues_border_noalpha
                                }
                                if(GACNow == "04-GAC.R382D"){
                                  colorValues_noalpha_all[7:8] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[7:8] -> colorValues_border_noalpha
                                }
                                if(GACNow == "01-Mock"){
                                  colorValues_noalpha_all[1:2] -> colorValues_noalpha
                                  colorValues_border_noalpha_all[1:2] -> colorValues_border_noalpha
                                }
                              }
                            }
                            
                            #####
                            
                            if(PositivityNow == "0-Negative"){
                              
                              colorValues_noalpha[1] -> colorValues_noalpha
                              colorValues_border_noalpha[1] -> colorValues_border_noalpha
                            } else {
                              colorValues_noalpha[2] -> colorValues_noalpha
                              colorValues_border_noalpha[2] -> colorValues_border_noalpha
                            }
                            
                            
                            
                            ####
                            TempTable_plyr_toLabel[which(
                              TempTable_plyr_toLabel$OverallInfo01 %in% GACNow &
                                TempTable_plyr_toLabel$OverallCellLine %in% CellNOW &
                                TempTable_plyr_toLabel$OverallInfo03 %in% MitophagyNow &
                                TempTable_plyr_toLabel$Classification %in% PositivityNow
                            ),] -> DataNow
                            
                            

                            # DEFINE AXIS LIMITS
                            yAxisLimits <- c(0, 
                                             floor(max(
                                               
                                               sapply(unique(DataNow$OverallUniqueConditionTag), function(w){
                                                 quantile(DataNow$FractionMitophagy[which(DataNow$OverallUniqueConditionTag == w)],
                                                          0.75, na.rm=TRUE)
                                               }), na.rm=TRUE)+10))
                            ####
                            
                            #### 
                            p <- ggplot(data = DataNow,
                                        aes(x = OverallLabeling, #OverallUniqueConditionTag
                                            y = FractionMitophagy,
                                            colour = interaction(OverallInfo01, MitoCat, sep=':'),
                                            fill= interaction(OverallInfo01, MitoCat, sep=':'))) + #, fill = Classification
                              geom_violin(trim=TRUE,
                                          position=position_dodge(width=0.8),
                                          lwd=0.35)+
                              # scale_fill_manual(values=colorValues_noalpha) +
                              scale_fill_manual(values = rep(colorValues_noalpha, 3),
                                                aesthetics = c("fill")) +
                              scale_colour_manual(values = rep(colorValues_border_noalpha, 3),
                                                  aesthetics = c("colour")) +
                              theme_classic() +
                              coord_cartesian(clip = "off",
                                              ylim = yAxisLimits) +
                              theme(plot.margin = unit(c(1,1,1,3), "lines"),
                                    axis.text.x=element_text(size=10, face="bold", color="black"),
                                    axis.text.y=element_text(size=8, color="black"),
                                    axis.title.y=element_text(size=8, color="black"),
                                    axis.title.x=element_text(size=5, color="black"),
                                    legend.position="none" ) +
                              ylab("Mitochondria in mitophagy\nper cell (%)") + 
                              ggnewscale::new_scale_colour()+
                              ggnewscale::new_scale_fill() +
                              ggforce::geom_sina(data = DataNow[DataNow$SampleNow,],
                                                 aes(x = OverallLabeling, #OverallUniqueConditionTag
                                                     y = FractionMitophagy,
                                                     colour=FractionMitophagy
                                                 ),
                                                 scale = "count", # area"count" "width"
                                                 method = "density", # density # counts
                                                 #maxwidth,
                                                 position=position_dodge(width=0.8),
                                                 alpha = 0.5,
                                                 shape=16,
                                                 size=0.07) +
                              scale_colour_steps(low=rgb(100,100,100, 255,max=255),
                                                 high="black",
                                                 breaks=c(0,0.1) ) +
                              ggnewscale::new_scale_colour()+
                              ggnewscale::new_scale_fill()+
                              scale_colour_manual(values = c(rep("black", 100)),
                                                  aesthetics = c("colour")) +
                              geom_boxplot(position=position_dodge(width=0.8),
                                           width=0.5,
                                           fill="transparent",
                                           color=rgb(50,50,50, 200, max=255),
                                           outlier.shape = NA,
                                           coef=0) +
                              xlab(paste0(ActualRData_name, "\n",
                                          "MinL_", minProteinLabel,
                                          "minMaxD_",minMaxDProteinLabel,
                                          "MinMax_", minProteinMaxLabel,
                                          "MinMedian_", minProteinMedianLabel,
                                          "\n",
                                          "MinStdDev_", minProteinStdDevLabel,
                                          "MinMitophagy_", minMitophagyLabel,
                                          "_ratio_", RatioMitophagyCutoff,
                                          "\n",
                                          "minP_", minPerimNow, "maxP_", maxPerimNow, " - ", PositivityNow))
                            
                            # 
                            png(paste0(rootFolder, "/",
                                       "09-MitophagyTubular", "/",
                                       ActualRData_name, "/",
                                       "MinL_", minProteinLabel,
                                       "minMaxD_",minMaxDProteinLabel,
                                       "MinMax_", minProteinMaxLabel,
                                       "MinMedian_", minProteinMedianLabel,
                                       "MinStdDev_", minProteinStdDevLabel,
                                       "MinMitophagy_", minMitophagyLabel,
                                       "/",
                                       "04-", paste0(GACNow, collapse="_"),
                                       "_", paste0(CellNOW, collapse="_"),
                                       "_My_", paste0(MitophagyNow, collapse="_"),
                                       "_R_", RatioMitophagyCutoff,
                                       "P", minPerimNow*10,
                                       "-",maxPerimNow*10,
                                       "_", PositivityNow,
                                       ".png"),
                                width=(760+(210* length(unique(as.character(DataNow$OverallLabeling)))  )), height=1600, res=600)
                            print(
                              p+geom_text(color="black",
                                          x = 0.7,
                                          y = -5,
                                          label= MyLabelsNow,
                                          hjust = 1,
                                          vjust = 1,
                                          size = 2.7)
                            )
                            dev.off()
                            
                            remove(p)
                            
                            
                            
                            
                            
                            
                            ########
                            
                            ####
                            StoredByMitoType_toBar[which(
                              StoredByMitoType_toBar$OverallInfo01 %in% GACNow &
                                StoredByMitoType_toBar$OverallCellLine %in% CellNOW &
                                StoredByMitoType_toBar$OverallInfo03 %in% MitophagyNow &
                                StoredByMitoType_toBar$Classification %in% PositivityNow
                            ),] -> DataNow
                            
                            
                            
                            df.summary2 <- DataNow %>%
                              group_by(OverallUniqueConditionTag,
                                       OverallLabeling,
                                       Classification,
                                       MitoCat,
                                       OverallInfo01,
                                       OverallUniqueConditionTag_Content) %>%
                              summarise(mean = mean(Fraction_mean, na.rm=TRUE),
                                        sd = sd(Fraction_mean, na.rm = TRUE),
                                        n = n())
                            df.summary2$sem <- (df.summary2$sd / sqrt(df.summary2$n))
                            
                            
                            
                            #
                            p <- ggplot(data = df.summary2,
                                        aes(x=OverallLabeling,
                                            y=mean,
                                            colour = interaction(OverallInfo01, MitoCat, sep=':'),
                                            fill = interaction(OverallInfo01, MitoCat, sep=':'))) +
                              geom_bar(stat = "identity",
                                       position=position_dodge()) +
                              geom_errorbar(
                                aes(ymin = mean-sem,
                                    ymax = mean+sem),
                                colour="black",
                                width = 0.5,
                                position = position_dodge(.9)) +
                              geom_jitter(data=DataNow,
                                          aes(x=OverallLabeling,
                                              y=Fraction_mean,
                                              fill = interaction(OverallInfo01, MitoCat, sep=':'),
                                              colour = interaction(OverallInfo01, MitoCat, sep=':')),
                                          colour="black",
                                          alpha=0.3,
                                          position = position_jitterdodge(0.1)) +
                              scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                              scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                              theme_classic() +
                              coord_cartesian(clip = "off") + #, ylim=c(0, 100)
                              theme(plot.margin = unit(c(1,1,1,2), "lines"),
                                    axis.text.x=element_text(angle = 0, size=10, face="bold", color="black"),
                                    axis.text.y=element_text(size=8, color="black"),
                                    axis.title.y=element_text(size=8, color="black"),
                                    axis.title.x=element_text(size=5, color="black"),
                                    legend.position="top" ) +
                              ylab(paste0("Mean cell ", strsplit(PositivityNow, "-")[[1]][2],"\nmitochondria\nnetwork category (%)")) +
                              xlab(paste0(ActualRData_name, "\n",
                                          "MinL_", minProteinLabel,
                                          "minMaxD_",minMaxDProteinLabel,
                                          "MinMax_", minProteinMaxLabel,
                                          "MinMedian_", minProteinMedianLabel,
                                          "\n",
                                          "MinStdDev_", minProteinStdDevLabel,
                                          "MinMitophagy_", minMitophagyLabel,
                                          "_ratio_", RatioMitophagyCutoff,
                                          "\n",
                                          "minP_", minPerimNow, "maxP_", maxPerimNow, " - ", PositivityNow))
                            
                            
                            
                            
                            # (2) Bar plots of means + individual jitter points + errors
                            png(file=paste0(rootFolder, "/",
                                            "09-MitophagyTubular", "/",
                                            ActualRData_name, "/",
                                            "MinL_", minProteinLabel,
                                            "minMaxD_",minMaxDProteinLabel,
                                            "MinMax_", minProteinMaxLabel,
                                            "MinMedian_", minProteinMedianLabel,
                                            "MinStdDev_", minProteinStdDevLabel,
                                            "MinMitophagy_", minMitophagyLabel,
                                            "/",
                                            "05-", paste0(GACNow, collapse="_"),
                                            "_", paste0(CellNOW, collapse="_"),
                                            "_My_", paste0(MitophagyNow, collapse="_"),
                                            "_R_", RatioMitophagyCutoff,
                                            "P", minPerimNow*10,
                                            "-",maxPerimNow*10,
                                            "_", PositivityNow,
                                            ".png"),
                                width = (630+(length(unique(df.summary2$OverallLabeling))*340)), #4000
                                height = 2600,
                                res=600)
                            
                            print(
                              p + geom_text(color="black",
                                            x = 0.4,
                                            y = -2,
                                            label= MyLabelsNow,
                                            hjust = 1,
                                            vjust = 1,
                                            size = 2.7)
                            )
                            dev.off()
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            df.summary2 <- DataNow %>%
                              group_by(OverallUniqueConditionTag,
                                       OverallLabeling,
                                       Classification,
                                       MitoCat,
                                       OverallInfo01,
                                       OverallUniqueConditionTag_Content) %>%
                              summarise(mean = mean(Fraction_median, na.rm=TRUE),
                                        sd = sd(Fraction_median, na.rm = TRUE),
                                        n = n())
                            df.summary2$sem <- (df.summary2$sd / sqrt(df.summary2$n))
                            
                            
                            
                            #
                            p <- ggplot(data = df.summary2,
                                        aes(x=OverallLabeling,
                                            y=mean,
                                            colour = interaction(OverallInfo01, MitoCat, sep=':'),
                                            fill = interaction(OverallInfo01, MitoCat, sep=':'))) +
                              geom_bar(stat = "identity",
                                       position=position_dodge()) +
                              geom_errorbar(
                                aes(ymin = mean-sem,
                                    ymax = mean+sem),
                                colour="black",
                                width = 0.5,
                                position = position_dodge(.9)) +
                              geom_jitter(data=DataNow,
                                          aes(x=OverallLabeling,
                                              y=Fraction_median,
                                              fill = interaction(OverallInfo01, MitoCat, sep=':'),
                                              colour = interaction(OverallInfo01, MitoCat, sep=':')),
                                          colour="black",
                                          alpha=0.3,
                                          position = position_jitterdodge(0.1)) +
                              scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                              scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                              theme_classic() +
                              coord_cartesian(clip = "off") + #, ylim=c(0, 100)
                              theme(plot.margin = unit(c(1,1,1,2), "lines"),
                                    axis.text.x=element_text(angle = 0, size=10, face="bold", color="black"),
                                    axis.text.y=element_text(size=8, color="black"),
                                    axis.title.y=element_text(size=8, color="black"),
                                    axis.title.x=element_text(size=5, color="black"),
                                    legend.position="top" ) +
                              ylab(paste0("Median cell ", strsplit(PositivityNow, "-")[[1]][2],"\nmitochondria\nnetwork category (%)")) +
                              xlab(paste0(ActualRData_name, "\n",
                                          "MinL_", minProteinLabel,
                                          "minMaxD_",minMaxDProteinLabel,
                                          "MinMax_", minProteinMaxLabel,
                                          "MinMedian_", minProteinMedianLabel,
                                          "\n",
                                          "MinStdDev_", minProteinStdDevLabel,
                                          "MinMitophagy_", minMitophagyLabel,
                                          "_ratio_", RatioMitophagyCutoff,
                                          "\n",
                                          "minP_", minPerimNow, "maxP_", maxPerimNow, " - ", PositivityNow))
                            
                            
                            
                            
                            # (2) Bar plots of means + individual jitter points + errors
                            png(file=paste0(rootFolder, "/",
                                            "09-MitophagyTubular", "/",
                                            ActualRData_name, "/",
                                            "MinL_", minProteinLabel,
                                            "minMaxD_",minMaxDProteinLabel,
                                            "MinMax_", minProteinMaxLabel,
                                            "MinMedian_", minProteinMedianLabel,
                                            "MinStdDev_", minProteinStdDevLabel,
                                            "MinMitophagy_", minMitophagyLabel,
                                            "/",
                                            "06-", paste0(GACNow, collapse="_"),
                                            "_", paste0(CellNOW, collapse="_"),
                                            "_My_", paste0(MitophagyNow, collapse="_"),
                                            "_R_", RatioMitophagyCutoff,
                                            "P", minPerimNow*10,
                                            "-",maxPerimNow*10,
                                            "_", PositivityNow,
                                            ".png"),
                                width = (630+(length(unique(df.summary2$OverallLabeling))*340)), #4000
                                height = 2600,
                                res=600)
                            
                            print(
                              p + geom_text(color="black",
                                            x = 0.4,
                                            y = -2,
                                            label= MyLabelsNow,
                                            hjust = 1,
                                            vjust = 1,
                                            size = 2.7)
                            )
                            dev.off()
                            
                            
                            
                            
                          }

                        }
                      }
                    }
                    
                    
                    
                    
                  }
                }
                
                
                
                
                
                
              
              }
            }
            
            
          }}}}}
  
  remove(RawData_size_mitol_list)
  remove(RawData_size_mitol_mitosig_list)
  remove(CellList)
  remove(TempTable_plyr_toLabel)
  remove(TemporaryData)
  remove(DataNow)
  

  
  
  
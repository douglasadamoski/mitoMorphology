#!/usr/bin/env /mnt/nfs/home/douglas/anaconda3/envs/R413/bin/Rscript --vanilla
#
#  02_HPC_Morphology_v1.R
#
#
#

#
args = commandArgs(trailingOnly=TRUE)

# Available folders
#ActualRData <- "A99.29_ectopica_parkin_DRP_minMitoLabel_300.RData"

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
dir.create(path = paste0(rootFolder, "/", "04-Morphology_Sandra"), showWarnings = FALSE)

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
dir.create(path = paste0(rootFolder, "/", "04-Morphology_Sandra", "/", ActualRData_name, "/"),
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
minProteinLabel <- 0
minProteinMaxLabel <- 1000
minProteinMedianLabel <- 0
minProteinStdDevLabel <- 0

# RawData_size_mitol_list 
for(minProteinStdDevLabel in c(0)){
  for(minProteinMedianLabel in c(0)){

      for(minProteinMaxLabel in c(0)){ #, 1000, 2000, 2500
        for(minProteinLabel in c(500, 700, 1000, 1500, 2000, 2500, 3000)){ #0, 300, 500, 
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
                                     "04-Morphology_Sandra", "/",
                                     ActualRData_name, "/",
                                     "MinL_", minProteinLabel,
                                     "minMaxD_",minMaxDProteinLabel,
                                     "MinMax_", minProteinMaxLabel,
                                     "MinMedian_", minProteinMedianLabel,
                                     "MinStdDev_", minProteinStdDevLabel,
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
                          "MinStdDev_", minProteinStdDevLabel))
            
            
            
            
            
            #
            png(paste0(rootFolder, "/",
                       "04-Morphology_Sandra", "/",
                       ActualRData_name, "/",
                       "MinL_", minProteinLabel,
                       "minMaxD_",minMaxDProteinLabel,
                       "MinMax_", minProteinMaxLabel,
                       "MinMedian_", minProteinMedianLabel,
                       "MinStdDev_", minProteinStdDevLabel,
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
                          "MinStdDev_", minProteinStdDevLabel))
            
            
            png(paste0(rootFolder, "/",
                       "04-Morphology_Sandra", "/",
                       ActualRData_name, "/",
                       "MinL_", minProteinLabel,
                       "minMaxD_",minMaxDProteinLabel,
                       "MinMax_", minProteinMaxLabel,
                       "MinMedian_", minProteinMedianLabel,
                       "MinStdDev_", minProteinStdDevLabel,
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
                                   "04-Morphology_Sandra", "/",
                                   ActualRData_name, "/",
                                   "MinL_", minProteinLabel,
                                   "minMaxD_",minMaxDProteinLabel,
                                   "MinMax_", minProteinMaxLabel,
                                   "MinMedian_", minProteinMedianLabel,
                                   "MinStdDev_", minProteinStdDevLabel,
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
                                   "04-Morphology_Sandra", "/",
                                   ActualRData_name, "/",
                                   "MinL_", minProteinLabel,
                                   "minMaxD_",minMaxDProteinLabel,
                                   "MinMax_", minProteinMaxLabel,
                                   "MinMedian_", minProteinMedianLabel,
                                   "MinStdDev_", minProteinStdDevLabel,
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
            
            
            
            
            
            
            
            
            
            # PERIMETER DEFINITION
            minPerimNow <- 4.3
            maxPerimNow <- 7.3
            
            # D DEFINITION
            minDNow <- 0.6
            maxDNow <- 0.9
            
            
            
            print(paste0("Iniciando loops de perimetro para esse corte intensidade"))
            
            ### ONLY PERIMETER
            for(minPerimNow in c(3.7, 4.0, 4.2, 4.3)   ){ #seq(4,6,by=0.1)
              for(maxPerimNow in c(6.5, 7.0, 7.1, 7.3) ){

                
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
                
                
                #
                TemporaryData <- (data.frame(OverallFileOrigin_base,
                                             OverallUniqueTag,
                                             OverallUniqueConditionTag,
                                             OverallMitoCat_Perim,
                                             OverallCellCat_Perim,
                                             
                                             
                                             OverallCellCat_Perim_negative,
                                             
                                             OverallCellLine,
                                             OverallInfo01,
                                             OverallInfo02,
                                             OverallInfo03,
                                             OverallInfo04,
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
                remove(OverallCellLine)
                remove(OverallInfo01)
                remove(OverallInfo02)
                remove(OverallInfo03)
                remove(OverallInfo04)
                remove(OverallLabeling)
                
                print(paste0("Salvando tabela")) 
                
                write.csv2(TemporaryData,
                           file=paste0(rootFolder, "/",
                                       "04-Morphology_Sandra", "/",
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
                
                
                MitophagyList <- list()
                for(i in 1:length( unique(TemporaryData$OverallInfo03)  ) ){
                  MitophagyList[[i]] <- unique(TemporaryData$OverallInfo03)[i]
                }
                
                
                TemporaryData -> BACKUP
                TemporaryData <- BACKUP
                
                #####
                for(MitophagyNow in MitophagyList ){
                  
                  ############## BY CELL - PERIMETER
                  print(paste0("Plotando positivo"))
                  
                  #
                  
                  
                  #
                  if("02-GAC.wt" %in% unique(TemporaryData$OverallInfo01) ){
                    
                    # "006010_04-MEF_DRP_03-GAC.K320A_01-withGLN_01-Ctl_01-DMSO" -> filterNow
                    
                    TemporaryData$TagToKeep <- TRUE
                    
                    
                   # TemporaryData
                   for(filterNow in c("006010_04-MEF_DRP_03-GAC.K320A_01-withGLN_01-Ctl_01-DMSO",
                                      "007010_04-MEF_DRP_03-GAC.K320A_01-withGLN_01-Ctl_02-CB839",
                                      "007011_04-MEF_DRP_03-GAC.K320A_01-withGLN_01-Ctl_02-CB839",
                                      "007012_04-MEF_DRP_03-GAC.K320A_01-withGLN_01-Ctl_02-CB839",
                                      "011010_04-MEF_DRP_03-GAC.K320A_02-woGLN_01-Ctl_02-CB839",
                                      "011011_04-MEF_DRP_03-GAC.K320A_02-woGLN_01-Ctl_02-CB839",
                                      "011012_04-MEF_DRP_03-GAC.K320A_02-woGLN_01-Ctl_02-CB839")){
                      
                     
                      #
                      MyTableNow <- TemporaryData[which(TemporaryData$OverallUniqueTag == filterNow &
                                                          TemporaryData$OverallCellCat_Perim == "3-Tubular"), ]
                      
                      
                      MyTableNow[(MyTableNow$OverallFileOrigin_base %in% names(table(MyTableNow[,c("OverallFileOrigin_base")]))[unname(table(MyTableNow[,c("OverallFileOrigin_base")]) > 60)]  ), "TagToKeep"] <- FALSE
                      
                      #
                      MyTableNow -> TemporaryData[which(TemporaryData$OverallUniqueTag == filterNow &
                                                          TemporaryData$OverallCellCat_Perim == "3-Tubular"), ]
                    
                    }
                    
                    
                    ####
                    for(filterNow in c(
                      #"010010_04-MEF_DRP_03-GAC.K320A_02-woGLN_01-Ctl_01-DMSO",
                      "010011_04-MEF_DRP_03-GAC.K320A_02-woGLN_01-Ctl_01-DMSO",
                      "010012_04-MEF_DRP_03-GAC.K320A_02-woGLN_01-Ctl_01-DMSO"
                    ) ){
                      
                      
                      #
                      MyTableNow <- TemporaryData[which(TemporaryData$OverallUniqueTag == filterNow &
                                                          TemporaryData$OverallCellCat_Perim == "1-Fragmented"), ]
                      
                      
                      MyTableNow[(MyTableNow$OverallFileOrigin_base %in% names(table(MyTableNow[,c("OverallFileOrigin_base")]))[unname(table(MyTableNow[,c("OverallFileOrigin_base")]) > 25 )]  ), "TagToKeep"] <- FALSE
                      
                      #
                      MyTableNow -> TemporaryData[which(TemporaryData$OverallUniqueTag == filterNow &
                                                          TemporaryData$OverallCellCat_Perim == "1-Fragmented"), ]
                      
                      
                      
                      #
                      MyTableNow <- TemporaryData[which(TemporaryData$OverallUniqueTag == filterNow &
                                                          TemporaryData$OverallCellCat_Perim == "2-Intermediate"), ]
                      
                      
                      
                      MyTableNow[(MyTableNow$OverallFileOrigin_base %in%   names(table(MyTableNow[,c("OverallFileOrigin_base")]))[unname(table(MyTableNow[,c("OverallFileOrigin_base")]) > 50)] ), "TagToKeep"] <- FALSE
                      
                      #
                      MyTableNow -> TemporaryData[which(TemporaryData$OverallUniqueTag == filterNow &
                                                          TemporaryData$OverallCellCat_Perim == "2-Intermediate"), ]
                      
                      
                      }
                    
                    
                    
                    #
                    TemporaryData <- TemporaryData[which(TemporaryData$TagToKeep),]
                    TemporaryData$TagToKeep <- NULL
                    
                    
                    
                    df.summary <- unique(TemporaryData[which(TemporaryData$OverallInfo01 != "01-Mock" &
                                                               TemporaryData$OverallInfo03 %in% MitophagyNow ), c("OverallUniqueTag",
                                                                                                                  "OverallLabeling",
                                                                                                                  "OverallUniqueConditionTag",
                                                                                                                  "OverallCellCat_Perim",
                                                                                                                  "OverallFileOrigin_base")]) %>%
                      group_by(OverallUniqueTag, OverallUniqueConditionTag, OverallLabeling) %>%
                      count(OverallCellCat_Perim)
                    


                  } else {
                    
                    
                    
                    
                    df.summary <- unique(TemporaryData[which( TemporaryData$OverallInfo03 %in% MitophagyNow ),c("OverallUniqueTag",
                                                                                                                "OverallLabeling",
                                                                                                                "OverallUniqueConditionTag",
                                                                                                                "OverallCellCat_Perim",
                                                                                                                "OverallFileOrigin_base")]) %>%
                      group_by(OverallUniqueTag, OverallUniqueConditionTag, OverallLabeling) %>%
                      count(OverallCellCat_Perim)
                    
                  }
                  
                  df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
                  df.summary$Percentage <- NA
                  MyCategoryNow <- "OverallCellCat_Perim"
                  df.summary <- df.summary[!is.na(df.summary[,MyCategoryNow]),]
                  
                  
                  # Add zeroes
                  for(condNow in unique(df.summary$OverallUniqueTag)){
                    if("1-Fragmented" %in% df.summary[df.summary$OverallUniqueTag==condNow, MyCategoryNow]){
                      #
                    } else {
                      MyLine <- df.summary[df.summary$OverallUniqueTag==condNow, ]
                      MyLine[,"n"] <- 0
                      MyLine[,MyCategoryNow] <- "1-Fragmented"
                      df.summary <- rbind(df.summary, MyLine[1,])
                    }
                    
                    if("2-Intermediate" %in% df.summary[df.summary$OverallUniqueTag==condNow, MyCategoryNow]){
                      #print("tem")
                    } else {
                      MyLine <- df.summary[df.summary$OverallUniqueTag==condNow, ]
                      MyLine[, "n"] <- 0
                      MyLine[, MyCategoryNow] <- "2-Intermediate"
                      df.summary <- rbind(df.summary, MyLine[1,])
                    }
                    
                    if("3-Tubular" %in% df.summary[df.summary$OverallUniqueTag==condNow, MyCategoryNow]){
                      #print("tem")
                    } else {
                      MyLine <- df.summary[df.summary$OverallUniqueTag==condNow, ]
                      MyLine[,"n"] <- 0
                      MyLine[, MyCategoryNow] <- "3-Tubular"
                      df.summary <- rbind(df.summary, MyLine[1,])
                    }
                    
                    # Add percentage
                    df.summary[df.summary$OverallUniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$OverallUniqueTag==condNow, "n"]/
                      sum(df.summary[df.summary$OverallUniqueTag==condNow, "n"])
                    
                  }
                  
                  #
                  df.summary[,MyCategoryNow] -> df.summary[,"Classification"] 
                  

                  
                  
                  df.summary2 <- df.summary %>%
                    group_by(OverallUniqueConditionTag, OverallLabeling, Classification) %>%
                    summarise(mean = mean(Percentage, na.rm=TRUE),
                              sd = sd(Percentage, na.rm = TRUE))
                  
                  
                  write.csv2(df.summary,
                             file=paste0(rootFolder, "/",
                                         "04-Morphology_Sandra", "/",
                                         ActualRData_name, "/",
                                         "MinL_", minProteinLabel,
                                         "minMaxD_",minMaxDProteinLabel,
                                         "MinMax_", minProteinMaxLabel,
                                         "MinMedian_", minProteinMedianLabel,
                                         "MinStdDev_", minProteinStdDevLabel,
                                         "/", "04_Table_Pos_",paste0(MitophagyNow, collapse="_"),"MinP_", minPerimNow,
                                         "_MaxP_", maxPerimNow,
                                         "_pontos",
                                         ".csv"))
                  write.csv2(df.summary2,
                             file=paste0(rootFolder, "/",
                                         "04-Morphology_Sandra", "/",
                                         ActualRData_name, "/",
                                         "MinL_", minProteinLabel,
                                         "minMaxD_",minMaxDProteinLabel,
                                         "MinMax_", minProteinMaxLabel,
                                         "MinMedian_", minProteinMedianLabel,
                                         "MinStdDev_", minProteinStdDevLabel,
                                         "/", "05_Table_Pos_",paste0(MitophagyNow, collapse="_"),"MinP_", minPerimNow,
                                         "_MaxP_", maxPerimNow,
                                         "_medias",
                                         ".csv"))
                  
                  #
                  p <- ggplot(data = df.summary2,
                              aes(x=OverallLabeling,
                                  y=mean,
                                  colour = Classification,
                                  fill = Classification)) +
                    geom_bar(stat = "identity",
                             position=position_dodge()) +
                    geom_errorbar(
                      aes(ymin = mean-sd,
                          ymax = mean+sd),
                      colour="black",
                      width = 0.5,
                      position = position_dodge(.9)) +
                    geom_jitter(data=df.summary,
                                aes(x=OverallLabeling,
                                    y=Percentage,
                                    fill = Classification,
                                    colour = Classification),
                                colour="black",
                                alpha=0.3,
                                position = position_jitterdodge(0.1)) +
                    scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                    scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                    theme_classic() +
                    coord_cartesian(clip = "off", ylim=c(0, 100)) +
                    theme(plot.margin = unit(c(1,1,1,2), "lines"),
                          axis.text.x=element_text(angle = 0, size=10, face="bold", color="black"),
                          axis.text.y=element_text(size=8, color="black"),
                          axis.title.y=element_text(size=8, color="black"),
                          axis.title.x=element_text(size=5, color="black"),
                          legend.position="top" ) +
                    ylab("Cell positive mitochondria\nnetwork category (%)") +
                    xlab(paste0(ActualRData_name, "\n",
                                "MinL_", minProteinLabel,
                                "minMaxD_",minMaxDProteinLabel,
                                "MinMax_", minProteinMaxLabel,
                                "MinMedian_", minProteinMedianLabel,
                                "\n",
                                "MinStdDev_", minProteinStdDevLabel,
                                "MinP_", minPerimNow,
                                "MaxP_", maxPerimNow))
                  
                  
                  
                  # (2) Bar plots of means + individual jitter points + errors
                  png(file=paste0(rootFolder, "/",
                                  "04-Morphology_Sandra", "/",
                                  ActualRData_name, "/",
                                  "MinL_", minProteinLabel,
                                  "minMaxD_",minMaxDProteinLabel,
                                  "MinMax_", minProteinMaxLabel,
                                  "MinMedian_", minProteinMedianLabel,
                                  "MinStdDev_", minProteinStdDevLabel,
                                  "/", "06_PositiveMitochondria_",paste0(MitophagyNow, collapse="_"),"_MinP_", minPerimNow*10,
                                  "_MaxP_", maxPerimNow*10,
                                  ".png"),
                      width = (630+(length(unique(df.summary2$OverallLabeling))*210)), #4000
                      height = 2600,
                      res=600)
                  
                  print(
                    p + geom_text(color="black",
                                  x = 0.4,
                                  y = -9,
                                  label= MyLabelsNow,
                                  hjust = 1,
                                  vjust = 1,
                                  size = 2.7)
                  )
                  dev.off()
                  
                  
                  
                  
                  
                  ############ NEGATIVE BY CELL - PERIMETER
                  
                  print(paste0("Plotando negativo"))
                  
                  
                  
                  df.summary <- unique(TemporaryData[ which( TemporaryData$OverallInfo03 %in% MitophagyNow ), c("OverallUniqueTag",
                                                                                                                "OverallLabeling",
                                                                                                                "OverallUniqueConditionTag",
                                                                                                                "OverallCellCat_Perim_negative",
                                                                                                                "OverallFileOrigin_base")]) %>%
                    group_by(OverallUniqueTag, OverallUniqueConditionTag, OverallLabeling) %>%
                    count(OverallCellCat_Perim_negative)
                  
                  df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
                  
                  df.summary$Percentage <- NA
                  
                  MyCategoryNow <- "OverallCellCat_Perim_negative"
                  
                  df.summary <- df.summary[!is.na(df.summary[,MyCategoryNow]),]
                  
                  # Add zeroes
                  for(condNow in unique(df.summary$OverallUniqueTag)){
                    
                    if("1-Fragmented" %in% df.summary[df.summary$OverallUniqueTag==condNow, MyCategoryNow]){
                      #
                    } else {
                      MyLine <- df.summary[df.summary$OverallUniqueTag==condNow, ]
                      MyLine[,"n"] <- 0
                      MyLine[,MyCategoryNow] <- "1-Fragmented"
                      df.summary <- rbind(df.summary, MyLine[1,])
                    }
                    
                    if("2-Intermediate" %in% df.summary[df.summary$OverallUniqueTag==condNow, MyCategoryNow]){
                      #print("tem")
                    } else {
                      MyLine <- df.summary[df.summary$OverallUniqueTag==condNow, ]
                      MyLine[, "n"] <- 0
                      MyLine[, MyCategoryNow] <- "2-Intermediate"
                      df.summary <- rbind(df.summary, MyLine[1,])
                    }
                    
                    if("3-Tubular" %in% df.summary[df.summary$OverallUniqueTag==condNow, MyCategoryNow]){
                      #print("tem")
                    } else {
                      MyLine <- df.summary[df.summary$OverallUniqueTag==condNow, ]
                      MyLine[,"n"] <- 0
                      MyLine[, MyCategoryNow] <- "3-Tubular"
                      df.summary <- rbind(df.summary, MyLine[1,])
                    }
                    
                    # Add percentage
                    df.summary[df.summary$OverallUniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$OverallUniqueTag==condNow, "n"]/
                      sum(df.summary[df.summary$OverallUniqueTag==condNow, "n"])
                    
                  }
                  
                  #
                  df.summary[,MyCategoryNow] -> df.summary[,"Classification"] 
                  df.summary2 <- df.summary %>%
                    group_by(OverallUniqueConditionTag, OverallLabeling, Classification) %>%
                    summarise(mean = mean(Percentage, na.rm=TRUE),
                              sd = sd(Percentage, na.rm = TRUE))
                  
                  #
                  write.csv2(df.summary,
                             file=paste0(rootFolder, "/",
                                         "04-Morphology_Sandra", "/",
                                         ActualRData_name, "/",
                                         "MinL_", minProteinLabel,
                                         "minMaxD_",minMaxDProteinLabel,
                                         "MinMax_", minProteinMaxLabel,
                                         "MinMedian_", minProteinMedianLabel,
                                         "MinStdDev_", minProteinStdDevLabel,
                                         "/", "04_Table_Neg_",paste0(MitophagyNow, collapse="_"),"MinP_", minPerimNow,
                                         "_MaxP_", maxPerimNow,
                                         "_pontos",
                                         ".csv"))
                  write.csv2(df.summary2,
                             file=paste0(rootFolder, "/",
                                         "04-Morphology_Sandra", "/",
                                         ActualRData_name, "/",
                                         "MinL_", minProteinLabel,
                                         "minMaxD_",minMaxDProteinLabel,
                                         "MinMax_", minProteinMaxLabel,
                                         "MinMedian_", minProteinMedianLabel,
                                         "MinStdDev_", minProteinStdDevLabel,
                                         "/", "05_Table_Neg_",paste0(MitophagyNow, collapse="_"),"MinP_", minPerimNow,
                                         "_MaxP_", maxPerimNow,
                                         "_medias",
                                         ".csv"))
                  
                  
                  p <- ggplot(data = df.summary2,
                              aes(x=OverallLabeling,
                                  y=mean,
                                  colour = Classification,
                                  fill = Classification)) +
                    
                    geom_bar(stat = "identity",
                             position=position_dodge()) +
                    geom_errorbar(
                      aes(ymin = mean-sd,
                          ymax = mean+sd),
                      colour="black",
                      width = 0.5,
                      position = position_dodge(.9)) +
                    geom_jitter(data=df.summary,
                                aes(x=OverallLabeling,
                                    y=Percentage,
                                    fill = Classification,
                                    colour = Classification),
                                colour="black",
                                alpha=0.3,
                                position = position_jitterdodge(0.1)) +
                    scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                    scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
                    theme_classic() +
                    coord_cartesian(clip = "off", ylim=c(0, 100)) +
                    theme(plot.margin = unit(c(1,1,1,2), "lines"),
                          axis.text.x=element_text(angle = 0, size=10, face="bold", color="black"),
                          axis.text.y=element_text(size=8, color="black"),
                          axis.title.y=element_text(size=8, color="black"),
                          axis.title.x=element_text(size=5, color="black"),
                          legend.position="top" ) +
                    ylab("Cell negative mitochondria\nnetwork category (%)") +
                    xlab(paste0(ActualRData_name, "\n",
                                "MinL_", minProteinLabel,
                                "minMaxD_",minMaxDProteinLabel,
                                "MinMax_", minProteinMaxLabel,
                                "MinMedian_", minProteinMedianLabel,
                                "\n",
                                "MinStdDev_", minProteinStdDevLabel,
                                "MinP_", minPerimNow,
                                "MaxP_", maxPerimNow))
                  # (2) Bar plots of means + individual jitter points + errors
                  png(file=paste0(rootFolder, "/",
                                  "04-Morphology_Sandra", "/",
                                  ActualRData_name, "/",
                                  "MinL_", minProteinLabel,
                                  "minMaxD_",minMaxDProteinLabel,
                                  "MinMax_", minProteinMaxLabel,
                                  "MinMedian_", minProteinMedianLabel,
                                  "MinStdDev_", minProteinStdDevLabel,
                                  "/", "07_NegativeMitochondria_",paste0(MitophagyNow, collapse="_"),"MinP_", minPerimNow*10,
                                  "_MaxP_", maxPerimNow*10,
                                  ".png"),
                      width = (630+(length(unique(df.summary2$OverallLabeling))*210)), #4000
                      height = 2600,
                      res=600)
                  
                  print(
                    p + geom_text(color="black",
                                  x = 0.4,
                                  y = -9,
                                  label= MyLabelsNow,
                                  hjust = 1,
                                  vjust = 1,
                                  size = 2.7)
                  )
                  dev.off()
                  
                }
               
                
                
                
                #Cleaning 
                remove(df.summary)
                remove(df.summary2)
                remove(TemporaryData)
                remove(condNow)
                remove(MyCategoryNow)
                remove(MyLine)
                remove(p)
                #
                
                
                
                ##########
                ########## END
                ##########
                ########## 06 - MORPHOLOGY PLOTS
                ##########
                
                
              }
              }
            
          }}}}}

  
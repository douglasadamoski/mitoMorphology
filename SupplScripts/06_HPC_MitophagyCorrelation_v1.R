#!/usr/bin/env /mnt/nfs/home/douglas/anaconda3/envs/R413/bin/Rscript --vanilla
#
#  06_HPC_MitophagyCorrelation_all2_v1.R
#
#
#

#
args = commandArgs(trailingOnly=TRUE)

# Available folders
#ActualRData <- "A99.28_ectopica_parkin_controle_minMitoLabel_300.RData"

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
library(ggpmisc)

# If loop to run everything at once
# utils::memory.limit(size=32000)

# Create main data
dir.create(path = paste0(rootFolder, "/", "08-MitophagyCorrelation_all2"), showWarnings = FALSE)

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
dir.create(path = paste0(rootFolder, "/", "08-MitophagyCorrelation_all2", "/", ActualRData_name, "/"),
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
minProteinLabel <- 500
minProteinMaxLabel <- 0
minProteinMedianLabel <- 0
minProteinStdDevLabel <- 0




# RawData_size_mitol_list 
for(minProteinStdDevLabel in c(0)){
  for(minProteinMedianLabel in c(0)){
    for(minMitophagyLabel in c(100, 500)){
      for(minProteinMaxLabel in c(0)){
        for(minProteinLabel in c(300, 500, 700, 1000, 1500, 2000, 2500)){
          #
          
          #Do no all zero
          if((sum(c(minProteinStdDevLabel,
                    minProteinMedianLabel,
                    minProteinMaxLabel,
                    minProteinLabel)) >0)){
            
            
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
                                     "08-MitophagyCorrelation_all2", "/",
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

            RawData_size_mitophagy_list <- RawData_size_mitol_mitosig_list
            remove(RawData_size_mitol_mitosig_list)
            gc()
            #

            #
            print("Coletando colunas")
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
            
            
            OverallProtein_Mean <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                     FUN=ColumnRetrieval,
                                                     columnNow="Protein_Mean")))
            
            OverallProtein_Max <- unname(unlist(pblapply(X=RawData_size_mitophagy_list,
                                                     FUN=ColumnRetrieval,
                                                     columnNow="Protein_Max")))
            
            
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
                                         OverallRatioMitophagy_Cellwise,
                                         OverallPositiveMitophagySensor,
                                         OverallPositiveCell,
                                         OverallPositiveMitochondria,
                                         OverallRatioMitophagy_Cellwise_MinMean,
                                         OverallProtein_Mean,
                                         OverallProtein_Max
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
            remove(OverallProtein_Mean)
            remove(OverallProtein_Max)
            
            # Block positive mitochondria in NEGATIVE cells
            # V0.3 addition
            TemporaryData[which(!TemporaryData$OverallPositiveCell), "OverallPositiveMitochondria"] <- FALSE
            
            
            # Remove a lista pra liberar memória
            remove(RawData_size_mitophagy_list)
            gc()
            gc()
            
            print("Ordenando a Temporary colunas")
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
            
            
            ###
            list()
            
            for(MitophagyNow in  unique(TemporaryData$OverallInfo03) ){
              for(CompoundNow in  unique(TemporaryData$OverallInfo04) ){

                myData <- TemporaryData[which(
                  TemporaryData$OverallInfo01 %in% c("02-GAC.wt", "03-GAC.K320A", "04-GAC.R382D") &
                    TemporaryData$OverallInfo02 %in% c("01-withGLN", "02-woGLN") &
                    TemporaryData$OverallInfo03 %in% MitophagyNow &
                    TemporaryData$OverallInfo04 %in% CompoundNow &
                    TemporaryData$OverallPositiveMitochondria
                ),]
                
                myData$Sample <- FALSE
                
                #
                myData <- myData[which(
                  !is.na(myData$OverallInfo01) &
                    !is.na(myData$OverallInfo02)
                ),]
                
                
                for(OverallInfo01NOW in c("02-GAC.wt", "03-GAC.K320A", "04-GAC.R382D") ){
                  for(OverallInfo02NOW in c("01-withGLN", "02-woGLN") ){
                    #
                    myData[sample(x = which(
                      myData$OverallInfo01 == OverallInfo01NOW &
                        myData$OverallInfo02 == OverallInfo02NOW &
                        myData$OverallPositiveMitochondria),
                      size=10000,
                      replace=FALSE), "Sample" ] <- TRUE
                  }
                }
                
                #
                myData <- myData[which(myData$Sample), ]
                
                
                #
                # 
                
                p <- ggplot(data=myData,
                            aes(x = (OverallRatioMitophagy_Cellwise_MinMean),
                                y = log2(OverallProtein_Max) )) +
                  geom_point(size=1,
                             shape=16,
                             alpha=0.1,
                             color=rgb(84, 39, 143, max=255) ) +
                  theme_classic() +
                  #geom_smooth(method = "lm",
                  #            se = FALSE) +
                  stat_poly_line() +
                  stat_poly_eq(aes(label = after_stat(eq.label))) +
                  coord_cartesian(xlim = c(0,3)) + 
                  ylab("Protein max signal") +
                  facet_grid(vars(OverallInfo02), vars(OverallInfo01))+
                  theme(plot.margin = unit(c(1,1,1,3), "lines"),
                        axis.text.x=element_text(size=10, face="bold", color="black"),
                        axis.text.y=element_text(size=8, color="black"),
                        axis.title.y=element_text(size=8, color="black"),
                        axis.title.x=element_text(size=6, color="black"),
                        legend.position="none" ) +
                  xlab(paste0("Mitophagy ratio", "\n",
                              ActualRData_name, "\n",
                              "MinL_", minProteinLabel,
                              "minMaxD_",minMaxDProteinLabel,
                              "MinMax_", minProteinMaxLabel,
                              "MinMedian_", minProteinMedianLabel,
                              "\n",
                              "MinStdDev_", minProteinStdDevLabel,
                              "MinMitophagy_", minMitophagyLabel) )+
                  ggtitle(paste0(MitophagyNow, " ", CompoundNow))
                
                
                
                
                
                
                png(paste0(rootFolder, "/",
                           "08-MitophagyCorrelation_all2", "/",
                           ActualRData_name, "/",
                           "MinL_", minProteinLabel,
                           "minMaxD_",minMaxDProteinLabel,
                           "MinMax_", minProteinMaxLabel,
                           "MinMedian_", minProteinMedianLabel,
                           "MinStdDev_", minProteinStdDevLabel,
                           "MinMitophagy_", minMitophagyLabel,
                           "/",
                           "01_PrtMax_Compound_", paste0(CompoundNow, collapse="_"),
                           "_Mitophagy_", paste0(MitophagyNow, collapse="_"),
                           ".png"),
                    width=1700, height=1400, res=300)
                print(p) 
                dev.off()
                
                #
                
                
                p <- ggplot(data=myData,
                            aes(x = (OverallRatioMitophagy_Cellwise_MinMean),
                                y = log2(OverallProtein_Mean) )) +
                  geom_point(size=1,
                             shape=16,
                             alpha=0.1,
                             color=rgb(84, 39, 143, max=255) ) +
                  theme_classic() +
                  #geom_smooth(method = "lm",
                  #            se = FALSE) +
                  stat_poly_line() +
                  stat_poly_eq(aes(label = after_stat(eq.label))) +
                  coord_cartesian(xlim = c(0,3)) + 
                  ylab("Protein MEAN signal") +
                  facet_grid(vars(OverallInfo02), vars(OverallInfo01))+
                  theme(plot.margin = unit(c(1,1,1,3), "lines"),
                        axis.text.x=element_text(size=10, face="bold", color="black"),
                        axis.text.y=element_text(size=8, color="black"),
                        axis.title.y=element_text(size=8, color="black"),
                        axis.title.x=element_text(size=6, color="black"),
                        legend.position="none" ) +
                  xlab(paste0("Mitophagy ratio", "\n",
                              ActualRData_name, "\n",
                              "MinL_", minProteinLabel,
                              "minMaxD_",minMaxDProteinLabel,
                              "MinMax_", minProteinMaxLabel,
                              "MinMedian_", minProteinMedianLabel,
                              "\n",
                              "MinStdDev_", minProteinStdDevLabel,
                              "MinMitophagy_", minMitophagyLabel) )+
                  ggtitle(paste0(MitophagyNow, " ", CompoundNow))
                
                
                
                
                
                
                png(paste0(rootFolder, "/",
                           "08-MitophagyCorrelation_all2", "/",
                           ActualRData_name, "/",
                           "MinL_", minProteinLabel,
                           "minMaxD_",minMaxDProteinLabel,
                           "MinMax_", minProteinMaxLabel,
                           "MinMedian_", minProteinMedianLabel,
                           "MinStdDev_", minProteinStdDevLabel,
                           "MinMitophagy_", minMitophagyLabel,
                           "/",
                           "01_PrtMean_Compound_", paste0(CompoundNow, collapse="_"),
                           "_Mitophagy_", paste0(MitophagyNow, collapse="_"),
                           ".png"),
                    width=1700, height=1400, res=300)
                print(p) 
                dev.off()
              }
            }
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
          }}}}}}
            
          
            
            





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
dir.create(path = paste0(rootFolder, "/", "10-Amounts"), showWarnings = FALSE)

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
dir.create(path = paste0(rootFolder, "/", "10-Amounts", "/", ActualRData_name, "/"),
           showWarnings = FALSE)


### Create function
ProteinSumPerCell <- function(tableNow){
  tableNow$ProteinSumPerCell <- NA
  for(originNow in unique(tableNow$FileOrigin_base)){
    tableNow[which(tableNow$FileOrigin_base == originNow), "ProteinSumPerCell"] <- sum(tableNow[which(tableNow$FileOrigin_base == originNow), "Protein_IntDen"])
  }
   return(tableNow)
}



### COLLECT RDATA
print(paste0("Coletando os dados novamente"))
load(paste0(rootFolder, "/", "02-PlatesRData", "/", ActualRData))
print(paste0("Dados coletados"))


#
# Trimm table
RawData_size_mitol_mitosig_list <- pblapply(X=RawData_size_mitol_list,
                                            FUN=ProteinSumPerCell)
remove(RawData_size_mitol_list)
gc()

#
RawData_size_mitol_mitosig_list[[1]]

#####
# Retrieve column with data

OverallUniqueConditionTag <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                    FUN=ColumnRetrieval,
                                                    columnNow="UniqueConditionTag")))
OverallCellMeanProtein <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                 FUN=ColumnRetrieval,
                                                 columnNow="CellMeanProtein"))) 
OverallProteinSumPerCell <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                                   FUN=ColumnRetrieval,
                                                   columnNow="ProteinSumPerCell")))
OverallUniqueTag <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                           FUN=ColumnRetrieval,
                                           columnNow="UniqueTag"))) 

OverallCellLine <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                          FUN=ColumnRetrieval,
                                          columnNow="CellLine"))) 

OverallInfo01 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                        FUN=ColumnRetrieval,
                                        columnNow="Info01"))) 
OverallInfo02 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                        FUN=ColumnRetrieval,
                                        columnNow="Info02")))
OverallInfo03 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                        FUN=ColumnRetrieval,
                                        columnNow="Info03")))
OverallInfo04 <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                        FUN=ColumnRetrieval,
                                        columnNow="Info04")))
OverallLabeling <- unname(unlist(pblapply(X=RawData_size_mitol_mitosig_list,
                                          FUN=ColumnRetrieval,
                                          columnNow="Labeling")))



TemporaryData <- (data.frame(OverallUniqueConditionTag,
                             
                             OverallCellMeanProtein,
                             OverallProteinSumPerCell,
                             OverallUniqueTag,
                             OverallCellLine,
                             OverallInfo01,
                             OverallInfo02,
                             OverallInfo03,
                             OverallInfo04,
                             OverallLabeling))

# Making unique
remove(OverallUniqueConditionTag)
remove(OverallCellMeanProtein)
remove(OverallProteinSumPerCell)
remove(OverallUniqueTag)
remove(OverallCellLine)
remove(OverallInfo01)
remove(OverallInfo02)
remove(OverallInfo03)
remove(OverallInfo04)
remove(OverallLabeling)

gc()



# Dplyr para gerar os valores
TempPositiveTable_plyr <- TemporaryData %>%
  group_by(OverallUniqueConditionTag,
           OverallCellMeanProtein,
           OverallProteinSumPerCell,
           OverallUniqueTag,
           OverallCellLine,
           OverallInfo01,
           OverallInfo02,
           OverallInfo03,
           OverallInfo04,
           OverallLabeling) %>%
  dplyr::summarize(
                   TotalMito = n()
  )
data.frame(TempPositiveTable_plyr) -> TemporaryData



print(paste0("Ordenando temporary data para quantidade mito"))
TemporaryData <- TemporaryData[with(TemporaryData, order(
  OverallInfo03,
  OverallCellLine,
  OverallInfo01,
  OverallInfo04,
  OverallInfo02
)), ]
TemporaryData$OverallLabeling <- factor(TemporaryData$OverallLabeling, levels=unique(TemporaryData$OverallLabeling))

remove(RawData_size_mitol_mitosig_list)
#
gc()

if(length(unique(TemporaryData$OverallInfo01)) > 1){
  Info01Now <- unique(TemporaryData$OverallInfo01)[unique(TemporaryData$OverallInfo01) != "01-Mock"]

  myData <- TemporaryData[which(TemporaryData$OverallInfo01 %in% Info01Now &
                                  TemporaryData$OverallInfo02 %in% "01-withGLN" &
                                  TemporaryData$OverallInfo03 %in% "01-Ctl" &
                                  TemporaryData$OverallInfo04 %in% "01-DMSO"),]
  
  
  #
  colorValues_noalpha <- c(rgb(255,255,255, max=255)) #, rgb(189,189,189, max=255)
  colorValues_border_noalpha <- c(rgb(0,0,0, max=255)) #, rgb(0,0,0, max=255)
  
  p <- ggplot(data = myData,
              aes(x = OverallLabeling, #OverallUniqueConditionTag
                  y = log2(OverallProteinSumPerCell),
                  fill=OverallLabeling)) + #, fill = Classification
    geom_violin(trim=TRUE,
                position=position_dodge(width=0.95),
                lwd=0.30,
                width=0.95) +
    # scale_fill_manual(values=colorValues_noalpha) +
    scale_fill_manual(values = rep(colorValues_noalpha, 8),
                      aesthetics = c("fill")) +
    scale_colour_manual(values = rep(colorValues_border_noalpha, 8),
                        aesthetics = c("colour")) +
    theme_classic() +
    coord_cartesian(clip = "off") + #c(0,100)
    theme(plot.margin = unit(c(1,1,1,3), "lines"),
          axis.text.x=element_text(size=10, face="bold", color="black"),
          axis.text.y=element_text(size=8, color="black"),
          axis.title.y=element_text(size=8, color="black"),
          axis.title.x=element_text(size=5, color="black"),
          legend.position="none" ) +
    ylab(paste0("Fluorescence intensity\nlog2(sum)")) + 
    #  scale_y_continuous(expand = c(0,0)) +
    #ggnewscale::new_scale_colour()+
    #ggnewscale::new_scale_fill() +
    geom_boxplot(position=position_dodge(width=0.95),
                 width=0.7,
                 fill="transparent",
                 color=rgb(50,50,50, 200, max=255),
                 outlier.shape = NA,
                 coef=0)
  
  write.csv2(myData, file=paste0(rootFolder, "/",
                    "10-Amounts", "/",
                    ActualRData_name, "/",
                    "Boxplot_",
                    "OnlyGln",
                    ".csv"))
  
  # (2) Bar plots of means + individual jitter points + errors
  png(file=paste0(rootFolder, "/",
                  "10-Amounts", "/",
                  ActualRData_name, "/",
                  "Boxplot_",
                  "OnlyGln",
                  ".png"),
      width = (630+(1.5*250)), #4000
      height = 1400,
      res=500)
  
  print(
    p + geom_text(color="black",
                  x = 0.5,
                  y = min(log2(myData$OverallProteinSumPerCell), na.rm=TRUE)-(0.1*(max(log2(myData$OverallProteinSumPerCell), na.rm=TRUE)-min(log2(myData$OverallProteinSumPerCell), na.rm=TRUE))),
                  label= MyLabelsNow,
                  hjust = 1,
                  vjust = 1,
                  size = 2.8)
  )
  dev.off()
  
  
  }

########## myData
for(Info01Now in unique(TemporaryData$OverallInfo01)  ){
  myData <- TemporaryData[which(TemporaryData$OverallInfo01 %in% Info01Now),]

#
  colorValues_noalpha <- c(rgb(255,255,255, max=255)) #, rgb(189,189,189, max=255)
  colorValues_border_noalpha <- c(rgb(0,0,0, max=255)) #, rgb(0,0,0, max=255)
  
  p <- ggplot(data = myData,
              aes(x = OverallLabeling, #OverallUniqueConditionTag
                  y = log2(OverallProteinSumPerCell),
                  fill=OverallLabeling)) + #, fill = Classification
    geom_violin(trim=TRUE,
                position=position_dodge(width=0.95),
                lwd=0.30,
                width=0.95) +
    # scale_fill_manual(values=colorValues_noalpha) +
    scale_fill_manual(values = rep(colorValues_noalpha, 8),
                      aesthetics = c("fill")) +
    scale_colour_manual(values = rep(colorValues_border_noalpha, 8),
                        aesthetics = c("colour")) +
    theme_classic() +
    coord_cartesian(clip = "off") + #c(0,100)
    theme(plot.margin = unit(c(1,1,1,3), "lines"),
          axis.text.x=element_text(size=10, face="bold", color="black"),
          axis.text.y=element_text(size=8, color="black"),
          axis.title.y=element_text(size=8, color="black"),
          axis.title.x=element_text(size=5, color="black"),
          legend.position="none" ) +
    ylab(paste0("Fluorescence intensity\nlog2(sum)")) + 
    #  scale_y_continuous(expand = c(0,0)) +
    #ggnewscale::new_scale_colour()+
    #ggnewscale::new_scale_fill() +
    geom_boxplot(position=position_dodge(width=0.95),
                 width=0.7,
                 fill="transparent",
                 color=rgb(50,50,50, 200, max=255),
                 outlier.shape = NA,
                 coef=0)
  
  # (2) Bar plots of means + individual jitter points + errors
  png(file=paste0(rootFolder, "/",
                  "10-Amounts", "/",
                  ActualRData_name, "/",
                  "Boxplot_AllConditions_",
                  Info01Now,
                  ".png"),
      width = (630+(4*250)), #4000
      height = 1400,
      res=500)
  
  print(
    p + geom_text(color="black",
                  x = 0.5,
                  y = min(log2(myData$OverallProteinSumPerCell), na.rm=TRUE)-(0.1*(max(log2(myData$OverallProteinSumPerCell), na.rm=TRUE)-min(log2(myData$OverallProteinSumPerCell), na.rm=TRUE))),
                  label= MyLabelsNow,
                  hjust = 1,
                  vjust = 1,
                  size = 2.8)
  )
  dev.off()

  
  write.csv2(myData, paste0(rootFolder, "/",
                    "10-Amounts", "/",
                    ActualRData_name, "/",
                    "Boxplot_AllConditions_",
                    Info01Now,
                    ".csv"))

}

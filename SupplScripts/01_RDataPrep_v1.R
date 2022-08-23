
# If loop to run everything at once
# utils::memory.limit(size=32000)

# Library load
library(dplyr)
library(data.table)
library(pbapply)
library(foreach)
library(doParallel)

n.cores <- 12

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()












# Create function
ReadMultipleFiles_function <- function(fileName){
  # print(fileName)
  temp <- fread(file=fileName,
                stringsAsFactors = FALSE,
                check.names = FALSE,
                header=TRUE,
                dec=",",
                showProgress=FALSE)
  temp$V1 <- NULL
  temp$ID <- NULL
  temp <- data.frame(temp, stringsAsFactors = FALSE, check.names=FALSE)
  return(temp)
}

# Create function
ColumnRetrieval <- function(tableNow, columnNow){
  # print(fileName)
  return(tableNow[,columnNow])
}

# Get rootFolder
rootFolder <- getwd()

# Create folder to receive RDatas
dir.create(path = paste0(rootFolder, "/", "02-PlatesRData"), showWarnings = FALSE)
dir.create(path = paste0(rootFolder, "/", "03-PlatesInfo"), showWarnings = FALSE)


# Available folders
AvailablePlates <- list.dirs(paste0(rootFolder, "/", "01-PlatesCSV"),
          full.names = TRUE, recursive = FALSE)

#

# list.dirs(paste0(rootFolder, "/", "01-PlatesCSV"),
#          full.names = TRUE, recursive = FALSE)[11] -> dirNow

foreach(
  dirNow = list.dirs(paste0(rootFolder, "/", "01-PlatesCSV"),
                     full.names = TRUE, recursive = FALSE)
) %dopar% {
  
  
library(dplyr)
library(data.table)
library(pbapply)
library(foreach)
library(doParallel)

  
  #
  splittedDirNow <- strsplit(dirNow, "/")[[1]]
  
  if(grepl("LC3", splittedDirNow[length(splittedDirNow)])){
    mitoInducerNow <- "CHQ"
  }
  if(grepl("parkin", splittedDirNow[length(splittedDirNow)])){
    mitoInducerNow <- "CCCP"
  }
  
  # Create file list
  fileList <- list.files(dirNow, full.names = TRUE)
  
  ## Perform file reading and list loading
  system.time({
    RawData_list <- pblapply(X=fileList,
                             FUN=ReadMultipleFiles_function)
  })
  
  # garbage collection
  gc()
  

  # Name list
  for(i in 1:length(RawData_list)){
    #
    names(RawData_list)[i] <- RawData_list[[i]][1, "UniqueTag"]
  }
  
  remove(fileList)
  remove(i)
  
  # garbage
  gc()
  
  
  
  for(minMitoLabel in c(300)){
    
    
    # Create info folder
    dir.create(path = paste0(rootFolder, "/", "03-PlatesInfo", "/", splittedDirNow[length(splittedDirNow)], "_minMitoLabel_", minMitoLabel, "/"), showWarnings = FALSE)
    
    ##########
    ########## 01 - HISTOGRAM MITOCHONDRIA PER CELL
    ##########
    ########## START
    ##########
    
    OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_list,
                                                     FUN=ColumnRetrieval,
                                                     columnNow="FileOrigin_base")))
    
    
    # Remove images with less than specific number of mitochondria
    png(paste0(rootFolder, "/", "03-PlatesInfo", "/", splittedDirNow[length(splittedDirNow)], "_minMitoLabel_", minMitoLabel, "/", "01_HistogramMitoPerCell.png"),
        width=1500, height=1200, res=300)
    hist(table(OverallFileOrigin_base),
         breaks=1000,
         xlim=c(0,300))
    dev.off()
    
    remove(OverallFileOrigin_base)
    
    #####
    #####
    # Create function
    GateNumeroMitocondrias_function <- function(tableNow){
      # print(fileName)
      SelectedTempCells <- rownames(table(tableNow$FileOrigin_base))[which(table(tableNow$FileOrigin_base) > 10 &
                                                                             table(tableNow$FileOrigin_base) < 500)]
      tableNow <- tableNow[which(tableNow$FileOrigin_base %in% SelectedTempCells),]
      #
      return(tableNow)
    }
    
    # Trimm table
    RawData_size_list <- pblapply(X=RawData_list,
                                  FUN=GateNumeroMitocondrias_function)
    # Remove old list
    # remove(RawData_list)
    #
    OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_size_list,
                                                     FUN=ColumnRetrieval,
                                                     columnNow="FileOrigin_base")))
    
    
    # Remove images with less than specific number of mitochondria
    png(paste0(rootFolder, "/", "03-PlatesInfo", "/", splittedDirNow[length(splittedDirNow)], "_minMitoLabel_", minMitoLabel, "/", "02_HistogramMitoPerCell_aftercut.png"),
        width=1500, height=1200, res=300)
    hist(table(OverallFileOrigin_base),
         breaks=1000,
         xlim=c(0,300))
    dev.off()
    
    remove(OverallFileOrigin_base)
    remove(GateNumeroMitocondrias_function)
    
    ##########
    ########## END
    ##########
    ########## 02 - HISTOGRAM MITOCHONDRIA PER CELL
    ##########
    
    
    
    
    
    
    
    ##########
    ########## 03 - HISTOGRAM MITOCHONDRIA PER LABELLING
    ##########
    ########## START
    ##########
    
    # Retrieve column with data
    OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_size_list,
                                                     FUN=ColumnRetrieval,
                                                     columnNow="FileOrigin_base")))
    OverallCellMeanMito <- unname(unlist(pblapply(X=RawData_size_list,
                                                  FUN=ColumnRetrieval,
                                                  columnNow="CellMeanMito"))) 
    
    TemporaryData <- unique(data.frame(OverallFileOrigin_base, OverallCellMeanMito))[,"OverallCellMeanMito"]
    
    # Use Mitochondrial probe intensity to filter mitochondria BY CELL MEAN
    # Thresholding
    
    maxMitoLabel <- quantile(x=TemporaryData, probs=1.00)
    
    #
    png(paste0(rootFolder, "/", "03-PlatesInfo", "/", splittedDirNow[length(splittedDirNow)], "_minMitoLabel_", minMitoLabel, "/", "03_HistogramMitochondrialLabelling.png"),
        width=1500, height=1200, res=300)
    
    hist(TemporaryData,
         breaks=100,
         main="Mitochondrial labelling",
         xlab="Mitochondrial labelling intensity",
         ylab="Amount",
         border="white", 
         col="#e31a1c")
    abline(v=minMitoLabel,
           lwd=2,
           col="#4a1486")
    abline(v=maxMitoLabel,
           lwd=2,
           col="#4a1486")
    
    dev.off()
    
    remove(OverallFileOrigin_base)
    remove(OverallCellMeanMito)
    remove(TemporaryData)
    
    ### Create function
    GateMarcacaoMitocondria_function <- function(tableNow){
      # print(fileName)
      tableNow <- tableNow[which(tableNow$CellMeanMito > minMitoLabel &
                                   tableNow$CellMeanMito < maxMitoLabel), ]
      
      return(tableNow)
    }
    
    # Trimm table
    RawData_size_mitol_list <- pblapply(X=RawData_size_list,
                                        FUN=GateMarcacaoMitocondria_function)
    
    # Remove old list
    remove(RawData_size_list)
    gc()
    
    
    # Retrieve column with data
    OverallFileOrigin_base <- unname(unlist(pblapply(X=RawData_size_mitol_list,
                                                     FUN=ColumnRetrieval,
                                                     columnNow="FileOrigin_base")))
    OverallCellMeanMito <- unname(unlist(pblapply(X=RawData_size_mitol_list,
                                                  FUN=ColumnRetrieval,
                                                  columnNow="CellMeanMito"))) 
    
    TemporaryData <- unique(data.frame(OverallFileOrigin_base, OverallCellMeanMito))[,"OverallCellMeanMito"]
    
    #
    png(paste0(rootFolder, "/", "03-PlatesInfo", "/", splittedDirNow[length(splittedDirNow)], "_minMitoLabel_", minMitoLabel, "/", "03_HistogramMitochondrialLabelling_aftercut.png"),
        width=1500, height=1200, res=300)
    
    hist(TemporaryData,
         breaks=100,
         main="Mitochondrial labelling",
         xlab="Mitochondrial labelling intensity",
         ylab="Amount",
         border="white", 
         col="#e31a1c")
    abline(v=minMitoLabel,
           lwd=2,
           col="#4a1486")
    abline(v=maxMitoLabel,
           lwd=2,
           col="#4a1486")
    
    dev.off()
    
    ###
    remove(GateMarcacaoMitocondria_function)
    remove(OverallFileOrigin_base)
    remove(OverallCellMeanMito)
    remove(TemporaryData)
    
    
    ##########
    ########## END
    ##########
    ########## 03 - HISTOGRAM MITOCHONDRIA PER LABELLING
    ##########
    
    
    
    
    #####
    ##### ADD -+ LABELLING
    ### Create function
    PositiveNegativeLabelling_function <- function(tableNow){
      
      # Add Glutaminases columns
      # GAC.wt
      tableNow[,"GAC.wt"] <- "-"
      tableNow[which(tableNow$Info01 == "02-GAC.wt"), "GAC.wt"] <- "+"
      # GAC.K320A
      tableNow[,"GAC.K320A"] <- "-"
      tableNow[which(tableNow$Info01 == "03-GAC.K320A"), "GAC.K320A"] <- "+"
      # GAC.R382D
      tableNow[,"GAC.R382D"] <- "-"
      tableNow[which(tableNow$Info01 == "04-GAC.R382D"), "GAC.R382D"] <- "+"
      # Gln
      tableNow[,"Gln"] <- "-"
      tableNow[which(tableNow$Info02 == "01-withGLN"), "Gln"] <- "+"
      # Mitophagy
      tableNow[,mitoInducerNow] <- "+"
      tableNow[which(tableNow$Info03 == "01-Ctl"), mitoInducerNow] <- "-"
      # CB839
      tableNow[,"CB839"] <- "-"
      tableNow[which(tableNow$Info04 == "02-CB839"), "CB839"] <- "+"
      # DRP
      tableNow[,"DRP"] <- "+"
      tableNow[which(grepl("_Ctl", tableNow$CellLine)), "DRP"] <- "-"
      
      #
      tableNow$Labeling <- sapply(1:dim(tableNow)[1], function(w){
        paste(tableNow[w, c("Gln",
                            "CB839",
                            "DRP",
                            mitoInducerNow,
                            "GAC.wt",
                            "GAC.K320A",
                            "GAC.R382D"
        ) ], collapse = "\n")
      })
      
      return(tableNow)
    }
    
    # Trimm table
    RawData_size_mitol_list2 <- pblapply(X=RawData_size_mitol_list,
                                         FUN=PositiveNegativeLabelling_function)
    
    RawData_size_mitol_list2 -> RawData_size_mitol_list
    remove(RawData_size_mitol_list2)
    
    
    
    # 
    MyLabelsNow <- paste(c("Gln",
                           "CB839",
                           "DRP",
                           mitoInducerNow,
                           "GAC.wt",
                           "GAC.K320A",
                           "GAC.R382D"), collapse = "\n")
    
    
    ######
    ###### SALVA O RDS
    ######
    print("Saving RData")
    save(RawData_size_mitol_list,
         MyLabelsNow,
         file = paste0(rootFolder, "/", "02-PlatesRData", "/",
                       splittedDirNow[length(splittedDirNow)],
                       "_minMitoLabel_", minMitoLabel, ".RData"))
    print("RData saved")
    remove(RawData_size_mitol_list)
    remove(MyLabelsNow)
    gc()
    gc()
  }
  remove(RawData_list)
  gc()
  gc()
}


# Stops cluster
parallel::stopCluster(cl = my.cluster)



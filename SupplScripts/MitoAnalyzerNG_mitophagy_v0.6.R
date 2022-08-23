
#
#
#
#
#
#
# CHANGELOG
#
#
# v0.3 - 20220322
#       Third version
#       Multiple cutoff testing
#
#
#
# v0.2 - 20220320
#       Second version
#       D gating
#
#
# v0.1 - 20220308
#       First version
#

# If loop to run everything at once

  
  #library(rJava)
  library(rChoiceDialogs)
  library(yarrr)
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
  library(data.table)
  library(readr)
  library(ggfortify)
  library(ggpubr)
  library(pbapply)
  
  #library(dplyr)
  #library(ggpubr)
  
  ###### ALPHA COLOR BLOCK
  ### START
  
  # Code from Markus Gesmann https://gist.github.com/mages/5339689
  
  ## Add an alpha value to a colour
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  ### END
  ###### ALPHA COLOR BLOCK
  
  
  
  # Define the folder to perform the analysis
  SelectedFolder <- rchoose.dir(default = getwd(),
                                caption = "Select Results Folder")
  
  
  
  # Folder Creation
  dir.create(path = paste0(SelectedFolder, "\\", "10-DescriptiveResults"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondria"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "12-PiratePlots"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "13-OtherPlots"), showWarnings = FALSE)
  
  
  
  
  
  
  
  
  ##################
  ################## FILE LOADING START
  ##################

  
  # Select MASK folder
  possibleMaskChoices <- list.dirs(SelectedFolder,
                                   full.names = FALSE,
                                   recursive = FALSE)
  
  # Remove non related folders
  possibleMaskChoices <- possibleMaskChoices[!(possibleMaskChoices %in% c("00-CellROIs",
                                                   "01-Images"))]
  
  #### Mitochondrial mask
  MaskSelection <- rselect.list(choices=possibleMaskChoices,
                                preselect = NULL,
                                multiple = FALSE,
                                title = "Mito",
                                graphics = getOption("menu.graphics"))
  
  
  possibleMaskChoices <- possibleMaskChoices[!(possibleMaskChoices %in% MaskSelection)]
  
  
  ### Fractal mask
  FRACTALSelection <- rselect.list(choices=possibleMaskChoices,
                                   preselect = NULL,
                                   multiple = FALSE,
                                   title = "Mito FRACTAL",
                                   graphics = getOption("menu.graphics"))
  
  
  possibleMaskChoices <- possibleMaskChoices[!(possibleMaskChoices %in% FRACTALSelection)]
  
  #### Protein mask
  ProteinMaskSelection <- rselect.list(choices=possibleMaskChoices,
                                preselect = NULL,
                                multiple = FALSE,
                                title = "Protein",
                                graphics = getOption("menu.graphics"))
  
  
  possibleMaskChoices <- possibleMaskChoices[!(possibleMaskChoices %in% ProteinMaskSelection)]
  
  
  ### Fractal mask
  ProteinFRACTALSelection <- rselect.list(choices=possibleMaskChoices,
                                   preselect = NULL,
                                   multiple = FALSE,
                                   title = "Protein FRACTAL",
                                   graphics = getOption("menu.graphics"))
  
  
  
  
  ####
  

  
  
  # Define the folder to perform the analysis
  PlateInfoLocation <- jchoose.files(default = getwd(),
                                     caption = "Select Plate Info",
                                     multi=FALSE)
  # Read plate info
  PlateInfo <- read.xlsx(
    xlsxFile=PlateInfoLocation,
    sheet="WellData",
    startRow = 1,
    colNames = TRUE,
    rowNames = FALSE,
    detectDates = FALSE,
    skipEmptyRows = TRUE,
    skipEmptyCols = TRUE,
    #rows = NULL,
    cols = 1:6,
    check.names = FALSE,
    sep.names = ".",
    namedRegion = NULL,
    na.strings = "NA",
    fillMergedCells = FALSE
  )
  
  
  
  
  # remove empty lines
  PlateInfo <- PlateInfo[which(PlateInfo$CellLine != ""),]
  
  # Create uniquetag
  PlateInfo$UniqueTag <- sapply(1:dim(PlateInfo)[1], function(y){paste0(PlateInfo[y, c("Well", "CellLine", "Info01", "Info02", "Info03", "Info04")], collapse="_")  })
  PlateInfo$UniqueConditionTag <- sapply(1:dim(PlateInfo)[1], function(y){paste0(PlateInfo[y, c("CellLine", "Info01", "Info02", "Info03", "Info04")], collapse="_")  })
  
  
  
  
  

  ## Collect all results
  
  # Create function
  ReadMultipleFiles_function <- function(fileName){
    # print(fileName)
    temp <- fread(file=fileName,
                  stringsAsFactors = FALSE,
                  check.names = FALSE,
                  header=TRUE)
    colnames(temp)[1] <- "ID"
    temp$FileOrigin <- strsplit(fileName, "/")[[1]][2]
    return(temp)
  }
  
  
  ######### MASK GENERAL RESULTS
  # Collect all 
  fileListMask <- grep("_Mito_Results.csv",
                       list.files(paste0(SelectedFolder, "\\", MaskSelection),
                                  full.names = TRUE),
                       value=TRUE)

  # Create the receiving list
  MaskResults_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    MaskResults_list <- pblapply(X=fileListMask,
                                     FUN=ReadMultipleFiles_function)
  })
  

  # Remove if there is any results for now
  if( exists("MaskResults") ){remove(MaskResults)}
  MaskResults <- bind_rows(MaskResults_list)
  remove(MaskResults_list)
  
  
  
  ######### MITOPHAGY GENERAL RESULTS
  # Collect all 
  fileListMitophagy <- grep("_Mitophagy_Results.csv",
                       list.files(paste0(SelectedFolder, "\\", MaskSelection),
                                  full.names = TRUE),
                       value=TRUE)
  

  # Create the receiving list
  MitophagyResults_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    MitophagyResults_list <- pblapply(X=fileListMitophagy,
                                 FUN=ReadMultipleFiles_function)
  })
  
  
  # Remove if there is any results for now
  if( exists("MitophagyResults") ){remove(MitophagyResults)}
  MitophagyResults <- bind_rows(MitophagyResults_list)
  remove(MitophagyResults_list)
  
  

  
  ######### PROTEIN GENERAL RESULTS
  # Collect all 
  fileListProtein <- grep("_Mito_Protein_Results.csv",
                            list.files(paste0(SelectedFolder, "\\", MaskSelection),
                                       full.names = TRUE),
                            value=TRUE)

  # Create the receiving list
  ProteinResults_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    ProteinResults_list <- pblapply(X=fileListProtein,
                                      FUN=ReadMultipleFiles_function)
  })
  
  
  # Remove if there is any results for now
  if( exists("ProteinResults") ){remove(ProteinResults)}
  ProteinResults <- bind_rows(ProteinResults_list)
  remove(ProteinResults_list)
  

  ####
  

  #### Create uniqueTag to remove missing mitochondrias
  MaskResults$CollapsedTag <- do.call(paste, c(MaskResults[, c("X", "Y", "FileOrigin")], sep="-"))
  MitophagyResults$CollapsedTag <- do.call(paste, c(MitophagyResults[, c("X", "Y", "FileOrigin")], sep="-"))
  ProteinResults$CollapsedTag <- do.call(paste, c(ProteinResults[, c("X", "Y", "FileOrigin")], sep="-"))
  
  MaskResults$CollapsedTag <- gsub("_Mito_Results.csv", "", MaskResults$CollapsedTag)
  MitophagyResults$CollapsedTag <- gsub("_Mitophagy_Results.csv", "", MitophagyResults$CollapsedTag)
  ProteinResults$CollapsedTag <- gsub("_Mito_Protein_Results.csv", "", ProteinResults$CollapsedTag)
  
  
  #
  CollapsedTagUnique <- unique(intersect(intersect(MaskResults$CollapsedTag, MitophagyResults$CollapsedTag), ProteinResults$CollapsedTag))
  
  #
  MaskResults <- MaskResults[which(MaskResults$CollapsedTag %in% CollapsedTagUnique), ]
  MitophagyResults <- MitophagyResults[which(MitophagyResults$CollapsedTag %in% CollapsedTagUnique), ]
  ProteinResults <- ProteinResults[which(ProteinResults$CollapsedTag %in% CollapsedTagUnique), ]
  
  
  
  ### Merge FRET and Mask results
  MaskResults <- MaskResults[order(MaskResults$CollapsedTag), ]
  MitophagyResults <- MitophagyResults[order(MitophagyResults$CollapsedTag), ]
  ProteinResults <- ProteinResults[order(ProteinResults$CollapsedTag), ]
  
  
  dim(MaskResults)
  dim(MitophagyResults)
  dim(ProteinResults)
  

  # Select relevant colnames
  MitophagyResults <- MitophagyResults[,c("Mean", "StdDev", "Mode", "Min", "Max", "IntDen", "Median", "RawIntDen")]
  colnames(MitophagyResults) <- paste("Mitophagy", colnames(MitophagyResults), sep="_")
  # Merge dataset
  FinalResults <- cbind(MaskResults, MitophagyResults)
  
  # Select relevant colnames
  ProteinResults <- ProteinResults[,c("Mean", "StdDev", "Mode", "Min", "Max", "IntDen", "Median", "RawIntDen")]
  colnames(ProteinResults) <- paste("Protein", colnames(ProteinResults), sep="_")
  # Merge dataset
  FinalResults <- cbind(FinalResults, ProteinResults)
  
  remove(MaskResults)
  remove(MitophagyResults)
  remove(ProteinResults)
  remove(fileListMask)
  remove(fileListMitophagy)
  remove(fileListProtein)
  remove(temp)
  remove(CollapsedTagUnique)
  remove(i)
  remove(j)
  remove(possibleMaskChoices)
  
  
  
  ######################
  #FinalResults$CollapsedTag <- NULL

  ###############
  
  ######### FRACTAL
  # Collect all 
  fileListFRACTAL <- grep(".csv",
                          list.files(paste0(SelectedFolder, "\\", FRACTALSelection),
                                     full.names = TRUE),
                          value=TRUE)

  
  # Create the receiving list
  FRACTALResults_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    FRACTALResults_list <- pblapply(X=fileListFRACTAL,
                                    FUN=ReadMultipleFiles_function)
  })
  
  
  # Remove if there is any results for now
  if( exists("FRACTALResults") ){remove(FRACTALResults)}
  FRACTALResults <- bind_rows(FRACTALResults_list)
  remove(FRACTALResults_list)
  remove(fileListFRACTAL)

  
###################
  FinalResults$WellNumbering <- sapply(FinalResults$FileOrigin, function(w){strsplit(w, "_")[[1]][1]})
  FinalResults$StackPos <- sapply(FinalResults$FileOrigin, function(w){strsplit(w, "_")[[1]][3]})
  FinalResults$FieldPos <- sapply(FinalResults$FileOrigin, function(w){strsplit(w, "_")[[1]][5]})
  FinalResults$Timepoint <- sapply(FinalResults$FileOrigin, function(w){strsplit(w, "_")[[1]][7]})
  FinalResults$Cell <- sapply(FinalResults$FileOrigin, function(w){strsplit(w, "_")[[1]][9]})
  
  
  # Fix to add FRACTAL results
  # Ordena por Origin e ID
  FinalResults <- FinalResults[order(FinalResults$FileOrigin, FinalResults$ID), ]
  
  
  # Cria o FileOrigin_base
  FRACTALResults$FileOrigin_base <- gsub("_Fractal_Results.csv", "", FRACTALResults$FileOrigin)
  FinalResults$FileOrigin_base <- gsub("_Mito_Results.csv", "", FinalResults$FileOrigin)
  
  # Faz a intersecção dos file origins
  FileOrigin_base_unique <- unique(intersect(FRACTALResults$FileOrigin_base, FinalResults$FileOrigin_base))
  
  #
  FinalResults <- FinalResults[which(FinalResults$FileOrigin_base %in% FileOrigin_base_unique), ]
  
  # go get my fractals
  # Add Base IDs
  FinalResults$FileOrigin_base_ID <- do.call(paste, c(FinalResults[, c("FileOrigin_base", "ID")], sep="-"))
  FRACTALResults$FileOrigin_base_ID <- do.call(paste, c(FRACTALResults[, c("FileOrigin_base", "ID")], sep="-"))
  
  # Add fractal value
  FinalResults <- cbind(FinalResults, FRACTALResults[match(FinalResults$FileOrigin_base_ID, FRACTALResults$FileOrigin_base_ID), "D"] )
  
  # Remove fractal info
  remove(FRACTALResults)
  
  ### ADD PLATE INFO TO TABLE
  FinalResults <- cbind(FinalResults, PlateInfo[match(FinalResults$WellNumbering, PlateInfo$Well), c("CellLine",
                                                                                                     "Info01",
                                                                                                     "Info02",
                                                                                                     "Info03",
                                                                                                     "Info04",
                                                                                                    "UniqueTag",
                                                                                                     "UniqueConditionTag")])
  
  
  
  #### ADD MITOCHONDRIAL MEAN STAINING
  FinalResults$CellMeanMito <- NA
  FinalResults$CellMeanProtein <- NA
  FinalResults$CellMeanMitophagy <- NA
  
  # Dplyt para gerar os valores
  intermediateForCellMean <- FinalResults %>%
    group_by(FileOrigin_base) %>%
    dplyr::summarize(CellMeanMito = mean(Mean, na.rm=TRUE),
                     CellMeanProtein = mean(Protein_Mean, na.rm=TRUE),
                     CellMeanMitophagy = mean(Mitophagy_Mean, na.rm=TRUE)
    )
  
  #
  intermediateForCellMean <- data.frame(intermediateForCellMean)
  intermediateForCellMean$FileOrigin_base -> rownames(intermediateForCellMean)
  
  #
  FinalResults$CellMeanMito <- as.numeric(intermediateForCellMean[FinalResults$FileOrigin_base , "CellMeanMito"])
  FinalResults$CellMeanProtein <- as.numeric(intermediateForCellMean[FinalResults$FileOrigin_base , "CellMeanProtein"])
  FinalResults$CellMeanMitophagy <- as.numeric(intermediateForCellMean[FinalResults$FileOrigin_base , "CellMeanMitophagy"])
  
  
  #### FINALIZANDO O FINALRESULTS
  
  remove(intermediateForCellMean)
 
  #
  write.csv2(FinalResults,
             file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "FinalResults_complete.csv"))
  
  
  
  
  ################## PROTEIN START

  ## Collect all results
  ######### MASK PROTEIN RESULTS
  # Collect all 
  fileListMaskProtein <- grep("_Protein_Results.csv",
                              list.files(paste0(SelectedFolder, "\\", ProteinMaskSelection),
                                         full.names = TRUE),
                              value=TRUE)

  # Create the receiving list
  MaskResultsProtein_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    MaskResultsProtein_list <- pblapply(X=fileListMaskProtein,
                                    FUN=ReadMultipleFiles_function)
  })
  
  
  # Remove if there is any results for now
  if( exists("MaskResultsProtein") ){remove(MaskResultsProtein)}
  MaskResultsProtein <- bind_rows(MaskResultsProtein_list)
  remove(MaskResultsProtein_list)
  
  
  
  
  
  
  ######### MITOPHAGY PROTEIN RESULTS
  # Collect all 
  fileListProteinMitophagy <- grep("_Mitophagy_Results.csv",
                                   list.files(paste0(SelectedFolder, "\\", ProteinMaskSelection),
                                              full.names = TRUE),
                                   value=TRUE)

  # Create the receiving list
  ProteinMitophagyResults_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    ProteinMitophagyResults_list <- pblapply(X=fileListProteinMitophagy,
                                    FUN=ReadMultipleFiles_function)
  })

  # Remove if there is any results for now
  if( exists("ProteinMitophagyResults") ){remove(ProteinMitophagyResults)}
  ProteinMitophagyResults <- bind_rows(ProteinMitophagyResults_list)
  remove(ProteinMitophagyResults_list)

  
  
  #### Create uniqueTag to remove missing mitochondrias
  MaskResultsProtein$CollapsedTag <- do.call(paste, c(MaskResultsProtein[, c("X", "Y", "FileOrigin")], sep="-"))
   ProteinMitophagyResults$CollapsedTag <- do.call(paste, c(ProteinMitophagyResults[, c("X", "Y", "FileOrigin")], sep="-"))
  
   MaskResultsProtein$CollapsedTag <- gsub("_Protein_Results.csv", "", MaskResultsProtein$CollapsedTag)
  ProteinMitophagyResults$CollapsedTag <- gsub("_Mitophagy_Results.csv", "", ProteinMitophagyResults$CollapsedTag)
  
  #
  CollapsedTagUnique <- unique((intersect(MaskResultsProtein$CollapsedTag, ProteinMitophagyResults$CollapsedTag)))
  
  #
  MaskResultsProtein <- MaskResultsProtein[which(MaskResultsProtein$CollapsedTag %in% CollapsedTagUnique), ]
  ProteinMitophagyResults <- ProteinMitophagyResults[which(ProteinMitophagyResults$CollapsedTag %in% CollapsedTagUnique), ]
  
  
  
  ### Merge and Mask results
  MaskResultsProtein <- MaskResultsProtein[order(MaskResultsProtein$CollapsedTag), ]
  ProteinMitophagyResults <- ProteinMitophagyResults[order(ProteinMitophagyResults$CollapsedTag), ]
  
  
  dim(MaskResultsProtein)
  dim(ProteinMitophagyResults)
  

  # Select relevant colnames
  ProteinMitophagyResults <- ProteinMitophagyResults[,c("Mean", "StdDev", "Mode", "Min", "Max", "IntDen", "Median", "RawIntDen")]
  colnames(ProteinMitophagyResults) <- paste("Mitophagy", colnames(ProteinMitophagyResults), sep="_")
  # Merge dataset
  FinalResultsProtein <- cbind(MaskResultsProtein, ProteinMitophagyResults)
  
  
  
  
  
  
  
  ######### FRACTAL
  # Collect all 
  fileListFRACTALProtein <- grep(".csv",
                          list.files(paste0(SelectedFolder, "\\", ProteinFRACTALSelection),
                                     full.names = TRUE),
                          value=TRUE)
  
  # Create the receiving list
  FRACTALResultsProtein_list <- list()
  
  # Create the dataframe list
  # with pblapply
  system.time({
    FRACTALResultsProtein_list <- pblapply(X = fileListFRACTALProtein,
                                         FUN = ReadMultipleFiles_function)
  })

  # Remove if there is any results for now
  if( exists("FRACTALResultsProtein") ){remove(FRACTALResultsProtein)}
  FRACTALResultsProtein <- bind_rows(FRACTALResultsProtein_list)
  remove(FRACTALResultsProtein_list)
  
  
  remove(fileListFRACTALProtein)

  
  ###################
  FinalResultsProtein$WellNumbering <- sapply(FinalResultsProtein$FileOrigin, function(w){strsplit(w, "_")[[1]][1]})
  FinalResultsProtein$StackPos <- sapply(FinalResultsProtein$FileOrigin, function(w){strsplit(w, "_")[[1]][3]})
  FinalResultsProtein$FieldPos <- sapply(FinalResultsProtein$FileOrigin, function(w){strsplit(w, "_")[[1]][5]})
  FinalResultsProtein$Timepoint <- sapply(FinalResultsProtein$FileOrigin, function(w){strsplit(w, "_")[[1]][7]})
  FinalResultsProtein$Cell <- sapply(FinalResultsProtein$FileOrigin, function(w){strsplit(w, "_")[[1]][9]})
  
  
  # Fix to add FRACTAL results
  # Ordena por Origin e ID
  FinalResultsProtein <- FinalResultsProtein[order(FinalResultsProtein$FileOrigin, FinalResultsProtein$ID), ]
  
  
  # Cria o FileOrigin_base
  FRACTALResultsProtein$FileOrigin_base <- gsub("_Fractal_Results.csv", "", FRACTALResultsProtein$FileOrigin)
  FinalResultsProtein$FileOrigin_base <- gsub("_Protein_Results.csv", "", FinalResultsProtein$FileOrigin)
  
  # Faz a intersecção dos file origins
  FileOrigin_base_unique <- unique(intersect(FRACTALResultsProtein$FileOrigin_base, FinalResultsProtein$FileOrigin_base))
  
  #
  FinalResultsProtein <- FinalResultsProtein[which(FinalResultsProtein$FileOrigin_base %in% FileOrigin_base_unique), ]
  
  # go get my fractals
  # Add Base IDs
  FinalResultsProtein$FileOrigin_base_ID <- do.call(paste, c(FinalResultsProtein[, c("FileOrigin_base", "ID")], sep="-"))
  FRACTALResultsProtein$FileOrigin_base_ID <- do.call(paste, c(FRACTALResultsProtein[, c("FileOrigin_base", "ID")], sep="-"))
  
  # Add fractal value
  FinalResultsProtein <- cbind(FinalResultsProtein, FRACTALResultsProtein[match(FinalResultsProtein$FileOrigin_base_ID, FRACTALResultsProtein$FileOrigin_base_ID), "D"] )
  
  # Remove fractal info
  remove(FRACTALResultsProtein)
  
  ### ADD PLATE INFO TO TABLE
  FinalResultsProtein <- cbind(FinalResultsProtein, PlateInfo[match(FinalResultsProtein$WellNumbering, PlateInfo$Well), c("CellLine",
                                                                                                     "Info01",
                                                                                                     "Info02",
                                                                                                     "Info03",
                                                                                                     "Info04",
                                                                                                     "UniqueTag",
                                                                                                     "UniqueConditionTag")])
  
  
  
  #### ADD PROTEIN MEAN STAINING
  FinalResultsProtein$CellMeanProtein <- NA
  FinalResultsProtein$CellMeanMitophagy <- NA
  
  # Dplyt para gerar os valores
  intermediateForCellMean <- FinalResultsProtein %>%
    group_by(FileOrigin_base) %>%
    dplyr::summarize(CellMeanProtein = mean(Mean, na.rm=TRUE),
                     CellMeanMitophagy = mean(Mitophagy_Mean, na.rm=TRUE)
    )
  
  #
  intermediateForCellMean <- data.frame(intermediateForCellMean)
  intermediateForCellMean$FileOrigin_base -> rownames(intermediateForCellMean)
  
  #
  FinalResultsProtein$CellMeanProtein <- as.numeric(intermediateForCellMean[FinalResultsProtein$FileOrigin_base , "CellMeanProtein"])
  FinalResultsProtein$CellMeanMitophagy <- as.numeric(intermediateForCellMean[FinalResultsProtein$FileOrigin_base , "CellMeanMitophagy"])
  
  
  #### FINALIZANDO O FINALRESULTS
  
  remove(intermediateForCellMean)

  #
  write.csv2(FinalResultsProtein,
             file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "FinalResultsProtein_complete.csv"))
  
  
  remove(FinalResults)
  remove(FinalResultsProtein)
  remove(FRACTALResults)
  remove(FRACTALSelection)
  remove(MaskResults)
  remove(MaskResultsProtein)
  remove(PlateInfo)
  remove(ProteinMitophagyResults)
  remove(fileListFRACTAL)
  remove(fileListProtein)
  remove(fileListFRACTALProtein)
  remove(fileListMaskProtein)
  remove(fileListMaskProtein)
  remove(fileListProteinMitophagy)
  remove(CollapsedTagUnique)
  remove(FileOrigin_base_unique)
  remove(MaskSelection)
  remove(PlateInfoLocation)
  remove(ProteinFRACTALSelection)
  remove(ProteinMaskSelection)
  remove(SelectedFolder)

  ################## PROTEIN END  
  
  
  
  ##################
  ################## FILE LOADING END
  ##################
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
############### START HERE
  
  #####
  FinalResults <- fread(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "FinalResults_complete.csv"),
        stringsAsFactors = FALSE,
        check.names = FALSE,
        header=TRUE,
        sep=";",
        dec=",")
  colnames(FinalResults)[1] <- "ID"
  
  
  #####
  FinalResultsProtein <- fread(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "FinalResultsProtein_complete.csv"),
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        header=TRUE,
                        sep=";",
                        dec=",")
  colnames(FinalResults)[1] <- "ID"

  
  
  FinalResults <- data.frame(FinalResults, check.names = FALSE, stringsAsFactors = FALSE)
  FinalResultsProtein <- data.frame(FinalResultsProtein, check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ########################
  
  ## BLACKLIST
  # Blacklist <- c("004018_01-PC3_01-Mock_02-woGLN_01-Ctl_.")
  # FinalResults <- FinalResults[!(FinalResults$UniqueTag %in% Blacklist),]
  # FinalResultsProtein <- FinalResultsProtein[!(FinalResultsProtein$UniqueTag %in% Blacklist),]

  #########################
  

  ####
  ####
  #### ANALYSIS 
  #### START
  ####
  
  
  
  
  
  ####
  ####
  #### MITOCHONDRIA 
  #### START
  ####
  
  # Correlation D and Intensity
  x_now <- FinalResults$Mean
  y_now <- FinalResults$D
  
  # Correlation plot
  png(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "01_D-Correlation_mitochondrialwise.png"),
      width=2000, height=2000, res=400)
  plot(x=x_now, y=y_now, pch=16,
       col=rgb(37,37,37, 10, max=255),
       cex=0.75, cex.main=0.75,
       main=paste0("pearson cor = ", cor(x=x_now, y=y_now), " p ", cor.test(x=x_now, y=y_now, method="pearson")[["p.value"]]),
       xlab="",
       ylab="Fractal D Value",
       cex.axis=1.3, yaxt="n", xaxt="n")
  axis(1, col.axis="black", las=2, lwd=2, cex.axis=1.25)
  axis(2, col.axis="black", las=2, lwd=2, cex.axis=1.25)
  box(lwd=2)
  abline(lm(y_now~x_now), lwd=2, lty=3, col=rgb(12,12,12,250, max=255))
  
  dev.off()
  
  
  
  
  ############################
  x_now <- FinalResults$CellMeanMito
  y_now <- FinalResults$D
  
  # Correlation plot
  png(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "01_D-Correlation_cellwise.png"),
      width=2000, height=2000, res=400)
  plot(x=x_now, y=y_now, pch=16,
       col=rgb(37,37,37, 10, max=255),
       cex=0.75, cex.main=0.75,
       main=paste0("pearson cor = ", cor(x=x_now, y=y_now), " p ", cor.test(x=x_now, y=y_now, method="pearson")[["p.value"]]),
       xlab="",
       ylab="Fractal D Value",
       cex.axis=1.3, yaxt="n", xaxt="n")
  axis(1, col.axis="black", las=2, lwd=2, cex.axis=1.25)
  axis(2, col.axis="black", las=2, lwd=2, cex.axis=1.25)
  box(lwd=2)
  abline(lm(y_now~x_now), lwd=2, lty=3, col=rgb(12,12,12,250, max=255))
  
  dev.off()
  
  
  
  
  
  
  
  ########## Cell filtering
  
  # Remove images with less than specific number of mitochondria
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "01_HistogramMitoPerCell.png"),
      width=1500, height=1200, res=300)
  hist(table(FinalResults$FileOrigin_base),
       breaks=1000,
       xlim=c(0,300))
  dev.off()

  ###
  #
  # GATE POR NUMERO DE MITOCONDRIAS
  # 20 / 500
  SelectedCells <- rownames(table(FinalResults$FileOrigin_base))[which(table(FinalResults$FileOrigin_base) > 10 &
                                                                         table(FinalResults$FileOrigin_base) < 500)]
  
  length(unique(FinalResults$FileOrigin_base))
  length(SelectedCells)
  
  FinalResults <- FinalResults[which(FinalResults$FileOrigin_base %in% SelectedCells),]
  
  
  # Use Mitochondrial probe intensity to filter mitochondria BY CELL MEAN
  # Thresholding
  minMitoLabel <- 0.10
  maxMitoLabel <- 1.00
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "01_HistogramYFPIntensity.png"),
      width=1500, height=1200, res=300)
  
  
  
  
  hist(as.data.frame(unique(FinalResults[,c("FileOrigin_base", "CellMeanMito")]))[,"CellMeanMito"],
       breaks=100,
       main="Mitochondrial labelling",
       xlab="YFP CellMeanMito intensity",
       ylab="Amount",
       border="white", 
       col="#e31a1c")
  abline(v=quantile(x=as.data.frame(unique(FinalResults[,c("FileOrigin_base", "CellMeanMito")]))[,"CellMeanMito"], probs=minMitoLabel),
         lwd=2,
         col="#4a1486")
  abline(v=quantile(x=as.data.frame(unique(FinalResults[,c("FileOrigin_base", "CellMeanMito")]))[,"CellMeanMito"], probs=maxMitoLabel),
         lwd=2,
         col="#4a1486")
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  # Cut
  FinalResults <- FinalResults[which(FinalResults$CellMeanMito < quantile(x=unique(FinalResults[,c("FileOrigin_base", "CellMeanMito")])[,"CellMeanMito"], probs=maxMitoLabel) &
                                       FinalResults$CellMeanMito > quantile(x=unique(FinalResults[,c("FileOrigin_base", "CellMeanMito")])[,"CellMeanMito"], probs=minMitoLabel)),]
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "02_HistogramYFPIntensity_filtered.png"),
      width=1500, height=1200, res=300)
  hist(as.data.frame(unique(FinalResults[,c("FileOrigin_base", "CellMeanMito")]))[,"CellMeanMito"],
       breaks=100,
       main="Mitochondrial labelling",
       xlab="YFP CellMeanMito intensity",
       ylab="Amount",
       border="white", 
       col="#e31a1c")
  dev.off()
  
  
  # Use YFP area to filter mitochondria
  # Thresholding
  minMitoSize <- 0.00
  maxMitoSize <- 1.00
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "03_HistogramMitoArea.png"),
      width=1500, height=1200, res=300)
  hist(FinalResults$Area,
       breaks=200,
       main="Area",
       xlab="Mitochondrial area",
       ylab="Amount",
       border="white", 
       col="#e31a1c")
  abline(v=quantile(x=FinalResults$Area, probs=minMitoSize),
         lwd=2,
         col="#4a1486")
  abline(v=quantile(x=FinalResults$Area, probs=maxMitoSize),
         lwd=2,
         col="#4a1486")
  dev.off()
  
  # Cut
  FinalResults <- FinalResults[which(FinalResults$Area < quantile(x=FinalResults$Area, probs=maxMitoSize) &
                                       FinalResults$Area > quantile(x=FinalResults$Area, probs=minMitoSize)),]
  
  
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "04_HistogramMitoArea_filtered.png"),
      width=1500, height=1200, res=300)
  hist(FinalResults$Area,
       breaks=200,
       main="Area",
       xlab="Mitochondrial area (um)",
       ylab="Amount",
       border="white", 
       col="#e31a1c")
  dev.off()
  
  
  
  
  
  ####
  #### Bounding shape descriptors
  ####
  
  #Define morphology columns
  morphologyColumms <- c("Area",
                         "Perim.",
                         "Major",
                         "Circ.",
                         #"Round",
                         #"Solidity",
                         "Feret",
                         "D")
  #"Minor",#"AR",
  
  # Remove NA lines
  dim(FinalResults)
  FinalResults <- FinalResults[(rowSums(is.na(FinalResults[,morphologyColumms])) == 0), ]
  dim(FinalResults)
  
  # Define seed for reproductibility
  set.seed(165454)
  
  # Perform PCA
  pca_res <- prcomp(FinalResults[,morphologyColumms],
                    scale. = TRUE)
  
  # PCA_Contributions
  write.csv2(pca_res[[2]],
             file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "PCA_Contributions.csv"))
  
  # 
  kmeans_res <- kmeans(FinalResults[,morphologyColumms],
                       centers = 10,
                       iter.max = 200,
                       nstart = 1,
                       algorithm = c("Hartigan-Wong", "Lloyd", "MacQueen")[1]
  )
  
  # ClusterMeans
  write.csv2(kmeans_res$centers,
             file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "Kmeans_centers.csv"))
  
  
  # kmeans_res$cluster
  
  
  #
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "05_Kmeans_clustering.png"),
      width=1500, height=1200, res=300)
  
  print(
    autoplot(kmeans_res,
             data = FinalResults[,morphologyColumms],
             loadings = TRUE,
             loadings.colour = "blue",
             loadings.label = TRUE,
             loadings.label.size = 3)
  )
  
  dev.off()
  
  
  
  
  
  
  
  ##### DEFINE POSITIVE MITOCHONDRIA
  # Descriptive plots
  
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "06_ProteinLabel_control.png"),
      width=1500, height=1000, res=300)
  hist(FinalResults[which(FinalResults$Info01 == "01-Mock"), "Protein_Mean"],
       breaks=200,
       main="Stain Negatives",
       xlab="Mitochondrial GAC stain",
       ylab="Amount",
       ylim=c(0,0.040),
       xlim=c(0, 800),
       freq=FALSE,
       border="white", 
       col="#e31a1c")
  
  dev.off()
  
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "07_ProteinLabel_positives.png"),
      width=1500, height=1000, res=300)
  hist(FinalResults[which(FinalResults$Info01 != "01-Mock"), "Protein_Mean"],
       breaks=2000,
       main="Stain Positives",
       xlab="Mitochondrial GAC stain",
       ylab="Amount",
       ylim=c(0,0.040),
       xlim=c(0, 800),
       freq=FALSE,
       border="white", 
       col="blue")
  
  dev.off()
  
  
  ######
  dir.create(path = paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "PossibleThresholds", "\\"), showWarnings = FALSE)
  ######
  dir.create(path = paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "PossibleThresholds_csv", "\\"), showWarnings = FALSE)
  
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  
  #
  possibleProteinThresholds <- seq(from=100, to=1500, by=10)  # c(100, 400)
  possibleMitoNumberThresholds <- seq(from=1, to=10, by=9)  # c(1, 2)
  
  myRowNames <- vector()
  for(NumberPositiveMito in possibleMitoNumberThresholds){
    for(minProteinLabel in possibleProteinThresholds){
      myRowNames <- c(myRowNames,
        paste0(NumberPositiveMito, "_", minProteinLabel)
      )
    }
  }
   
  possibleProteinThresholdsDF <- data.frame(CutNames=myRowNames)

  possibleProteinThresholdsDF$MitoNumber <- as.numeric(unname(sapply(possibleProteinThresholdsDF$CutNames, function(h){
    unlist(strsplit(h, "_"))[1]
  })))
  possibleProteinThresholdsDF$Threshold <- as.numeric(unname(sapply(possibleProteinThresholdsDF$CutNames, function(h){
    unlist(strsplit(h, "_"))[2]
  })))
  
  possibleProteinThresholdsDF$MockMeanPositivity <- NA
  row.names(possibleProteinThresholdsDF) <- myRowNames
  remove(myRowNames)

  
  #
  for(NumberPositiveMito in possibleMitoNumberThresholds){
    for(minProteinLabel in possibleProteinThresholds){
      
      #      
      FinalResults$PositiveCell <- FALSE
      FinalResults$PositiveMitochondria <- FALSE

      # Use protein stain to define positive mitochondria
      # Thresholding
      # minProteinLabel <- 100
      
      # Marca as mitocôndrias positivas
      FinalResults$PositiveMitochondria[FinalResults$Protein_Mean > minProteinLabel] <- TRUE
  
      #
      #unique(FinalResults$FileOrigin_base)
      
      #       # Marca as CÉLULAS positivas
      # Tentativa com DPLYR para agilizar
      TempPositiveTable_plyr_toLabel <- FinalResults %>%
        group_by(FileOrigin_base) %>%
        summarize(positiveMito = sum(PositiveMitochondria, na.rm=TRUE),
                  totalMito = n()
        )
      FinalResults[(FinalResults$FileOrigin_base %in% TempPositiveTable_plyr_toLabel$FileOrigin_base[TempPositiveTable_plyr_toLabel$positiveMito > NumberPositiveMito]), "PositiveCell"] <- TRUE
      remove(TempPositiveTable_plyr_toLabel)
      # Marca as CÉLULAS positivas
      #for(i in unique(FinalResults$FileOrigin_base)){
      #  if(sum(FinalResults[FinalResults$FileOrigin_base == i, "PositiveMitochondria"]) > NumberPositiveMito){
      #    FinalResults[FinalResults$FileOrigin_base == i, "PositiveCell"] <- TRUE
      #  } else {
      #  }
      #}

      # Plot conditions
      TempPositiveTable <- unique(FinalResults[,c("UniqueConditionTag", "UniqueTag", "FileOrigin_base", "PositiveCell")])
      TempPositiveTable$PositiveCell[TempPositiveTable$PositiveCell] <- 1
      TempPositiveTable$PositiveCell[TempPositiveTable$PositiveCell == FALSE] <- 0
      
      
      
      
      # Dplyr para gerar os valores
      TempPositiveTable_plyr <- TempPositiveTable %>%
        group_by(UniqueTag, UniqueConditionTag) %>%
        dplyr::summarize(PositiveSum = sum(PositiveCell, na.rm=TRUE),
                         TotalCell = n()
        )
      
      TempPositiveTable_plyr$PositiveSum/TempPositiveTable_plyr$TotalCell -> TempPositiveTable_plyr$PositivePercent
      #
      TempPositiveTable_plyr$PositivePercent * 100 -> TempPositiveTable_plyr$PositivePercent
      
      
      myComparisons <- list()
      k <- 0
      for(i in 1:length(unique(TempPositiveTable_plyr$UniqueConditionTag))){
        for(j in 1:length(unique(TempPositiveTable_plyr$UniqueConditionTag))){
          if(j>i){
            k<- k+1
            myComparisons[[k]] <- c(unique(TempPositiveTable_plyr$UniqueConditionTag)[i],
                                    unique(TempPositiveTable_plyr$UniqueConditionTag)[j])
          }
        }
      }
      
      
      png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "PossibleThresholds", "\\", "MinMitos_", NumberPositiveMito, "_Thr_", stringr::str_pad(minProteinLabel, nchar(as.character(max(possibleProteinThresholds))), pad = "0"), ".png"),
          width=1500, height=1500, res=200)
      
      #bottom, left, top and right
      # par(mar=c(5.1,20.1,4.1,2.1))
      print(
        ggbarplot(TempPositiveTable_plyr,
                  x = "UniqueConditionTag",
                  y = "PositivePercent", 
                  add = c("mean_se", "jitter")) +
          rotate_x_text(90) +
          
          stat_compare_means() # comparisons=myComparisons
      )
      
      dev.off()
      
      #bottom, left, top and right
      #par(mar=c(5.1,4.1,4.1,2.1))
      TempPositiveTable_plyr <- data.frame(TempPositiveTable_plyr)
      
      
      # Save
      write.csv2(TempPositiveTable_plyr, 
        file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "PossibleThresholds_csv", "\\", "MinMitos_", NumberPositiveMito, "_Thr_", stringr::str_pad(minProteinLabel, nchar(as.character(max(possibleProteinThresholds))), pad = "0"), ".csv")
      )
      
      
      ### Add to the table
      possibleProteinThresholdsDF[which(possibleProteinThresholdsDF$MitoNumber == NumberPositiveMito &
                                          possibleProteinThresholdsDF$Threshold == minProteinLabel), "MockMeanPositivity"] <- mean((TempPositiveTable_plyr[grep("_01-Mock_", TempPositiveTable_plyr$UniqueConditionTag), "PositivePercent" ]))
      
     }
  }
  
  
  # Add suggestion column 
  possibleProteinThresholdsDF$Suggestion <- "NO"
  #
  CuttedTemp <- possibleProteinThresholdsDF[which(possibleProteinThresholdsDF$MockMeanPositivity < 5), ]
  #  
  possibleProteinThresholdsDF[rownames(CuttedTemp)[which.min(CuttedTemp$Threshold)[1]], "Suggestion"] <- "YES"
  #
  remove(CuttedTemp)
  #
  if("-Inf" %in% row.names(possibleProteinThresholdsDF) ){
    possibleProteinThresholdsDF <- possibleProteinThresholdsDF[row.names(possibleProteinThresholdsDF) != "-Inf",]
  }
  
  #
  possibleProteinThresholdsDF[is.na(possibleProteinThresholdsDF$MockMeanPositivity), "MockMeanPositivity"] <- 0
  
  # Salva a tabela
  write.csv2(possibleProteinThresholdsDF,
             file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "PossibleThresholds", "\\", "SummaryTable", ".csv"))
  
  ######

  # Use protein stain to define positive mitochondria
  # Thresholding
  if(
    length(possibleProteinThresholdsDF$Threshold[possibleProteinThresholdsDF$Suggestion == "YES"]) == 1
  ){
    # AQUI EU PUXO A MARCAÇÃO SUGERIDA, ALTERAR SE NAO CONCORDAR
    minProteinLabel <- possibleProteinThresholdsDF$Threshold[possibleProteinThresholdsDF$Suggestion == "YES"]
    NumberPositiveMito <- possibleProteinThresholdsDF$MitoNumber[possibleProteinThresholdsDF$Suggestion == "YES"]
  } else {
    #
    minProteinLabel <- 310
    # Numer mito test
    NumberPositiveMito <- 10
  }

  
  
  ### MANUAL DEFINITION
  minProteinLabel <- 600
  # Numer mito test
  NumberPositiveMito <- 1
  ###
  
  
  
  
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "06_ProteinLabel_superpose.png"),
      width=1500, height=1000, res=300)
  hist(FinalResults[which(FinalResults$Info01 == "01-Mock"), "Protein_Mean"],
       breaks=200,
       main="Stain Negatives",
       xlab="Mitochondrial GAC stain",
       ylab="Amount",
       ylim=c(0,0.01),
       xlim=c(0, 1200),
       freq=FALSE,
       border="white", 
       col=rgb(1,0,0,1/4))
  abline(v=minProteinLabel, lwd=2)
  
  hist(FinalResults[which(FinalResults$Info01 != "01-Mock"), "Protein_Mean"],
       breaks=2000,
       main="Stain Positives",
       xlab="Mitochondrial GAC stain",
       ylab="Amount",
       ylim=c(0,0.01),
       xlim=c(0, 1200),
       freq=FALSE,
       border="white", 
       col=rgb(0,0,1,1/4),
       add=TRUE)
  
  dev.off()
  
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "07_ProteinLabel_positives.png"),
      width=1500, height=1000, res=300)
  hist()
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Limpa a marcação atual
  FinalResults$PositiveCell <- FALSE
  FinalResults$PositiveMitochondria <- FALSE
  
  # Marca as mitocôndrias positivas
  FinalResults$PositiveMitochondria[FinalResults$Protein_Mean > minProteinLabel] <- TRUE
  
  # Marca as CÉLULAS positivas
  # Tentativa com DPLYR para agilizar
  TempPositiveTable_plyr_toLabel <- FinalResults %>%
    group_by(FileOrigin_base) %>%
    summarize(positiveMito = sum(PositiveMitochondria, na.rm=TRUE),
              totalMito = n()
    )
  FinalResults[(FinalResults$FileOrigin_base %in% TempPositiveTable_plyr_toLabel$FileOrigin_base[TempPositiveTable_plyr_toLabel$positiveMito > NumberPositiveMito]), "PositiveCell"] <- TRUE
  remove(TempPositiveTable_plyr_toLabel)
  # Marca as CÉLULAS positivas
  #for(i in unique(FinalResults$FileOrigin_base)){
  #  if(sum(FinalResults[FinalResults$FileOrigin_base == i, "PositiveMitochondria"]) > NumberPositiveMito){
  #    FinalResults[FinalResults$FileOrigin_base == i, "PositiveCell"] <- TRUE
  #  } else {
  #  }
  #}
  
  
  # Plot conditions
  TempPositiveTable <- unique(FinalResults[,c("UniqueConditionTag", "UniqueTag", "FileOrigin_base", "PositiveCell")])
  TempPositiveTable$PositiveCell[TempPositiveTable$PositiveCell] <- 1
  TempPositiveTable$PositiveCell[TempPositiveTable$PositiveCell == FALSE] <- 0
  
  # Dplyr para gerar os valores
  TempPositiveTable_plyr <- TempPositiveTable %>%
    group_by(UniqueTag, UniqueConditionTag) %>%
    dplyr::summarize(PositiveSum = sum(PositiveCell, na.rm=TRUE),
                     TotalCell = n()
    )
  
  TempPositiveTable_plyr$PositiveSum/TempPositiveTable_plyr$TotalCell -> TempPositiveTable_plyr$PositivePercent
  #
  TempPositiveTable_plyr$PositivePercent * 100 -> TempPositiveTable_plyr$PositivePercent
  
  head(TempPositiveTable_plyr)
  myComparisons <- list()
  k <- 0
  for(i in 1:length(unique(TempPositiveTable_plyr$UniqueConditionTag))){
    for(j in 1:length(unique(TempPositiveTable_plyr$UniqueConditionTag))){
      if(j>i){
        k<- k+1
        myComparisons[[k]] <- c(unique(TempPositiveTable_plyr$UniqueConditionTag)[i],
                                unique(TempPositiveTable_plyr$UniqueConditionTag)[j])
      }
    }
  }
  
  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "08_PercentPositivity.png"),
      width=1500, height=1500, res=200)
  
  #bottom, left, top and right
  # par(mar=c(5.1,20.1,4.1,2.1))
  print(
    ggbarplot(TempPositiveTable_plyr,
              x = "UniqueConditionTag",
              y = "PositivePercent", 
              add = c("mean_se", "jitter")) +
      rotate_x_text(90) +
      
      stat_compare_means() # comparisons=myComparisons
  )
  dev.off()
  
  #bottom, left, top and right
  #par(mar=c(5.1,4.1,4.1,2.1))
  
  # + # Add pairwise comparisons p-value comparisons = my_comparisons
  #stat_compare_means(label.y = 50)
  
  
  
  
  
  #### FRACTION POSITIVE MITOCHONDRIA WITHIN POSITIVE CELLS
  FinalResults$FractionMitoPositive <- NA
  
  ##### Find positive mitochondria
  TempTable_plyr_toLabel <- FinalResults[which(FinalResults$PositiveCell), ] %>%
    group_by(FileOrigin_base, UniqueConditionTag) %>%
    count(PositiveMitochondria)
  TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)

  #  
  TempTable_plyr_toLabel$PositiveMitochondria[TempTable_plyr_toLabel$PositiveMitochondria == "TRUE"] <- "Positive"
  #
  TempTable_plyr_toLabel$PositiveMitochondria[TempTable_plyr_toLabel$PositiveMitochondria == "FALSE"] <- "Negative"

  #
  for(BaseNow in unique(FinalResults$FileOrigin_base)){
  #
    if(length(which(TempTable_plyr_toLabel$FileOrigin_base == BaseNow)) > 0){
      
      TempLineNow <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == BaseNow), ]
      
      if(sum(TempLineNow$PositiveMitochondria == "Positive")>0){
        #
        a <- TempLineNow[which(TempLineNow$PositiveMitochondria == "Positive"), "n"]
      } else {
        a <- 0
      }
      
      if(sum(TempLineNow$PositiveMitochondria == "Negative")>0){
        #
        b <- (a+(TempLineNow[which(TempLineNow$PositiveMitochondria == "Negative"), "n"]))
      } else {
        b <- a
      }
      #      
      (a/b)*100 -> FinalResults[which(FinalResults$FileOrigin_base == BaseNow), "FractionMitoPositive"]
    }
    #
  }
  
  remove(TempTable_plyr_toLabel)
  
  
  
  #######
  
  #
  write.csv2(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "12_PositiveMitochondriaPerCondition", ".csv"),
             unique(FinalResults[!is.na(FinalResults$FractionMitoPositive), c("UniqueConditionTag", "FractionMitoPositive")]) )
  
  
  
  p <- ggplot(data=unique(FinalResults[!is.na(FinalResults$FractionMitoPositive), c("UniqueConditionTag", "FractionMitoPositive")]),
              aes(x=UniqueConditionTag, y=FractionMitoPositive)) + 
    geom_boxplot() + theme_classic() + geom_jitter(shape=16,
                  position=position_jitter(0.2),
                  alpha=0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

  png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "12_PositiveMitochondriaPerCondition", ".png"),
      width=2300, height=1500, res=250)
  print(p)
  dev.off()
  #
  
  
  #
  
  

  
  
  #### D MITO BY POSITIVITY
  
  head(FinalResults)
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  
  for(i in unique(FinalResults$UniqueConditionTag) ){
    
    temp <- FinalResults[FinalResults$UniqueConditionTag == i,]
    
    
    
    p <- ggboxplot(temp, x = "PositiveMitochondria", y = "D",
                   # color = "supp", palette = "jco",
                   add = "jitter",
                   color="gray50",
                   title=paste0(i,"\nRatio positive/negative ",
                                mean(temp[temp$PositiveMitochondria, "D"])/
                                  mean(temp[!temp$PositiveMitochondria, "D"]) )) +
      stat_compare_means()
    
    png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "09_D_label_", i, ".png"),
        width=1000, height=1500, res=200)
    print(p)
    dev.off()
    
    
    
    
    
    
    
    
    
    p <- ggboxplot(temp, x = "PositiveMitochondria", y = "Perim.",
                   # color = "supp", palette = "jco",
                   add = "jitter",
                   color="gray50",
                   title=paste0(i,"\nRatio positive/negative ",
                                mean(temp[temp$PositiveMitochondria, "Perim."])/
                                  mean(temp[!temp$PositiveMitochondria, "Perim."]) )) +
      stat_compare_means()
    
    png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "09_Perim_label_", i, ".png"),
        width=1000, height=1500, res=200)
    print(p)
    dev.off()
    
  }
  
  
  
  
  
  
  
  
  
  #### D MITO BY GLN
  
  for(i in unique(FinalResults$Info01) ){
    
    temp <- FinalResults[(FinalResults$Info01 == i &
                            FinalResults$PositiveMitochondria),]
    
    
    if(dim(temp)[1] > 1){
      
      
      temp$Division <- sapply(1:dim(temp)[1], function(y){paste0(temp[y, c("Info02", "Info03")], collapse="_")  })
      
      
      #
      myComparisons <- list()
      k <- 0
      for(l in 1:length(unique(temp$Division))){
        for(j in 1:length(unique(temp$Division))){
          if(j>l){
            k<- k+1
            myComparisons[[k]] <- c(unique(temp$Division)[l],
                                    unique(temp$Division)[j])
          }
        }
      }
      
      p <- ggboxplot(temp, x = "Division", y = "D",
                     # color = "supp", palette = "jco",
                     add = "jitter",
                     color="gray50",
                     title=paste0(i,"\nOnlyPositiveMito "
                     )) +
        stat_compare_means(comparisons=myComparisons) +
        rotate_x_text(90)
      
      png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "10_D_Conditions_", i, ".png"),
          width=1000, height=1500, res=200)
      print(p)
      dev.off()
      
      
      
      
      
      
      
      p <- ggboxplot(temp, x = "Division", y = "Perim.",
                     # color = "supp", palette = "jco",
                     add = "jitter",
                     color="gray50",
                     title=paste0(i,"\nOnlyPositiveMito "
                     )) +
        stat_compare_means(comparisons=myComparisons) +
        rotate_x_text(90)
      
      png(paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "10_Perim_Conditions_", i, ".png"),
          width=1000, height=1500, res=200)
      print(p)
      dev.off()
      
    }
    
    
  }
  

  
  #### Tamanho mito vs intensidade GLS
  for(i in unique(FinalResults$UniqueConditionTag) ){
    
    temp <- FinalResults[(FinalResults$UniqueConditionTag == i &
                            FinalResults$PositiveMitochondria),]
    
    
    if(dim(temp)[1] > 10){
      
      
      
      x_now <- temp$Protein_Mean
      y_now <- temp$D
      
      # Correlation plot
      png(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "11_Correlation_D_", i, ".png"),
          width=2000, height=2000, res=400)
      plot(x=x_now, y=y_now, pch=16,
           col=rgb(37,37,37, 10, max=255),
           cex=0.75, cex.main=0.75,
           main=paste0("pearson cor = ", cor(x=x_now, y=y_now), " p ", cor.test(x=x_now, y=y_now, method="pearson")[["p.value"]]),
           xlab=paste0(i),
           ylab="Fractal D Value",
           cex.axis=1.3, yaxt="n", xaxt="n")
      axis(1, col.axis="black", las=2, lwd=2, cex.axis=1.25)
      axis(2, col.axis="black", las=2, lwd=2, cex.axis=1.25)
      box(lwd=2)
      abline(lm(y_now~x_now), lwd=2, lty=3, col=rgb(12,12,12,250, max=255))
      
      dev.off()
      
      
      
      
      
      
      
      
      x_now <- temp$Protein_Mean
      y_now <- temp$Perim.
      
      # Correlation plot
      png(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "11_Correlation_Perim_", i, ".png"),
          width=2000, height=2000, res=400)
      plot(x=x_now, y=y_now, pch=16,
           col=rgb(37,37,37, 10, max=255),
           cex=0.75, cex.main=0.75,
           main=paste0("pearson cor = ", cor(x=x_now, y=y_now), " p ", cor.test(x=x_now, y=y_now, method="pearson")[["p.value"]]),
           xlab=paste0(i),
           ylab="Fractal D Value",
           cex.axis=1.3, yaxt="n", xaxt="n")
      axis(1, col.axis="black", las=2, lwd=2, cex.axis=1.25)
      axis(2, col.axis="black", las=2, lwd=2, cex.axis=1.25)
      box(lwd=2)
      abline(lm(y_now~x_now), lwd=2, lty=3, col=rgb(12,12,12,250, max=255))
      
      dev.off()
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####
  #### D-thresholding test ONLY POSITIVE MITOCONDRIA
  #### START
  ####
  

  
  
  minimumDs <- seq(from=0.67,
                   to=0.70,
                   by=0.01)
  
  maximumDs <- seq(from=0.95,
                   to=0.98,    #round(x=max(FinalResults$D[FinalResults$PositiveCell], na.rm=TRUE)-0.05, digits = 2),
                   by=0.01)
  
  
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByD", "\\"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByPerimeter", "\\"), showWarnings = FALSE)
  
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByD", "\\", "ByMito", "\\"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByD", "\\", "ByCell", "\\"), showWarnings = FALSE)
  
  # Ordena a tabela
  FinalResults <- FinalResults[order(FinalResults$FileOrigin_base), ]
  

  
  ### Create multiple single tables in a LIST
  # Clean mitoCats
  # Add columns for categories
  FinalResults$MitoCat <- NA
  FinalResults$CellCat <- NA

  usefullColumns <- c("FileOrigin_base", "D", "Perim.", "PositiveMitochondria", "MitoCat", "FractionMitoPositive")
  
  if(exists("TukeyTableForEvaluation_byCellOnly_DScores")){
    remove(TukeyTableForEvaluation_byCellOnly_DScores)
  }
  if(exists("TukeyTableForEvaluation_byMitoOnly_DScores")){
    remove(TukeyTableForEvaluation_byMitoOnly_DScores)
  }
  
  
  # Loop para os diversos cortes de D
  for(minDNow in minimumDs){
    for(maxDNow in maximumDs){
      if(maxDNow > minDNow){
        
        
        ####
        # Clean mitoCats
        # Add columns for categories
        FinalResults$MitoCat <- NA
        FinalResults$CellCat <- NA
        #FinalResults$ <- NA  
        
        # Add intermediate
        FinalResults[which(
          !is.na(FinalResults$D)
        ),"MitoCat"] <- "2-Intermediate"

        # Add fragmented
        FinalResults[which(
          !is.na(FinalResults$D) &
            FinalResults$D < minDNow
        ),"MitoCat"] <- "1-Fragmented"
        
        # Add tubular
        FinalResults[which(
          !is.na(FinalResults$D) &
            FinalResults$D > maxDNow
        ),"MitoCat"] <- "3-Tubular"
        
        
        #####
        
        #####
        
        
        
        
        
        ##### Add CellCat for:
        # Cells labelled with "01-Mock" will be FULLY SELECTED
        # Cells labelled otherwise:
        # Only FinalResults$PositiveMitochondria TRUE
        #
        TempTable_plyr_toLabel <- FinalResults[which(
          FinalResults$Info01 == "01-Mock" | FinalResults$PositiveMitochondria
        ), ] %>%
          group_by(FileOrigin_base) %>%
          count(MitoCat)
        TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)

        ###
        FinalResultsSplitted <- split(FinalResults[,usefullColumns],
                                      f = FinalResults[,"FileOrigin_base"])
        
        
        
        
        #FinalResultsSplitted[6102] -> myTableNow
        # define function
        CellularClassificationFunction <- function(myTableNow){
          #print(myTableNow$FileOrigin_base[1])
          if(length(which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]))>0){
            myTableNow$CellCat <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), "n"])]
          } else {
          }
          return(myTableNow)
        }
        
        # Create the dataframe list
        system.time({
          myFinalResultsSplitted <- lapply(X=FinalResultsSplitted,
                                           FUN=CellularClassificationFunction)
        })
        
        #
        myFinalResultsSplittedjoined <- bind_rows(myFinalResultsSplitted) #, .id = "column_label"
        
        #
        # Ordena a tabela
        myFinalResultsSplittedjoined <- myFinalResultsSplittedjoined[order(myFinalResultsSplittedjoined$FileOrigin_base), ]
        
        # Return data
        FinalResults[,c("MitoCat", "CellCat")] <- myFinalResultsSplittedjoined[,c("MitoCat", "CellCat")]
        
        #
        remove(FinalResultsSplitted)
        remove(myFinalResultsSplitted)
        remove(myFinalResultsSplittedjoined)
        remove(TempTable_plyr_toLabel)
        
        
        #######
        
        
        #system.time({
        
        #p <-0
        #for(CellNow in unique(TempTable_plyr_toLabel$FileOrigin_base)){
        #  p<-p+1
        #  print(p/length(unique(TempTable_plyr_toLabel$FileOrigin_base)))
        #  FinalResults[which(FinalResults$FileOrigin_base == CellNow), "CellCat"] <-
        #    TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == CellNow), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == CellNow), "n"])]
        #}
        
        # })
        #FinalResults[, "CellCat"] <- unname( unlist(sapply(unique(TempTable_plyr_toLabel$FileOrigin_base),
        #                                                   function(w){
        #                                                     rep(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == w), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == w), "n"])], length(which(FinalResults$FileOrigin_base == w)))
        #                                                   })))
        #
        # Marca as CÉLULAS positivas
        
        #for(CellNow in unique(FinalResults$FileOrigin_base)){
        #  FinalResults[which(FinalResults$FileOrigin_base == CellNow), "CellCat"] <-
        #    names(which.max(table(FinalResults[which(FinalResults$FileOrigin_base == CellNow), "MitoCat"])))
        #}
        
        
        
        # table(FinalResults$CellCat)
        
        #
        # CreateTempTable
        # Only transfected cells for non-mock conditions
        #
        TempTable <- FinalResults[which(
          FinalResults$Info01 == "01-Mock" | FinalResults$PositiveMitochondria
        ), ]
        
        
        
        #
        theme_set(theme_pubclean())
        
        # TempTable$MitoCat <- as.factor(TempTable$MitoCat)
        # TempTable$CellCat <- as.factor(TempTable$CellCat)
        # unique(TempTable[,c("UniqueTag", "UniqueConditionTag", "CellCat")])
        
        ############## BY CELL
        df.summary <- unique(TempTable[,c("UniqueTag", "UniqueConditionTag", "CellCat", "FileOrigin_base")]) %>%
          group_by(UniqueTag, UniqueConditionTag) %>%
          count(CellCat)
        df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
        
        
        df.summary$Percentage <- NA
        
        # Add zeroes
        for(condNow in unique(df.summary$UniqueTag)){

          if("1-Fragmented" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"CellCat"] <- "1-Fragmented"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          if("2-Intermediate" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"CellCat"] <- "2-Intermediate"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          
          if("3-Tubular" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"CellCat"] <- "3-Tubular"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          
          
          #
          # Add percentage
          df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
            sum(df.summary[df.summary$UniqueTag==condNow, "n"])
          
        }
        
        
        df.summary2 <- df.summary %>%
          group_by(UniqueConditionTag, CellCat) %>%
          summarise(mean = mean(Percentage, na.rm=TRUE),
                    sd = sd(Percentage, na.rm = TRUE))
        
        #df.summary2
        
        # (2) Bar plots of means + individual jitter points + errors
        png(file=paste0(SelectedFolder, "\\",
                        "11-SelectedMitochondriaOnly",
                        "\\",
                        "ByD",
                        "\\", "ByCell", "\\", "ByCellCat",
                        "_MinD_", stringr::str_pad(minDNow, nchar(as.character(max(minimumDs))), pad = "0"), "_MaxD_", stringr::str_pad(maxDNow, nchar(as.character(maximumDs)), pad = "0"), ".png"),
            width = 4600, height = 2200,
            res=300)
        
        print(
          ggplot(data = df.summary2,
                 aes(x=UniqueConditionTag,
                     y=mean,
                     colour = CellCat,
                     fill = CellCat)) +
            
            geom_bar(stat = "identity",
                     position=position_dodge()) +
            geom_errorbar(
              aes(ymin = mean-sd,
                  ymax = mean+sd),
              colour="black",
              width = 0.5,
              position = position_dodge(.9)) +
            geom_jitter(data=df.summary,
                        aes(x=UniqueConditionTag,
                            y=Percentage,
                            fill = CellCat,
                            colour = CellCat),
                        colour="black",
                        alpha=0.3,
                        position = position_jitterdodge(0.1)) +
            scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
        )
        dev.off()
        
        
        ### STAT TEST
        
        
        #
        for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
          
          # Make anova
          aov_now <- aov(formula=Percentage~UniqueConditionTag,
                         data=df.summary[df.summary[,"CellCat"] == mitoCatNow, ])
          
          tukey_now <- TukeyHSD(x=aov_now)
          
          
          
          
          if(exists("TukeyTableForEvaluation_byCellOnly_DScores")){
            
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinD_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumDs))), pad = "0"),
                                                                "_MaxD_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumDs))), pad = "0"))
            # add the column
            TukeyTableForEvaluation_byCellOnly_DScores <- cbind(TukeyTableForEvaluation_byCellOnly_DScores, tukey_now$UniqueConditionTag[rownames(TukeyTableForEvaluation_byCellOnly_DScores), 4,drop=FALSE])
            
          } else {
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinD_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumDs))), pad = "0"),
                                                                "_MaxD_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumDs))), pad = "0"))
            # Just copy
            TukeyTableForEvaluation_byCellOnly_DScores <- tukey_now$UniqueConditionTag[, 4,drop=FALSE]
          }
          
          
          
        }
        
        
        
        
        ########### BY MITOCHONDRIA
        
        
        
        df.summary <- TempTable[!is.na(TempTable$MitoCat),] %>%
          group_by(UniqueTag, UniqueConditionTag) %>%
          count(MitoCat)
        df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
        
        
        df.summary$Percentage <- NA
        
        # Add zeroes
        for(condNow in unique(df.summary$UniqueTag)){
          
          
          if("1-Fragmented" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"MitoCat"] <- "1-Fragmented"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          if("2-Intermediate" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"MitoCat"] <- "2-Intermediate"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          
          if("3-Tubular" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"MitoCat"] <- "3-Tubular"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          #
          # Add percentage
          df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
            sum(df.summary[df.summary$UniqueTag==condNow, "n"])
          
        }
        
        
        df.summary2 <- df.summary %>%
          group_by(UniqueConditionTag, MitoCat) %>%
          summarise(mean = mean(Percentage, na.rm=TRUE),
                    sd = sd(Percentage, na.rm = TRUE))
        
        #df.summary2
        
        #####
        
        # (2) Bar plots of means + individual jitter points + errors
        png(file=paste0(SelectedFolder, "\\",
                        "11-SelectedMitochondriaOnly",
                        "\\",
                        "ByD",
                        "\\", "ByMito", "\\", "ByMitoCat",
                        "_MinD_", stringr::str_pad(minDNow, nchar(as.character(max(minimumDs))), pad = "0"), "_MaxD_", stringr::str_pad(maxDNow, nchar(as.character(max(maximumDs))), pad = "0"), ".png"),
            width = 4600, height = 2200,
            res=300)
        
        print(
          ggplot(data = df.summary2,
                 aes(x=UniqueConditionTag,
                     y=mean,
                     colour = MitoCat,
                     fill = MitoCat)) +
            
            geom_bar(stat = "identity",
                     position=position_dodge()) +
            geom_errorbar(
              aes(ymin = mean-sd,
                  ymax = mean+sd),
              colour="black",
              width = 0.5,
              position = position_dodge(.9)) +
            geom_jitter(data=df.summary,
                        aes(x=UniqueConditionTag,
                            y=Percentage,
                            fill = MitoCat,
                            colour = MitoCat),
                        colour="black",
                        alpha=0.3,
                        position = position_jitterdodge(0.1)) +
            scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
        )
        dev.off()
        
        
        ### STAT TEST
        
        
        #
        for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
          
          # Make anova
          aov_now <- aov(formula=Percentage~UniqueConditionTag,
                         data=df.summary[df.summary[,"MitoCat"] == mitoCatNow, ])
          
          tukey_now <- TukeyHSD(x=aov_now)
          
          
          
          
          if(exists("TukeyTableForEvaluation_byMitoOnly_DScores")){
            
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinD_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumDs))), pad = "0"),
                                                                "_MaxD_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumDs))), pad = "0"))
            # add the column
            TukeyTableForEvaluation_byMitoOnly_DScores <- cbind(TukeyTableForEvaluation_byMitoOnly_DScores, tukey_now$UniqueConditionTag[rownames(TukeyTableForEvaluation_byMitoOnly_DScores), 4,drop=FALSE])
            
          } else {
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinD_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumDs))), pad = "0"),
                                                                "_MaxD_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumDs))), pad = "0"))
            # Just copy
            TukeyTableForEvaluation_byMitoOnly_DScores <- tukey_now$UniqueConditionTag[, 4,drop=FALSE]
          }
          
          
          
        }
        
      }
    }
  }
  
  
  # Save byCell
  
  #
  for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
    
    write.csv2(file=paste0(SelectedFolder, "\\",
                           "11-SelectedMitochondriaOnly",
                           "\\",
                           "ByD",
                           "\\", "ByCellOnly_",
                           mitoCatNow,".csv"),
               TukeyTableForEvaluation_byCellOnly_DScores[ , grepl(mitoCatNow, colnames(TukeyTableForEvaluation_byCellOnly_DScores)), drop=FALSE] )
    
    write.csv2(file=paste0(SelectedFolder, "\\",
                           "11-SelectedMitochondriaOnly",
                           "\\",
                           "ByD",
                           "\\", "ByMitoOnly_",
                           mitoCatNow,".csv"),
               TukeyTableForEvaluation_byMitoOnly_DScores[ , grepl(mitoCatNow, colnames(TukeyTableForEvaluation_byMitoOnly_DScores)), drop=FALSE] )
    
  }
  
  
  
  ####
  #### D-thresholding test ONLY POSITIVE MITOCONDRIA
  #### END
  ####
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####
  #### Perimeter test ONLY POSITIVE MITOCONDRIA
  #### START
  ####
  
  
  
  
  minimumPerimeter <- seq(from=5.3,
                          to=5.5,
                          by=0.1)
  
  maximumPerimeter <- seq(from=8.6,
                          to=9.1,
                          by=0.1)
  
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByPerimeter", "\\"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByPerimeter", "\\", "ByMito", "\\"), showWarnings = FALSE)
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\", "ByPerimeter", "\\", "ByCell", "\\"), showWarnings = FALSE)
  
  # Ordena a tabela
  FinalResults <- FinalResults[order(FinalResults$FileOrigin_base), ]
  
  
  
  ### Create multiple single tables in a LIST
  # Clean mitoCats
  # Add columns for categories
  FinalResults$MitoCat <- NA
  FinalResults$CellCat <- NA
  
  usefullColumns <- c("FileOrigin_base", "D", "Perim.", "PositiveMitochondria", "MitoCat", "FractionMitoPositive")
  
  if(exists("TukeyTableForEvaluation_byCellOnly_Perimeter")){
    remove(TukeyTableForEvaluation_byCellOnly_Perimeter)
  }
  if(exists("TukeyTableForEvaluation_byMitoOnly_Perimeter")){
    remove(TukeyTableForEvaluation_byMitoOnly_Perimeter)
  }
  
  
  # Loop para os diversos cortes de D
  for(minDNow in minimumPerimeter){
    for(maxDNow in maximumPerimeter){
      if(maxDNow > minDNow){
        
        
        ####
        # Clean mitoCats
        # Add columns for categories
        FinalResults$MitoCat <- NA
        FinalResults$CellCat <- NA
        #FinalResults$ <- NA  
        
        
        # Add intermediate
        FinalResults[which(
          !is.na(FinalResults$Perim.)
        ),"MitoCat"] <- "2-Intermediate"
        
        
        
        # Add fragmented
        FinalResults[which(
          !is.na(FinalResults$Perim.) &
            FinalResults$Perim. < minDNow
        ),"MitoCat"] <- "1-Fragmented"
        
        # Add tubular
        FinalResults[which(
          !is.na(FinalResults$Perim.) &
            FinalResults$Perim. > maxDNow
        ),"MitoCat"] <- "3-Tubular"
        
        
        
        
        #####
        
        #####
        
        
        
        
        
        ##### Add CellCat for:
        # Cells labelled with "01-Mock" will be FULLY SELECTED
        # Cells labelled otherwise:
        # Only FinalResults$PositiveMitochondria TRUE
        #
        TempTable_plyr_toLabel <- FinalResults[which(
          FinalResults$Info01 == "01-Mock" | FinalResults$PositiveMitochondria & FinalResults$PositiveCell
        ), ] %>%
          group_by(FileOrigin_base) %>%
          count(MitoCat)
        TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
        
        ###
        FinalResultsSplitted <- split(FinalResults[,usefullColumns],
                                      f = FinalResults[,"FileOrigin_base"])
        
        
        
        
        #FinalResultsSplitted[6102] -> myTableNow
        # define function
        CellularClassificationFunction <- function(myTableNow){
          #print(myTableNow$FileOrigin_base[1])
          if(length(which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]))>0){
            myTableNow$CellCat <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), "n"])]
          } else {
          }
          return(myTableNow)
        }
        
        # Create the dataframe list
        system.time({
          myFinalResultsSplitted <- lapply(X=FinalResultsSplitted,
                                           FUN=CellularClassificationFunction)
        })
        
        #
        myFinalResultsSplittedjoined <- bind_rows(myFinalResultsSplitted) #, .id = "column_label"
        
        #
        # Ordena a tabela
        myFinalResultsSplittedjoined <- myFinalResultsSplittedjoined[order(myFinalResultsSplittedjoined$FileOrigin_base), ]
        
        # Return data
        FinalResults[,c("MitoCat", "CellCat")] <- myFinalResultsSplittedjoined[,c("MitoCat", "CellCat")]
        
        #
        remove(FinalResultsSplitted)
        remove(myFinalResultsSplitted)
        remove(myFinalResultsSplittedjoined)
        remove(TempTable_plyr_toLabel)
        
        
        #######
        
        
        #system.time({
        
        #p <-0
        #for(CellNow in unique(TempTable_plyr_toLabel$FileOrigin_base)){
        #  p<-p+1
        #  print(p/length(unique(TempTable_plyr_toLabel$FileOrigin_base)))
        #  FinalResults[which(FinalResults$FileOrigin_base == CellNow), "CellCat"] <-
        #    TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == CellNow), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == CellNow), "n"])]
        #}
        
        # })
        #FinalResults[, "CellCat"] <- unname( unlist(sapply(unique(TempTable_plyr_toLabel$FileOrigin_base),
        #                                                   function(w){
        #                                                     rep(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == w), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == w), "n"])], length(which(FinalResults$FileOrigin_base == w)))
        #                                                   })))
        #
        # Marca as CÉLULAS positivas
        
        #for(CellNow in unique(FinalResults$FileOrigin_base)){
        #  FinalResults[which(FinalResults$FileOrigin_base == CellNow), "CellCat"] <-
        #    names(which.max(table(FinalResults[which(FinalResults$FileOrigin_base == CellNow), "MitoCat"])))
        #}
        
        
        
        # table(FinalResults$CellCat)
        
        #
        # CreateTempTable
        # Only transfected cells for non-mock conditions
        #
        TempTable <- FinalResults[which(
          FinalResults$Info01 == "01-Mock" | FinalResults$PositiveMitochondria & FinalResults$PositiveCell
        ), ]
        
        
        
        #
        theme_set(theme_pubclean())
        
        # TempTable$MitoCat <- as.factor(TempTable$MitoCat)
        # TempTable$CellCat <- as.factor(TempTable$CellCat)
        # unique(TempTable[,c("UniqueTag", "UniqueConditionTag", "CellCat")])
        
        ############## BY CELL
        df.summary <- unique(TempTable[,c("UniqueTag", "UniqueConditionTag", "CellCat", "FileOrigin_base")]) %>%
          group_by(UniqueTag, UniqueConditionTag) %>%
          count(CellCat)
        df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
        
        
        df.summary$Percentage <- NA
        
        # Add zeroes
        for(condNow in unique(df.summary$UniqueTag)){
          
          if("1-Fragmented" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"CellCat"] <- "1-Fragmented"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          if("2-Intermediate" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"CellCat"] <- "2-Intermediate"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          
          if("3-Tubular" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"CellCat"] <- "3-Tubular"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          
          
          #
          # Add percentage
          df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
            sum(df.summary[df.summary$UniqueTag==condNow, "n"])
          
        }
        
        
        df.summary2 <- df.summary %>%
          group_by(UniqueConditionTag, CellCat) %>%
          summarise(mean = mean(Percentage, na.rm=TRUE),
                    sd = sd(Percentage, na.rm = TRUE))
        
        #df.summary2
        
        # (2) Bar plots of means + individual jitter points + errors
        png(file=paste0(SelectedFolder, "\\",
                        "11-SelectedMitochondriaOnly",
                        "\\",
                        "ByPerimeter",
                        "\\", "ByCell", "\\", "ByCellCat",
                        "_MinPerimeter_", stringr::str_pad(minDNow, nchar(as.character(max(minimumPerimeter))), pad = "0"), "_MaxPerimeter_", stringr::str_pad(maxDNow, nchar(as.character(maximumPerimeter)), pad = "0"), ".png"),
            width = 4600, height = 2200,
            res=300)
        
        print(
          ggplot(data = df.summary2,
                 aes(x=UniqueConditionTag,
                     y=mean,
                     colour = CellCat,
                     fill = CellCat)) +
            
            geom_bar(stat = "identity",
                     position=position_dodge()) +
            geom_errorbar(
              aes(ymin = mean-sd,
                  ymax = mean+sd),
              colour="black",
              width = 0.5,
              position = position_dodge(.9)) +
            geom_jitter(data=df.summary,
                        aes(x=UniqueConditionTag,
                            y=Percentage,
                            fill = CellCat,
                            colour = CellCat),
                        colour="black",
                        alpha=0.3,
                        position = position_jitterdodge(0.1)) +
            scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
        )
        dev.off()
        
        
        ### STAT TEST
        
        
        #
        for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
          
          # Make anova
          aov_now <- aov(formula=Percentage~UniqueConditionTag,
                         data=df.summary[df.summary[,"CellCat"] == mitoCatNow, ])
          
          tukey_now <- TukeyHSD(x=aov_now)
          
          
          
          
          if(exists("TukeyTableForEvaluation_byCellOnly_Perimeter")){
            
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinPerimeter_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumPerimeter))), pad = "0"),
                                                                "_MaxPerimeter_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumPerimeter))), pad = "0"))
            # add the column
            TukeyTableForEvaluation_byCellOnly_Perimeter <- cbind(TukeyTableForEvaluation_byCellOnly_Perimeter, tukey_now$UniqueConditionTag[rownames(TukeyTableForEvaluation_byCellOnly_Perimeter), 4,drop=FALSE])
            
          } else {
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinPerimeter_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumPerimeter))), pad = "0"),
                                                                "_MaxPerimeter_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumPerimeter))), pad = "0"))
            # Just copy
            TukeyTableForEvaluation_byCellOnly_Perimeter <- tukey_now$UniqueConditionTag[, 4,drop=FALSE]
          }
          
          
          
        }
        
        
        
        
        ########### BY MITOCHONDRIA
        
        
        
        df.summary <- TempTable[!is.na(TempTable$MitoCat),] %>%
          group_by(UniqueTag, UniqueConditionTag) %>%
          count(MitoCat)
        df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
        
        
        df.summary$Percentage <- NA
        
        # Add zeroes
        for(condNow in unique(df.summary$UniqueTag)){
          
          
          if("1-Fragmented" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"MitoCat"] <- "1-Fragmented"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          if("2-Intermediate" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"MitoCat"] <- "2-Intermediate"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          
          
          if("3-Tubular" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
            #print("tem")
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,"MitoCat"] <- "3-Tubular"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          
          #
          # Add percentage
          df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
            sum(df.summary[df.summary$UniqueTag==condNow, "n"])
          
        }
        
        
        df.summary2 <- df.summary %>%
          group_by(UniqueConditionTag, MitoCat) %>%
          summarise(mean = mean(Percentage, na.rm=TRUE),
                    sd = sd(Percentage, na.rm = TRUE))
        
        #df.summary2
        
        #####
        
        # (2) Bar plots of means + individual jitter points + errors
        png(file=paste0(SelectedFolder, "\\",
                        "11-SelectedMitochondriaOnly",
                        "\\",
                        "ByPerimeter",
                        "\\", "ByMito", "\\", "ByMitoCat",
                        "_MinPerimeter_", stringr::str_pad(minDNow, nchar(as.character(max(minimumPerimeter))), pad = "0"), "_MaxPerimeter_", stringr::str_pad(maxDNow, nchar(as.character(max(maximumPerimeter))), pad = "0"), ".png"),
            width = 4600, height = 2200,
            res=300)
        
        print(
          ggplot(data = df.summary2,
                 aes(x=UniqueConditionTag,
                     y=mean,
                     colour = MitoCat,
                     fill = MitoCat)) +
            
            geom_bar(stat = "identity",
                     position=position_dodge()) +
            geom_errorbar(
              aes(ymin = mean-sd,
                  ymax = mean+sd),
              colour="black",
              width = 0.5,
              position = position_dodge(.9)) +
            geom_jitter(data=df.summary,
                        aes(x=UniqueConditionTag,
                            y=Percentage,
                            fill = MitoCat,
                            colour = MitoCat),
                        colour="black",
                        alpha=0.3,
                        position = position_jitterdodge(0.1)) +
            scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
        )
        dev.off()
        
        
        ### STAT TEST
        
        
        #
        for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
          
          # Make anova
          aov_now <- aov(formula=Percentage~UniqueConditionTag,
                         data=df.summary[df.summary[,"MitoCat"] == mitoCatNow, ])
          
          tukey_now <- TukeyHSD(x=aov_now)
          
          
          
          
          if(exists("TukeyTableForEvaluation_byMitoOnly_Perimeter")){
            
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinPerimeter_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumPerimeter))), pad = "0"),
                                                                "_MaxPerimeter_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumPerimeter))), pad = "0"))
            # add the column
            TukeyTableForEvaluation_byMitoOnly_Perimeter <- cbind(TukeyTableForEvaluation_byMitoOnly_Perimeter, tukey_now$UniqueConditionTag[rownames(TukeyTableForEvaluation_byMitoOnly_Perimeter), 4,drop=FALSE])
            
          } else {
            #Change colname for adjusted pvalues
            colnames(tukey_now$UniqueConditionTag)[4] <- paste0(mitoCatNow, "_MinPerimeter_",
                                                                stringr::str_pad(minDNow, nchar(as.character(max(minimumPerimeter))), pad = "0"),
                                                                "_MaxPerimeter_",
                                                                stringr::str_pad(maxDNow, nchar(as.character(max(maximumPerimeter))), pad = "0"))
            # Just copy
            TukeyTableForEvaluation_byMitoOnly_Perimeter <- tukey_now$UniqueConditionTag[, 4,drop=FALSE]
          }
          
          
          
        }
        
      }
    }
  }
  
  
  # Save byCell
  
  #
  for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
    
    write.csv2(file=paste0(SelectedFolder, "\\",
                           "11-SelectedMitochondriaOnly",
                           "\\",
                           "ByPerimeter",
                           "\\", "ByCellOnly_",
                           mitoCatNow,".csv"),
               TukeyTableForEvaluation_byCellOnly_Perimeter[ , grepl(mitoCatNow, colnames(TukeyTableForEvaluation_byCellOnly_Perimeter)), drop=FALSE] )
    
    write.csv2(file=paste0(SelectedFolder, "\\",
                           "11-SelectedMitochondriaOnly",
                           "\\",
                           "ByPerimeter",
                           "\\", "ByMitoOnly_",
                           mitoCatNow,".csv"),
               TukeyTableForEvaluation_byMitoOnly_Perimeter[ , grepl(mitoCatNow, colnames(TukeyTableForEvaluation_byMitoOnly_Perimeter)), drop=FALSE] )
    
  }
  
  
  
  ####
  #### Perimeter test ONLY POSITIVE MITOCONDRIA
  #### END
  ####
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####
  #### Perimeter test ONLY POSITIVE MITOCONDRIA
  #### START
  ####
  
  
  
  ####
  #### Perimeter DEFINITION ONLY POSITIVE MITOCONDRIA
  #### END
  ####
  
  
  # PERIMETER DEFINITION
  minDNow <- 5.3
  maxDNow <- 8.7
  usefullColumns <- c("FileOrigin_base", "D", "Perim.", "PositiveMitochondria", "MitoCat", "FractionMitoPositive")
  
  
  
  # Ordena a tabela
  FinalResults <- FinalResults[order(FinalResults$FileOrigin_base), ]
  
  ####
  # Clean mitoCats
  # Add columns for categories
  FinalResults$MitoCat <- NA
  FinalResults$CellCat <- NA
  #FinalResults$ <- NA  
  
  
  # Add intermediate
  FinalResults[which(
    !is.na(FinalResults$Perim.)
  ),"MitoCat"] <- "2-Intermediate"
  
  
  
  # Add fragmented
  FinalResults[which(
    !is.na(FinalResults$Perim.) &
      FinalResults$Perim. < minDNow
  ),"MitoCat"] <- "1-Fragmented"
  
  # Add tubular
  FinalResults[which(
    !is.na(FinalResults$Perim.) &
      FinalResults$Perim. > maxDNow
  ),"MitoCat"] <- "3-Tubular"
  
  
  
  
  #####
  
  
  ##### Add CellCat for:
  # Cells labelled with "01-Mock" will be FULLY SELECTED
  # Cells labelled otherwise:
  # Only FinalResults$PositiveMitochondria TRUE
  # 
  TempTable_plyr_toLabel <- FinalResults[which(
    FinalResults$Info01 == "01-Mock" | FinalResults$PositiveMitochondria & FinalResults$PositiveCell
  ), ] %>%
    group_by(FileOrigin_base) %>%
    count(MitoCat)
  TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)


  ###
  FinalResultsSplitted <- split(FinalResults[,usefullColumns],
                                f = FinalResults[,"FileOrigin_base"])
  
  
  # myTableNow <- FinalResultsSplitted[20000]
  
  #FinalResultsSplitted[6102] -> myTableNow
  # define function
  CellularClassificationFunction <- function(myTableNow){
    #print(myTableNow$FileOrigin_base[1])
    if(length(which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]))>0){
      myTableNow$CellCat <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), "MitoCat"][which.max(TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), "n"])]
    } else {
    }
    return(myTableNow)
  }
  
  # Create the dataframe list
  system.time({
    myFinalResultsSplitted <- lapply(X=FinalResultsSplitted,
                                     FUN=CellularClassificationFunction)
  })
  
  #
  myFinalResultsSplittedjoined <- bind_rows(myFinalResultsSplitted) #, .id = "column_label"
  
  #
  # Ordena a tabela
  myFinalResultsSplittedjoined <- myFinalResultsSplittedjoined[order(myFinalResultsSplittedjoined$FileOrigin_base), ]
  
  # Return data
  FinalResults[,c("MitoCat", "CellCat")] <- myFinalResultsSplittedjoined[,c("MitoCat", "CellCat")]
  

  #
  remove(FinalResultsSplitted)
  remove(myFinalResultsSplitted)
  remove(myFinalResultsSplittedjoined)
  remove(TempTable_plyr_toLabel)
  
  #
  # CreateTempTable
  # Only transfected cells for non-mock conditions
  #    #  & !is.na(FinalResults$CellCat
  TempTable <- FinalResults[which(
    FinalResults$Info01 == "01-Mock" | FinalResults$PositiveMitochondria & FinalResults$PositiveCell
  ), ]


  #  
  theme_set(theme_pubclean())
  
  # TempTable$MitoCat <- as.factor(TempTable$MitoCat)
  # TempTable$CellCat <- as.factor(TempTable$CellCat)
  # unique(TempTable[,c("UniqueTag", "UniqueConditionTag", "CellCat")])
  
  ############## BY CELL
  df.summary <- unique(TempTable[,c("UniqueTag", "UniqueConditionTag", "CellCat", "FileOrigin_base")]) %>%
    group_by(UniqueTag, UniqueConditionTag) %>%
    count(CellCat)
  df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
  
  
  df.summary$Percentage <- NA
  
  # Add zeroes
  for(condNow in unique(df.summary$UniqueTag)){
    
    if("1-Fragmented" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
      #
    } else {
      MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
      MyLine[,"n"] <- 0
      MyLine[,"CellCat"] <- "1-Fragmented"
      df.summary <- rbind(df.summary, MyLine[1,])
    }
    
    
    if("2-Intermediate" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
      #print("tem")
    } else {
      MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
      MyLine[,"n"] <- 0
      MyLine[,"CellCat"] <- "2-Intermediate"
      df.summary <- rbind(df.summary, MyLine[1,])
    }
    
    
    
    if("3-Tubular" %in% df.summary[df.summary$UniqueTag==condNow, "CellCat"]){
      #print("tem")
    } else {
      MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
      MyLine[,"n"] <- 0
      MyLine[,"CellCat"] <- "3-Tubular"
      df.summary <- rbind(df.summary, MyLine[1,])
    }
    
    
    
    
    #
    # Add percentage
    df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
      sum(df.summary[df.summary$UniqueTag==condNow, "n"])
    
  }
  
  
  df.summary2 <- df.summary %>%
    group_by(UniqueConditionTag, CellCat) %>%
    summarise(mean = mean(Percentage, na.rm=TRUE),
              sd = sd(Percentage, na.rm = TRUE))
  
  #df.summary2
  
  
  
  dir.create(path = paste0(SelectedFolder, "\\", "11-SelectedMitochondriaOnly", "\\"), showWarnings = FALSE)
  
  # (2) Bar plots of means + individual jitter points + errors
  png(file=paste0(SelectedFolder, "\\",
                  "11-SelectedMitochondriaOnly",
                  "\\", 
                  "ByCell_MinPerimeter_", minDNow, "_MaxPerimeter_", maxDNow, ".png"),
      width = 4600, height = 2200,
      res=300)
  
  print(
    ggplot(data = df.summary2,
           aes(x=UniqueConditionTag,
               y=mean,
               colour = CellCat,
               fill = CellCat)) +
      
      geom_bar(stat = "identity",
               position=position_dodge()) +
      geom_errorbar(
        aes(ymin = mean-sd,
            ymax = mean+sd),
        colour="black",
        width = 0.5,
        position = position_dodge(.9)) +
      geom_jitter(data=df.summary,
                  aes(x=UniqueConditionTag,
                      y=Percentage,
                      fill = CellCat,
                      colour = CellCat),
                  colour="black",
                  alpha=0.3,
                  position = position_jitterdodge(0.1)) +
      scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
      scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
  )
  dev.off()
  
  
  
  #
  for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
    
    # Make anova
    aov_now <- aov(formula=Percentage~UniqueConditionTag,
                   data=df.summary[df.summary[,"CellCat"] == mitoCatNow, ])
    
    tukey_now <- TukeyHSD(x=aov_now)
    
    write.csv2(tukey_now$UniqueConditionTag,
               file=paste0(SelectedFolder, "\\",
                           "11-SelectedMitochondriaOnly",
                           "\\", mitoCatNow, "_",
                           "ByCell_MinPerimeter_", minDNow, "_MaxPerimeter_", maxDNow, ".csv")  )
  }
  
  
  
  
  ########### BY MITOCHONDRIA
  df.summary <- TempTable[!is.na(TempTable$MitoCat),] %>%
    group_by(UniqueTag, UniqueConditionTag) %>%
    count(MitoCat)
  df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
  
  
  df.summary$Percentage <- NA
  
  # Add zeroes
  for(condNow in unique(df.summary$UniqueTag)){
    
    
    if("1-Fragmented" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
      #
    } else {
      MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
      MyLine[,"n"] <- 0
      MyLine[,"MitoCat"] <- "1-Fragmented"
      df.summary <- rbind(df.summary, MyLine[1,])
    }
    
    
    if("2-Intermediate" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
      #print("tem")
    } else {
      MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
      MyLine[,"n"] <- 0
      MyLine[,"MitoCat"] <- "2-Intermediate"
      df.summary <- rbind(df.summary, MyLine[1,])
    }
    
    
    
    if("3-Tubular" %in% df.summary[df.summary$UniqueTag==condNow, "MitoCat"]){
      #print("tem")
    } else {
      MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
      MyLine[,"n"] <- 0
      MyLine[,"MitoCat"] <- "3-Tubular"
      df.summary <- rbind(df.summary, MyLine[1,])
    }
    
    #
    # Add percentage
    df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
      sum(df.summary[df.summary$UniqueTag==condNow, "n"])
    
  }
  
  
  df.summary2 <- df.summary %>%
    group_by(UniqueConditionTag, MitoCat) %>%
    summarise(mean = mean(Percentage, na.rm=TRUE),
              sd = sd(Percentage, na.rm = TRUE))
  
  #df.summary2
  
  #####
  
  # (2) Bar plots of means + individual jitter points + errors
  png(file=paste0(SelectedFolder, "\\",
                  "11-SelectedMitochondriaOnly",
                  "\\", 
                  "ByMito_MinPerimeter_", minDNow, "_MaxPerimeter_", maxDNow, ".png"),
      width = 4600, height = 2200,
      res=300)
  
  print(
    ggplot(data = df.summary2,
           aes(x=UniqueConditionTag,
               y=mean,
               colour = MitoCat,
               fill = MitoCat)) +
      
      geom_bar(stat = "identity",
               position=position_dodge()) +
      geom_errorbar(
        aes(ymin = mean-sd,
            ymax = mean+sd),
        colour="black",
        width = 0.5,
        position = position_dodge(.9)) +
      geom_jitter(data=df.summary,
                  aes(x=UniqueConditionTag,
                      y=Percentage,
                      fill = MitoCat,
                      colour = MitoCat),
                  colour="black",
                  alpha=0.3,
                  position = position_jitterdodge(0.1)) +
      scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
      scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
  )
  dev.off()
  
  
  ### STAT TEST
  
  
  
  #
  for(mitoCatNow in c("1-Fragmented", "2-Intermediate", "3-Tubular")){
    
    # Make anova
    aov_now <- aov(formula=Percentage~UniqueConditionTag,
                   data=df.summary[df.summary[,"MitoCat"] == mitoCatNow, ])
    
    tukey_now <- TukeyHSD(x=aov_now)
    
    write.csv2(tukey_now$UniqueConditionTag,
               file=paste0(SelectedFolder, "\\",
                           "11-SelectedMitochondriaOnly",
                           "\\", mitoCatNow, "_",
                           "ByMito_MinPerimeter_", minDNow, "_MaxPerimeter_", maxDNow, ".csv")  )
  }
  
  
  ####
  ####
  #### FINAL MITOCHONDRIAL CLASSES DEFINITION 
  #### END
  ####
  
  
  
  

  
  
  
  
  
  
  ####### SAVE
  write.csv2(FinalResults,
             file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "FinalResults_complete_mitoClass.csv"))
  
  
  
  ####
  ####
  #### MITOCHONDRIA 
  #### END
  ####
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####
  ####
  #### MITOPHAGY 
  #### START
  ####
  
  
  
  ############### START HERE
  
  #####
  FinalResults <- fread(file=paste0(SelectedFolder, "\\", "10-DescriptiveResults", "\\", "FinalResults_complete_mitoClass.csv"),
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        header=TRUE,
                        sep=";",
                        dec=",")
  colnames(FinalResults)[1] <- "ID"
  
  
  
  FinalResults <- data.frame(FinalResults, check.names = FALSE, stringsAsFactors = FALSE)
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  FinalResults$ID <- NULL
  
  
  
  ##### DEFINE POSITIVE MITOPHAGY
  FinalResults$PositiveMitophagySensor <- FALSE
  FinalResults$RatioMitophagy_Cellwise <- (FinalResults$Mitophagy_Mean / FinalResults$CellMeanMitophagy)
  # Calcula as razões na célula entre média 

  
  
  # Descriptive plots
  # CellMeanMitophagy
  
  dir.create(path = paste0(SelectedFolder, "\\", "12-Mitophagy", "\\"), showWarnings = FALSE)
  
  
  png(paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "20_MitophagySensor_levels.png"),
      width=1500, height=1200, res=300)
  hist(FinalResults[, "CellMeanMitophagy"],
       breaks=200,
       main="Mitophagy Sensor Label",
       xlab="Mitophagy mean cell level",
       ylab="Amount",
       border="white", 
       col="#e31a1c")
  
  dev.off()
  
  
  
  
  
  ########### MITOPHAGY LOOPS
  
  

  ###########
  ###########
  ########### MITOPHAGY LOOP: CATEGORIES 1/2/3+
  ########### START
  ###########
  
  minMitophagyPossibilities <- c(500, 1000, 1500, 3000)
  RatioMitophagyPossibilities <- seq(from=1.1, to=2.1, by=0.1)

  for(minMitophagyLabel in minMitophagyPossibilities){
    
    print(paste0("Estou no corte de fluorescência de mitofagia ", minMitophagyLabel))
    
    dir.create(path = paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\"), showWarnings = FALSE)
    
    
    
    ##### DEFINE POSITIVE MITOPHAGY
    FinalResults$PositiveMitophagySensor <- FALSE
    
    # Marca as células positivas
    FinalResults$PositiveMitophagySensor[FinalResults$CellMeanMitophagy > minMitophagyLabel] <- TRUE

    
    
    temp <- FinalResults[ FinalResults$Info03 == "01-Ctl" ,  ]
    png(paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\", "01-MitophagySensor_levels_hist_control.png"),
        width=1500, height=700, res=300)
    hist(temp$RatioMitophagy_Cellwise, breaks=200,
         main=mean(temp$RatioMitophagy_Cellwise)
         #,
         #xlim=c(5,20)
    )
    dev.off()
    temp <- FinalResults[ FinalResults$Info03 != "01-Ctl" ,  ]
    png(paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\", "01-MitophagySensor_levels_hist_treatment.png"),
        width=1500, height=700, res=300)
    hist(temp$RatioMitophagy_Cellwise, breaks=200,
         main=mean(temp$RatioMitophagy_Cellwise)
         #,
         #xlim=c(5,20)
    )
    dev.off()
    
    
    
    
    
    for(RatioMitophagyCutoff in RatioMitophagyPossibilities){

      print(paste0("Estou no corte de razão de mitofagia ", RatioMitophagyCutoff))
      
      FinalResults$PositiveMitophagy <- FALSE
      FinalResults$PositiveMitophagyDots_posMito <- FALSE
      FinalResults$PositiveMitophagyDots_negMito <- FALSE
      
      FinalResults$PosMito_number <- FALSE
      FinalResults$NegMito_number <- FALSE
      
      FinalResults$PositiveMitophagyCell <- FALSE
      
      FinalResults$PositiveMitophagyDots_TotalMito <- FALSE
      FinalResults$TotalMito_number <- FALSE
      
    
      # Define positive mitochondrias using Cutoff
      FinalResults$PositiveMitophagy[(
        FinalResults$PositiveMitophagySensor &
          FinalResults$RatioMitophagy_Cellwise > RatioMitophagyCutoff
      )] <- TRUE
      #
      
      usefullColumns <- c("FileOrigin_base",
                          "D", "PositiveMitochondria", "MitoCat", "CellCat",
                          "PositiveMitophagy", "PositiveMitophagySensor",
                          "PositiveMitophagyDots_posMito",
                          "PositiveMitophagyDots_negMito",
                          "PosMito_number",
                          "NegMito_number",
                          "PositiveMitophagyCell")
      
      
      #####
      ##### Add Number of Mitophagic dots for POSITIVE MITOCHONDRIA
      # 
      TempTable_plyr_toLabel <- FinalResults[which(FinalResults$PositiveMitochondria & FinalResults$PositiveCell
      ), ] %>%
        group_by(FileOrigin_base) %>%
        count(PositiveMitophagy)
      TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
      ###
      FinalResultsSplitted <- split(FinalResults[,usefullColumns],
                                    f = FinalResults[,"FileOrigin_base"])
      
      
      
      

      # define function
      CellularClassificationFunction <- function(myTableNow){
        #print(myTableNow$FileOrigin_base[1])
        if(length(which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]))>0){
          
          InternalTemp  <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), ]
          
          # If there is positive mitochondria, else...
          if(TRUE %in% InternalTemp[,"PositiveMitophagy"]){
            myTableNow[,"PositiveMitophagyDots_posMito"] <- InternalTemp[InternalTemp[,"PositiveMitophagy"], "n"]
          } else {
            myTableNow[,"PositiveMitophagyDots_posMito"] <- 0
          }
          myTableNow[,"PosMito_number"] <- sum(InternalTemp[,"n"])
         } else {
           myTableNow[,"PositiveMitophagyDots_posMito"] <- 0
           myTableNow[,"PosMito_number"] <- 0
        }
        return(myTableNow)
      }
      
      # Create the dataframe list
      system.time({
        myFinalResultsSplitted <- lapply(X=FinalResultsSplitted,
                                         FUN=CellularClassificationFunction)
      })
      
      #
      myFinalResultsSplittedjoined <- bind_rows(myFinalResultsSplitted) #, .id = "column_label"

      # Ordena a tabela
      myFinalResultsSplittedjoined <- myFinalResultsSplittedjoined[order(myFinalResultsSplittedjoined$FileOrigin_base), ]
      
      # Return data
      FinalResults[,c("PositiveMitophagyDots_posMito", "PosMito_number")] <- myFinalResultsSplittedjoined[,c("PositiveMitophagyDots_posMito", "PosMito_number")]
      
      
      
      #length(myFinalResultsSplitted)
      #sum(myFinalResultsSplitted[[49869]]$PositiveMitochondria)
      #length(myFinalResultsSplitted[[49869]]$PositiveMitochondria)
      #length(myFinalResultsSplitted[[49869]]$PositiveMitochondria) - sum(myFinalResultsSplitted[[49869]]$PositiveMitochondria)
      
      for(i in c(      "FinalResultsSplitted",
                       "TempTable_plyr_toLabel",
                       "InternalTemp",
                       "myTableNow",
                       "myFinalResultsSplitted",
                       "myFinalResultsSplittedjoined")){
       
        if(exists(i)){
          remove(list=i)
        } 
      }
      
      #####
      ##### Add Number of Mitophagic dots for NEGATIVE MITOCHONDRIA
      # 
      TempTable_plyr_toLabel <- FinalResults[which(!FinalResults$PositiveMitochondria
      ), ] %>%
        group_by(FileOrigin_base) %>%
        count(PositiveMitophagy)
      TempTable_plyr_toLabel <- data.frame(TempTable_plyr_toLabel)
      ###
      FinalResultsSplitted <- split(FinalResults[,usefullColumns],
                                    f = FinalResults[,"FileOrigin_base"])
      
      
      
      
      
      # define function
      CellularClassificationFunction <- function(myTableNow){
        #print(myTableNow$FileOrigin_base[1])
        if(length(which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]))>0){
          
          InternalTemp  <- TempTable_plyr_toLabel[which(TempTable_plyr_toLabel$FileOrigin_base == myTableNow$FileOrigin_base[1]), ]
          
          # If there is positive mitochondria, else...
          if(TRUE %in% InternalTemp[,"PositiveMitophagy"]){
            myTableNow[,"PositiveMitophagyDots_negMito"] <- InternalTemp[InternalTemp[,"PositiveMitophagy"], "n"]
          } else {
            myTableNow[,"PositiveMitophagyDots_negMito"] <- 0
          }
          myTableNow[,"NegMito_number"] <- sum(InternalTemp[,"n"])
        } else {
          myTableNow[,"PositiveMitophagyDots_negMito"] <- 0
          myTableNow[,"NegMito_number"] <- 0
        }
        return(myTableNow)
      }
      
      # Create the dataframe list
      system.time({
        myFinalResultsSplitted <- lapply(X=FinalResultsSplitted,
                                         FUN=CellularClassificationFunction)
      })
      
      #
      myFinalResultsSplittedjoined <- bind_rows(myFinalResultsSplitted) #, .id = "column_label"
      
      # Ordena a tabela
      myFinalResultsSplittedjoined <- myFinalResultsSplittedjoined[order(myFinalResultsSplittedjoined$FileOrigin_base), ]
      
      # Return data
      FinalResults[,c("PositiveMitophagyDots_negMito", "NegMito_number")] <- myFinalResultsSplittedjoined[,c("PositiveMitophagyDots_negMito", "NegMito_number")]
      
      
      #myFinalResultsSplitted[[49869]]
      #21/(21+33)
      #length(myFinalResultsSplitted)
      #sum(myFinalResultsSplitted[[49869]]$PositiveMitochondria)
      #length(myFinalResultsSplitted[[49869]]$PositiveMitochondria)
      #length(myFinalResultsSplitted[[49869]]$PositiveMitochondria) - sum(myFinalResultsSplitted[[49869]]$PositiveMitochondria)
      
      
      for(i in c(      "FinalResultsSplitted",
                       "TempTable_plyr_toLabel",
                       "InternalTemp",
                       "myTableNow",
                       "myFinalResultsSplitted",
                       "myFinalResultsSplittedjoined")){
        
        if(exists(i)){
          remove(list=i)
        } 
      }
      
      
      
      
      #
######### EVALUATIONS
      # FinalResults$PositiveMitophagyDots_TotalMito
      # FinalResults$TotalMito_number <- NA
      # Total Mitochondria
      FinalResults[,c("PositiveMitophagyDots_TotalMito")] <- rowSums(FinalResults[,c("PositiveMitophagyDots_negMito", "PositiveMitophagyDots_posMito")])
      FinalResults[,c("TotalMito_number")] <- rowSums(FinalResults[,c("NegMito_number", "PosMito_number")])
      
      # Create per cell table
      FinalResults_PerCell <- unique(FinalResults[,c("UniqueTag",
                                  "UniqueConditionTag",
                                  "CellCat",
                                  "PositiveMitophagySensor",
                                  "FileOrigin_base",
                                  "PositiveMitophagyDots_negMito",
                                  "NegMito_number",
                                  "PositiveMitophagyDots_posMito",
                                  "PosMito_number",
                                  "PositiveMitophagyDots_TotalMito",
                                  "TotalMito_number")])
      
      
      ####
      
      
      
      ###
      ### ALL MITCOCHONDRIAS
      ###
      SELECTEDCOLUMN <- "PositiveMitophagyDots_TotalMito"
      
      for(SELECTEDCOLUMN in c("PositiveMitophagyDots_TotalMito",
                              "PositiveMitophagyDots_posMito",
                              "PositiveMitophagyDots_negMito")){
        
        
        #### Seleciona os dados
        TempTable <- FinalResults_PerCell[which(
          FinalResults_PerCell$PositiveMitophagySensor
        ),]
        
        # Renomeia as categorias
        TempTable[(TempTable[,SELECTEDCOLUMN] >= 3), SELECTEDCOLUMN] <- "3+"
        
        
        #
        df.summary <- TempTable[,] %>%
          group_by(UniqueTag, UniqueConditionTag) %>%
          #### AQUI ESCOLHE O TIPO
          count(get(SELECTEDCOLUMN))
        df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
        colnames(df.summary)[3] <- SELECTEDCOLUMN
        df.summary$Percentage <- NA
        
        # Add zeroes
        for(condNow in unique(df.summary$UniqueTag)){
          
          
          if("0" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "0"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          if("1" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "1"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          if("2" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "2"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          if("3+" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "3+"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          #
          # Add percentage
          df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
            sum(df.summary[df.summary$UniqueTag==condNow, "n"])
          
        }
        
        
        df.summary2 <- df.summary %>%
          group_by(UniqueConditionTag, get(SELECTEDCOLUMN)) %>%
          summarise(mean = mean(Percentage, na.rm=TRUE),
                    sd = sd(Percentage, na.rm = TRUE))
        df.summary2 <- data.frame(df.summary2, stringsAsFactors = FALSE)
        
        colnames(df.summary2)[2] <- SELECTEDCOLUMN
        #df.summary2
        
        
        # REMOVE ZEROES FROM PLOT
        df.summary <- df.summary[!(df.summary[,SELECTEDCOLUMN] == "0"),]
        df.summary2 <- df.summary2[!(df.summary2[,SELECTEDCOLUMN] == "0"),]
        
        #####
        
        
        
        # (2) Bar plots of means + individual jitter points + errors
        png(file=paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\", "02-",
                        SELECTEDCOLUMN,
                        "_CutOff_", RatioMitophagyCutoff, ".png"),
            width = 4600, height = 2200,
            res=300)
        
        print(
          ggplot(data = df.summary2,
                 aes_string(x="UniqueConditionTag",
                            y="mean",
                            colour = SELECTEDCOLUMN ,
                            fill = SELECTEDCOLUMN )) +
            
            geom_bar(stat = "identity",
                     position=position_dodge()) +
            geom_errorbar(
              aes(ymin = mean-sd,
                  ymax = mean+sd),
              colour="black",
              width = 0.5,
              position = position_dodge(.9)) +
            geom_jitter(data=df.summary,
                        aes_string(x="UniqueConditionTag",
                                   y="Percentage",
                                   fill = SELECTEDCOLUMN,
                                   colour = SELECTEDCOLUMN),
                        colour="black",
                        alpha=0.3,
                        position = position_jitterdodge(0.1)) +
            scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
        )
        dev.off()
        
        
        #
        #
        for(mitoCatNow in c("1", "2", "3+")){
          
          # Make anova
          aov_now <- aov(formula=Percentage~UniqueConditionTag,
                         data=df.summary[df.summary[,SELECTEDCOLUMN] == mitoCatNow, ])
          
          tukey_now <- TukeyHSD(x=aov_now)
          
          write.csv2(tukey_now$UniqueConditionTag,
                     file=paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\", "02-",
                                 SELECTEDCOLUMN,
                                 "Anova_", mitoCatNow, "_CutOff_", RatioMitophagyCutoff, ".csv")  )
        }
        
        
        
        
      }
      
      
      
      
      
      
      
      
      
      ###
      ### TUBULAR CELLS
      ###
      SELECTEDCOLUMN <- "PositiveMitophagyDots_TotalMito"
      
      for(SELECTEDCOLUMN in c("PositiveMitophagyDots_TotalMito",
                              "PositiveMitophagyDots_posMito",
                              "PositiveMitophagyDots_negMito")){
        
        
        #### Seleciona os dados
        TempTable <- FinalResults_PerCell[which(
          FinalResults_PerCell$PositiveMitophagySensor &
            FinalResults_PerCell$CellCat == "3-Tubular"
        ),]
        
        # Renomeia as categorias
        TempTable[(TempTable[,SELECTEDCOLUMN] >= 3), SELECTEDCOLUMN] <- "3+"
        
        
        #
        df.summary <- TempTable[,] %>%
          group_by(UniqueTag, UniqueConditionTag) %>%
          #### AQUI ESCOLHE O TIPO
          count(get(SELECTEDCOLUMN))
        df.summary <- data.frame(df.summary, stringsAsFactors = FALSE)
        colnames(df.summary)[3] <- SELECTEDCOLUMN
        df.summary$Percentage <- NA
        
        # Add zeroes
        for(condNow in unique(df.summary$UniqueTag)){
          
          
          if("0" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "0"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          if("1" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "1"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          if("2" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "2"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          if("3+" %in% df.summary[df.summary$UniqueTag==condNow, SELECTEDCOLUMN]){
            #
          } else {
            MyLine <- df.summary[df.summary$UniqueTag==condNow, ]
            MyLine[,"n"] <- 0
            MyLine[,SELECTEDCOLUMN] <- "3+"
            df.summary <- rbind(df.summary, MyLine[1,])
          }
          #
          # Add percentage
          df.summary[df.summary$UniqueTag==condNow, "Percentage"] <- 100*df.summary[df.summary$UniqueTag==condNow, "n"]/
            sum(df.summary[df.summary$UniqueTag==condNow, "n"])
          
        }
        
        
        df.summary2 <- df.summary %>%
          group_by(UniqueConditionTag, get(SELECTEDCOLUMN)) %>%
          summarise(mean = mean(Percentage, na.rm=TRUE),
                    sd = sd(Percentage, na.rm = TRUE))
        df.summary2 <- data.frame(df.summary2, stringsAsFactors = FALSE)
        
        colnames(df.summary2)[2] <- SELECTEDCOLUMN
        #df.summary2
        
        
        # REMOVE ZEROES FROM PLOT
        df.summary <- df.summary[!(df.summary[,SELECTEDCOLUMN] == "0"),]
        df.summary2 <- df.summary2[!(df.summary2[,SELECTEDCOLUMN] == "0"),]
        
        #####
        
        
        
        # (2) Bar plots of means + individual jitter points + errors
        png(file=paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\", "02-Tubular_",
                        SELECTEDCOLUMN,
                        "_CutOff_", RatioMitophagyCutoff, ".png"),
            width = 4600, height = 2200,
            res=300)
        
        print(
          ggplot(data = df.summary2,
                 aes_string(x="UniqueConditionTag",
                            y="mean",
                            colour = SELECTEDCOLUMN ,
                            fill = SELECTEDCOLUMN )) +
            
            geom_bar(stat = "identity",
                     position=position_dodge()) +
            geom_errorbar(
              aes(ymin = mean-sd,
                  ymax = mean+sd),
              colour="black",
              width = 0.5,
              position = position_dodge(.9)) +
            geom_jitter(data=df.summary,
                        aes_string(x="UniqueConditionTag",
                                   y="Percentage",
                                   fill = SELECTEDCOLUMN,
                                   colour = SELECTEDCOLUMN),
                        colour="black",
                        alpha=0.3,
                        position = position_jitterdodge(0.1)) +
            scale_colour_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            scale_fill_manual(values=c("#7fc97f", "#fdc086", "#beaed4")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
        )
        dev.off()
        
        
        #
        #
        for(mitoCatNow in c("1", "2", "3+")){
          
          # Make anova
          aov_now <- aov(formula=Percentage~UniqueConditionTag,
                         data=df.summary[df.summary[,SELECTEDCOLUMN] == mitoCatNow, ])
          
          tukey_now <- TukeyHSD(x=aov_now)
          
          write.csv2(tukey_now$UniqueConditionTag,
                     file=paste0(SelectedFolder, "\\", "12-Mitophagy", "\\", "Cut_", minMitophagyLabel, "\\", "02-Tubular_",
                                 SELECTEDCOLUMN,
                                 "Anova_", mitoCatNow, "_CutOff_", RatioMitophagyCutoff, ".csv")  )
        }
        
        
        
        
      }
      
      
      
      
      
      
      
      
      
      
      
    }
  }
  
  ###########
  ########### END
  ########### MITOPHAGY LOOP: CATEGORIES 1/2/3+
  ########### 
  ###########
  
  
  
  
  ####
  ####
  #### MITOPHAGY 
  #### END
  ####
  
  
  
  
  
  
  
  
  
  
  ####
  ####
  #### ANALYSIS 
  #### END
  ####
  
  
  
  
  
  
  
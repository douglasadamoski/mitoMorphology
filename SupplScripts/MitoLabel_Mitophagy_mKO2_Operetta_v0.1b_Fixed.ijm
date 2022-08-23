/*
Mitochondrial Label analysis using Columbus Export

Script:
MitoLabel_Mitophagy_mKO2_Operetta_v0.1_Fixed.ijm

Douglas Adamoski
douglas.adamoski@gmail.com

Input:

Files sorted by ColumbusFileSorter_v0.2 in Well-Stack-Field-Timepoint mode.

ALERT 1:
The experiment should be performed on Operetta equipement using 40X magnification objetive

ALERT 2:
NO other files/folders than images should be in same folder.

ALERT 3:
Remove []'s from path!

CHANGELOG
v0.1 - First version, feb 09th 2022
		First version
*/


// Create funcion ArrayUnique
// objective of this function is clean one vector in order to remove duplicates
// gently provided by Richard Wheeler:
// http://www.richardwheeler.net/contentpages/textgallery.php?gallery=ImageJ_Macros
function ArrayUnique(array) {
	array 	= Array.sort(array);
	array 	= Array.concat(array, 999999);
	uniqueA = newArray();
	i = 0;	
   	while (i<(array.length)-1) {
		if (array[i] == array[(i)+1]) {
			//print("found: "+array[i]);			
		} else {
			uniqueA = Array.concat(uniqueA, array[i]);
		}
   		i++;
   	}
	return uniqueA;
}


// run("Memory & Threads...", "maximum=30000 parallel=4 run");

// Asks for main dir
  dir = getDirectory("Choose a Directory ");
// Request variables from user
// changed from original example: https://imagej.nih.gov/ij/macros/DialogDemo.txt

// Maximum scratch find
// Valor1=30000;
// Valor2="002.tif";
Valor5=0.248;

// Create dialog window
  Dialog.create("Analysis parameters (Pre-entered suggestions)");
//  Dialog.addMessage("Minimum spheroid size");
//  Dialog.addNumber("Value:", Valor1);
//  Dialog.addMessage("Brightfield channel");
//  Dialog.addString("Value:", Valor2);
  Dialog.addMessage("Pixel size (um) - Operetta's 40X Objective:");
  Dialog.addNumber("Value:", Valor5);
  Dialog.addCheckbox("Check this box to run in batch mode", false);
  Dialog.show();
  
//  Valor1 = Dialog.getNumber();
//  Valor2 = Dialog.getString();
  Valor5 = Dialog.getNumber();
  BatchStatus = Dialog.getCheckbox();

// Prints ImageJ Version
print("Fiji/ImageJ Version:", getVersion());
print("MitoLabel_Mitophagy_mKO2_Operetta_v0.1_Fixed.ijm");

// Report the values  
  print("Main dir:", dir);
  print("Pixel size:", Valor5);
//  print("Brightfield channel:", Valor2);
    
print("Script started!");

// Prints the starttime
print("Script started at:");
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     print(TimeString+second);

// Creates an array with all folders in dir
WellsToEvaluate = getFileList(dir);





// Create results dir
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator);
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator + "00-CellROIs" + File.separator);
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator + "01-Images" + File.separator);
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator + "02-Mito_Results" + File.separator);
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator + "03-Fractal_Results" + File.separator);
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator + "04-Protein_Results" + File.separator);
File.makeDirectory(dir + File.separator + ".." + File.separator + "Results" + File.separator + "05-Protein_Fractal_Results" + File.separator);



// Enters BatchMode
if (BatchStatus==true) setBatchMode(true);
     // LOOP 01: there are no multiple collections for this experiment   
    
    // LOOP 02: WELLS
    // START
    for (Quiabo=0; Quiabo<WellsToEvaluate.length; Quiabo++){

	    // Cria o nome atual da pasta
	    WellNow = dir + WellsToEvaluate[Quiabo] + File.separator;


	    // Creates an array with all folders in dir
        StacksToEvaluate = getFileList(WellNow);




  // LOOP 03: STACKS
    // START
    for (Chuchu=0; Chuchu<StacksToEvaluate.length; Chuchu++){

	    // Cria o nome atual da pasta
	    StackNow = dir + WellsToEvaluate[Quiabo] + File.separator + StacksToEvaluate[Chuchu] + File.separator;


	    // Creates an array with all folders in dir
        FieldsToEvaluate = getFileList(StackNow);




    // LOOP 04: FIELDS
    // START
    for (Melao=0; Melao<FieldsToEvaluate.length; Melao++){

	    // Cria o nome atual da pasta
	    FieldNow = dir + WellsToEvaluate[Quiabo] + File.separator + StacksToEvaluate[Chuchu] + File.separator + FieldsToEvaluate[Melao] + File.separator;


	    // Creates an array with all folders in dir
        TimepointsToEvaluate = getFileList(FieldNow);



    // LOOP 05: TIMEPOINTS - IMAGE ANALYSIS
    // START
    for (Morango=0; Morango<TimepointsToEvaluate.length; Morango++){


	    // Cria o nome atual da pasta
	    TimepointNow = dir + WellsToEvaluate[Quiabo] + File.separator + StacksToEvaluate[Chuchu] + File.separator + FieldsToEvaluate[Melao] + File.separator + TimepointsToEvaluate[Morango] + File.separator;

	    // Creates an array with all Images in dir
        ImagesToEvaluate = getFileList(TimepointNow);

	   // Image opening loop
	   for (Mamao=0; Mamao<ImagesToEvaluate.length; Mamao++){
		ImageNow = dir + WellsToEvaluate[Quiabo] + File.separator + StacksToEvaluate[Chuchu] + File.separator + FieldsToEvaluate[Melao] + File.separator + TimepointsToEvaluate[Morango] + ImagesToEvaluate[Mamao];
	    // Open the selected channel
	   open(ImageNow);
		}


		// Correct image details
		// Valor5=0.248;
		setOption("ScaleConversions", true);
		// run("Properties...", "unit=um pixel_width="+ Valor5 + " pixel_height=" + Valor5);
   		 run("Set Scale...", "distance=1 known=" + Valor5 + " unit=um global");
		// 

		//
		// CHANNEL RENAMING START
		//
		// DAPI Channel
		selectWindow("003.tif");
		rename("Nuclei");

		// Mitochondria Channel
		selectWindow("002.tif");
		rename("Mitochondria");

		// Protein Channel
		selectWindow("004.tif");
		rename("Protein");

		// Mitophagy Channel
		selectWindow("001.tif");
		rename("Mitophagy");
		
		//
		// CHANNEL RENAMING END
		//

		//
		// ANALYSIS START
		//



		//
		// PROTEIN ROI DEFINITION
		// START
		//

		// Select adequated channel
		selectWindow("Protein");
		// Create a copy of the image
        run("Duplicate...", "title=Protein_adjust");

		// Select adequated channel
		selectWindow("Protein_adjust");
	    // Background adjustments
	    // Threshold image to remove ONLY BACKGROUND SIGNAL
	    // empty well area
	    setMinAndMax(350, 5000);
	    run("Apply LUT");
		
		// block size
        // the size of the local region around a pixel for which the histogram is equalized. This size should be larger than the size of features to be preserved.

        // histogram bins
        // the number of histogram bins used for histogram equalization. The implementation internally works with byte resolution, so values larger than 256 are not meaningful. This value also limits the quantification of the output when processing 8bit gray or 24bit RGB images. The number of histogram bins should be smaller than the number of pixels in a block.

        // max slope
        // limits the contrast stretch in the intensity transfer function. Very large values will let the histogram equalization do whatever it wants to do, that is result in maximal local contrast. The value 1 will result in the original image.

		// Rodar 4 vezes
		// run("Subtract Background...", "rolling=5 sliding");
		run("Enhance Local Contrast (CLAHE)", "blocksize=200 histogram=80 maximum=3 mask=*None* fast_(less_accurate)"); // 

        //
		run("Despeckle", "stack");
		
        // eles usaram 1 um
		// run("Subtract Background...", "rolling=3 stack");
		// run("Subtract Background...", "rolling=0.2 sliding");
		run("Subtract Background...", "rolling=10");

        //
		run("Despeckle", "stack");

		// run("Sigma Filter Plus", "radius=1 use=3 minimum=0.9 outlier stack");

		//run("Maximum...", "radius=1");
		run("Gamma...", "value=0.95 stack");

		// Create a copy of the image
        run("Duplicate...", "title=Protein_mask");

		//Select the original image
		selectWindow("Protein_mask");

		// Change to 8-bit
		run("8-bit");
		// Perform adaptative threshold to make the actual mask
		run("adaptiveThr ", "using=[Weighted mean] from=60 then=-15");

		// Run multiple Despeckle to remove isolated points
		for (i = 0; i < 10; i++) {
			run("Despeckle", "stack");
		}

		// Remove outliers
		run("Remove Outliers...", "radius=1.5 threshold=500 which=Bright stack");

		// unlink mitochondria
        run("Watershed Irregular Features", "erosion=3 convexity_threshold=0 separator_size=0-Infinity");
        // run("Watershed Irregular Features", "erosion=5 convexity_threshold=0 separator_size=0-Infinity");

		// setOption("BlackBackground", false);
		run("Colors...", "foreground=black background=white selection=yellow");
        run("Convert to Mask");

		//
		// PROTEIN ROI DEFINITION
		// END
		//



		//
		// MITOCHONDRIAL ROI DEFINITION
		// START
		// 

		// Select adequated channel
		selectWindow("Mitochondria");
		// Create a copy of the image
        run("Duplicate...", "title=Mitochondria_adjust");

		// Select adequated channel
		selectWindow("Mitochondria_adjust");
	    // Background adjustments
	    // Threshold image to remove ONLY BACKGROUND SIGNAL
	    // empty well area
	    setMinAndMax(100, 4000);
	    run("Apply LUT");

		// block size
        // the size of the local region around a pixel for which the histogram is equalized. This size should be larger than the size of features to be preserved.

        // histogram bins
        // the number of histogram bins used for histogram equalization. The implementation internally works with byte resolution, so values larger than 256 are not meaningful. This value also limits the quantification of the output when processing 8bit gray or 24bit RGB images. The number of histogram bins should be smaller than the number of pixels in a block.

        // max slope
        // limits the contrast stretch in the intensity transfer function. Very large values will let the histogram equalization do whatever it wants to do, that is result in maximal local contrast. The value 1 will result in the original image.

		// Rodar 4 vezes
		// run("Subtract Background...", "rolling=5 sliding");
		run("Enhance Local Contrast (CLAHE)", "blocksize=200 histogram=80 maximum=3 mask=*None* fast_(less_accurate)"); // 

        //
		run("Despeckle", "stack");
		
        // eles usaram 1 um
		// run("Subtract Background...", "rolling=3 stack");
		// run("Subtract Background...", "rolling=0.2 sliding");
		run("Subtract Background...", "rolling=3");

        //
		run("Despeckle", "stack");

		// run("Sigma Filter Plus", "radius=1 use=3 minimum=0.9 outlier stack");

		//run("Maximum...", "radius=1");
		run("Gamma...", "value=0.95 stack");

		// Create a copy of the image
        run("Duplicate...", "title=Mitochondria_mask");

		//Select the original image
		selectWindow("Mitochondria_mask");

		// Change to 8-bit
		run("8-bit");
		// Perform adaptative threshold to make the actual mask
		run("adaptiveThr ", "using=[Weighted mean] from=60 then=-15");

		// Run multiple Despeckle to remove isolated points
		for (i = 0; i < 10; i++) {
			run("Despeckle", "stack");
		}

		// Remove outliers
		run("Remove Outliers...", "radius=1.5 threshold=500 which=Bright stack");

		// unlink mitochondria
        run("Watershed Irregular Features", "erosion=3 convexity_threshold=0 separator_size=0-Infinity");
        // run("Watershed Irregular Features", "erosion=5 convexity_threshold=0 separator_size=0-Infinity");

		// setOption("BlackBackground", false);
		run("Colors...", "foreground=black background=white selection=yellow");
        run("Convert to Mask");


	    // FIND CELL START
		// Select nuclei image
		selectWindow("Nuclei");
	    // Background adjustments
	    // Threshold image to remove ONLY BACKGROUND SIGNAL
	    // empty well area
	    setMinAndMax(500, 4500);
	    run("Apply LUT");

		// Select mitochondria image
		selectWindow("Mitochondria");
		// Create a copy of the image
        run("Duplicate...", "title=Mitochondria_wholecell");
		//Select the original image
		selectWindow("Mitochondria_wholecell");
	    // Background adjustments
	    // Threshold image to make the best blurred cell area
	    // empty well area
	    setMinAndMax(50, 6000);
	    run("Apply LUT");
	    
	    // Apply a sum between nuclei and mitochondria channel
		imageCalculator("Max create", "Nuclei", "Mitochondria_wholecell");
		// Close mitochondrial cell mask
		selectWindow("Mitochondria_wholecell");
		close();
		// Select result image
		selectWindow("Result of Nuclei");
		// Enhance contrast
        run("8-bit");	
	    run("Enhance Local Contrast (CLAHE)", "blocksize=200 histogram=80 maximum=2 mask=*None* fast_(less_accurate)"); // 
		
		// Blur to remove mitochondrial aspect
		run("Gaussian Blur...", "sigma=3 scaled");
		// Find maxima
		run("Find Maxima...", "prominence=25 exclude output=[Segmented Particles]");
		// Select result window
		selectWindow("Result of Nuclei Segmented");
		
		// Analyze cell masks and add as ROIs
		run("Analyze Particles...", "add");
		
		// Close masks
		selectWindow("Result of Nuclei");
		run("Close");
		selectWindow("Result of Nuclei Segmented");
		run("Close");





		//
		// MITOPHAGY DEFINITION
		// START
		// 

		// Select adequated channel
		selectWindow("Mitophagy");
		// Create a copy of the image
        run("Duplicate...", "title=Mitophagy_adjust");

		// Select adequated channel
		selectWindow("Mitophagy_adjust");
	    // Background adjustments
	    // Threshold image to remove ONLY BACKGROUND SIGNAL
	    // empty well area
	    setMinAndMax(0, 16000);
	    run("Apply LUT");
	    
		//
		// MITOPHAGY DEFINITION
		// END
		// 





	    // Sanity check to evaluate if there are any ROI
		nROIs_Cell = roiManager("count");
		print(nROIs_Cell);
		if (nROIs_Cell > 0){

			// Salva as ROIs das células, caso existam
			roiManager("Save", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "00-CellROIs" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_" + "CellROI" + ".zip");

			// Limpa as ROIs para processar o restante
			// A ideia é sempre reabrir o arquivo com todas as ROIs no início do loop, por isso limpar aqui
			roiManager("Delete");

			// Reporta quantas ROIS a celula tem
			// nROIs_Cell = 10;
			print("Este campo tem:" + nROIs_Cell + " células.");
			// Iterate within each Cell ROI
			for (ROINow = 0; ROINow < nROIs_Cell; ROINow++) {
				// Reporta a roi atual
				print("Estou na célula " + ROINow + " de um total de " + nROIs_Cell + " células.");

 				//
				// MITOCHONDRIAL MASK MEASUREMENTS 
				// START
				// 
				
				// Seleciona a máscara mitocondrial
      		  	selectWindow("Mitochondria_mask");
				// duplica com um titulo generico descartavel para agora
        		run("Duplicate...", "title=masks");
       			 // Collect ROIs from file
				roiManager("Open", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "00-CellROIs" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_" + "CellROI" + ".zip");
				// Desselect any roi
       			roiManager("deselect");
       			// Seleciona a cópia genérica
    			selectWindow("masks");
				// Pega uma única roi
    			roiManager("Select", ROINow);
    			// roiManager("Select", 17);
    			run("Colors...", "foreground=white background=black selection=yellow");
    			// Deleta toda a informação externa a célula
    			run("Clear Outside");
    			// Clear ROI manager
				roiManager("deselect")
				roiManager("Delete");
				run("Select None");

				// Define as características a serem analisadas
				run("Set Measurements...", "add redirect=None decimal=1");
				// encontra as mitocondrias
				run("Analyze Particles...", "  show=[Overlay Masks] exclude add");

				// Sanity check to evaluate if there are any ROI
				nROIs = roiManager("count");
				if (nROIs>1){
					// Reporta a roi atual
					print("A célula " + ROINow + " de um total de " + nROIs + " mitocôndrias.");

					//
					// MITOCHONDRIAL MEASUREMENTS 
					// START
				    // 
					run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack add redirect=None decimal=6");
					
					// MEASURE MITO IMAGE
					// Seleciona as mitocondrias ajustadas, não a máscara binária
					selectWindow("Mitochondria");
					roiManager("deselect");
					roiManager("Measure");

					// Se os resultados existirem, salva eles
					if (isOpen("Results")){
					selectWindow("Results"); 
					saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "02-Mito_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Mito_Results" + ".csv");
					run("Close");
 					}

					// MEASURE PROTEIN IMAGE
					// Seleciona as proteinas ajustadas, não a máscara binária
					selectWindow("Protein");
					roiManager("deselect");
					roiManager("Measure");

					// Se os resultados existirem, salva eles
					if (isOpen("Results")){
					selectWindow("Results"); 
					saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "02-Mito_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Mito_Protein_Results" + ".csv");
					run("Close");
 					}


					// MEASURE MITOPHAGY IMAGE
					// Seleciona as proteinas ajustadas, não a máscara binária
					selectWindow("Mitophagy");
					roiManager("deselect");
					roiManager("Measure");

					// Se os resultados existirem, salva eles
					if (isOpen("Results")){
					selectWindow("Results"); 
					saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "02-Mito_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Mitophagy_Results" + ".csv");
					run("Close");
 					}

 					//
					// MITOCHONDRIAL MEASUREMENTS 
					// END
				    // 


 					
					//
					// FRACTAL MEASUREMENTS 
					// START
				    // 
    				selectWindow("masks");
 				  	roiManager("Show All with labels");
   					setOption("BlackBackground", true);
    				// Loop para cortar cada mitocondria
					for (u = 0; u < nROIs; u++) {
    					// Desseleciona todas as ROIs
    					roiManager("deselect");
    					// Tira a marcação também
    					run("Select None");
   						// Pega a imagem da máscara binária
   						selectWindow("Mitochondria_mask");
   						// Duplica ela e chama de crop
   						run("Duplicate...", "title=crop");
   						// Seleciona essa imagem
    					selectWindow("crop");
						// Seleciona a mitocondria atual
						roiManager("Select", u);
						//roiManager("Select", 22);
						//
   						run("Colors...", "foreground=white background=black selection=yellow");
   					 	// Corta o resto da imagem
   					 	run("Clear Outside");
   						// Se cortar a figura toda altera o resultado do fractal, não fazer
   						//  run("Crop");
    					//  run("Canvas Size...", "width=" + getWidth*5 + " height=" + getHeight*5 + " position=Center zero");
    					run("Fractal Box Count...", "box=2,3,4,6,8,12,16,32,64 black");
    					close();
    					close("crop");
					    // Next round!
    					// Fecha o loop do fractal
						}
    
					// Save Fractal image results caso existam
					if (isOpen("Results")){
						selectWindow("Results"); 
						saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "03-Fractal_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "")  + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Fractal_Results" + ".csv");
						// print(dir  + File.separator + ".." + File.separator + "Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + "YFP_results" + ".csv");
 						run("Close");
 					}

 					//
					// FRACTAL MEASUREMENTS 
					// END
					//
 					
					// fecha a imagem com as máscaras mitocondriais dessa célula
					selectWindow("masks");
    				close();


					// Clear ROI manager
					// Jogando fora todas as mitocondrias daquela célula
					roiManager("deselect");
					roiManager("Delete");
					
					// Fecha o if check para ver se existem mitocondrias
					}

					// Fecha a janela "masks" caso ela exista
					if (isOpen("masks")){
					selectWindow("masks"); 
					run("Close");
 					}

 				//
				// MITOCHONDRIAL MASK MEASUREMENTS 
				// END
				// 


				
 				//
				// PROTEIN MASK MEASUREMENTS 
				// START
				// 

				// Seleciona a máscara mitocondrial
      		  	selectWindow("Protein_mask");
				// duplica com um titulo generico descartavel para agora
        		run("Duplicate...", "title=masks");
       			// Collect ROIs from file
				roiManager("Open", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "00-CellROIs" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_" + "CellROI" + ".zip");
				// Desselect any roi
       			roiManager("deselect");
       			// Seleciona a cópia genérica
    			selectWindow("masks");
				// Pega uma única roi
    			roiManager("Select", ROINow);
    			// roiManager("Select", 17);
    			run("Colors...", "foreground=white background=black selection=yellow");
    			// Deleta toda a informação externa a célula
    			run("Clear Outside");
    			// Clear ROI manager
				roiManager("deselect")
				roiManager("Delete");
				run("Select None");

				// Define as características a serem analisadas
				run("Set Measurements...", "add redirect=None decimal=1");
				// encontra as mitocondrias
				run("Analyze Particles...", "  show=[Overlay Masks] exclude add");

				// Sanity check to evaluate if there are any ROI
				nROIs = roiManager("count");
				if (nROIs>1){
					// Reporta a roi atual
					print("A célula " + ROINow + " de um total de " + nROIs + " proteínas.");

					//
					// PROTEIN MEASUREMENTS 
					// START
				    // 
					run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack add redirect=None decimal=6");
					
					// MEASURE PROTEIN IMAGE
					// Seleciona as mitocondrias ajustadas, não a máscara binária
					selectWindow("Protein");
					roiManager("deselect");
					roiManager("Measure");

					// Se os resultados existirem, salva eles
					if (isOpen("Results")){
					selectWindow("Results"); 
					saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "04-Protein_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Protein_Results" + ".csv");
					run("Close");
 					}


					// MEASURE MITOPHAGY IMAGE
					// Seleciona as mitocondrias ajustadas, não a máscara binária
					selectWindow("Mitophagy");
					roiManager("deselect");
					roiManager("Measure");

					// Se os resultados existirem, salva eles
					if (isOpen("Results")){
					selectWindow("Results"); 
					saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "04-Protein_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Mitophagy_Results" + ".csv");
					run("Close");
 					}

 					//
					// PROTEIN MEASUREMENTS 
					// END
				    // 


 					
					//
					// FRACTAL MEASUREMENTS 
					// START
				    // 
    				selectWindow("masks");
 				  	roiManager("Show All with labels");
   					setOption("BlackBackground", true);
    				// Loop para cortar cada proteína
					for (u = 0; u < nROIs; ++u) {
    					// Desseleciona todas as ROIs
    					roiManager("deselect");
    					// Tira a marcação também
    					run("Select None");
   						// Pega a imagem da máscara binária
   						selectWindow("Protein_mask");
   						// Duplica ela e chama de crop
   						run("Duplicate...", "title=crop");
   						// Seleciona essa imagem
    					selectWindow("crop");
						// Seleciona a proteína atual
						roiManager("Select", u);
						//roiManager("Select", 22);
						//
   						run("Colors...", "foreground=white background=black selection=yellow");
   					 	// Corta o resto da imagem
   					 	run("Clear Outside");
   						// Se cortar a figura toda altera o resultado do fractal, não fazer
   						//  run("Crop");
    					//  run("Canvas Size...", "width=" + getWidth*5 + " height=" + getHeight*5 + " position=Center zero");
    					run("Fractal Box Count...", "box=2,3,4,6,8,12,16,32,64 black");
    					close();
    					close("crop");
					    // Next round!
    					// Fecha o loop do fractal
						}
    
					// Save Fractal image results caso existam
					if (isOpen("Results")){
						selectWindow("Results"); 
						saveAs("Results", dir  + File.separator + ".." + File.separator + "Results" + File.separator + "05-Protein_Fractal_Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "")  + "_" + replace(TimepointsToEvaluate[Morango], "/*$", "") + "_Cell_" + ROINow + "_Fractal_Results" + ".csv");
						// print(dir  + File.separator + ".." + File.separator + "Results" + File.separator + replace(WellsToEvaluate[Quiabo], "/*$", "") + "_" + replace(StacksToEvaluate[Chuchu], "/*$", "") + "_" + replace(FieldsToEvaluate[Melao], "/*$", "") + "_" + "YFP_results" + ".csv");
 						run("Close");
 					}

 					//
					// FRACTAL MEASUREMENTS 
					// END
					//
 					
					// fecha a imagem com as máscaras proteina dessa célula
					selectWindow("masks");
    				close();


					// Clear ROI manager
					// Jogando fora todas as proteina daquela célula
					roiManager("deselect");
					roiManager("Delete");
					
					// Fecha o if check para ver se existem proteina
					}

				// Fecha a janela "masks" caso ela exista
				if (isOpen("masks")){
				selectWindow("masks"); 
				run("Close");
 				}


 				//
				// PROTEIN MASK MEASUREMENTS 
				// END
				// 

					
			// Fecha o Loop para cada célula
			}
		// Fecha o if check para ver se existiam células
		}

    // Close open images
	 while (nImages>0) { 
      selectImage(nImages); 
       close(); 
      }


   	 // Try to release memory to the system
		run("Collect Garbage");
   		run("Close All");
  	 	call("java.lang.System.gc"); 


//
// ANALYSIS END
//




    // LOOP 05: TIMEPOINTS - IMAGE ANALYSIS
    // END
    }

    // LOOP 04: FIELDS
    // END
    }
    	
    // LOOP 03: STACKS
    // END
    }
    	
    // LOOP 02: WELLS
    // END
    }


// 
print("Results Analyzed!");


// Exits BatchMode
setBatchMode(false);

// Prints the endtime
print("Script finished at:");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     print(TimeString+second);

// Saves the Summary
// SummaryName = dir  + File.separator + ".." + File.separator + "Results" + File.separator + "Summary_" + Valor1 + ".csv";
// selectWindow("Summary");  //select Log-window 
// saveAs("Text", SummaryName);
// run("Close");


// Saves the log
LogName = dir  + File.separator + ".." + File.separator + "Results" + File.separator + "Log_" + ".txt";
selectWindow("Log");  //select Log-window 
saveAs("Text", LogName);
run("Close");


// Open an dialog screen on script finish asking about to quit imagej
  Dialog.create("Script finished");
  Dialog.addMessage("Good luck!")
  Dialog.addCheckbox("Check this box to quit ImageJ", false);
  Dialog.show();
  CloseStatus = Dialog.getCheckbox();
  if (CloseStatus==true) run("Quit");


//////////////////////////////////////////////




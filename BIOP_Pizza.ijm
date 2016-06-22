/*
 * BIOP PTS Measurement Tool (Pizza Target Slice)
 * By Olivier Burri, EPFL - SV - PTECH BIOP
 * Last Edit: Dec 17th 2015
 * 
 * DESCRIPTION:
 * ------------
 * This tool performs analyses on a "Siemens Star"-like pattern image developed at the 
 * BioImaging and Optics Core Facility (BIOP) at the Ecole Polytechnique Fédérale de Lausanne
 * 
 * The approach is based on the analysis of concentric circles whose intensity profiles are
 * subjected to a 1D local maxima-minima detection.
 * 
 * To fully install this macro, you also require the BIOP Library 
 *    Install the BIOP Library by going to FIJI > UPDATE... > MANAGE UPDATE SITES
 *    And selected the PTBIOP Update site
 * 
 * Finally, this tool makes use of ActionBar from Jerôme Mutterer
 *    Install ActionBar by going to FIJI > UPDATE... > MANAGE UPDATE SITES
 *    And selected the IBMP-CNRS Update site
 * 
 *  For convenience of installation, we recommend that you use FIJI instead of the vanilla  ImageJ
 *  
 */



 
// Install the BIOP Library
call("BIOP_LibInstaller.installLibrary", "BIOP"+File.separator+"BIOPLib.ijm");


runFrom = "jar:file:BIOP/BIOP_PTS.jar!/BIOP_Pizza.ijm";

//////////////////////////////////////////////////////////////////////////////////////////////
// The line below is for debugging. Place this VSI file in the ActionBar folder within Plugins
//////////////////////////////////////////////////////////////////////////////////////////////
//runFrom = "/plugins/ActionBar/BIOP_Pizza.ijm"";


run("Action Bar",runFrom);
exit();

<codeLibrary>

function toolName() {
	return "BIOP Pizza Target Detection Settings";
}

/*
 * Detection settings to use for the analysis
 */
function detectionSettings(){
	names = newArray("Oversample Image", "Values in Microns", "From radius", "Till radius", "Increment", "Line Width", "Distance Allowed Variation Percent", "Dark Center", "Center Noise Tolerance", "User Draws ROI", "Fine-Tune Center", "Measure on Flattened Image", "Summary Options", "Max Allowed Delta Peak-to Peak", "Min Consecutive Good Distances", "Min Good Peaks Percent");
	types = newArray("c", "m", "n", "n", "n", "n", "n", "c", "n", "c", "c", "c", "m", "n", "n", "n");
	defaults = newArray(false, "", 20, 50, 1, 20, 25, false, true, 10, true, false, "", 0.010, 50, 50);
	promptParameters(names, types, defaults);
	
}

/*
 * Main function to process the image
 */
function processImage(){
	ori = getTitle();
	close("\\Others");

	idx = findRoisWithName("R:.*");
	if (idx.length >0) {
		roiManager("Select", idx);
		roiManager("Delete");
	}
	isOverSampled = getBool("Oversample Image");
	lineWidth= parseInt(getData("Line Width"));
	
	// Oversample the image to help the detection of local maxima
	if (isOverSampled) {
			name = getTitle();
			run("Scale...", "x=2 y=2 interpolation=Bilinear average create title=[Oversampled - "+name+"]");
	}

	// Gather image data
	name = getTitle();

	getDimensions(w,h,c,z,t);
	getVoxelSize(vx,vy,vz,U);
	setData("Unit", U);

	// Measure tolerance percentage to use for peak to peak detection, 
	// it should be more than 2% of the max-min intensity, as estimated manually
	getStatistics(area, mean, min, max, std);
	tolerance = 0.02 * (max-min);
	// Make sure it makes sense (Not smaller than 2 grey levels)
	if (tolerance < 2) tolerance = 2;
	print("Maxima Detection Tolerance: "+tolerance);


	// Detect center of pattern
	isUserROI = getBool("User Draws ROI");
	isTuneCenter = getBool("Fine-Tune Center");
	
	// OPTION Flatten large variations in intensity
	isMeasuredOnFlattened = getBool("Measure on Flattened Image");
	if (isMeasuredOnFlattened) {
		run("Select None");
		//Flatten low frequency intensity changes
		run("Duplicate...", "title=FlatField");
		run("32-bit");
		run("Gaussian Blur...", "sigma=4.0 scaled");
		getStatistics(area, mean, min, max);
		run("Divide...", "value="+max);
		run("Enhance Contrast", "saturated=0.35");
		
		imageCalculator("Divide create 32-bit", name,"FlatField");
		rename("Flattened - "+name);
		name = "Flattened - "+name;
	}



	// Start working on center
	selectImage(name);

	idx = findRoiWithName("Center");

	selectImage(name);
	if (isUserROI) {
		drawAnalysisRegion(false);
	} else if (idx == -1) {
		makePoint(getWidth()/2, getHeight()/2 );			
	}

	if (isTuneCenter) {
		findCenter();
	}


	
	/* 
	 *  Start Analysis from the center of the target, 
	 *  represented bz a point ROI
	 */
	//There should only be one
	getSelectionCoordinates(xpoints, ypoints);
	centerX = xpoints[0];
	centerY = ypoints[0];	
	
	// From here make concentric circles
	start    = parseFloat(getData("From radius"));
	end      = parseFloat(getData("Till radius"));
	step     = parseFloat(getData("Increment"));
	pctRange = parseFloat(getData("Distance Allowed Variation Percent"))/100;
	
	nPos  = round((end-start) / step)+1;
	
	for (i=0; i<nPos; i++) {
		selectImage(name);
		rad = ((i)*step+start );
		//makeCircleLine(centerX, centerY, diam/vx);
		print("Circle", centerX, centerY, rad/vx);
		cleanCircle(centerX, centerY, rad/vx);
		rName = "R: "+rad;
		Roi.setName(rName);
		roiManager("Add");
		data = measureSelection(ori, pctRange, lineWidth, tolerance);
		//Plot.add("line", xValues, yValues)

	}


}

function saveData() {
	// Save the image and ROI set
	saveRois("Open");

}

/*
 * Helper function to find the center of the target automatically
 */
function findCenter() {
	// Center point
	idx = findRoiWithName("Center");
	if (idx == -1) {
		cx = getWidth()/2;
		cy = getHeight()/2;
	} else {
		roiManager("Select", idx);
		getSelectionCoordinates(xcord,ycord);
		cx = xcord[0];
		cy = ycord[0];
		roiManager("Delete");
	}
	
	noise = getData("Center Noise Tolerance");
	
	// Make a centered rectangle around cx,cy to limit maxima finding
	wr = getWidth()/10;
	hr = getHeight()/10;
	makeRectangle(cx - wr/2, cy - hr/2, wr, hr);

	// Blur
	run("Gaussian Blur...", "sigma=3");
	
	// Use find maxima
	isDark = getBool("Dark Center");
	if (!isDark) {
		dark = "";
	} else {
		dark = "light";
	}
	
	run("Find Maxima...", "noise="+noise+" "+dark+" output=[Point Selection]");
	run("Undo");
	// Make sure we found something
	if (selectionType != 10) {
		makePoint(cx,cy);
	}
	Roi.setName("Center");
	roiManager("Add");
	// Cancel Gaussian Blur
	
	
	
}

/* 
 *  1D Local maxima detection
 */
function measureDistances(goodBadRangePct, circleRad) {
	run("Remove Overlay");
	getSelectionCoordinates(peakPosx, peakPosy);
	profile = getProfile();
	
	//Get dimensions of image to overlay good and bad peaks
	getDimensions(x,y,c,z,t);
	maxAllowed = (1+goodBadRangePct)*circleRad*sin(2*PI/360);
	goodDist = circleRad*sin(2*PI/360);
	minAllowed = maxOf(0,(1-goodBadRangePct)*circleRad*sin(2*PI/360));

	
	theY = y/6;
	ybad = 5*y/6;

	height=2;

	makeRectangle(5,ybad, minAllowed, height);
	run("Properties... ", "fill=blue");
	run("Add Selection...");
	
	makeRectangle(5+minAllowed+5,ybad, goodDist, height);
	run("Properties... ", "fill=green");
	run("Add Selection...");
	
	makeRectangle(5+minAllowed+goodDist+10,ybad, maxAllowed, height);
	run("Properties... ", "fill=red");
	run("Add Selection...");
	
	d = newArray(0);
	goodD = 0;
	badD = 0;
	avg = 0;
	// Count max consecutive good peaks
	consecutiveGood = 0;
	lastGood = 0;
	for(i=1; i<peakPosx.length;i++) {
		//print(consecutiveGood);
		distance=dist(peakPosx[i-1], peakPosy[i-1], peakPosx[i], peakPosy[i]);
		if (distance <= maxAllowed && distance >= minAllowed) {
			d = Array.concat(d,distance);
			goodD++;
			consecutiveGood++;
			makeRectangle(peakPosx[i-1]+0.5,theY,peakPosx[i]-peakPosx[i-1]-0.5,height );
			run("Properties... ", "fill=green");
			if(peakPosx[i-1] < profile.length) 
				avg+=profile[peakPosx[i-1]];
			
		} else {
			badD++;
			makeRectangle(peakPosx[i-1]+0.5,theY,peakPosx[i]-peakPosx[i-1]-0.5,height );
			
			if(distance > maxAllowed) {
				run("Properties... ", "fill=red");
			} else {
				run("Properties... ", "fill=blue");
			}
			if (consecutiveGood > lastGood) {
				lastGood = consecutiveGood;
			}
			consecutiveGood = 0;
		}
		run("Add Selection...");
	}

	if (consecutiveGood > lastGood)
		lastGood = consecutiveGood;
				
	Array.getStatistics(d, min, max, mean, stdDev);
	res = newArray(mean, stdDev, goodD, badD, lastGood, avg/goodD);
	return res;
}

function dist(x1,y1,x2,y2) {
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

/*
 * Parses typical acquisition names to get all the metadata
 */
function parseImageName(n) {
	//Example 20150330_IN2_63x1.4_white_damage 1.lsm
	data = newArray(-1,-1,-1,-1,-1,-1);
	pattern = "\\d{8}_.*_\\d{1,3}x\\d\\.\\d{1,2}_.*_.*\\s?\\d{0,5}\\..*";
	if(matches(n,pattern)) {
		// Pick up the pieces
		date = substring(n,0,8);

		n = substring(n,9);
		mic = substring(n,0, indexOf(n,"_"));
		
		
		n = substring(n,indexOf(n,"_")+1);
		obj =  substring(n,0, indexOf(n,"x")+1);

		n = substring(n,indexOf(n,"x")+1);
		na = substring(n,0,indexOf(n,"_"));
		
		n = substring(n,indexOf(n,"_")+1);
		wave = substring(n,0,indexOf(n,"_"));
		
		n = substring(n,indexOf(n,"_")+1);
		if (matches(n, ".*\\s.*")) {
			note = substring(n,0,lastIndexOf(n," "));
		} else {
			note = substring(n,0,lastIndexOf(n,"."));
		}

		data = newArray(date, mic, obj, na, wave, note);
	}

	return data;
	
}
function measureSelection(original_name, pctRange, lineWidth, tolerance) {
	getVoxelSize(vx,vy,vz,u);
	name = getTitle();
	isOversampled = "No";
	if (startsWith(name, "Oversampled"))	isOversampled = "Yes";
	
	radiusRaw = Roi.getName;

	// Try to parse the original name to use as columns
	metadata = parseImageName(original_name);

	radius = parseFloat(substring(radiusRaw, 3, lengthOf(radiusRaw)));

	// Find bright peaks and dark peacks
	findPeaks(true, lineWidth, tolerance);
	distanceBright = measureDistances(pctRange, radius/vx);
	close();
	selectImage(name);
	findPeaks(false, lineWidth, tolerance);
	distanceDark = measureDistances(pctRange, radius/vx);
	close();
	
	res = distanceBright;
	
	/*
	 * 0: mean
	 * 1: stdDec
	 * 2: nGoodPeaks
	 * 3: nBadPeaks
	 * 4: max number of consecutive good detections
	 * 5: avgIntensity
	 */
	prepareTable("Pattern Detector Detailed Results");
	n=nResults;
	setResult("Circle Radius ["+u+"]", n, radius);


	maxGoodDelta = 0.01;
	maxGoodSd = 0.01;
	//quality = maxGoodDelta / pow(res[0]*vx,2)+maxGoodDelta) * maxGoodSd /(pow(res[1],2)+maxGoodSd) * ;
	ep = radius*sin(2*PI/360);
	setResult("Label", n, name);
	setResult("Oversampled", n, isOversampled);
	if (metadata[0] != -1) {
		setResult("Date", n, metadata[0]);
		setResult("Machine", n, metadata[1]);
		setResult("Objective", n, metadata[2]);
		setResult("NA Objective", n, metadata[3]);
		setResult("Wavelength", n, metadata[4]);
		setResult("Note", n, metadata[5]);
	}
	setResult("Expected Mean Peak to Peak Distance ["+u+"]", n, ep);
	setResult("Mean Peak to Peak Distance ["+u+"]", n, res[0]*vx);
	setResult("StDev", n, res[1]);
	setResult("Expected FWHM ["+u+"]", n, 2*ep/3);
	setResult("FWHM ["+u+"]", n, 2*res[0]*vx/3);
	setResult("StDev", n, res[1]*vx/3);
	setResult("Total Peaks", n, (res[3]+res[2]));
	setResult("Good Distance [%]", n, (res[2]/(res[3]+res[2])*100));
	setResult("Max Consecutive Good Distances", n, res[4]);
	setResult("Michelson Contrast", n, (distanceBright[5]-distanceDark[5]) / (distanceBright[5]+distanceDark[5]));
	//setResult("Quality", n, quality);
	
	closeTable("Pattern Detector Detailed Results");
	selectWindow("Results");
	run("Close");

	res = newArray(ep, res[0]*vx, res[1], 2*ep/3, 2*res[0]*vx/3, res[1]*vx/3, (res[3]+res[2]), (res[2]/(res[3]+res[2])*100), res[4]);
	return res;
}

/*
 * Based on the given criteria, return a single line per image
 * containing the expected best fit
 */
function summarizeTable() {
	maxP2PDelta = parseFloat(getData("Max Allowed Delta Peak-to Peak"));
	minGoodCnt  = parseInt(getData("Min Consecutive Good Distances"));
	minDistPct  = parseFloat(getData("Min Good Peaks Percent"));
	// Check result table exists
	winds = getList("window.titles");
	for (i=0; i<winds.length; i++) {
		if (winds[i] == "Pattern Detector Detailed Results") {
			prepareTable("Pattern Detector Detailed Results");
		} // Otherwise assume that it's the current result table.
	}
	
	n=nResults;
	lastName = getResultLabel(0);
	lastOversampled = getResultString("Oversampled",0);
	//getVoxelSize(vx,vy,vz,u);

	// Get units somehow...
	u = getData("Unit");

	// Summary data
	SNameA 		= newArray(0);
	SContrastA 	= newArray(0);
	SradiusA	= newArray(0);
	SepA		= newArray(0);
	Sp2pAvgA	= newArray(0);
	Sp2pStdDevA	= newArray(0);
	SepFWHMA	= newArray(0);
	SFWHMA 		= newArray(0);
	SFWHMStdDevA= newArray(0);
	StotPeakA	= newArray(0);
	SgoodDPctA  = newArray(0);
	//SQA		= newArray(0);
	SmaxGoodCntA= newArray(0);
	SOversampledA= newArray(0);
	SContrastA  = newArray(0);
	dateA       = newArray(0);
	machA 		= newArray(0);
	objA 		= newArray(0);
	naA 		= newArray(0);
	waveA 		= newArray(0);
	noteA 		= newArray(0);

	isFullData = false;
	if(!isNaN(getResult("Date", 0))) {
		isFullData = true;
		lastDate = getResultString("Date", 0);
		lastMach = getResultString("Machine", 0);
		lastObj  = getResultString("Objective", 0);
		lastNa   = getResultString("NA Objective", 0);
		lastWave = getResultString("Wavelength", 0);
		lastNote = getResultString("Note", 0);
	}
	
	for(i=0; i<n; i++) {
		// Read name and results
		currentName = getResultLabel(i);
		radius 		= getResult("Circle Radius ["+u+"]", i);
		ep 			= getResult("Expected Mean Peak to Peak Distance ["+u+"]", i);
		p2pAvg		= getResult("Mean Peak to Peak Distance ["+u+"]", i);
		p2pStdDev	= getResult("StDev", i);
		epFWHM		= getResult("Expected FWHM ["+u+"]", i);
		FWHM 		= getResult("FWHM ["+u+"]", i);
		FWHMStdDev	= getResult("StDev", i);
		totPeak		= getResult("Total Peaks", i);
		goodDPct    = getResult("Good Distance [%]", i);
		//Q			= getResult("Quality", i);
		maxGoodCnt	= getResult("Max Consecutive Good Distances", i);
		oversampled = getResultString("Oversampled", i);
		currentContrast = getResult("Michelson Contrast", i);
		
		if (isFullData) {
			date = getResultString("Date", i);
			mach = getResultString("Machine", i);
			obj  = getResultString("Objective", i);
			na   = getResultString("NA Objective", i);
			wave = getResultString("Wavelength", i);
			note = getResultString("Note", i);
			
		}

		if (currentName != lastName) {

			//Parse current results and reset
			// Criteria: Expected Delta Peak to Peak <= 10nm
			// AND Max good detections > 50
			// AND goodpct > 75
			
			for(k=0; k<radiusA.length; k++) {
				//print(abs(epA[k] - p2pAvgA[k]), maxGoodCntA[k], goodDPct);
				
				if(abs(epA[k] - p2pAvgA[k]) < maxP2PDelta && maxGoodCntA[k] > minGoodCnt && goodDPctA[k] > minDistPct) {
					
					SNameA	    = Array.concat(SNameA, lastName);
					SOversampledA = Array.concat(SOversampledA, lastOversampled);
					SradiusA	= Array.concat(SradiusA, radiusA[k]);
					SepA		= Array.concat(SepA,epA[k]);
					Sp2pAvgA	= Array.concat(Sp2pAvgA,p2pAvgA[k]);
					Sp2pStdDevA	= Array.concat(Sp2pStdDevA,p2pStdDevA[k]);
					SepFWHMA	= Array.concat(SepFWHMA,epFWHMA[k]);
					SFWHMA 		= Array.concat(SFWHMA,FWHMA[k]);
					SFWHMStdDevA= Array.concat(SFWHMStdDevA,FWHMStdDevA[k]);
					StotPeakA	= Array.concat(StotPeakA,totPeakA[k]);
					SgoodDPctA  = Array.concat(SgoodDPctA,goodDPctA[k]);
					//SQA			= Array.concat(SQA,QA[k]);
					SmaxGoodCntA= Array.concat(SmaxGoodCntA,maxGoodCntA[k]);
					SContrastA = Array.concat(SContrastA, currentContrast);
					k=radiusA.length;

					if (isFullData) {
						dateA       = Array.concat(dateA,lastDate);
						machA 		= Array.concat(machA,lastMach);
						objA 		= Array.concat(objA,lastObj);
						naA 		= Array.concat(naA,lastNa);
						waveA 		= Array.concat(waveA,lastWave);
						noteA 		= Array.concat(noteA,lastNote);
					}
				}
			}
			
			// Reset
			lastName = currentName;
			lastOversampled = oversampled;
			if (isFullData) {
				lastDate = date;
				lastOversampled = oversampled;
				lastMach = mach;
				lastObj = obj;
				lastNa = na;
				lastWave = wave;
				lastNote = note;
			}
			
			radiusA		= newArray(0);
			epA			= newArray(0);
			p2pAvgA		= newArray(0);
			p2pStdDevA	= newArray(0);
			epFWHMA		= newArray(0);
			FWHMA 		= newArray(0);
			FWHMStdDevA	= newArray(0);
			totPeakA	= newArray(0);
			goodDPctA   = newArray(0);
			//QA			= newArray(0);
			maxGoodCntA	= newArray(0);
		}

		// Concat Results
			radiusA		= Array.concat(radiusA, radius);
			epA			= Array.concat(epA, ep);
			p2pAvgA		= Array.concat(p2pAvgA, p2pAvg);
			p2pStdDevA	= Array.concat(p2pStdDevA, p2pStdDev);
			epFWHMA		= Array.concat(epFWHMA, epFWHM);
			FWHMA 		= Array.concat(FWHMA, FWHM);
			FWHMStdDevA	= Array.concat(FWHMStdDevA, FWHMStdDev);
			totPeakA	= Array.concat(totPeakA, totPeak);
			goodDPctA   = Array.concat(goodDPctA, goodDPct);
			//QA			= Array.concat(QA, Q);
			maxGoodCntA	= Array.concat(maxGoodCntA, maxGoodCnt);

		
	}
	
	
	//Last Result
	for(k=0; k<radiusA.length; k++) {
		if(abs(epA[k] - p2pAvgA[k]) < maxP2PDelta && maxGoodCntA[k] > minGoodCnt && goodDPctA[k] > minDistPct) {
			SOversampledA = Array.concat(SOversampledA, lastOversampled);
			SNameA		= Array.concat(SNameA, lastName);
			SradiusA	= Array.concat(SradiusA, radiusA[k]);
			SepA		= Array.concat(SepA,epA[k]);
			Sp2pAvgA	= Array.concat(Sp2pAvgA,p2pAvgA[k]);
			Sp2pStdDevA	= Array.concat(Sp2pStdDevA,p2pStdDevA[k]);
			SepFWHMA	= Array.concat(SepFWHMA,epFWHMA[k]);
			SFWHMA 		= Array.concat(SFWHMA,FWHMA[k]);
			SFWHMStdDevA= Array.concat(SFWHMStdDevA,FWHMStdDevA[k]);
			StotPeakA	= Array.concat(StotPeakA,totPeakA[k]);
			SgoodDPctA  = Array.concat(SgoodDPctA,goodDPctA[k]);
			//SQA			= Array.concat(SQA,QA[k]);
			SmaxGoodCntA= Array.concat(SmaxGoodCntA,maxGoodCntA[k]);
		    SContrastA = Array.concat(SContrastA, currentContrast);
			k=radiusA.length;

			if (isFullData) {
				dateA       = Array.concat(dateA,lastDate);
				machA 		= Array.concat(machA,lastMach);
				objA 		= Array.concat(objA,lastObj);
				naA 		= Array.concat(naA,lastNa);
				waveA 		= Array.concat(waveA,lastWave);
				noteA 		= Array.concat(noteA,lastNote);
			}
					
		}
	}
	closeTable("Pattern Detector Detailed Results");



	//Make the summary table

	prepareTable("Pattern Detector Summary");
	f=nResults;
	for(k=0; k<SNameA.length; k++) {
		
		setResult("Label", f+k, SNameA[k]);
		setResult("Oversampled", f+k, SOversampledA[k]);
		if (isFullData) {
			setResult("Date", f+k, dateA[k]);
			setResult("Machine", f+k, machA[k]);
			setResult("Objective", f+k, objA[k]);
			setResult("NA Objective", f+k, naA[k]);
			setResult("Wavelength", f+k, waveA[k]);
			setResult("Note", f+k, noteA[k]);
		}
		
		setResult("Circle Radius ["+u+"]", f+k, SradiusA[k]);
		
		setResult("Expected Mean Peak to Peak Distance ["+u+"]", f+k, SepA[k]);
		setResult("Mean Peak to Peak Distance ["+u+"]", f+k, Sp2pAvgA[k]);
		setResult("StDev", f+k, Sp2pStdDevA[k]);
		setResult("Expected FWHM ["+u+"]", f+k, SepFWHMA[k]);
		setResult("FWHM ["+u+"]", f+k, SFWHMA[k]);
		setResult("StDev", f+k, SFWHMStdDevA[k]);
		setResult("Total Peaks", f+k, StotPeakA[k]);
		setResult("Good Distance [%]", f+k, SgoodDPctA[k]);
		//setResult("Quality", f+k, QA[k]);
		setResult("Max Consecutive Good Distances", f+k, SmaxGoodCntA[k]);
		setResult("Michelson Contrast", f+k, SContrastA[k]);
	}
		closeTable("Pattern Detector Summary");
		selectWindow("Results");
		run("Close");	
}

/*
 * Peak finder
 */
function findPeaks(darkBackground, lineWidth, tolerance){
	// Get the profile
	run("Straighten...", "line="+lineWidth);
	if(darkBackground) {
		rename("Bright");
	} else {
		rename("Dark");
	}
	
	getDimensions(x,y,c,z,t);
	getVoxelSize(vx,vy,vz,u);
	makeLine(0,y/2,x,y/2);
	profile = getProfile();


	
	x1 = 0;
	y1 = y/2;
	x2 = x;
	y2 = y/2;

	//Taken and adapted from 
	/*
	 * PeakFinder Tool
	 * N.Vischer, 10.11.13 23:29
     * http://simon.bio.uva.nl/objectj/examples/PeakFinder/peakfinder.html
	*/
	dMin = 1; //minimum separation distance (pixels)
	//tolerance = 1; //minimum separation amplitude
	priority = "Amplitude";
	prior = "A"; //choose"A". "L" or "R" (Amplitude, Left, Right)
	includeEnds = false; //endpoints are included though they are not peaks

	
	dx = x2 - x1;
	dy = y2 - y1;
	sin2 = sin(atan2(dy, dx));
	cos2 = cos(atan2(dy, dx));


	len = profile.length;
	if (darkBackground) {
		peakArr = Array.findMaxima(profile, tolerance);
	}
	else {
		peakArr = Array.findMinima(profile, tolerance);
	}
	nMaxima = peakArr.length;
	qualifiedArr = newArray(len);
	Array.fill(qualifiedArr, 1);
	nQualified = 0;
	if (prior != "A"){
		Array.sort(peakArr);
	}
	if (prior == "R"){
		Array.invert(peakArr);
	}
	for (jj = 0; jj < nMaxima; jj++){
		pos = peakArr[jj];
		if (qualifiedArr[pos] == 1){
			nQualified ++;
			if (dMin > 1)
			for (kk = pos - (dMin -1); kk <= pos + (dMin - 1); kk++){
				if (kk >= 0 && kk < len){
					qualifiedArr[kk] = 0;
				}
			}
		}
		else
			peakArr[jj] = -1;
	}
	Array.sort(peakArr);
	Array.invert(peakArr);
	peakArr = Array.trim(peakArr, nQualified);
	if (prior != "R")
		Array.invert(peakArr);
	nVertices = nQualified + 2;//include end points for now
	arr4x = newArray(nVertices);
	arr4y = newArray(nVertices);
	for (jj = 0; jj < nQualified; jj++){
		arr4x[jj + 1] = x1 + cos2 * peakArr[jj];
		arr4y[jj + 1] = y1 + sin2 * peakArr[jj];
	}
	arr4x[0] = x1;
	arr4y[0] = y1;
	arr4x[nVertices -1] = x2;
	arr4y[nVertices -1] = y2;
	
	peakPosx = Array.slice(arr4x, 1, nVertices-2);
	peakPosy = Array.slice(arr4y, 1, nVertices-2);

	
	makeSelection("polyline", peakPosx, peakPosy);
	
}

/* 
 *  It is not sufficient to make a circle ROI and straighten it. To ensure optimal sampling,
 *  we generate the circle ourselves as a highly sampled line profile.
 */
function cleanCircle(cx,cy,r) {
	nSteps = 360*2;
	increment = 2*PI/nSteps;
	px = newArray(nSteps+1);
	py = newArray(nSteps+1);
	
	for(i=0; i<=nSteps; i++) {
		px[i] = r*cos(i*increment)+cx;
		py[i] = r*sin(i*increment)+cy;
	}

	makeSelection("line", px,py);
}


function drawAnalysisRegion(showConfirm) {
	getVoxelSize(vx,vy,vz,U);
	// Center point
	
	idx = findRoiWithName("Center");
	if (idx == -1) {
		cx = getWidth()/2;
		cy = getHeight()/2;
	} else {
		roiManager("Select", idx);
		getSelectionCoordinates(xcord,ycord);
		cx = xcord[0];
		cy = ycord[0];
		roiManager("Delete");
	}

	
	// Radii Start and End
	r1 = parseFloat(getData("From radius"))/vx;
	r2 = parseFloat(getData("Till radius"))/vx;
	print(r1);
	// Be smart estimating the initial radius, min peak 2 peak 
	// p2p = radius1*sin(2*PI/360)
	// radius = p2p / sin(2*PI/360)
	if (isNaN(r1)) r1 = 2 / sin(2*PI/360) * 0.75;
	if (isNaN(r2)) r2 = r1*3;
	
	params = newArray(cx,cy,r1,r2, 0, vx, U);
	
	shift=1;
	ctrl=2; 
	rightButton=4;
	alt=8;
	leftButton=16;
	
	x2=-1; y2=-1; z2=-1; flags2=-1;
	getCursorLoc(x, y, z, flags);
	Overlay.remove;
	makeOneObject(params);
	Overlay.show();
	setOption("DisablePopupMenu", true);
	setTool("hand");
	
	while (flags&rightButton==0) {
		getCursorLoc(x, y, z, flags);
		if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
			params = identifyAndDrawPosition(x,y, params);
		}
		x2=x; y2=y; z2=z; flags2=flags;
		wait(10);
	}
	
	// Save the resulting values for cx,cy,r1,r2
	makePoint(params[0],params[1]);
	Roi.setName("Center");
	roiManager("Add");
	setData("From radius", d2s(params[2]*vx,0));
	setData("Till radius", d2s(params[3]*vx,0));

	if (showConfirm)
		showMessage("Extents Updated");
		wait(300);
}


function identifyAndDrawPosition(x,y, params) {
	zoom = getZoom();
	
	nearPx = 40/zoom;
	transX = false;
	transY = false;
	r1change = false;
	r2change = false;
	

	cx = params[0];
	cy = params[1];
	r1 = params[2];
	r2 = params[3];
	highlight = params[4];
	cal = params[5];
	u  = params[6];
	
	d = dist(x,y,cx,cy);
	if(x-nearPx/2 <= cx && x+nearPx/2 >= cx) {
		transX = true;
	}
		
	if(y-nearPx/2 <= cy && y+nearPx/2 >= cy)
		transY = true;

	if(d >= r1-nearPx/2 && d<= r1+nearPx/2)
		r1change = true;

	if(d >= r2-nearPx/2 && d<= r2+nearPx/2)
		r2change = true;


	res = 0;
	
	if(transX || transY) res =  1;
	if(r1change) res = 2;
	if(r2change) res = 3;

	if (flags&leftButton !=0 && res ==1) { cx = x; cy = y; }
	
	if (flags&leftButton !=0 && res ==2) r1 = d;
	
	if (flags&leftButton !=0 && res ==3) r2 = d;

	params = newArray(cx,cy,r1,r2, res, cal, u);
	
	if (res != 0) { Overlay.remove; makeOneObject(params); Overlay.show;}
	return params;
}

function dist(x1,y1,x2,y2) {
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

function makeOneObject(paramsArray) {
	
	// Draw a rectangle based on an angle
	centerX = paramsArray[0];
	centerY = paramsArray[1];
	radius1 = paramsArray[2];
	radius2 = paramsArray[3];
	h       = paramsArray[4];
	cal     = paramsArray[5];
	unit    = paramsArray[6];

	lw = 10;
	txtr1 = "P2P "+d2s(radius1*sin(2*PI/360)*cal,1)+" "+unit;
	txtr2 = "P2P "+d2s(radius2*sin(2*PI/360)*cal,1)+" "+unit;

	
	setColor("Cyan");
	if (h==1) {setColor("Orange");}
	setLineWidth(2);
	Overlay.drawLine(centerX-20, centerY, centerX+20, centerY);
	Overlay.drawLine(centerX, centerY-20, centerX, centerY+20);
	Overlay.show
	setColor("Cyan");
	
	if (h==2) {setColor("Orange");}
	setLineWidth(lw);
	Overlay.drawEllipse(centerX-radius1, centerY-radius1, 2*radius1, 2*radius1);
	// Draw text over the first ellipse
	// Get width of text
	setFont("SanSerif", 24);
	w1 = getStringWidth(txtr1);
    Overlay.drawString(txtr1, centerX-w1/2,centerY - radius1-5);
	setColor("Cyan");

	if (h==3) {setColor("Orange");}
	setLineWidth(lw);
	Overlay.drawEllipse(centerX-radius2, centerY-radius2, 2*radius2, 2*radius2);
	// Draw text over the first ellipse
	// Get width of text
	setFont("SanSerif", 24);
	w2 = getStringWidth(txtr2);
    Overlay.drawString(txtr2, centerX-w2/2,centerY - radius2-5);
	setColor("Cyan");
	
}
</codeLibrary>

/*
 * ActionBar Buttons to make the interface pretty
 */
<text><html><font size=FS color=#0C2981>Parameters
<line>
<button>
label=Save Parameters
arg=<macro>
saveParameters();
</macro>

<button>
label=Load Parameters
arg=<macro>
loadParameters();
</macro>
</line>


<line>
<button>
label=Select Folder
icon=noicon
arg=<macro>
//Open the file and parse the data
openParamsIfNeeded();
setImageFolder("Select Working Folder");
getSaveFolder();
</macro>
</line>

<line>
<button>
label=Detection Settings
icon=noicon
arg=<macro>
detectionSettings();
</macro>
</line>

<text><html><font size=3 color=#0C2981> Image
<line>
<button>
label= Select to Open
icon=noicon
arg=<macro>
selectImageDialog();
</macro>

<button>
label=Close Others 
icon=noicon
arg=<macro>
close("\\Others");
</macro>

</line>

<text><html><font size=3 color=#0C2981> Process

<line>
<button>
label= Calibrate Image
icon=noicon
arg=<macro>
	Dialog.create("Pixel Size Calculator");
	Dialog.addNumber("Magnification", 0, 0, 3, "x");
	Dialog.addNumber("C-mount", 1.0, 1, 3, "x");
	Dialog.addNumber("CCD pixel size", 6.45, 2, 4, "um");
	Dialog.addNumber("Binning", 1, 0, 3, "");	
	Dialog.show();
	
	mag=Dialog.getNumber();
	cm=Dialog.getNumber();
	ccd=Dialog.getNumber();
	bin=Dialog.getNumber();
	
	pixelsize=(ccd*bin)/(mag*cm);
	
	print ("The pixel size is "+pixelsize+" microns");
	
	run("Properties...", "   unit=micron pixel_width="+pixelsize+" pixel_height="+pixelsize);
</macro>

</line>


<line>
<button>
label= Current image
icon=noicon
arg=<macro>
processImage();
saveData();
</macro>

<button>
label= Choose Folder
icon=noicon
arg=<macro>
dir = getDirectory("Folder with Files");
setData("Image Folder", dir);
nI = getNumberImages();
images = getImagesList();
images = Array.sort(images);
for (i=0; i< nI; i++) {
	open(dir+images[i]);
	openRoiSet(dir+images[i]);
	processImage();
	saveData();
	run("Close All");
}
</macro>

</line>
<line>
<button>
label= Check Peaks Current ROI
icon=noicon
arg=<macro>
getVoxelSize(vx,vy,vz,u);
pctRange = parseFloat(getData("Distance Allowed Variation Percent"))/100;
lineWidth= parseInt(getData("Line Width"));

// Measure tolerance percentage to use, it should be 2% of the max-min intensity, as estimated manually
getStatistics(area, mean, min, max, std);
tolerance = 0.02 * max-min;
// Make sure it makes sense (Not smaller than 2 grey levels)
if (tolerance < 2) tolerance = 2;
	

name = getTitle();
radRaw = Roi.getName;
radius = parseFloat(substring(radRaw, 3, lengthOf(radRaw)));

findPeaks(true, lineWidth);
distanceBright = measureDistances(pctRange, radius/vx);
	

</macro>
</line>

<line>
<button>
label= Summarize Table
icon=noicon
arg=<macro>
summarizeTable();
</macro>
</line>

<line>
<button>
label= Set Extents
icon=noicon
arg=<macro>
drawAnalysisRegion(true);
</macro>
</line>
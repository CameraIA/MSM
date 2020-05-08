/*
 * Postprocessing after segmentation and before metrics as in the Journal of Sync Radia 2017
 * 
 * Sep 19th 2017
 * CAMERA Image Processing
 */
 
var inputpath = "/Users/dani/Dropbox/prog/ImageJ/ALS/groundtruth/";
var inputfile = "";
var outputwatershed = "";
var outputvoronoi = "";

macro "getMask"{
	run("Close All");

	getParameters();
	//setBatchMode(true);
	getWatershed();
	//setBatchMode(false);
}

function getParameters()
{
	waitForUser( "Setting up the experiment","Select a binary image in exp folder using the next dialog");
	open();
	inputpath = getInfo("image.directory");
	inputfile = getInfo("image.filename");
	inputfile = split(inputfile, "_");  //assumes that the prefix comes before _
	inputfile = inputfile[0];
	print(inputfile);
	run("Close All");

	Dialog.create("Get particles watershed from images:");
	Dialog.addString("Input path:",inputpath,80);
	Dialog.addString("Output path for watershed:",inputpath+"watershed/",80);
	Dialog.addString("Output path for voronoi:",inputpath+"voronoi/",80);
	Dialog.show();

	inputpath = Dialog.getString();
	outputwatershed = Dialog.getString();
	outputvoronoi = Dialog.getString();
	if (!File.exists(outputwatershed))
		File.makeDirectory(outputwatershed);
	if (!File.exists(outputvoronoi))
		File.makeDirectory(outputvoronoi);
}

function getWatershed(){
	print("Processing data in "+inputpath);
	print("Output for watershed is in "+outputwatershed);
	print("Output for voronoi is in "+outputvoronoi);
	run("Image Sequence...", "open=["+inputpath+"] increment=1 scale=100 file="+inputfile+" or=[] sort");
    
	//run("Image Sequence...", "open="+inputpath+inputfile+" sort");
	rename("input");

	for (i=0; i<nSlices; i++) {
		setSlice(i+1);
		run("Duplicate...", "title=slice");
		run("Watershed");
		wait(100); 
		saveAs("TIFF", outputwatershed+"watershed_"+IJ.pad(i,4)+".tif");
		//run("Invert");
		run("Voronoi");
		wait(100); 
		setThreshold(1,255);
		run("Convert to Mask", "  black");
		saveAs("TIFF", outputvoronoi+"voronoi_"+IJ.pad(i,4)+".tif");
		close();
		
	}
	run("Close All");
}

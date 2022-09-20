title=getTitle;
subtitle=substring(title,5,lastIndexOf(title,"."));
input=getInfo("image.Directory");

run("Duplicate...", "duplicate");

run("Subtract Background...", "rolling=200 light sliding");
run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=1024 maximum=3 mask=*None*");
setTool("line");
	beep();  
	
	waitForUser("Draw a horizontal line between the start and end positions, then click OK.");
//	roiManager("Add");
	getLine(x1, y1, x2, y2, lineWidth);
	angle=floor(-180*atan((y2-y1)/(x2-x1))/PI);
run("Rotate... ", "angle="+toString(angle)+" grid=1 interpolation=Bicubic fill enlarge");	
setTool("freehand");
	beep();  
	waitForUser("Select a bounding line for the region of boundaries to be analyzed, then click OK.");

run("Duplicate...", "duplicate");
run("8-bit");
run("Select All");
run("Enhance Contrast...", "saturated=0.3 normalize equalize");
run("Auto Local Threshold", "method=Phansalkar radius=10 parameter_1=0 parameter_2=0 white");
run("Restore Selection");
run("Make Inverse");
run("Colors...", "foreground=white background=black selection=yellow");
run("Fill", "slice");
run("Select All");
run("Duplicate...", " ");
run("Invert");

saveAs("Tiff", input+"sample"+subtitle+".tif");
selectWindow("sample"+subtitle+".tif");

run("Simple Neurite Tracer", " ");

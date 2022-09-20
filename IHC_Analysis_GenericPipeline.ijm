input = getDirectory("Input directory"); //select the folder including nd2 images of the data to be analyzed
list = getFileList(input);
ll=list.length;
for (i = 0; i < ll; i++) {
	open(list[i]);
}
for (im=1; im<ll+1; im++){
	if (im==1) {
	abodyarr=newArray("Nuclei");
	Dialog.create("How many stainings in this experiments?");
	Dialog.addNumber("# of antibodies:", 2);
	beep();  
	Dialog.show();
	abodyc=Dialog.getNumber();
	Dialog.create("Which channel is the nuclear marker?");
	Dialog.addNumber("Nuclei Channel #:", 1);
	beep();  
	Dialog.show();
	charr=newArray(Dialog.getNumber());
	mipchannels=newArray(abodyc);
	zchannels=newArray(abodyc);
	for (i = 0; i < abodyc; i++) {
		strray=newArray("Antibody/protein #",toString(i+1),":");
		diastr=String.join(strray,"");
		Dialog.create(diastr);
		Dialog.addString("Name:", "ppERK");
		Dialog.addNumber("Channel:", 2); 
		Dialog.addChoice("Analysis:", newArray("MIP","Single Layer"));
		beep();  
		Dialog.show();	
		abodyarr=Array.concat(abodyarr,Dialog.getString());
		charr=Array.concat(charr,Dialog.getNumber());
		analysis = Dialog.getChoice();
		if (analysis=="MIP") {
			mipchannels[i]=1;
			zchannels[i]=0;}
		if (analysis=="Single Layer") {
			mipchannels[i]=0;
			zchannels[i]=1;}
	}
	}
		
  
input=getInfo("image.Directory");
title=getTitle;
extract=".nd2";
subtitle=substring(title,0,indexOf(title,extract));
folder=input+subtitle+"/";
File.makeDirectory(folder);

run("Duplicate...", "title=NewStack.tif duplicate"); 
Stack.getDimensions(width, height, channels, slices, frames); 
run("Despeckle");

run("Duplicate...", "duplicate channels=1");
run("32-bit");
run("ROF Denoise", "theta=25");
run("16-bit");
for (k=1;k<=slices;k++){
	setSlice(k);
run("Subtract Background...", "sliding=20");
run("Enhance Local Contrast (CLAHE)", "blocksize=31 histogram=512 maximum=3 mask=*None*");
}
setOption("BlackBackground", true);
run("Make Binary", "method=IJ_IsoData background=Dark calculate black");
run("16-bit");
a=4095/255;
run("Multiply...", "value=a stack");
saveAs("Tiff", folder+"Mask.tif");

run("Duplicate...", "duplicate");
run("Invert", "stack");
run("Subtract...", "value=61440 stack");
saveAs("Tiff", folder+"Cyto_Mask.tif");

for (c = 0; c < abodyc; c++) {
folder1=folder+"/"+abodyarr[c+1]+"/";
File.makeDirectory(folder1);

selectWindow("NewStack.tif"); 
dupstra=newArray("duplicate channels=",toString(charr[c+1]));
dupstr=String.join(dupstra, "");
run("Duplicate...", dupstr);
run("Remove Outliers...", "radius=2 threshold=50 which=Bright stack");
saveAs("Tiff", folder1+abodyarr[c+1]+".tif");
imageCalculator("Min create stack", abodyarr[c+1]+".tif","Mask.tif");
selectWindow("Result of "+abodyarr[c+1]+".tif");
saveAs("Tiff", folder1+abodyarr[c+1]+"_Masked.tif");

imageCalculator("Min create stack", abodyarr[c+1]+".tif","Cyto_Mask.tif");
selectWindow("Result of "+abodyarr[c+1]+".tif");
saveAs("Tiff", folder1+abodyarr[c+1]+"_Cyto_Masked.tif");
}
	
a0=0;
side=0;
selectWindow(abodyarr[1]+"_Cyto_Masked.tif");
  Dialog.create("How is that sample mounted?");
  Dialog.addChoice("Mounting:", newArray("Lateral","Flat","Skip"));
  beep();  
  Dialog.show();
  mount = Dialog.getChoice();
  if (mount=="Flat") {
  sides=2;
  }
  if (mount=="Lateral") {
  sides=1;
  }
  if (mount=="Skip") {
  sides=0;
  }
  while (side<sides) {
  	setTool("polyline"); 
  	beep();  
  	waitForUser("Select Line To Analyze, then click OK.");
  	roiManager("Add");
  	
  	for (c = 0; c < abodyc; c++) {
  		folder1=folder+"/"+abodyarr[c+1]+"/";  		
  		if (zchannels[c]==1){
  			for (j=1; j<slices+1; j++) {
			selectWindow(abodyarr[c+1]+"_Masked.tif");	
			setTool("polyline");  
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_"+toString(j)+toString(side)+".csv");

			selectWindow(abodyarr[c+1]+"_Cyto_Masked.tif");	
			setTool("polyline");  
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_Cyto_"+toString(j)+toString(side)+".csv");
			

			selectWindow(abodyarr[c+1]+".tif");	
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"unmasked_"+toString(j)+toString(side)+".csv");

			
			selectWindow("Mask.tif"); 
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Mask_"+toString(j)+toString(side)+".csv");

			selectWindow("Cyto_Mask.tif"); 
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(j)+toString(side)+".csv");
  		}
  		}
  		if (mipchannels[c]==1){			
  		a0=0;
  		for (j=1; j<slices+1; j++) {
  			a=floor((j-3)/3)+1;	
  			if (a>a0) {
  				str=a*3-1;
  				stp=a*3+1;
		if (stp<slices+1) {
			selectWindow(abodyarr[c+1]+"_Masked.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abodyarr[c+1]+"_Masked.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_"+toString(j+1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_"+toString(j-1)+toString(side)+".csv");
			selectWindow("MAX_"+abodyarr[c+1]+"_Masked.tif"); 
			close();

			selectWindow(abodyarr[c+1]+"_Cyto_Masked.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abodyarr[c+1]+"_Cyto_Masked.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_Cyto_"+toString(j+1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_Cyto_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_Cyto_"+toString(j-1)+toString(side)+".csv");			
			selectWindow("MAX_"+abodyarr[c+1]+"_Cyto_Masked.tif"); 
			close();
			
			selectWindow("Mask.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_Mask.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Mask_"+toString(j+1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_Mask_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_Mask_"+toString(j-1)+toString(side)+".csv");
			selectWindow("MAX_Mask.tif"); 
			close();

			selectWindow("Cyto_Mask.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_Cyto_Mask.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(j+1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(j-1)+toString(side)+".csv");			
			selectWindow("MAX_Cyto_Mask.tif"); 
			close();
			

			selectWindow(abodyarr[c+1]+".tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abodyarr[c+1]+".tif");  
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_unmasked_"+toString(j+1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_unmasked_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_unmasked_"+toString(j-1)+toString(side)+".csv");
			selectWindow("MAX_"+abodyarr[c+1]+".tif"); 
			close();	
			
			}
		}
		a0=floor((j-3)/3)+1;
  		}
  		saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_1"+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_"+toString(slices)+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_Cyto_"+"_1"+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_Cyto_"+"_"+toString(slices)+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_Mask_1"+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_Mask_"+toString(slices)+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_Cyto_Mask_1"+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(slices)+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_unmasked_1"+toString(side)+".csv");
  		saveAs("Measurements", folder1+"Values_"+abodyarr[c+1]+"_unmasked_"+toString(slices)+toString(side)+".csv");  		
		}
		}
		roiManager("Delete");	
		side=side+1;
	}
	
	selectWindow(title);
	close();
	selectWindow("Mask.tif");
	close();
	selectWindow("Cyto_Mask.tif");
	close();
	for (c = 0; c < abodyc; c++) {
		selectWindow(abodyarr[c+1]+".tif");
		close();
		selectWindow(abodyarr[c+1]+"_Masked.tif");
		close();
		selectWindow(abodyarr[c+1]+"_Cyto_Masked.tif");
		close();
	}
	selectWindow("NewStack.tif"); 
	close();
}
	
wlist = getList("window.titles"); 
for (i=0; i<wlist.length; i++){      
	window = wlist[i]; 
	selectWindow(window); 
	run("Close"); 
	} 

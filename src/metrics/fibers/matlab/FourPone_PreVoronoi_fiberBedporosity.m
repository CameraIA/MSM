clear

slice = 1;
srcFiles = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\*.tif'); % *Original Image* Gives the location of the images stored as the variable "srcFiles"
filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\',srcFiles(slice).name); % *Original Image* this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
Imginfo = imfinfo(filename); % all the images have the same hieght and width so in lines 8-9 I store these values as h and w.
I = imread(filename); 
w = Imginfo.Width;
h = Imginfo.Height;
z = length(srcFiles); % number of images in the folder

z=40;
%The purpose of this code is to create an image for the voronoi tesselation
%that I will do in FIJI

for slice = 1 : z
    %Get watershedded fiber images (NOT THE COMPARED ONES BC THOSE MERGE
    %THE FIBERS BACK TOGETHER TOO MUCH)
    srcFiles = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Fiber Finding\2_Fiber_Compare_1 Output\2p05_bw_watershed_of_inverse_inversed_inFIJI\*.tif'); 
    filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Fiber Finding\2_Fiber_Compare_1 Output\2p05_bw_watershed_of_inverse_inversed_inFIJI\',srcFiles(slice).name); 
    Imginfo = imfinfo(filename); % all the images have the same hieght and width so in lines 8-9 I store these values as h and w.
    fibers = imread(filename); 
    fibers = im2double(fibers); %converts image to double precision
    
    %Get NotCMC images
    srcFiles2 = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4_notCMC_compare_1\bw\*.tif'); % *BW output of Two_fiber_compare_1.m* Gives the location of the images stored as the variable "srcFiles"
    filename2 = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4_notCMC_compare_1\bw\',srcFiles2(slice).name); % *BW output of Two_fiber_compare_1.m* this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
    NotCMC = imread(filename2); 
    NotCMC = im2double(NotCMC);
    
    %Next five lines output the image I will voronoi in FIJI
    fibersWnotCMC = fibers;
    fibersWnotCMC(NotCMC==1) = 1;
    
    fibersWnotCMCI = imcomplement(fibersWnotCMC);
    
    f = sprintf('%04d',slice); %assigns index for naming images in sequence       
    imwrite(fibersWnotCMCI,['G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p1_fibersWnotCMC_invertedForVoronoi\4p1_fibWnotCMCI_', f,'.tif'],'tif'); 
    
    %Remainder of code calculates the fiber bed porosity
    
    %Make vector of number of fiber pixels in each slice
    fibers_fixed = fibers;
    fibers_fixed(NotCMC==1) = 0;
    fibers1 = find(fibers_fixed==1);
    fibers2 = length(fibers1);
    slice_fibers(slice) = fibers2;
    
    f = sprintf('%04d',slice); %assigns index for naming images in sequence       
    imwrite(fibers_fixed,['G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p1_watershedded_fibers_fixed\4p1_fibers_fixed_', f,'.tif'],'tif');
    
    %Make vector of number of matrix + pore pixels in each slice
    matrixPore1 = find(fibersWnotCMCI==1);
    matrixPore2 = length(matrixPore1);
    slice_matrixPore(slice) = matrixPore2;


    slice
end

Table_slice = table(slice_fibers,slice_matrixPore);
writetable(Table_slice,'G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p1_FiberBedPorosityAnalysis\slice_pixelvalues.xlsx');

%Calculate area fractions on a per slice basis (A indicates area)
A_fiber_bed_porosity = (slice_matrixPore)./(slice_matrixPore + slice_fibers); %This is epsilon in all of my calculations

Table_areaCalcs = table(A_fiber_bed_porosity);
writetable(Table_areaCalcs,'G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p1_FiberBedPorosityAnalysis\AreaAnalysis_calcs.xlsx');

%Calculate number of pixels for entire 3D volume
pixels_fibers = sum(slice_fibers);
pixels_matrixPore = sum(slice_matrixPore);

Table_volCalcs = table(pixels_fibers,pixels_matrixPore);
writetable(Table_volCalcs,'G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p1_FiberBedPorosityAnalysis\VolPixel_Calcs.xlsx');

%Calculate various volume fractions
fiber_bed_porosity = (pixels_matrixPore)/(pixels_matrixPore + pixels_fibers); %This is epsilon in all of my calculations

Table_volCalcs2 = table(fiber_bed_porosity);
writetable(Table_volCalcs2,'G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p1_FiberBedPorosityAnalysis\VolAnalysis_Calcs.xlsx');

clear

slice = 1;
srcFiles = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\*.tif'); % *Original Image* Gives the location of the images stored as the variable "srcFiles"
filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\',srcFiles(slice).name); % *Original Image* this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
Imginfo = imfinfo(filename); % all the images have the same hieght and width so in lines 8-9 I store these values as h and w.
I = imread(filename); 
w = Imginfo.Width;
h = Imginfo.Height;
z = length(srcFiles); % number of images in the folder
NHOOD = ones(7,7); %Size of neighbrhood
z=40;

for slice = 1 : z
    % Not needed for this code: srcFiles2 = dir('G:\March 2015\Image Processing\121_2_novac_1p2cm_800ms_8bit_enhcont_crop\Fiber Finding\Fiber_Compare_1 Output\bw\*.tif'); % Gives the location of the images stored as the variable "srcFiles"
    % Not needed for this code: filename2 = strcat('G:\March 2015\Image Processing\121_2_novac_1p2cm_800ms_8bit_enhcont_crop\Fiber Finding\Fiber_Compare_1 Output\bw\',srcFiles2(slice).name); % this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
    % Not needed for this code: Ibw = imread(filename2); 
    
    filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\',srcFiles(slice).name); % *Original Image*
    I = imread(filename);
    I = im2double(I); %converts image to double precision
    % Not needed for this code: I(Ibw==255)=1;
    
    thresh = multithresh(I,5);
    seg_I = imquantize(I,thresh);
    RGB = label2rgb(seg_I); %RGB image of multithresh
    
    bw = zeros(h,w);
    bw(seg_I ==4) = 1; %Creates bw image of regions with mainly glass tube white
    bw = bwareaopen(bw,25,4); %Gets rid of small white things in above bw image
    
    B = ordfilt2(bw,25,NHOOD);
    B = bwareaopen(B,600,4);
    
    B =imcomplement(B); %This line and the next two lines fill the holes in the tube. This line inverts (complements) the bw image (now we want to get rid of small white objects)
    B =bwareaopen(B,400,4); %PLAY HERE, This line gets rid of small white objects
    B =imcomplement(B);
    
    %Not needed for this code:C = B;
    %Not needed for this code:C = imfill(B,'holes');
    %Not needed for this code:B(C==0) = 1;
    %Not needed for this code:B = bwareaopen(B,10000);
    
    f = sprintf('%04d',slice); %assigns index for naming images in sequence       
    imwrite(B,['G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\3_multithresh_1\3_bw_', f,'.tif'],'tif');
    slice    % tells the slice number in the command window so you know how far along the run is

end
    
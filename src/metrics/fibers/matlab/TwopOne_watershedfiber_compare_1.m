clear

slice = 1;
srcFiles = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\*.tif'); % *Here input BW images from 1_circles.m output*. Gives the location of the images stored as the variable "srcFiles"
filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\',srcFiles(slice).name); % *Here input BW images from 1_circles.m output*. this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
Imginfo = imfinfo(filename); % all the images have the same hieght and width so in lines 8-9 I store these values as h and w.
I = imread(filename); 
w = Imginfo.Width;
h = Imginfo.Height;
z = length(srcFiles); % number of images in the folder
r=11; %PLAY HERE,this is the total number of images that are compared in a single comparison. Except for endpoints of the stack, the middle image in the stack of r images is rewritten based on comparison with the r/2-0.5 images surrounding it. 
z=40;

for a=1:z %(not needed code: a=r/2+.5:z-r/2+.5)
    b=a-r/2+.5; %This and the next three lines deal with the first r/2+0.5 images that do not have any images to compare before them. this sets all of the negative slices equal to 1 for comparison. 
    if b<1
        b=1;
    end
    c=a+r/2-.5; %This and the next three lines deal with the last r/2+0.5 images that do not have any images to compare before them. this sets all of the slices>z equal to z for comparison.
    if c>z
        c=z;
    end
    srcFiles2 = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\*.tif'); % *Original image* Grabbing original images for use in the RGB image
    filename2 = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\',srcFiles2(a).name); % *Original image* Grabbing original images for use in the RGB image
    I4=imread(filename2); %I4=original image
    I4 = im2double(I4);
    RGB = zeros(h,w,3);
    RGB(:,:,1) = I4; %This and next two lines placing original image I4 as all three RGB layers to keep it greyscale
    RGB(:,:,2) = I4;
    RGB(:,:,3) = I4;
    I2=zeros(h,w);
    for slice = b:c %This for loop creates a sum of r images centered around image a. We then pick regions of high enough value (meaning a fiber was identified in that location many times) to make a new bw image for slice a
        srcFiles = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Fiber Finding\2_Fiber_Compare_1 Output\2p05_bw_watershed_of_inverse_inversed_inFIJI\*.tif'); % *Here input BW images from One_circles.m output.* Gives the location of the images stored as the variable "srcFiles"
        filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Fiber Finding\2_Fiber_Compare_1 Output\2p05_bw_watershed_of_inverse_inversed_inFIJI\',srcFiles(slice).name); % *Here input BW images from One_circles.m output.* this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder          
        I = imread(filename);
        I = im2double(I);
        I2=I2+I;
    end
    t=c-b+1; %this is the number of images in the image sum
    I2=I2/t; %normalizes the sum of image so we can pick a percent certainty that a fiber is there
    bw =zeros(h,w);
    bw(I2>=.5)=1; %PLAY HERE, If the sum of images has > the specifed fraction of normalized white summed up, then that region is identified as fiber. This will overall tend to erode the fibers if too large an r is chosen. 
    
    bw2 = bwperim(bw); %creates RGB image that makes circles around fibers. 
    for i =1:h
        for j = 1:w
            if bw2(i,j) == 1;
                RGB(i,j,1) = 0;
                RGB(i,j,2) = 0;
                RGB(i,j,3) = 1;
            end
        end
    end
    bw1 = uint8(zeros(h,w));
    bw1(bw==1) = 255;
    
    f = sprintf('%04d',a); 
    imwrite(RGB,['G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Fiber Finding\2_Fiber_Compare_1 Output\2p1_RGB_watershed_compare\2_RGB_watershed_', f,'.tif'],'tif'); % *Output folder for 2_fiber_compare_1.m*
    imwrite(bw1,['G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Fiber Finding\2_Fiber_Compare_1 Output\2p1_bw_watershed_compare\2_bw_watershed_', f,'.tif'],'tif');% *Output folder for 2_fiber_compare_1.m*
    a  
end
    

    
    
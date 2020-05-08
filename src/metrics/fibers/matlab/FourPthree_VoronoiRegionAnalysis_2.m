%clear

slice = 1;
%srcFiles = dir('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\*.tif'); % *Original Image* Gives the location of the images stored as the variable "srcFiles"
%filename = strcat('G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\',srcFiles(slice).name); % *Original Image* this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
%Imginfo = imfinfo(filename); % all the images have the same hieght and width so in lines 8-9 I store these values as h and w.
%I = imread(filename); 
%w = 2082;%Imginfo.Width;
%h = 2073;%Imginfo.Height;
%z = length(srcFiles); % number of images in the folder
z = 2160;

slice_Mean_regions_area = zeros(z,1);
slice_Std_regions_area = zeros(z,1);
                   
slice_Mean_Frac_unocc = zeros(z,1);
slice_Std_Frac_unocc = zeros(z,1);
slice_GEV_k_Frac_unocc = zeros(z,1);

slice_Mean_Area_unocc = zeros(z,1);
slice_Std_Area_unocc = zeros(z,1);
slice_GEV_k_Area_unocc = zeros(z,1);

slice_Mean_area_weighted_normalized_frac_unocc_all = zeros(z,1);
slice_Std_area_weighted_normalized_frac_unocc_all = zeros(z,1);
slice_GEV_shape_factor_k = zeros(z,1);
                   
slice_alpha = zeros(z,1);
slice_alpha_unocc = zeros(z,1);

slice_Mean_porosity = zeros(z,1);
slice_Mean_porosity_term_forF = zeros(z,1);
slice_Mean_porosity_term_forK0 = zeros(z,1);
                   



for slice = 1 : z
    %Get Voronoi Images
    srcFiles = dir(name1); %*Original Image*  Gives the location of the images stored as the variable "srcFiles"
    filename = strcat(name2,srcFiles(slice).name); %*Original Image*  this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
    Imginfo = imfinfo(filename); % all the images have the same hieght and width so in lines 8-9 I store these values as h and w.
    I = imread(filename); 
    I = im2double(I); %converts image to double precision
    
                   
   %% [Mike Comment: Pulls in the Voronoi images, inverts it. Then creates a second grey-scale one.]
                   Changes to double to scale from 0-1.
    Voronoi = zeros(h,w);
    Voronoi(I==0)= 1;
    Voronoi2 = uint8(zeros(h,w));
    Voronoi2(Voronoi==1) = 255;
    
    %Get binary fiber watershedded compare image
    srcFiles2 = dir(name3); % *BW output of Two_fiber_compare_1.m* Gives the location of the images stored as the variable "srcFiles"
    filename2 = strcat(name4,srcFiles2(slice).name); % *BW output of Two_fiber_compare_1.m* this combines the folder that the images are in with the name of the individual images so that I can call up certain images. In line 3, slice = 1, so srcFiles(slice).name is the first image in the folder  
    If = imread(filename2); 
    If = im2double(If);
   
    %Make image showing voronoi overlapping fiber bed
    Voronoi_fibers = zeros(h,w);
    for i =1:h
        for j = 1:w
            if If(i,j) == 1;
                Voronoi_fibers(i,j) = 1;
            end
            if Voronoi(i,j) == 0;
                Voronoi_fibers(i,j) = 0.5;
            end
        end
    end
    
    
    % Determine connected regions of voronoi
    CC = bwconncomp(Voronoi,4); %This determines the connected components (regions) of the img
    L=bwlabel(Voronoi,4); %Labels connected components of matrix/pore
    regions = regionprops(CC,If,'Area','PixelValues','PixelIdxList'); %Useful for regionprops output and 
    %the next few lines: http://blogs.mathworks.com/steve/2009/02/27/using-ismember-with-the-output-of-regionprops/
    
    number_of_regions1 = [regions.Area];
    number_of_regions = length(number_of_regions1);
    
    %This next set of code dilates each individual region by 1, calculates
    %the perimeter length, half of which I will add to each entry in area
    %and as unoccupied area in frac_unocc
    %extra = zeros(1,number_of_regions);
    %SE = strel('square', 3);
    %for m = 1:number_of_regions
    %    Individ_regions = zeros(h,w);
    %    Individ_regions(regions(m).PixelIdxList) = 1;
    %    dilated = imdilate(Individ_regions,SE);
    %    perim = bwperim(dilated,4);
    %    sum_perim = sum(sum(perim));
    %    extra(m) = sum_perim; 
    %    
    %end
    
    %This section determines the areas and then average and stdev of the voronoi regions
    regions_area_all = transpose([regions.Area]); %vector of all regions in proper order
    regions_area_all = regions_area_all(2:number_of_regions); %region labeled #1 in list is the large region identifying the notCMC region. We do not want this. 
    regions_area_all = transpose(regions_area_all);
    
    Mean_regions_area = mean(regions_area_all);
    Std_regions_area = std(regions_area_all);
    
    slice_Mean_regions_area(slice)= Mean_regions_area;
    slice_Std_regions_area(slice)= Std_regions_area;
    
    %This section determines the number of occupied pixels in each voronoi
    %region
    Occupied = zeros(1,number_of_regions);
    for g = 1: number_of_regions
        regions_PixelValues = [regions(g).PixelValues];
        a1 = sum(regions_PixelValues);
        Occupied(g) = a1;
    end
    
    Occupied = Occupied(2:number_of_regions);
    
    %This section calculates the fraction unoccupied for each cell, and the
    %mean and stdev of those and the k value of the GEV fit. Then we save
    %the entry for each slice
    Frac_unocc_all = (regions_area_all-Occupied)./regions_area_all;
    Mean_Frac_unocc = mean(Frac_unocc_all);
    Std_Frac_unocc = std(Frac_unocc_all);
    
    Frac_unocc_all2 = transpose(Frac_unocc_all);
    pd1=fitdist(Frac_unocc_all2,'Generalized Extreme Value');
    
    slice_Mean_Frac_unocc(slice)= Mean_Frac_unocc;
    slice_Std_Frac_unocc(slice) = Std_Frac_unocc;
    slice_GEV_k_Frac_unocc(slice)=pd1.k;
                   
                   %% [Mike Comment:: Have translated up to here. 12:37pm Sept. 14.]
    
    %This section calculates the area weighted fraction unoccupied for each
    %cell (equivalently the unoccupied area of each cell), and the mean,
    %stdev and k value of the GEV fit
    Area_unocc_all = regions_area_all.*Frac_unocc_all;
    Mean_Area_unocc = mean(Area_unocc_all);
    Std_Area_unocc = std(Area_unocc_all);
    
    Area_unocc_all2=transpose(Area_unocc_all);
    pd2=fitdist(Area_unocc_all2,'Generalized Extreme Value');
    
    slice_Mean_Area_unocc(slice)= Mean_Area_unocc;
    slice_Std_Area_unocc(slice) = Std_Area_unocc;
    slice_GEV_k_Area_unocc(slice)=pd2.k;
    
    %This section calculates the area weighted normalized (by total area)
    %fraction unnocupied for each cell AND then calculates the shape
    %parameter (k) of a Generalized Extreme Value fit to the distribution
    %of area_weighted_normalized_frac_unocc_all. 
    %Websites:
    %http://www.mathworks.com/help/stats/distributionfitting-app.html?searchHighlight=histogram%20fit,
    %http://www.mathworks.com/help/stats/fitdist.html
    totalArea = sum(regions_area_all);
    area_weighted_normalized_frac_unocc_all = (regions_area_all.*Frac_unocc_all)./totalArea;
    slice_Mean_area_weighted_normalized_frac_unocc_all(slice) = mean(area_weighted_normalized_frac_unocc_all);
    slice_Std_area_weighted_normalized_frac_unocc_all(slice) = std(area_weighted_normalized_frac_unocc_all);
    
    area_weighted_normalized_frac_unocc_all2 = transpose(area_weighted_normalized_frac_unocc_all);
    pd = fitdist(area_weighted_normalized_frac_unocc_all2,'Generalized Extreme Value');
    slice_GEV_shape_factor_k(slice) = pd.k;
    
    %Transpose the per slice data that contains information about EACH
    %voronoi cell
    regions_area_all = transpose(regions_area_all);
    Frac_unocc_all = transpose(Frac_unocc_all);
    area_weighted_normalized_frac_unocc_all = transpose(area_weighted_normalized_frac_unocc_all);
    Area_unocc_all = transpose(Area_unocc_all);
    
    %This section calculates the alpha parameter used in source 209
    slice_alpha(slice) = mean(regions_area_all.*regions_area_all)/(mean(regions_area_all).*mean(regions_area_all));
    
    %This section calculates the alpha parameter used in source 209, but
    %instead of using the entire cell area, uses only the unoccupied cell
    %area to account for the fact that the fibers are different sizes
    
    slice_alpha_unocc(slice) = mean(Area_unocc_all.*Area_unocc_all)/(mean(Area_unocc_all).*mean(Area_unocc_all));
    
    %Calculate the average porosities for each slice using area weighted
    %sums:
    Frac_unocc_all2=Frac_unocc_all;
    Frac_unocc_all2(Frac_unocc_all2>=0.98) = NaN;
    
    Mean_porosity = nansum(regions_area_all.*Frac_unocc_all2)./totalArea;
    slice_Mean_porosity(slice) = Mean_porosity;
    
    one1 = ones(number_of_regions,1);
    one = one1(2:number_of_regions);
    Mean_porosity_term_forF = nansum(regions_area_all.*Frac_unocc_all2./(one-Frac_unocc_all2))./totalArea;
    slice_Mean_porosity_term_forF(slice) = Mean_porosity_term_forF;
    
    Mean_porosity_term_forK0 = nansum(regions_area_all.*Frac_unocc_all2.^2./((one-Frac_unocc_all2).^2))./totalArea;
    slice_Mean_porosity_term_forK0(slice) = Mean_porosity_term_forK0;
    
    %% [Mike Comment: Translated to here]
                   
    %Print images
    f = sprintf('%04d',slice); %assigns index for naming images in sequence       
    imwrite(Voronoi2,[name5, f,'.tif'],'tif');
    imwrite(Voronoi_fibers,[name6, f,'.tif'],'tif');
    
    %export regions_area_all, Frac_unocc_all, and
    %area_weighted_normalized_frac_unocc_all for each slice
%    g = sprintf('%04d',slice); %assigns index for naming images in sequence
%    Table_voronoidist_lastslice = table(regions_area_all,Frac_unocc_all,area_weighted_normalized_frac_unocc_all);
%    writetable(Table_voronoidist_lastslice,['G:\March_2016_2bunch\recon_20160325_180509_235p2_wet_0p8cm_cont_4097im_1500ex_17keV_17_\Not CMC Finding\4p3_voronoi_region_analysis\4p3_voronoi_analysis\VoronoiCellAnalysis_',g,'.xlsx']);
    
    slice    % tells the slice number in the command window so you know how far along the run is

end

%Print table of values per slice
Table_slice = table(slice_Mean_regions_area,slice_Std_regions_area,slice_Mean_Frac_unocc,slice_Std_Frac_unocc,slice_GEV_k_Frac_unocc,slice_Mean_Area_unocc,slice_Std_Area_unocc,slice_GEV_k_Area_unocc,slice_alpha,slice_alpha_unocc,slice_Mean_porosity,slice_Mean_porosity_term_forF,slice_Mean_porosity_term_forK0,slice_Mean_area_weighted_normalized_frac_unocc_all,slice_Std_area_weighted_normalized_frac_unocc_all,slice_GEV_shape_factor_k);
writetable(Table_slice,name7);

%Calculate stack averages of avg and stdev for regions_area and Frac_unocc

stack_Mean_regions_area = mean(slice_Mean_regions_area);
stack_Std_regions_area = mean(slice_Std_regions_area);
stack_FracStd_regions_area = stack_Std_regions_area/stack_Mean_regions_area; %fractional standard deviation

stack_Mean_Frac_unocc = mean(slice_Mean_Frac_unocc);
stack_Std_Frac_unocc = mean(slice_Std_Frac_unocc);
stack_FracStd_Frac_unocc = stack_Std_Frac_unocc/stack_Mean_Frac_unocc; %fractional standard deviation
stack_GEV_k_Frac_unocc = mean(slice_GEV_k_Frac_unocc); %Mean shape factor value for GEV fit to the distribution of unoccupied fractions

stack_Mean_Area_unocc = mean(slice_Mean_Area_unocc);
stack_Std_Area_unocc = mean(slice_Std_Area_unocc);
stack_FracStd_Area_unocc = stack_Std_Area_unocc/stack_Mean_Area_unocc; %fractional standard deviation
stack_GEV_k_Area_unocc = mean(slice_GEV_k_Area_unocc); %Mean shape factor value for GEV fit to the distribution of unoccupied areas

stack_alpha = mean(slice_alpha);
stack_alpha_unocc = mean(slice_alpha_unocc);

stack_Mean_porosity = mean(slice_Mean_porosity); %Use this value of porosity to calculate permeability
stack_Std_Mean_porosity = std(slice_Mean_porosity); %Use this as the standard deviation of the above value 

stack_Mean_porosity_term_forF = mean(slice_Mean_porosity_term_forF);%Use this value of porosity term to calculate F
stack_Std_Mean_porosity_term_forF = std(slice_Mean_porosity_term_forF);%Use this as the standard deviation of the above value

stack_Mean_porosity_term_forK0 = mean(slice_Mean_porosity_term_forK0);%Use this value of porosity term to calculate K0
stack_Std_Mean_porosity_term_forK0 = std(slice_Mean_porosity_term_forK0);%Use this as the standard deviation of the above value

stack_Mean_area_weighted_normalized_frac_unocc_all = mean(slice_Mean_area_weighted_normalized_frac_unocc_all); %Mean of area weight fraction unoccupied
stack_Std_area_weighted_normalized_frac_unocc_all = mean(slice_Std_area_weighted_normalized_frac_unocc_all); %This is the stdev of frac unoccupied, BUT WEIGHTED BY AREA. Another possible measure of variation to plot agains K0 and F (we used 
stack_FracStd_area_weighted_normalized_frac_unocc_all = stack_Std_area_weighted_normalized_frac_unocc_all/stack_Mean_area_weighted_normalized_frac_unocc_all; %Fractional standard deviation

stack_GEV_shape_factor_k = mean(slice_GEV_shape_factor_k); %Mean shape factor value

Table_stackCalcs2 = table(stack_Mean_regions_area,stack_Std_regions_area,stack_FracStd_regions_area,stack_Mean_Frac_unocc,stack_Std_Frac_unocc,stack_FracStd_Frac_unocc,stack_GEV_k_Frac_unocc,stack_Mean_Area_unocc,stack_Std_Area_unocc,stack_FracStd_Area_unocc,stack_GEV_k_Area_unocc,stack_alpha,stack_alpha_unocc,stack_Mean_porosity,stack_Std_Mean_porosity,stack_Mean_porosity_term_forF,stack_Std_Mean_porosity_term_forF,stack_Mean_porosity_term_forK0,stack_Std_Mean_porosity_term_forK0,stack_Mean_area_weighted_normalized_frac_unocc_all,stack_Std_area_weighted_normalized_frac_unocc_all,stack_FracStd_area_weighted_normalized_frac_unocc_all,stack_GEV_shape_factor_k);
writetable(Table_stackCalcs2,name8);








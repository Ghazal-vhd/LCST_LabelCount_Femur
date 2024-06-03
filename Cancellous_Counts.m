%Created by C Heveran and G Vahidi- Montana State University

%load images
%last modified Aug 17, 2022

clear all; clc; close all;

display('select the folder where the images are located')
x=uigetdir;
addpath(x)
cd(x)


file1 = uigetfile({'*.Tiff'},'Select the file for the reflectance image with PO and EC borders indicated'); 
%file2 = uigetfile({'*.Tiff'},'Select the file for reflectance image with lamellar/woven regions indicated');
file3 = uigetfile({'*.Tiff'},'Select the file for the GREEN lacunae');
file4 = uigetfile({'*.Tiff'},'Select the file for the RED lacunae');


IRF  = imread(file1); %read the reflectance image with PO and EC borders indicated
%IRF2 = imread(file2); % read the reflectance image with the lamellar/woven regions indicafted
IG = imread(file3); % read the green image
IR = imread(file4); % read the red image

%% Read the cancellous boundaries from the reflectance image

% note: outline the trabecular regions in blue

BlueImage   = IRF(:,:,3);   %cancellous 'islands'
RedImage    = IRF(:,:,1);
GreenImage  = IRF(:,:,2); 

BlueOnly  = BlueImage - RedImage - GreenImage;

ROI = imfill(BlueOnly,'holes');

ROIprops = regionprops(ROI,'all');

for i = 1:length(ROIprops)
    CancArea(i) = ROIprops(i).Area;
end

TotalArea = sum(CancArea); %pixels^2

%% Plot the regions and boundaries

figure(1); clf(1)
subplot(1,2,1)
imshow(BlueImage)
title('Blue channel of image')

subplot(1,2,2)
imshow(ROI)
title('Cancellous regions to be analyzed')

%% GREEN (2d) FLUORESCENCE IMAGES: Blue = labeled, Red = unlabeled

[rows columns channels] = size(IG);

FindBlue  = IG(:,:,3);  %find the blue pixels - these are labeled lacunae
FindRed   = IG(:,:,1);  %find the red pixels - these are unlabeled lacunae
FindGreen = IG(:,:,2);  %find the green pixels

% Create the mask for finding blue marked lacunae

BlueMask = FindBlue - FindRed - FindGreen; 
BlueMask_filled = imfill(BlueMask,'holes');
BlueMask_filled = logical(BlueMask_filled);

% Create the mask for finding red marked lacunae

RedMask = FindRed - FindBlue - FindGreen;
RedMask_filled = imfill(RedMask,'holes');
RedMask_filled = logical(RedMask_filled);

% Remove extraneous pixels outside of the periosteal or endocortical
% boundaries

%BlueMask_filled_clean = logical(BlueMask_filled.*ROI);
%RedMask_filled_clean  = logical(RedMask_filled.*ROI);

BlueMask_filled_clean = BlueMask_filled;
RedMask_filled_clean = RedMask_filled;

% Analyze the regions for the Blue and Red Masks

s1 = regionprops(BlueMask_filled_clean,'all');
s2 = regionprops(RedMask_filled_clean,'all');


for i = 1:length(s1)
    bluearea(i) = s1(i).Area;
    if s1(i).Area > 10
    bluecentroid = s1(i).Centroid;
    bluecentx(i) = bluecentroid(1);
    bluecenty(i) = bluecentroid(2);
    end
end

for i = 1:length(s2)
    redarea(i) = s2(i).Area;
    if s2(i).Area > 10
    redcentroid = s2(i).Centroid;
    redcentx(i) = redcentroid(1);
    redcenty(i) = redcentroid(2);
    end
end

% Plot a marker on the original green image for each of the filled lacunae

figure(2); clf(2)
imshow(FindBlue);
hold on
plot(bluecentx, bluecenty,'b*',redcentx,redcenty,'r*')
title('Lacunae colored blue (calcein label) or red (no label)')


for i = 1:length(s1)
    Blue_area(i)     = s1(i).Area;
    Blue_majaxis(i)  = s1(i).MajorAxisLength;
    Blue_minaxis(i)  = s1(i).MinorAxisLength;
    Blue_circ(i)     = Blue_minaxis(i) / Blue_majaxis(i);
    Blue_centroid    = s1(i).Centroid;
    Blue_centx(i)    = Blue_centroid(1);
    Blue_centy(i)    = Blue_centroid(2);
end

if length(s2) > 0 

for i = 1:length(s2)
    Red_area(i)     = s2(i).Area;
    Red_majaxis(i)  = s2(i).MajorAxisLength;
    Red_minaxis(i)  = s2(i).MinorAxisLength;
    Red_circ(i)     = Red_minaxis(i) / Red_majaxis(i);
end
end

GREEN_Labeled_Num   = length(Blue_area);
GREEN_Unlabeled_Num = length(Red_area);

% calculate number densities in lamellar and woven regions
GREEN_Labeled_Num_dens   = GREEN_Labeled_Num / TotalArea;
GREEN_Unlabeled_Num_dens = GREEN_Unlabeled_Num / TotalArea;


%% RED FLUORESCENCE IMAGES: Blue = labeled, Green = unlabeled

[rows columns channels] = size(IR);

FindBlue2  = IR(:,:,3);  %find the blue pixels - these are labeled lacunae
FindRed2   = IR(:,:,1);  %find the red pixels 
FindGreen2 = IR(:,:,2);  %find the green pixels - these are unlabeled lacunae

% Create the mask for finding blue marked lacunae

BlueMask2 = FindBlue2 - FindRed2 - FindGreen2; 
BlueMask_filled2 = imfill(BlueMask2,'holes');
BlueMask_filled2 = logical(BlueMask_filled2);

% Create the mask for finding green marked lacunae

GreenMask2 = FindGreen2 - FindBlue2 - FindRed2;
GreenMask_filled2 = imfill(GreenMask2,'holes');
GreenMask_filled2 = logical(GreenMask_filled2);

% Remove extraneous pixels outside of the periosteal or endocortical
% boundaries

%BlueMask_filled_clean2 = logical(BlueMask_filled2.*ROI);
%GreenMask_filled_clean2  = logical(GreenMask_filled2.*ROI);

BlueMask_filled_clean2 = BlueMask_filled2;
GreenMask_filled_clean2  = GreenMask_filled2;


% Analyze the regions for the Blue and Red Masks

s3 = regionprops(BlueMask_filled_clean2,'all');
s4 = regionprops(GreenMask_filled_clean2,'all');

 
for i = 1:length(s3)
    bluearea2(i) = s3(i).Area;
    if s3(i).Area > 10
    bluecentroid2 = s3(i).Centroid;
    bluecentx2(i) = bluecentroid2(1);
    bluecenty2(i) = bluecentroid2(2);
    end
end

for i = 1:length(s4)
    greenarea2(i) = s4(i).Area;
    if s4(i).Area > 10
    greencentroid2 = s4(i).Centroid;
    greencentx2(i) = greencentroid2(1);
    greencenty2(i) = greencentroid2(2);
    end
end
   

% Plot a marker on the original RED image for each of the filled lacunae

figure(3); clf(3)
imshow(FindRed2);
hold on
plot(bluecentx2, bluecenty2,'b*',greencentx2,greencenty2,'g*')
title('Lacunae colored blue (calcein label) or green (no label)')

if length(s3) > 0
for i = 1:length(s3)
    Blue_area2(i)     = s3(i).Area;
    Blue_majaxis2(i)  = s3(i).MajorAxisLength;
    Blue_minaxis2(i)  = s3(i).MinorAxisLength;
    Blue_circ2(i)     = Blue_minaxis2(i) / Blue_majaxis2(i);
    Blue_centroid2    = s3(i).Centroid;
    Blue_centx2(i)    = Blue_centroid2(1);
    Blue_centy2(i)    = Blue_centroid2(2);
end
end

if length(s4) > 0
for i = 1:length(s4)
    Green_area2(i)     = s4(i).Area;
    Green_majaxis2(i)  = s4(i).MajorAxisLength;
    Green_minaxis2(i)  = s4(i).MinorAxisLength;
    Green_circ2(i)     = Green_minaxis2(i) / Green_majaxis2(i);
end
end

RED_Labeled_Num       = length(Blue_area2);
RED_Unlabeled_Num     = length(Green_area2);

% calculate number densities 
RED_Labeled_Num_dens  = RED_Labeled_Num / TotalArea;
RED_Labeled_Num_dens  = RED_Labeled_Num / TotalArea;

%% Some additional comparisons about GREEN (2d) and RED (4d, 6d, or 8d) lacunae

%Total lacunae
TotalLacunae = GREEN_Labeled_Num + GREEN_Unlabeled_Num;

%percent of GREEN IMAGE lacunae that are labeled
TotalGreenLacunae = GREEN_Labeled_Num + GREEN_Unlabeled_Num;
GREEN_LabeledPct = GREEN_Labeled_Num / TotalGreenLacunae * 100;

%percent of RED IMAGE lacunae that are labeled
TotalRedLacunae = RED_Labeled_Num + RED_Unlabeled_Num;
RED_LabeledPct = RED_Labeled_Num / TotalRedLacunae * 100;

%% double labeled

% find which lacunae are double labeled 

LabelOverlap = zeros(rows,columns);
for i = 1:rows
    for j = 1:columns
        if BlueMask_filled(i,j) == 1 && BlueMask_filled(i,j) == BlueMask_filled2(i,j)
           LabelOverlap(i,j) = 1;
        end
    end
end

LabelOverlap = logical(LabelOverlap);
s5 = regionprops(LabelOverlap,'all');
Overlap_Num = length(s5);

LabelOverlapPct   = Overlap_Num / TotalLacunae * 100;


%% Reporting
names = {'GREEN total % labeled','RED total % labeled','Total Doub Lab Lac'};
T = table(GREEN_LabeledPct',RED_LabeledPct', LabelOverlapPct','VariableNames',names);




%Created by C Heveran and G Vahidi- Montana State University


%load images
%last modified Aug 17, 2022

clear all; clc; close all;

display('select the folder where the images are located')
x=uigetdir;
addpath(x)
cd(x)


file1 = uigetfile({'*.Tiff'},'Select the file for the reflectance image with PO and EC borders indicated'); 
file2 = uigetfile({'*.Tiff'},'Select the file for reflectance image with lamellar/woven regions indicated');
file3 = uigetfile({'*.Tiff'},'Select the file for the GREEN lacunae');
file4 = uigetfile({'*.Tiff'},'Select the file for the RED lacunae');


IRF  = imread(file1); %read the reflectance image with PO and EC borders indicated
IRF2 = imread(file2); % read the reflectance image with the lamellar/woven regions indicafted
IG = imread(file3); % read the green image
IR = imread(file4); % read the red image

%% Read the Periosteal and Endocortical boundaries from the reflectance image

RedBoundary   = IRF(:,:,1);   %endocortical
GreenBoundary = IRF(:,:,2);   %other borders
BlueBoundary  = IRF(:,:,3);   %periosteal

PeriostealBoundary_gray = BlueBoundary - RedBoundary - GreenBoundary;
PeriostealBoundary = PeriostealBoundary_gray > 100;

EndocorticalBoundary_gray = RedBoundary - GreenBoundary-BlueBoundary;
EndocorticalBoundary = EndocorticalBoundary_gray > 100;

OtherBoundaries_gray = GreenBoundary-RedBoundary-BlueBoundary;
OtherBoundaries = OtherBoundaries_gray > 100;

TotalBoundary = PeriostealBoundary + EndocorticalBoundary + OtherBoundaries; 

ROI = imfill(TotalBoundary,'holes');

ROIprops = regionprops(ROI,'all');
TotalArea = ROIprops(1).Area; %pixels^2

%% Find the lamellar and woven regions from the other reflectance image
% Blue = lamellar, red = woven. There may be multiple regions of each type.

RedBoundaryLW   = IRF2(:,:,1);   %woven
BlueBoundaryLW  = IRF2(:,:,3);   %lamellar
GreenBoundaryLW = IRF2(:,:,2);

WovenBoundaries_gray = RedBoundaryLW - BlueBoundaryLW - GreenBoundaryLW;
WovenBoundaries = WovenBoundaries_gray > 100;
WovenFilled = imfill(WovenBoundaries,'holes');

LamellarBoundaries_gray = BlueBoundaryLW - GreenBoundaryLW - RedBoundaryLW;
LamellarBoundaries = LamellarBoundaries_gray > 100;
LamellarFilled = imfill(LamellarBoundaries,'holes');

LamWovenBoundaries = LamellarBoundaries + WovenBoundaries;

WovenROIs = regionprops(WovenFilled,'all');

if length(WovenROIs) > 1% If there is a WovenROI that is larger than an accidental speck
for i = 1:length(WovenROIs)
    WovenAreaVec(i) = WovenROIs(i).Area;
end

WovenArea = sum(WovenAreaVec);

else
    WovenArea = 0;

end

LamellarROIs = regionprops(LamellarFilled,'all');

if length(LamellarROIs) > 0
for i = 1:length(LamellarROIs)
    LamellarAreaVec(i) = LamellarROIs(i).Area;
end

LamellarArea = sum(LamellarAreaVec);

else
    LamellarArea = [];
end

%% Plot the regions and boundaries

figure(1); clf(1)
subplot(1,3,1)
imshowpair(ROI,LamWovenBoundaries)
title('Total ROI')

subplot(1,3,2)
imshow(WovenFilled)
title('Woven regions')

subplot(1,3,3)
imshow(LamellarFilled)
title('Lamellar regions')

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

%Numbers in Woven Bone
BlueMask_filled_woven = logical(BlueMask_filled.*WovenFilled);
RedMask_filled_woven  = logical(RedMask_filled.*WovenFilled);

%Number in Lamellar Bone
BlueMask_filled_lamellar = logical(BlueMask_filled.*LamellarFilled);
RedMask_filled_lamellar  = logical(RedMask_filled.*LamellarFilled);

s3 = regionprops(BlueMask_filled_woven,'all');
s4 = regionprops(RedMask_filled_woven,'all');

s5 = regionprops(BlueMask_filled_lamellar,'all');
s6 = regionprops(RedMask_filled_lamellar,'all');

if length(s3) > 0   % if there is at least one woven region and it has lacunae in it
for i = 1:length(s3)
    Blue_woven_area(i)     = s3(i).Area;
    Blue_woven_majaxis(i)  = s3(i).MajorAxisLength;
    Blue_woven_minaxis(i)  = s3(i).MinorAxisLength;
    Blue_woven_circ(i)     = Blue_woven_minaxis(i) / Blue_woven_majaxis(i);
    Blue_woven_centroid    = s3(i).Centroid;
    Blue_woven_centx(i)    = Blue_woven_centroid(1);
    Blue_woven_centy(i)    = Blue_woven_centroid(2);
end

else
    Blue_woven_area     = [];
    Blue_woven_majaxis  = [];
    Blue_woven_minaxis  = [];
    Blue_woven_circ     = [];
    Blue_woven_centx    = [];
    Blue_woven_centy    = [];
    
    
end

if length(s4) > 0 

for i = 1:length(s4)
    Red_woven_area(i)     = s4(i).Area;
    Red_woven_majaxis(i)  = s4(i).MajorAxisLength;
    Red_woven_minaxis(i)  = s4(i).MinorAxisLength;
    Red_woven_circ(i)     = Red_woven_minaxis(i) / Red_woven_majaxis(i);
end

else  % if there is not a woven region
       
    Red_woven_area     = [];
    Red_woven_majaxis  = [];
    Red_woven_minaxis  = [];
    Red_woven_circ     = [];
    
end

if length(s5) > 0

for i = 1:length(s5)
    Blue_lamellar_area(i)     = s5(i).Area;
    Blue_lamellar_majaxis(i)  = s5(i).MajorAxisLength;
    Blue_lamellar_minaxis(i)  = s5(i).MinorAxisLength;
    Blue_lamellar_circ(i)     = Blue_lamellar_minaxis(i) / Blue_lamellar_majaxis(i);
    Blue_lamellar_centroid    = s5(i).Centroid;
    Blue_lamellar_centx(i)    = Blue_lamellar_centroid(1);
    Blue_lamellar_centy(i)    = Blue_lamellar_centroid(2);
end

for i = 1:length(s6)
    Red_lamellar_area(i)     = s6(i).Area;
    Red_lamellar_majaxis(i)  = s6(i).MajorAxisLength;
    Red_lamellar_minaxis(i)  = s6(i).MinorAxisLength;
    Red_lamellar_circ(i)     = Red_lamellar_minaxis(i) / Red_lamellar_majaxis(i);
end

else 
     Blue_lamellar_area     = [];
    Blue_lamellar_majaxis  = [];
    Blue_lamellar_minaxis  = [];
    Blue_lamellar_circ    = [];
    Blue_lamellar_centroid    = [];
    Blue_lamellar_centx   = [];
    Blue_lamellar_centy    = [];
    Red_lamellar_area     = [];
    Red_lamellar_majaxis  = [];
    Red_lamellar_minaxis  = [];
    Red_lamellar_circ     = [];
    
end

GREEN_Labeled_lamellar_Num   = length(Blue_lamellar_area);
GREEN_Labeled_woven_Num      = length(Blue_woven_area);
GREEN_Unlabeled_lamellar_Num = length(Red_lamellar_area);
GREEN_Unlabeled_woven_Num    = length(Red_woven_area);

% calculate number densities in lamellar and woven regions

if LamellarArea > 0
GREEN_Labeled_lamellar_Num_dens   = GREEN_Labeled_lamellar_Num / LamellarArea;
GREEN_Unlabeled_lamellar_Num_dens = GREEN_Unlabeled_lamellar_Num / LamellarArea;
GREEN_TotalLabeled_Num_dens       = (GREEN_Labeled_lamellar_Num + GREEN_Labeled_woven_Num) / (LamellarArea + WovenArea);
GREEN_TotalUnlabeled_Num_dens     = (GREEN_Unlabeled_lamellar_Num + GREEN_Unlabeled_woven_Num) / (LamellarArea + WovenArea);
end


GREEN_Labeled_woven_Num_dens      = GREEN_Labeled_woven_Num / WovenArea;
GREEN_Unlabeled_woven_Num_dens    = GREEN_Unlabeled_woven_Num / WovenArea;



figure(3)
subplot(1,3,1)
imshow(GreenBoundaryLW)
title('Boundaries')

subplot(1,3,2)
imshowpair(BlueMask_filled_lamellar,RedMask_filled_lamellar)
title('Lacunae in Lamellar Bone')

subplot(1,3,3)
imshowpair(BlueMask_filled_woven,RedMask_filled_woven)
title('Lacunae in Woven Bone')

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

s7 = regionprops(BlueMask_filled_clean2,'all');
s8 = regionprops(GreenMask_filled_clean2,'all');

 
for i = 1:length(s7)
    bluearea2(i) = s7(i).Area;
    if s7(i).Area > 10
    bluecentroid2 = s7(i).Centroid;
    bluecentx2(i) = bluecentroid2(1);
    bluecenty2(i) = bluecentroid2(2);
    end
end

for i = 1:length(s8)
    greenarea2(i) = s8(i).Area;
    if s8(i).Area > 10
    greencentroid2 = s8(i).Centroid;
    greencentx2(i) = greencentroid2(1);
    greencenty2(i) = greencentroid2(2);
    end
end
   

% Plot a marker on the original RED image for each of the filled lacunae

figure(4); clf(4)
imshow(FindRed2);
hold on
plot(bluecentx2, bluecenty2,'b*',greencentx2,greencenty2,'g*')
title('Lacunae colored blue (calcein label) or green (no label)')

%Numbers in Woven Bone
BlueMask_filled_woven2 = logical(BlueMask_filled2.*WovenFilled);
GreenMask_filled_woven2  = logical(GreenMask_filled2.*WovenFilled);

%Number in Lamellar Bone
BlueMask_filled_lamellar2 = logical(BlueMask_filled2.*LamellarFilled);
GreenMask_filled_lamellar2  = logical(GreenMask_filled2.*LamellarFilled);

s9 = regionprops(BlueMask_filled_woven2,'all');
s10 = regionprops(GreenMask_filled_woven2,'all');

s11 = regionprops(BlueMask_filled_lamellar2,'all');
s12 = regionprops(GreenMask_filled_lamellar2,'all');

if length(s9) > 0
for i = 1:length(s9)
    Blue_woven_area2(i)     = s9(i).Area;
    Blue_woven_majaxis2(i)  = s9(i).MajorAxisLength;
    Blue_woven_minaxis2(i)  = s9(i).MinorAxisLength;
    Blue_woven_circ2(i)     = Blue_woven_minaxis2(i) / Blue_woven_majaxis2(i);
    Blue_woven_centroid2    = s9(i).Centroid;
    Blue_woven_centx2(i)    = Blue_woven_centroid2(1);
    Blue_woven_centy2(i)    = Blue_woven_centroid2(2);
end

else
    Blue_woven_area2     = [];
    Blue_woven_majaxis2  = [];
    Blue_woven_minaxis2  = [];
    Blue_woven_circ2     = [];
    Blue_woven_centx2    = [];
    Blue_woven_centy2    = [];
end

if length(s10) > 0
for i = 1:length(s10)
    Green_woven_area2(i)     = s10(i).Area;
    Green_woven_majaxis2(i)  = s10(i).MajorAxisLength;
    Green_woven_minaxis2(i)  = s10(i).MinorAxisLength;
    Green_woven_circ2(i)     = Green_woven_minaxis2(i) / Green_woven_majaxis2(i);
end

else       
    
    Green_woven_area2     = [];
    Green_woven_majaxis2  = [];
    Green_woven_minaxis2  = [];
    Green_woven_circ2     = [];
   
end

if length(s11) > 1

for i = 1:length(s11)
    Blue_lamellar_area2(i)     = s11(i).Area;
    Blue_lamellar_majaxis2(i)  = s11(i).MajorAxisLength;
    Blue_lamellar_minaxis2(i)  = s11(i).MinorAxisLength;
    Blue_lamellar_circ2(i)     = Blue_lamellar_minaxis2(i) / Blue_lamellar_majaxis2(i);
    Blue_lamellar_centroid2    = s11(i).Centroid;
    Blue_lamellar_centx2(i)    = Blue_lamellar_centroid2(1);
    Blue_lamellar_centy2(i)    = Blue_lamellar_centroid2(2);
end

for i = 1:length(s12)
    Green_lamellar_area2(i)     = s12(i).Area;
    Green_lamellar_majaxis2(i)  = s12(i).MajorAxisLength;
    Green_lamellar_minaxis2(i)  = s12(i).MinorAxisLength;
    Green_lamellar_circ2(i)     = Green_lamellar_minaxis2(i) / Green_lamellar_majaxis2(i);
end

else
    Blue_lamellar_area2     = [];
    Blue_lamellar_majaxis2  = [];
    Blue_lamellar_minaxis2  = [];
    Blue_lamellar_circ2     = [];
    Blue_lamellar_centroid2    = [];
    Blue_lamellar_centx2    = [];
    Blue_lamellar_centy2    = [];
    Green_lamellar_area2     = [];
    Green_lamellar_majaxis2  = [];
    Green_lamellar_minaxis2  = [];
    Green_lamellar_circ2     = [];
    

end

RED_Labeled_lamellar_Num    = length(Blue_lamellar_area2);
RED_Labeled_woven_Num       = length(Blue_woven_area2);
RED_Unlabeled_lamellar_Num  = length(Green_lamellar_area2);
RED_Unlabeled_woven_Num     = length(Green_woven_area2);

% calculate number densities in lamellar and woven regions
if LamellarArea > 0
RED_Labeled_lamellar_Num_dens   = RED_Labeled_lamellar_Num / LamellarArea;
RED_Unlabeled_lamellar_Num_dens = RED_Unlabeled_lamellar_Num / LamellarArea;
RED_TotalLabeled_Num_dens       = (RED_Labeled_lamellar_Num + RED_Labeled_woven_Num) / (LamellarArea + WovenArea);
RED_TotalUnlabeled_Num_dens     = (RED_Unlabeled_lamellar_Num + RED_Unlabeled_woven_Num) / (LamellarArea + WovenArea);
end

RED_Labeled_woven_Num_dens      = RED_Labeled_woven_Num / WovenArea;
RED_Unlabeled_woven_Num_dens    = RED_Unlabeled_woven_Num / WovenArea;



figure(5)
subplot(1,3,1)
imshow(GreenBoundaryLW)
title('Boundaries')

subplot(1,3,2)
imshowpair(BlueMask_filled_lamellar2,GreenMask_filled_lamellar2)
title('Lacunae in Lamellar Bone')

subplot(1,3,3)
imshowpair(BlueMask_filled_woven2,GreenMask_filled_woven2)
title('Lacunae in Woven Bone')


%% Some additional comparisons about GREEN (2d) and RED (4d, 6d, or 8d) lacunae

%Total lacunae
TotalLamellarLacunae = GREEN_Labeled_lamellar_Num + GREEN_Unlabeled_lamellar_Num;
TotalWovenLacunae    = GREEN_Labeled_woven_Num + GREEN_Unlabeled_woven_Num;

%percent of GREEN IMAGE lacunae that are labeled
TotalGreenLacunae = GREEN_Labeled_lamellar_Num + GREEN_Labeled_woven_Num + GREEN_Unlabeled_lamellar_Num + GREEN_Unlabeled_woven_Num; 
GREEN_LabeledPct = (GREEN_Labeled_lamellar_Num + GREEN_Labeled_woven_Num) / TotalGreenLacunae * 100;

GREEN_LabeledPct_Lam = GREEN_Labeled_lamellar_Num / (GREEN_Labeled_lamellar_Num + GREEN_Unlabeled_lamellar_Num)*100;
GREEN_LabeledPct_Wov = GREEN_Labeled_woven_Num / (GREEN_Labeled_woven_Num + GREEN_Unlabeled_woven_Num)*100;

%percent of RED IMAGE lacunae that are labeled
TotalRedLacunae = RED_Labeled_lamellar_Num + RED_Labeled_woven_Num + RED_Unlabeled_lamellar_Num + RED_Unlabeled_woven_Num; 
RED_LabeledPct = (RED_Labeled_lamellar_Num + RED_Labeled_woven_Num) / TotalRedLacunae * 100;

RED_LabeledPct_Lam = RED_Labeled_lamellar_Num / (RED_Labeled_lamellar_Num + RED_Unlabeled_lamellar_Num)*100;
RED_LabeledPct_Wov = RED_Labeled_woven_Num / (RED_Labeled_woven_Num + RED_Unlabeled_woven_Num)*100;

%% double labeled

% find which lacunae are double labeled in lamellar bone

LamellarOverlap = zeros(rows,columns);
for i = 1:rows
    for j = 1:columns
        if BlueMask_filled_lamellar(i,j) == 1 && BlueMask_filled_lamellar(i,j) == BlueMask_filled_lamellar2(i,j)
           LamellarOverlap(i,j) = 1;
        end
    end
end

LamellarOverlap = logical(LamellarOverlap);
s10 = regionprops(LamellarOverlap,'all');
Overlap_lamellar_Num = length(s10);

LamellarOverlapPct   = Overlap_lamellar_Num / TotalLamellarLacunae * 100;


if isempty(WovenROIs) == 0  % if there is a woven region 

WovenOverlap = zeros(rows,columns);
for i = 1:rows
    for j = 1:columns
        if BlueMask_filled_woven(i,j) == 1 && BlueMask_filled_woven(i,j) == BlueMask_filled_woven2(i,j)
           WovenOverlap(i,j) = 1;
        end
    end
end

WovenOverlap = logical(WovenOverlap);
s11 = regionprops(WovenOverlap,'all');

Overlap_woven_Num = length(s11);

WovenOverlapPct   = Overlap_woven_Num / TotalWovenLacunae * 100;

else   
 WovenOverlapPct = NaN;
end



%% calculate ALL double labels
TotalLacunae = TotalLamellarLacunae + TotalWovenLacunae;

if WovenOverlapPct > 1  %asking if its not NaN

TotalOverlapPct = (Overlap_woven_Num + Overlap_lamellar_Num) / TotalLacunae * 100;

else
    TotalOverlapPct = NaN;
end


%% distance of lacunae from surfaces

%find minimum distance from each lacuna centroid to the endocortical surface
 
% s13 = regionprops(EndocorticalBoundary,'all');  %lets us find the pixels on this boundary
% s14 = regionprops(PeriostealBoundary,'all');
% 
% % find the biggest object, which should correspond to the boundary
% % (usually, but not always, object #1)
% 
% for i = 1:length(s13)
%     ECcelllength(i) = length(s13(i).PixelList);
%     EC_pixelcell{i} = s13(i).PixelList;
% end
% 
% [a b] = max(ECcelllength);
% EC_pixels = EC_pixelcell{b};
% 
% for i = 1:length(s14)
%  POcelllength(i) = length(s14(i).PixelList);
%  PO_pixelcell{i} = s14(i).PixelList;
% 
% end
% [a b] = max(POcelllength);
% PO_pixels = PO_pixelcell{b};
% 
% ECx = EC_pixels(:,1);
% ECy = EC_pixels(:,2);
% 
% POx = PO_pixels(:,1);
% POy = PO_pixels(:,2);
% 
% ECpixels = length(ECx);
% POpixels = length(POx);
% 
% figure(6);clf(6)
% imshow(BlueBoundary)
% hold on
% plot(ECx,ECy,'b*',POx,POy,'r*')
% legend('pixelized EC boundary','pixelized PO boundary')
% 
% % GREEN LACUNAE
% 
% % find the distances from the endocortical and periosteal surface of all
% % the labeled lacunae
% 
% % GREEN, LAMELLAR
%  for i = 1:length(Blue_lamellar_centx)
%      for j = 1:ECpixels     % endocortical surface, green labeled, lamellar
%          EC_green_lam_dist{i}(j) = sqrt((Blue_lamellar_centx(i) - ECx(j))^2 + (Blue_lamellar_centy(i) - ECy(j))^2);
%      end
%      
%      for j = 1:POpixels    % periosteal surface, green labeled, lamellar
%          PO_green_lam_dist{i}(j) = sqrt((Blue_lamellar_centx(i) - POx(j))^2 + (Blue_lamellar_centy(i) - POy(j))^2);
%      end
%      
%  EC_green_lam_dist_min(i) = min(EC_green_lam_dist{i});
%  PO_green_lam_dist_min(i) = min(PO_green_lam_dist{i});
%  EC_green_lam_rel_distance_labeled(i) = EC_green_lam_dist_min(i)/(EC_green_lam_dist_min(i) +  PO_green_lam_dist_min(i));
%  
%  end
%  
%  % GREEN, WOVEN
%  if length(Blue_woven_centx) > 0
%      
%  for i = 1:length(Blue_woven_centx)
%      for j = 1:ECpixels     % endocortical surface, green labeled, woven
%          EC_green_wov_dist{i}(j) = sqrt((Blue_woven_centx(i) - ECx(j))^2 + (Blue_woven_centy(i) - ECy(j))^2);
%      end
%      
%      for j = 1:POpixels    % periosteal surface, green labeled, woven
%          PO_green_wov_dist{i}(j) = sqrt((Blue_woven_centx(i) - POx(j))^2 + (Blue_woven_centy(i) - POy(j))^2);
%      end
%      
%  EC_green_wov_dist_min(i) = min(EC_green_wov_dist{i});
%  PO_green_wov_dist_min(i) = min(PO_green_wov_dist{i});
%  EC_green_wov_rel_distance_labeled(i) = EC_green_wov_dist_min(i)/(EC_green_wov_dist_min(i) +  PO_green_wov_dist_min(i));
%  
%  end  
%      
%  end
%  
%  
%  % RED LACUNAE
%  
%  % RED LAMELLAR
%  
%   for i = 1:length(Blue_lamellar_centx2)
%      for j = 1:ECpixels     % endocortical surface, red labeled, lamellar
%          EC_red_lam_dist{i}(j) = sqrt((Blue_lamellar_centx2(i) - ECx(j))^2 + (Blue_lamellar_centy2(i) - ECy(j))^2);
%      end
%      
%      for j = 1:POpixels    % periosteal surface, red labeled, lamellar
%          PO_red_lam_dist{i}(j) = sqrt((Blue_lamellar_centx2(i) - POx(j))^2 + (Blue_lamellar_centy2(i) - POy(j))^2);
%      end
%      
%  EC_red_lam_dist_min(i) = min(EC_red_lam_dist{i});
%  PO_red_lam_dist_min(i) = min(PO_red_lam_dist{i});
%  EC_red_lam_rel_distance_labeled(i) = EC_red_lam_dist_min(i)/(EC_red_lam_dist_min(i) +  PO_red_lam_dist_min(i));
%  
%   end
%  
%   
%  % RED WOVEN
%  
%  if length(Blue_woven_centx2) > 0
%      
%  for i = 1:length(Blue_woven_centx2)
%      for j = 1:ECpixels     % endocortical surface, red labeled, woven
%          EC_red_wov_dist{i}(j) = sqrt((Blue_woven_centx2(i) - ECx(j))^2 + (Blue_woven_centy2(i) - ECy(j))^2);
%      end
%      
%      for j = 1:POpixels    % periosteal surface, red labeled, woven
%          PO_red_wov_dist{i}(j) = sqrt((Blue_woven_centx2(i) - POx(j))^2 + (Blue_woven_centy2(i) - POy(j))^2);
%      end
%      
%  EC_red_wov_dist_min(i) = min(EC_red_wov_dist{i});
%  PO_red_wov_dist_min(i) = min(PO_red_wov_dist{i});
%  EC_red_wov_rel_distance_labeled(i) = EC_red_wov_dist_min(i)/(EC_red_wov_dist_min(i) +  PO_red_wov_dist_min(i));
%  
%  end  
%      
%  end
%  
% 
% 
% figure(7);clf(7)
% subplot(2,1,1)
% h = histogram(EC_green_lam_rel_distance_labeled);
% h.BinWidth = 0.1; %change to modify bind width
% h.FaceColor = 'green'; %change to modify color; either 'r', 'g','b'... or color triplets [a b c]
% h.FaceAlpha = 0.6; %change to modify transparency
% 
% hold on
% 
% if  length(Blue_woven_centx) > 0  % if there are green labelled lacunae in woven bone
%     
% h2 = histogram(EC_green_wov_rel_distance_labeled);
% h2.BinWidth = 0.1; %change to modify bind width
% h2.FaceColor = 'cyan'; %change to modify color; either 'r', 'g','b'... or color triplets [a b c]
% h2.FaceAlpha = 0.6; %change to modify transparency   
%     
% end
% 
% legend('green, lamellar', 'green, woven')
% title('Histograms of rel dist of GREEN lacunae 0 = EC surface, 1 = PO surface')
% 
% subplot(2,1,2)
% 
% h3 = histogram(EC_red_lam_rel_distance_labeled);
% h3.BinWidth = 0.1; %change to modify bind width
% h3.FaceColor = 'r'; %change to modify color; either 'r', 'g','b'... or color triplets [a b c]
% h3.FaceAlpha = 0.6; %change to modify transparency  
% 
% hold on
% 
% if  length(Blue_woven_centx2) > 0
% h4 = histogram(EC_red_wov_rel_distance_labeled);
% h4.BinWidth = 0.1; %change to modify bind width
% h4.FaceColor = 'y'; %change to modify color; either 'r', 'g','b'... or color triplets [a b c]
% h4.FaceAlpha = 0.7; %change to modify transparency 
% end
% 
% legend('red, lamellar', 'red, woven')
% title('Histograms of rel dist of RED lacunae 0 = EC surface, 1 = PO surface')
%  % calculate distances, where 'i' keeps count of individual lacunae
%  

%% Reporting
names = {'GREEN total % labeled','GREEN Lam % lab','GREEN Wov % lab','RED total % labeled','RED Lam % lab','RED wov % lab','Lamellar R & G % with both', 'Woven R & G % with both','Total Doub Lab Lac'};
T = table(GREEN_LabeledPct',GREEN_LabeledPct_Lam', GREEN_LabeledPct_Wov', RED_LabeledPct', RED_LabeledPct_Lam', RED_LabeledPct_Wov',LamellarOverlapPct',WovenOverlapPct',TotalOverlapPct','VariableNames',names);




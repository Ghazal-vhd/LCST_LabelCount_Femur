%Created by Connor Boone- Montana State University

% HOW TO
%
% To run the code, select the folder where the images are located 
    % ex: double click the 'L' folder in C-4-3
% the code will then identify the border, red,green and type files as long
% as they are named as:
    % FILE NAMING SCHEME 
    % red: -r, -R, -Red, -red
    % green: -g, -G, -Green, -green
    % border: -border, -Border
    % type: -type, -Type
% CODE OUTPUTS:
    % 1. a table called 'T' with the percent remodeling cells in each region 
        %  (0-10%,10-20%,20-30%,...) for both the green and red files
    % 2. Two data structures, called 'cell_data_green' and 'cell_data_red'
        % which house data for each of the marked cells specifiying their,
        % coordinates, percentage distance, distance to EC and PO edges, 
        % whether the cell is remodeling or not, (0 or 1), distance from 
        % the closest EC and closest PO points (used for calculating 
        % the percentage distance, and the percentage distance which is
        % calculated by dividing the distance to the EC border divided 
        % by distance between the closest EC and PO points
% OPTIONS:
    % hard_code_images: is a option to hard code the path of each of the 
    % files, by changing the variable "Hard code images" to true, and 
    % inputting the file paths in.
    %
    % create_excel_file: creates an excel file in the current folder of the
    % data in the T table
    %
    % Testing: set to true to display 5 figures showing the EC and PO 
    % borders and the marked cells (both positive and negative) for 
    % both red and green files. To ensure the code is properly 
    % identifying the EC and PO borders and the positive and 
    % negative remodeling cells
    %
    % Figures: set to true to show 2 figures, one for the green and one for
    % the red marked files, showing the percent remodeling cells in each 
    % region (0-10%,10-20%,20-30%,...) 
% WARNINGS:
    % A warning will be thrown when a percentage distance is calculated to
    % be greater than 100%, this might not be an error depending on the
    % size of the shape
    % A warning will be thrown if the files aren't selected properly
clc
clear variables
close all

testing = true;
hard_code_images = false;
show_figures = true;
create_excel_file = true;

%% file selection
if hard_code_images

    
    % INPUT FILE PATHS HERE
    BorderImage = imread("/Users/Connor/Desktop/College/Work/Distance/C-4-3/L/4-3-l-border.tiff");
    im_green = imread("/Users/Connor/Desktop/College/Work/Distance/C-4-3/L/4-3-l-g.tiff");
    im_red = imread("/Users/Connor/Desktop/College/Work/Distance/C-4-3/L/4-3-l-r.tiff");
    im_type = imread("/Users/Connor/Desktop/College/Work/Distance/C-4-3/L/4-3-l-Type.tiff");

else

    disp('select the folder with the quadrent where the images are located')
    x=string(uigetdir);
    file_path = strsplit(x,"/");

    main_folder = file_path(end-1);
    file_path = file_path(end);

    
    if ~file_path.contains(["A","L","M","P"])
        warning("Warning, folder not selected is not named A, L, M, P " + ...
            "autoselection might still work")
    end
    
    % checks files in the selected folder
    folder_contents = string(ls(x+"/"));
    disp("The folder selected has the path of "+x)

    % splits it into the individual files
    folder_contents = folder_contents.split;
    % removes the blank elements 
    folder_contents(folder_contents=="") = [];
    
    % the different naming schemes for each of the file types
    border_naming_scheme = ["border","Border"];
    % border_naming_scheme = ["adsbf"];
    green_naming_scheme = ["-G","-g","green","Green"];
    red_naming_scheme = ["-R","-r","red","Red"];
    type_naming_scheme = ["Type","type"];
    
    % locates the files in the folder
    border_file = folder_contents(contains(folder_contents,border_naming_scheme));
    green_file = folder_contents(contains(folder_contents,green_naming_scheme));
    red_file = folder_contents(contains(folder_contents,red_naming_scheme));
    type_file = folder_contents(contains(folder_contents,type_naming_scheme));
    
    
    % checks to make sure files were identified correctly, if not, a file will
    % be empty
    if isempty(border_file) || isempty(green_file) ...
        || isempty(red_file) || isempty(type_file)
        % figure out which file is empty
    
        if isempty(border_file)
            warning("Border file does not match naming scheme")
        end
    
        if isempty(green_file)
            warning("Green file does not match naming scheme")
        end
    
        if isempty(red_file)
            warning("Red file does not match naming scheme")   
        end
    
        if isempty(type_file)
            warning("Type file does not match naming scheme")   
        end
        % displays the folder selected and the files within the folder
        disp("These are the files found in selected folder")
        disp(x)
        disp(folder_contents)

        % code won't run with this error, so code quits
        disp("Quitting execution")
        return
    
    end
    
% loads the images in
BorderImage = imread(x+"/"+border_file);
im_green = imread(x+"/"+green_file);
im_red = imread(x+"/"+red_file);
im_type = imread(x+"/"+type_file);

% GV's old code to select the cells, this still works as long as you
% comment out everything above to the 'else' call
%
%   disp('select the folder where the images are located')
%     x=uigetdir;
%     addpath(x)
%     cd(x)
%     disp("Select the file for the reflectance image with PO and EC borders indicated")
%     file1 = uigetfile({'*.tiff'}); 
%     disp('Select the file for reflectance image with lamellar/woven regions indicated')
%     file2 = uigetfile({'*.tiff'});
%     disp('Select the file for the GREEN lacunae')
%     file3 = uigetfile({'*.tiff'});
%     disp('Select the file for the RED lacunae')
%     file4 = uigetfile({'*.tiff'});
%     
%     BorderImage = imread(file1);
%     im_type = imread(file2);
%     im_green = imread(file3);
%     im_red = imread(file4);

end



%%
% Sperate into RGB
RedImage_border = BorderImage(:,:,1);
BlueImage_border = BorderImage(:,:,3);

% Create images of only RGB
OnlyRed_border = RedImage_border - BlueImage_border;
OnlyBlue_border = BlueImage_border - RedImage_border;

% make into black and white

red_bw_border = imbinarize(OnlyRed_border);
red_mask_border = imfill(red_bw_border,"holes");

blue_bw_border = imbinarize(OnlyBlue_border);
blue_mask_border = imfill(blue_bw_border,"holes");

% type image
RedImage_type = im_type(:,:,1);
BlueImage_type = im_type(:,:,3);

% Create images of only RGB
OnlyRed_type = RedImage_type - BlueImage_type;
OnlyBlue_type = BlueImage_type - RedImage_type;

% make into black and white

red_bw_type = imbinarize(OnlyRed_type);
blue_bw_type = imbinarize(OnlyBlue_type);

% makes the outline into a closed shape
red_mask_type = imfill(red_bw_type,"holes");
blue_mask_type = imfill(blue_bw_type,"holes");
    
% creates a mask with only the outer points
periosteal_edge = edge(red_mask_border,"sobel");
endocortical_edge = edge(blue_mask_border,"sobel");

%% Finds the points of the endocortial and periosteal lines
% creates a region and gets data from mask with outer points 
% creates structure 

periosteal_struct = regionprops(periosteal_edge,"PixelList");
endocortical_struct = regionprops(endocortical_edge,"PixelList");

% collects all of the periosteal line points

% sometimes, region props above returns a structure of size >1
% this consolidates all of the data we want into one array

allPixelsPerio = [nan,nan]; % placeholder, initiates the array
for i = 1:length(periosteal_struct)
    % collects first element
    a = periosteal_struct(i).PixelList;
    % adds that to the total list
    allPixelsPerio = cat(1,a,allPixelsPerio);
end

% collects all endocortial line points
allPixelsEndo = [nan,nan]; % placeholder, initiates the array

for i = 1:length(endocortical_struct)
    
    % collects nth element
    a = endocortical_struct(i).PixelList;

    % adds that to the total list
    allPixelsEndo = cat(1,a,allPixelsEndo);
end

% the last row will be NaN do to placeholder so remove the placeholder
allPixelsEndo(end,:) = [];
allPixelsPerio(end,:) = [];

% breaks it up to the x and y components
perio_x = allPixelsPerio(:,1);
perio_y = allPixelsPerio(:,2);

endo_x = allPixelsEndo(:,1);
endo_y = allPixelsEndo(:,2);

% this test checks if the endocortial and periosteal lines are being
% identified correctly
if testing

    figure(10)
    clf(10)
    imshow(BorderImage(:,:,1:3))
    hold on
    plot(perio_x,perio_y,".","Color","magenta")
    plot(endo_x,endo_y,".","Color","yellow")
    title("Check to make sure endo (yellow) and perio (pink) lines are being identified correctly")
end

%% finds the marked points for green file

% Seperates images into R and B form
OnlyRed_Seperated_im_green = im_green(:,:,1);
OnlyBlue_Seperated_im_green = im_green(:,:,3);

% Create images of only RGB
OnlyRed_im_green = OnlyRed_Seperated_im_green - OnlyBlue_Seperated_im_green;
OnlyBlue_im_green = OnlyBlue_Seperated_im_green - OnlyRed_Seperated_im_green;

% binarizes
red_bw_im_green = imbinarize(OnlyRed_im_green);
blue_bw_im_green = imbinarize(OnlyBlue_im_green);

% creates a mask
red_mask_im_green = imfill(red_bw_im_green,"holes");
blue_mask_im_green = imfill(blue_bw_im_green,"holes");

% finds the centroid of the marked cells
marked_cells_im_green = centroidFinder(im_green,"g");

% unpacks data
marked_x_im_green = marked_cells_im_green{1};
marked_y_im_green = marked_cells_im_green{2};
marked_remodeling_im_green = marked_cells_im_green{3};

% creates data structure for marked cells
cell_data_green = struct("Centroid_coord","",...
    "closest_endo_point_coord","",...
    "closest_perio_point_coord","",...
    "distance_to_endo", "",...
    "distance_to_perio","", ...    
    "remodeling", "",...
    "total_dist_endo","",...
    "percent","");

for i = 1:length(marked_remodeling_im_green)
    % adds the remodleing and centroid coordinates to the cell_data
    % structure

    cell_data_green(i).remodeling = marked_remodeling_im_green(i);
    cell_data_green(i).Centroid_coord = [marked_x_im_green(i),marked_y_im_green(i)];
    
    % find the minimum distance to endocortial 

    % sets the minimum distance
    min_dist = inf;
    for j = 1:length(endo_x)
        % calculates the distance between the cell and the wall point
        distance = norm([endo_x(j),endo_y(j)] - [marked_x_im_green(i),marked_y_im_green(i)]);
        if min_dist > distance
           % if the min distance to a wall is found, it is the new min dist
            min_dist = distance;
            % save the index of min dist for later use
            loc_j = j;
        end

    end

    cell_data_green(i).distance_to_endo = min_dist;
    cell_data_green(i).closest_endo_point_coord = [endo_x(loc_j),endo_y(loc_j)];

    % find the minimum distance to periosteal 
    
    % sets the minimum distance
    min_dist = inf;
    for k = 1:length(perio_x)
        % calculates the distance between the cell and the wall point
        distance = norm([perio_x(k),perio_y(k)] - [marked_x_im_green(i),marked_y_im_green(i)]);
        if min_dist > distance
           % if the min distance to a wall is found, it is the new min dist
            min_dist = distance;
            % save the index of min dist for later use
            loc_k = k;
        end

    end
    
    cell_data_green(i).distance_to_perio = min_dist;
    cell_data_green(i).closest_perio_point_coord = [perio_x(loc_k),perio_y(loc_k)];

    % calculates the total distance between the closest endocortial and
    % periostal point

    cell_data_green(i).total_dist_endo = norm(cell_data_green(i).closest_endo_point_coord - cell_data_green(i).closest_perio_point_coord );

    % calculates the percent distance from endocortial based on total
    % distance calculated above
    cell_data_green(i).percent = cell_data_green(i).distance_to_endo/cell_data_green(i).total_dist_endo;
    
    % a warning is thrown if the percent is higher than calculated 
    if cell_data_green(i).percent > 1
        warning("Warning, cell data green is recording a percentage above 100% ")
        disp("This might not be an error, depending on the shape of the bone")
        disp("Percent is " + cell_data_green(i).percent)
        disp("Cell distance to endo " + cell_data_green(i).distance_to_endo)
        disp("The cell is located at " + cell_data_green(i).Centroid_coord)
        disp("Total distance is " + cell_data_green(i).total_dist_endo)
        disp("Perio point is " + cell_data_green(i).closest_endo_point_coord)
        disp("Endo point is " + cell_data_green(i).closest_perio_point_coord)
        disp(norm(cell_data_green(i).closest_endo_point_coord - cell_data_green(i).closest_perio_point_coord ))
        
    end
end



%% finds the marked points for red file


% Seperates images into R and B form
OnlyGreen_Seperated_im_red = im_red(:,:,2);
OnlyBlue_Seperated_im_red = im_red(:,:,3);
% Create images of only RGB
OnlyRed_im_red = OnlyGreen_Seperated_im_red - OnlyBlue_Seperated_im_red;
OnlyBlue_im_red = OnlyBlue_Seperated_im_red - OnlyGreen_Seperated_im_red;

% binarizes
red_bw_im_green = imbinarize(OnlyRed_im_red);
blue_bw_im_green = imbinarize(OnlyBlue_im_red);

% creates a mask
red_mask_im_red = imfill(red_bw_im_green,"holes");
blue_mask_im_red = imfill(blue_bw_im_green,"holes");

% finds the centroid of each marked cell
marked_cells_im_red = centroidFinder(im_red,"r");

% unpacks data
marked_x_im_red = marked_cells_im_red{1};
marked_y_im_red = marked_cells_im_red{2};
marked_remodeling_im_red = marked_cells_im_red{3};

% creates data structure for marked cells
cell_data_red = struct("Centroid_coord","",...
        "closest_endo_point_coord","",...
        "closest_perio_point_coord","",...
        "distance_to_endo", "",...
        "distance_to_perio","", ...    
        "remodeling", "",...
        "total_dist_endo","",...
        "percent","");

for i = 1:length(marked_remodeling_im_red)
    % adds the remodleing and centroid coordinates to the cell_data
    % structure

    cell_data_red(i).remodeling = marked_remodeling_im_red(i);
    cell_data_red(i).Centroid_coord = [marked_x_im_red(i),marked_y_im_red(i)];
    
    % find the minimum distance to endocortial 

    % sets the minimum distance
    min_dist = inf;
    for j = 1:length(endo_x)
        % calculates the distance between the cell and the wall point

        distance = norm([marked_x_im_red(i),marked_y_im_red(i)] - [endo_x(j),endo_y(j)]);
        if min_dist > distance
           % if the min distance to a wall is found, it is the new min dist
            min_dist = distance;
            % save the index of min dist for later use
            loc_j = j;
        end

    end

    cell_data_red(i).distance_to_endo = min_dist;
    cell_data_red(i).closest_endo_point_coord = [endo_x(loc_j),endo_y(loc_j)];

    % find the minimum distance to periosteal 
    
    % sets the minimum distance
    min_dist = inf;
    for k = 1:length(perio_x)
        % calculates the distance between the cell and the wall point

        distance = norm([marked_x_im_red(i),marked_y_im_red(i)] - [perio_x(k),perio_y(k)]);
        if min_dist > distance
           % if the min distance to a wall is found, it is the new min dist
            min_dist = distance;
            % save the index of min dist for later use
            loc_k = k;
        end

    end
    
    cell_data_red(i).distance_to_perio = min_dist;
    cell_data_red(i).closest_perio_point_coord = [perio_x(loc_k),perio_y(loc_k)];

    % calculates the total distance between the closest endocortial and
    % periostal point

    cell_data_red(i).total_dist_endo = norm(cell_data_red(i).closest_endo_point_coord - cell_data_red(i).closest_perio_point_coord );

    % calculates the percent distance from endocortial based on total
    % distance calculated above
    cell_data_red(i).percent = cell_data_red(i).distance_to_endo/cell_data_red(i).total_dist_endo;

    % a warning is thrown if the percent is higher than calculated 
    if cell_data_red(i).percent > 1
        warning("Warning, cell data red is recording a percentage above 100% ")
        disp("This might not be an error, depending on the shape of the bone")
        disp("Percent is " + cell_data_red(i).percent)
        disp("Cell distance to endo " + cell_data_red(i).distance_to_endo)
        disp("The cell is located at " + cell_data_red(i).Centroid_coord)
        disp("Total distance is " + cell_data_red(i).total_dist_endo)
        disp("Perio point is " + cell_data_red(i).closest_endo_point_coord)
        disp("Endo point is " + cell_data_red(i).closest_perio_point_coord)
        disp(norm(cell_data_red(i).closest_endo_point_coord - cell_data_red(i).closest_perio_point_coord ))
        
    end
end

%% UNPACKING DATA

for i = 1:length(cell_data_red)

    red_data(i,1) = cell_data_red(i).percent;
    red_data(i,2) = cell_data_red(i).remodeling;
    
end

for i = 1:length(cell_data_green)

    green_data(i,1) = cell_data_green(i).percent;
    green_data(i,2) = cell_data_green(i).remodeling;
    
end

% Analysing data
edges = 0:0.1:1;

% the data is broken up into bins based on the % distance 
red_data_disc = discretize(red_data(:,1),edges);
red_data_disc(:,2) = red_data(:,2);

green_data_disc = discretize(green_data(:,1),edges);
green_data_disc(:,2) = green_data(:,2);

for i = 1:length(edges)-1
    % finds the number for each percentile bin, i
    
    a = find(green_data_disc(:,1) == i);
    
    % gets the remodeling data from function above
    b = green_data_disc(a,:);
    
    % counts number of positive remodeling cells in region
    num_of_positive = length(find(b(:,2)==1));
    num_of_negative = length(find(b(:,2)==0));
    
    % calculates the percent remodeling in that region, and adds it to the
    % array
    percent_remodeling_green(i) = num_of_positive / (num_of_negative + num_of_positive);

end

for i = 1:length(edges)-1
    % finds the number for each percentile bin, i
    
    a = find(red_data_disc(:,1) == i);
    
    % gets the remodeling data from function above
    b = red_data_disc(a,:);
    
    % counts number of positive remodeling cells in region
    num_of_positive = length(find(b(:,2)==1));
    num_of_negative = length(find(b(:,2)==0));
    
    % calculates the percent remodeling in that region, and adds it to the
    % array
    percent_remodeling_red(i) = num_of_positive / (num_of_negative + num_of_positive);

end

% disp(percent_remodeling_green)
% disp(percent_remodeling_red)
% total_remodeling_data(1,:) = percent_remodeling_green;
% total_remodeling_data(2,:) = percent_remodeling_red;

% Labels for the table
row_names = ["0-10%","10-20%","20-30%","30-40%",...
    "40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"];
col_names = ["Green";"Red"];
% creating the table
T = table(percent_remodeling_green',percent_remodeling_red');
T.Properties;
% updating the row and col names
T.Properties.VariableNames = col_names;
T.Properties.RowNames = row_names;
T

if create_excel_file
    excel_file_name = main_folder + "-" +file_path +".xlsx";
    writetable(T,excel_file_name,"WriteRowNames",true)
end

%% FIGURES
if show_figures

figure(1)
clf(1)

plot(percent_remodeling_red,".","MarkerSize",20)
title("Percent remodeling vs normalized distance from endocrotial")
subtitle("-r file")
xlabel("Distance from endocortial (10% intervals)")
ylabel("Percent remodeling ")
ylim([-0.01,1.01])

figure(2)
clf(2)

plot(percent_remodeling_green,".","MarkerSize",20)
title("Percent remodeling vs normalized distance from endocrotial")
subtitle("-g file")
xlabel("Distance from endocortial (10% intervals)")
ylabel("Percent remodeling")
ylim([-0.01, 1.01])
end

% a visual check to make sure all cells are idenitified correctly

if testing

    % this is for the -g file



    figure(11)
    clf(11)
    imshow(BorderImage(:,:,1:3))
    hold on
    for i = 1:length(cell_data_green)
        if cell_data_green(i).remodeling
            plot(cell_data_green(i).Centroid_coord(1),cell_data_green(i).Centroid_coord(2),".b","MarkerSize",25)
        else
            plot(cell_data_green(i).Centroid_coord(1),cell_data_green(i).Centroid_coord(2),".r","MarkerSize",25)
    
        end
    end
    title("For -g file, dots show marked cells, remodeling (blue) or not (red)")
    disp("the marked cells on figure 11 should match figure 12")    

    figure(12)
    clf(12)
    imshow(im_green(:,:,1:3))


    % this is for the -r file
    figure(13)
    clf(13)
    imshow(BorderImage(:,:,1:3))
    hold on
    for i = 1:length(cell_data_red)
        if cell_data_red(i).remodeling
            plot(cell_data_red(i).Centroid_coord(1),cell_data_red(i).Centroid_coord(2),".b","MarkerSize",25)
        else
            plot(cell_data_red(i).Centroid_coord(1),cell_data_red(i).Centroid_coord(2),".g","MarkerSize",25)
    
        end
    end
    title("For -r file, dots show marked cells, remodeling (blue) or not (green)","FontSize",15)
    subtitle("Marked file Should match -r file directly (next figure)")
    
    figure(14)
    clf(14)
    hold on

    imshow(im_red(:,:,1:3))
    disp("the marked cells on figure 13 should match figure 14")    

end
%% FUNCTIONS

function output = centroidFinder(image,g_or_r)
    
    % function finds the centroids of marked points
    % inputs: marked file (*-*-*-g or -r.tiff)
    % outputs: centroid location of all marked cells, and whether the cell
    % is remodeling or not 

    
    RedImage   = image(:,:,1);
    GreenImage = image(:,:,2);
    BlueImage   = image(:,:,3);
    
    
    % Create images of only RGB
    OnlyRed = RedImage - BlueImage;
    OnlyBlue = BlueImage - RedImage;
    OnlyGreen = GreenImage - RedImage;

    % creates a binary image
    OnlyBlue_BW = imbinarize(OnlyBlue);
    OnlyRed_BW = imbinarize(OnlyRed);
    OnlyGreen_BW = imbinarize(OnlyGreen);

    %% for positive marked points
    s1_positive = regionprops(OnlyBlue_BW,"all");

    s1_length_positive = length(s1_positive);
    centroid_x_positive = size(1,s1_length_positive);
    centroid_y_positive = size(1,s1_length_positive);
    
    
    for i = 1:s1_length_positive 
        % stores x and y values of the centroid
        centroid_temp = s1_positive(i).Centroid;
        centroid_x_positive(i) = centroid_temp(1);
        centroid_y_positive(i) = centroid_temp(2);
    end


    %% For negative marked points
    
    if g_or_r == "r"
        % if the file ends in -r, then green marked dots are negative 
        s1_negative = regionprops(OnlyGreen_BW,"all");

        s1_length_negative = length(s1_negative);
        centroid_x_negative = size(1,s1_length_negative);
        centroid_y_negative = size(1,s1_length_negative);
        
        
        for i = 1:s1_length_negative 
            % stores x and y values of the centroid
            centroid_temp = s1_negative(i).Centroid;
            centroid_x_negative(i) = centroid_temp(1);
            centroid_y_negative(i) = centroid_temp(2);
        end
    elseif g_or_r == "g"
        % if the file ends in -g, then red marked dots are negative 
        s1_negative = regionprops(OnlyRed_BW,"all");

        s1_length_negative = length(s1_negative);
        centroid_x_negative = size(1,s1_length_negative);
        centroid_y_negative = size(1,s1_length_negative);
        
        
        for i = 1:s1_length_negative 
            % stores x and y values of the centroid
            centroid_temp = s1_negative(i).Centroid;
            centroid_x_negative(i) = centroid_temp(1);
            centroid_y_negative(i) = centroid_temp(2);
        end

    end
    %% combines the negative and posative marked arrays
    
    centroid_x = cat(1,centroid_x_negative', centroid_x_positive');
    centroid_y = cat(1,centroid_y_negative', centroid_y_positive');
    
    
    remodeling = ones(length(centroid_y),1);
    
    % marks all negative remodeled cells as negative in the array
    for i = 1:length(centroid_x_negative')
        remodeling(i) = false;
    end

    output = {centroid_x,centroid_y,remodeling};


end




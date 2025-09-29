%% Analyze Filler Orientation from Microscope Images
% Author: GEUNUNG YOO
% Date: 2025-09-29
%
% --- DESCRIPTION ---
% This script analyzes a series of microscope images to quantify the
% orientation of rod-like fillers. It performs the following steps:
%   1. Reads multiple images from a specified subfolder.
%   2. Applies image processing techniques (grayscale conversion, filtering,
%      edge detection) to identify fillers.
%   3. Measures the orientation and area of each detected filler using
%      regionprops.
%   4. Filters out noise based on a minimum area threshold.
%   5. Generates an intermediate plot showing the detected filler outlines
%      overlaid in red on the original images for verification.
%   6. Calculates the overall angle distribution using Kernel Density
%      Estimation and visualizes it in a final, publication-quality plot.

%% --- INITIAL SETUP ---
clear; clc; close all;

%% --- USER PARAMETERS ---
% --- File Path Settings ---
params.image_folder = 'Sample_image'; % Name of the subfolder containing images.

% --- List of image files to be analyzed ---
params.image_files = {'horizontal_1.jpg'};

% --- Image Processing Settings ---
params.noise_filter_size = [3 3];   % Size of the median filter for noise reduction.
params.min_filler_area   = 10;      % Minimum pixel area to be considered a valid filler.
params.outline_thickness = 3;       % Thickness of the red outlines. Increase for thicker lines (e.g., 2, 3).

% --- Final Plot Settings ---
params.kde_bandwidth = 5;         % Bandwidth for the Kernel Density Estimation (controls smoothness).


%% --- SCRIPT START ---
[script_folder, ~, ~] = fileparts(mfilename('fullpath'));
all_angles = []; % Initialize an array to store angles from all images

% Create a figure for intermediate visualization results
figure('Name', 'Image Processing Steps', 'Color', 'w');
num_images = length(params.image_files);

% Loop through each image file for analysis
for k = 1:num_images
    % Construct the full path to the image file
    full_path_to_image = fullfile(script_folder, params.image_folder, params.image_files{k});
    fprintf('Processing image: %s\n', full_path_to_image);
    
    % 1. Image Loading and Pre-processing
    img_rgb = imread(full_path_to_image);
    img_gray = rgb2gray(img_rgb);
    img_filtered = medfilt2(img_gray, params.noise_filter_size);
    
    % 2. Edge Detection
    edges = edge(img_filtered, 'sobel');
    
    % 3. Identify Regions and Extract Properties
    % 'PixelIdxList' is used to reconstruct the shape of valid fillers later.
    regions = regionprops(edges, 'Orientation', 'Area', 'PixelIdxList');
    
    % 4. Filter out small regions (noise)
    valid_regions = regions([regions.Area] > params.min_filler_area);
    
    % Store the orientation angles of valid fillers
    valid_angles = [valid_regions.Orientation];
    all_angles = [all_angles, valid_angles];
    
    % 5. Visualize the detected outlines for the current image
    % Reconstruct a binary image containing only the valid fillers
    labeled_fillers = false(size(img_gray));
    for i = 1:length(valid_regions)
        labeled_fillers(valid_regions(i).PixelIdxList) = true;
    end
    
    % Get the perimeter (outline) of the valid fillers
    filler_outlines = bwperim(labeled_fillers);
    
    % Make the outlines thicker for better visibility
    se = strel('disk', params.outline_thickness);
    thicker_outlines = imdilate(filler_outlines, se);
    
    % --- Display Original vs. Outlined Image ---
    % Left subplot: Original grayscale image
    subplot(num_images, 2, 2*k - 1);
    imshow(img_gray);
    title(['Original: ' params.image_files{k}]);
    
    % Right subplot: Original image with red outlines overlaid
    subplot(num_images, 2, 2*k);
    % Create an RGB version of the grayscale image
    overlay_display = repmat(img_gray, [1, 1, 3]);
    % Get linear indices of the outline pixels
    outline_indices = find(thicker_outlines);

    % Modify each color channel individually to create the red overlay
    red_channel = overlay_display(:,:,1);
    green_channel = overlay_display(:,:,2);
    blue_channel = overlay_display(:,:,3);
    
    red_channel(outline_indices) = 255;   % Set Red to max
    green_channel(outline_indices) = 0;     % Set Green to 0
    blue_channel(outline_indices) = 0;      % Set Blue to 0
    
    % Recombine the channels
    overlay_display(:,:,1) = red_channel;
    overlay_display(:,:,2) = green_channel;
    overlay_display(:,:,3) = blue_channel;

    imshow(overlay_display);
    title(['Detected Outlines: ' params.image_files{k}]);
end

%% --- Final Analysis and Plotting ---
disp('All images processed. Generating final distribution plot...');

% Create a new figure for the final distribution plot
figure('Name', 'Final Angle Distribution Plot', 'Color', 'w');

% Use Kernel Density Estimation to get a smooth distribution curve
kde_fit = fitdist(all_angles', 'Kernel', 'BandWidth', params.kde_bandwidth);

% Define the x-axis range for the plot
x_values = linspace(-90, 90, 1000);

% Calculate the probability density function values
pdf_values = pdf(kde_fit, x_values);

% Normalize the PDF so the peak is at 1.0 for clear visualization
pdf_values = pdf_values / max(pdf_values);

% Create the final plot with professional styling
plot(x_values, pdf_values, 'LineWidth', 3);
grid on;
box on;
xlabel('Alignment Angle (degrees)');
ylabel('Normalized Probability Density');
title('Overall Filler Alignment Distribution');
set(gca, 'XTick', -90:30:90, 'XLim', [-90 90], 'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

disp('--- Analysis Complete ---');
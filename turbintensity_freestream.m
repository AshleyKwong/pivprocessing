% Ashley Kwong
% 02/10/2026
% This code takes the calibration coefficients (polynomial) and performs
% multicamera planar stitching for PIV. Assumes that you have the
% calibration struct from calib_tools and the vector fields from pivtools.
clear; clc; close all;
set(groot, 'defaultAxesFontName', 'Cambria Math');
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontName', 'Cambira Math');
set(groot, 'defaultTextFontSize', 12);
set(groot, 'DefaultFigureVisible', 'off');  % Turn off figure visibility globally

cameraList  = ["Cam1", "Cam2", "Cam3", "Cam4", "Cam5"];  
savePath = 'E:\Case 6 PIV pivtools Processed\'; % path where data lives and will be saved.

%% USER OPTIONS
SKIP_LOOP_PROCESSING = true;  % Set to true to skip loop processing (if already done)
FORCE_RECOMPUTE_AVG = false;   % Set to true to recompute average even if it exists

%%
% the size of the image should correspond to the physical location
% import the calibration
addpath(genpath('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\calib_tools\'));
load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\calib_tools\poly_2d2c\Step_01_Calibration\y235_aoan11_aoafn11\calib_AK.mat');% calib struct

% frame is the size of the raw PIV image --> need this for dewarping !
% case 6
frame = {zeros(3528,5312, 'single'),zeros(3472,5312,'single'),zeros(3512,5312,'single'),zeros(3536,5312,'single'), zeros(3520,5312,'single')}; % size of each camera from the calibration grid !

%% DEWARPING
cx_img = est_dewarp_size_2d(frame,calib);
[~,cx_img,px_img] = dewarp_image_2d(frame,calib, cx_img);

% -------- MAKE THE WORLD GRID. -----------------------------------
x1_min = cellfun(@(c) min(c.x1(1,:)), cx_img, 'UniformOutput', true);
x1_max = cellfun(@(c) max(c.x1(1,:)), cx_img, 'UniformOutput', true);
x2_min = cellfun(@(c) min(c.x2(:,1)), cx_img, 'UniformOutput', true);
x2_max = cellfun(@(c) max(c.x2(:,1)), cx_img, 'UniformOutput', true);

dt = 49.6 * 10^-6; % bw frames
magFactor = zeros(2,5);
for i = 1:5
    magFactor(1, i) = calib.var{i}.mmppix(1);
    magFactor(2, i) = calib.var{i}.mmppix(2);
end

mean_del_x1 = mean(magFactor(1,:));
mean_del_x2 = mean(magFactor(2, :));
windowsize = [16,16];
overlap_percentage = 20;
dx1 = mean_del_x1*windowsize(2)*(1-overlap_percentage/100);
dx2 = mean_del_x1*windowsize(1)*(1-overlap_percentage/100);
[worldx1, worldx2]= meshgrid(min(x1_min):dx1:max(x1_max), min(x2_max):-dx2:max(x2_min));
worldx1 = (single(worldx1));
worldx2= single(worldx2);

%% MAIN PROCESSING LOOP
% Get all directories and filter for "loop=XX" format only
d = dir(savePath); 
d = d([d.isdir]);  % Only directories
d = d(~ismember({d.name}, {'.', '..'}));  % Remove . and ..

% Filter for folders matching "loop=XX" pattern (case-insensitive)
loopPattern = '^loop\s*=\s*\d+$';  % Matches "loop=0", "loop = 1", "loop=99", etc.
validLoops = false(size(d));
for i = 1:length(d)
    validLoops(i) = ~isempty(regexpi(d(i).name, loopPattern));
end
totalLoops = d(validLoops);

if isempty(totalLoops)
    warning('No folders matching "loop=XX" pattern found in %s', savePath);
else
    fprintf('Found %d valid loop folders:\n', length(totalLoops));
    for i = 1:length(totalLoops)
        fprintf('  - %s\n', totalLoops(i).name);
    end
end
% SKIP LOOP PROCESSING IF ALREADY DONE
if SKIP_LOOP_PROCESSING
    fprintf('\n⚠ SKIP_LOOP_PROCESSING = true. Skipping loop processing...\n');
    fprintf('Loading existing data for averaging and visualization.\n\n');
else
    tic;
    final_pass_index = 5;  % Renamed from numOfWindows for clarity
    skip_frames = []; % Define frames to skip if needed

    for loopNo = 1:length(totalLoops)
        fprintf("....Currently processing %s (loop %d/%d)\n", totalLoops(loopNo).name, loopNo, length(totalLoops));
        saveFolder = fullfile(savePath, totalLoops(loopNo).name);
        
        % Check if physical window centers already exist
        physicalLocationFile = fullfile(saveFolder, 'windowCenterCameras_mm.mat');
        if isfile(physicalLocationFile)
            load(physicalLocationFile, 'windowCenterCameras_mm');
            fprintf("\t✓ Loaded existing physical window centers (mm) for all cameras\n");
            skipCalibration = true;
        else
            skipCalibration = false;
            % Initialize struct with descriptive field names
            % Each camera has unique physical coordinates based on its calibration
            windowCenterCameras_mm = struct();
            windowCenterCameras_mm.x1_mm = cell(1, 5);  % Physical x1 coordinate in mm
            windowCenterCameras_mm.x2_mm = cell(1, 5);  % Physical x2 coordinate in mm
            windowCenterCameras_mm.description = 'Physical window center locations in mm for each camera (post-calibration)';
            windowCenterCameras_mm.cameras = cameraList;
        end
        
        % Get number of frames
        d = dir(fullfile(saveFolder, 'uncalibrated_piv'));
        numImages = d(~ismember({d.name}, {'.', '..'}));
        
        % Get first camera's file list to determine total frames
        base_dir_sample = dir(fullfile(saveFolder, 'uncalibrated_piv', numImages(1).name, ...
            'Cam1', 'instantaneous', '*.mat'));
        nTotalFrames = length(base_dir_sample) - 1;
        nValidFrames = nTotalFrames - length(skip_frames);
        
        % Preallocate allCameras struct
        allCameras = struct();
        allCameras.u = cell(nValidFrames, 5);
        allCameras.v = cell(nValidFrames, 5);
        
        % Process each camera
        for a = 1:length(cameraList)
            fprintf("\t Processing %s (loop %s)\n", cameraList(a), totalLoops(loopNo).name);
            
            % Get file list for this camera
            base_dir = dir(fullfile(saveFolder, 'uncalibrated_piv', numImages(1).name, ...
                ['Cam' num2str(a)], 'instantaneous', '*.mat'));
            
            % Load first frame to get window centers (if not already loaded)
            if ~skipCalibration
                firstFrameData = load(fullfile(base_dir(1).folder, base_dir(1).name));
                window_centers_y = single(firstFrameData.piv_result(final_pass_index).win_ctrs_y);
                window_centers_x = single(firstFrameData.piv_result(final_pass_index).win_ctrs_x);
                
                % Compute window center grid ONCE per camera (in pixels)
                [wincenter_px1, wincenter_px2] = meshgrid(window_centers_x, window_centers_y);
                
                % Convert pixel centers to mm ONCE per camera
                % Each camera has unique calibration transformation
                cx_windowcenterx1_mm = polyvaln(calib.fit_cx1{a}.p, [wincenter_px1(:), wincenter_px2(:)]);
                cx_windowcenterx2_mm = polyvaln(calib.fit_cx2{a}.p, [wincenter_px1(:), wincenter_px2(:)]);
                cx_windowcenterx1_mm = reshape(cx_windowcenterx1_mm, size(wincenter_px1));
                cx_windowcenterx2_mm = reshape(cx_windowcenterx2_mm, size(wincenter_px2));
                
                % Store physical locations (unique per camera, no truncation)
                windowCenterCameras_mm.x1_mm{a} = cx_windowcenterx1_mm;
                windowCenterCameras_mm.x2_mm{a} = cx_windowcenterx2_mm;
                
                fprintf("\t  ✓ Computed physical locations for %s: [%d × %d] grid\n", ...
                    cameraList(a), size(cx_windowcenterx1_mm, 1), size(cx_windowcenterx1_mm, 2));
                
                clear firstFrameData;
            else
                % Load from saved file (camera-specific physical locations)
                cx_windowcenterx1_mm = windowCenterCameras_mm.x1_mm{a};
                cx_windowcenterx2_mm = windowCenterCameras_mm.x2_mm{a};
                
                fprintf("\t  ✓ Using saved physical locations for %s: [%d × %d] grid\n", ...
                    cameraList(a), size(cx_windowcenterx1_mm, 1), size(cx_windowcenterx1_mm, 2));
                
                % Reconstruct pixel grid for velocity transformation
                window_centers_y = linspace(1, size(cx_windowcenterx1_mm, 1), size(cx_windowcenterx1_mm, 1));
                window_centers_x = linspace(1, size(cx_windowcenterx1_mm, 2), size(cx_windowcenterx1_mm, 2));
                [wincenter_px1, wincenter_px2] = meshgrid(window_centers_x, window_centers_y);
            end
            
            % Process all frames for this camera
            counter = 0;
            for b = 1:length(base_dir)-1
                if ismember(b, skip_frames)
                    continue;
                end
                counter = counter + 1;

                frameFile = fullfile(base_dir(b).folder, base_dir(b).name);
                frameData = load(frameFile, 'piv_result');  % Load only piv_result field

                % Extract velocity in pixels
                velocity_uy = single(frameData.piv_result(final_pass_index).uy);
                velocity_ux = single(frameData.piv_result(final_pass_index).ux);
                clear frameData;

                % Compute pixel positions after displacement
                px1_grid = wincenter_px1 + velocity_uy;
                px2_grid = wincenter_px2 + velocity_ux;
                
                % Convert displaced positions to mm using camera-specific calibration
                cx_velgrid_x1_mm = polyvaln(calib.fit_cx1{a}.p, [px1_grid(:), px2_grid(:)]);
                cx_velgrid_x2_mm = polyvaln(calib.fit_cx2{a}.p, [px1_grid(:), px2_grid(:)]);
                cx_velgrid_x1_mm = reshape(cx_velgrid_x1_mm, size(velocity_uy));
                cx_velgrid_x2_mm = reshape(cx_velgrid_x2_mm, size(velocity_uy));
                
                % Compute velocity in m/s (no truncation)
                % Uses camera-specific physical window centers
                x1_ms = (cx_velgrid_x1_mm - cx_windowcenterx1_mm) * 1e-3 / dt;
                x2_ms = (cx_velgrid_x2_mm - cx_windowcenterx2_mm) * 1e-3 / dt;
                
                % Store without truncation
                allCameras.u{counter, a} = x2_ms;
                allCameras.v{counter, a} = x1_ms;
            end
            
            fprintf("\t  ✓ Processed %d frames for camera %d\n", counter, a);
        end
        
        % Save physical window centers if newly computed
        if ~skipCalibration
            save(physicalLocationFile, 'windowCenterCameras_mm', '-v7.3');
            fprintf("\t✓ Saved physical window centers (mm) for all cameras - unique per camera\n");
        end
        
        % Save velocity fields
        S = whos('allCameras');
        fprintf("\t Memory: %.2f GB\n", S.bytes/1e9);
        
        velocityFile = fullfile(saveFolder, 'processedvelocityfields.mat');
        save(velocityFile, 'allCameras', '-v7.3');
        fprintf("\t✓ Saved velocity fields for loop %d\n", loopNo);
        
        clear allCameras;
    end

    endTime = toc;
    fprintf("\n✓ Total processing time: %.2f min\n", endTime/60);
end

%% VISUALIZATION: Check single frame
% Commented out GIF generation - uncomment only if needed
% generateCorrelationGifs(saveFolder, cameraList, final_pass_index);

%% AVERAGE VELOCITY FIELDS
% Check if average already exists

avgFileName = sprintf('averagedvelfields_uv_%d.mat', length(totalLoops)*150); 
avgFilePath = fullfile(savePath, avgFileName);

if isfile(avgFilePath) && ~FORCE_RECOMPUTE_AVG
    fprintf('\n✓ Found existing averaged velocity fields: %s\n', avgFileName);
    fprintf('  Loading from file...\n');
    avgVelocityField= load(avgFilePath).meanCameras; % this will load in the 
    fprintf('  ✓ Loaded averaged fields (set FORCE_RECOMPUTE_AVG=true to recompute)\n\n');
else
    if FORCE_RECOMPUTE_AVG
        fprintf('\n⚠ FORCE_RECOMPUTE_AVG = true. Recomputing average...\n\n');
    else
        fprintf('\n⚠ Averaged velocity fields not found. Computing...\n\n');
    end
    avgVelocityField = averageVelocityFields(savePath, totalLoops);
end

function meanCameras = averageVelocityFields(savePath, totalLoops)
    tic;
    fprintf('\n=== Averaging Velocity Fields ===\n');
    
    % Get number of frames and cameras from first file
    firstPath = fullfile(savePath, totalLoops(1).name, "processedvelocityfields.mat");
    firstData = load(firstPath, 'allCameras');
    [nFrames, nCameras] = size(firstData.allCameras.u);
    
    % Initialize accumulators
    sumU = firstData.allCameras.u;
    sumV = firstData.allCameras.v;
    clear firstData;
    
    fprintf('\tInitialized with loop 1/%d\n', length(totalLoops));
    
    % Accumulate remaining loops
    for k = 2:length(totalLoops)
        nextPath = fullfile(savePath, totalLoops(k).name, "processedvelocityfields.mat");
        nextData = load(nextPath, 'allCameras');
        
        % Add element-wise
        for frame = 1:nFrames
            for cam = 1:nCameras
                sumU{frame, cam} = sumU{frame, cam} + nextData.allCameras.u{frame, cam};
                sumV{frame, cam} = sumV{frame, cam} + nextData.allCameras.v{frame, cam};
            end
        end
        
        clear nextData;
        fprintf('\tAdded loop %d/%d\n', k, length(totalLoops));
    end
    
    % Average across loops
    numLoops = length(totalLoops);
    for frame = 1:nFrames
        for cam = 1:nCameras
            sumU{frame, cam} = sumU{frame, cam} / numLoops;
            sumV{frame, cam} = sumV{frame, cam} / numLoops;
        end
    end
    
    % Average across frames (150 frames -> 1 per camera)
    meanU = cell(1, nCameras);
    meanV = cell(1, nCameras);
    
    for cam = 1:nCameras
        % Stack all frames for this camera
        stackU = cat(3, sumU{:, cam});
        stackV = cat(3, sumV{:, cam});
        
        % Average across 3rd dimension
        meanU{cam} = mean(stackU, 3);
        meanV{cam} = mean(stackV, 3);
    end
    
    meanCameras.u = meanU;
    meanCameras.v = meanV;
    
    endTime = toc;
    fprintf('✓ Averaged %d loops × %d frames to 1×%d fields in %.1f sec\n', ...
        numLoops, nFrames, nCameras, endTime);
    
    % Save
    S = whos('meanCameras');
    fprintf("\tMemory: %.2f GB\n", S.bytes/1e9);
    
    avgFramesNo = sprintf('averagedvelfields_uv_%d', length(totalLoops)*150);
    filename = fullfile(savePath, [avgFramesNo, '.mat']);
    save(filename, 'meanCameras', '-v7.3');
    fprintf("✓ Saved %s\n", avgFramesNo);
end

%% VISUALIZATION: Mean velocity fields
figure('Position', [100, 100, 1200, 800], 'Visible', 'on');
meanU = avgVelocityField.u(1, :);

for cam = 1:5
    subplot(2, 3, cam);
    imagesc(meanU{cam});
    colormap(redblue(10));
    colorbar();
    clim([0 30]);
    title(sprintf('Camera %d: Mean U', cam));
    xlabel('X'); ylabel('Y');
    axis image;
end
sgtitle('Raw Mean U Velocity Fields Across 5 Cameras', 'FontSize', 14);

%% CAMERA MERGING - Median with Tukey window
velocityU = avgVelocityField.u(1, :);
velocityV = avgVelocityField.v(1, :);
masks = {};
% addpath for the functions
addpath(genpath('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV')); 


% Load physical locations for merging (from first loop folder)
load(fullfile(savePath, totalLoops(1).name, 'windowCenterCameras_mm.mat'), 'windowCenterCameras_mm');

[worldX, worldY, U_tukey, ~] = merge_cameras_python_style_median(...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'tukey', 0.5);

[~, ~, U_hann, ~] = merge_cameras_python_style_median(...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'hann', []);

% Compare
figure('Visible', 'on');
subplot(2,1,1); 
imagesc(worldX(1,:), flipud(worldY(:, 1)), U_tukey); 
title('Tukey'); colorbar; colormap(jet); 
axis image;

subplot(2,1,2); 
imagesc(worldX(1,:), flipud(worldY(:, 1)), U_hann); 
title('Hann'); colorbar;colormap(jet); 
axis image;

sgtitle("Median Merging: Tukey vs Hann");

%% CAMERA MERGING - Mean with Tukey window
[worldX, worldY, U_tukey_mean, ~] = merge_cameras_python_style_mean(...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'tukey', 0.5);

[~, ~, U_hann_mean, ~] = merge_cameras_python_style_mean(...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'hann', []);

% Compare
figure('Visible', 'on');
subplot(2,1,1); 
imagesc(worldX(1,:), flipud(worldY(:, 1)), U_tukey_mean); 
title('Tukey'); colorbar;
axis image;

subplot(2,1,2); 
imagesc(worldX(1,:), flipud(worldY(:, 1)), U_hann_mean); 
title('Hann'); colorbar;
axis image;

colormap(parula);
sgtitle("Mean Merging: Tukey vs Hann");

%% OPTIONAL: Function to generate correlation GIFs (commented out by default)
function generateCorrelationGifs(saveFolder, cameraList, final_pass_index)
    % Uncomment and call this function only when you need correlation diagnostics
    d = dir(fullfile(saveFolder, 'uncalibrated_piv'));
    numImages = d(~ismember({d.name}, {'.', '..'}));
    
    for a = 1:length(cameraList)
        fprintf("\tGenerating GIF for %s\n", cameraList(a));
        gifName = fullfile(saveFolder, sprintf('peakMag_Camera%d.gif', a));
        delayTime = 0.05;
        
        camDir = dir(fullfile(saveFolder, 'uncalibrated_piv', numImages(1).name, ...
            ['Cam' num2str(a)], 'instantaneous', '*.mat'));
        nFrames = numel(camDir);
        
        fig = figure('Visible', 'off');
        for f = 1:nFrames
            cameraData = load(fullfile(camDir(f).folder, camDir(f).name));
            pm = squeeze(cameraData.piv_result(final_pass_index).peak_mag);
            imagesc(pm);
            axis image off;
            colorbar;
            title(sprintf('Camera %d - Frame %d', a, f));
            drawnow;
            
            frame = getframe(fig);
            [imind, cm] = rgb2ind(frame2im(frame), 256);
            
            if f == 1
                imwrite(imind, cm, gifName, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
            else
                imwrite(imind, cm, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
            end
            
            clear cameraData pm imind cm frame;
        end
        close(fig);
    end
end

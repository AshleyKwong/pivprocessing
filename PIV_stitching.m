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
cameraList  = ["Cam1", "Cam2", "Cam3", "Cam4", "Cam5"];  

savePath = 'E:\Case 6 PIV pivtools Processed\'; % path where data lives and will be saved.


%%
% the size of the image should correspond to the physical location
% import the calibration
addpath(genpath('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\calib_tools\'));

load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\calib_tools\poly_2d2c\Step_01_Calibration\y235_aoan11_aoafn11\calib_AK.mat');% calib struct
% load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\calib_tools\poly_2d2c\Step_01_Calibration\PIV_calibration_emptytunnel\calib_AK.mat');% calib struct

% for each camera.
% polyval the calib.fit_cx1{1,1}.p.Coefficients --> in the x dir and y dir
% frame is the size of the raw PIV image --> need this for dewarping !
% case 6
frame = {zeros(3528,5312, 'single'),zeros(3472,5312,'single'),zeros(3512,5312,'single'),zeros(3536,5312,'single'), zeros(3520,5312,'single')}; % size of each camera from the calibration grid !
% empty tunnel
% frame = {zeros(4600,5312, 'single'),zeros(4600,5312,'single'),zeros(4600,5312,'single'),zeros(4600,5312,'single'), zeros(4600,5312,'single')}; % size of each camera from the calibration grid !

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

%% checking the correlation to see if need to drop any frames.
d = dir(savePath); d= d(~ismember({d.name}, {'.', '..'}));
totalLoops = d;
tic;
numOfWindows = 5;

for loopNo = 1:length(totalLoops) % might want to change the num of loops to just ensemble running

    fprintf("....Currently doing %s\n", totalLoops(loopNo).name);
    saveFolder = sprintf('E:\\Case 6 PIV pivtools Processed\\%s',totalLoops(loopNo).name);

    correlation_data= struct();
    correlation_data.droppedFrames = cell(1, 5);

    % checking instantaneous frames - pixel displacement and GIF of peak_mag
    d = dir(fullfile(saveFolder, '\uncalibrated_piv\'));
    numImages = d(~ismember({d.name}, {'.', '..'}));          % remove . and ..

    pixelDisplacementFrames          = cell(1, numel(cameraList)); % check pixel displacement
    peak_magFrames       = cell(1, numel(cameraList)); % gets correlation value


    for a = 1:length(cameraList) % loop through each camera
        fprintf("\t On %s\n", cameraList(a));

        gifName   = fullfile(saveFolder, sprintf('peakMag_Camera%d.gif', a));
        delayTime = 0.05;  % seconds per frame

        % Get list of MAT files for this camera once
        camDir = dir(fullfile(saveFolder, 'uncalibrated_piv', numImages(1).name, ...
            ['Cam' num2str(a)], 'instantaneous', '*.mat'));
        nFrames = numel(camDir);  % number of frames

        % Create invisible figure once
        fig = figure('Visible','off');

        for f = 1:nFrames
            fprintf("On frame %i\n", f);

            cameraData = load(fullfile(camDir(f).folder, camDir(f).name));
            pm = squeeze(cameraData.piv_result(numOfWindows).peak_mag);

            imagesc(pm);
            axis image off;
            colorbar;
            title(sprintf('Camera %d - Frame %d', a, f));
            drawnow;

            frame = getframe(fig);
            [imind, cm] = rgb2ind(frame2im(frame), 256);

            if f == 1
                imwrite(imind, cm, gifName, 'gif', ...
                    'LoopCount', inf, 'DelayTime', delayTime);
            else
                imwrite(imind, cm, gifName, 'gif', ...
                    'WriteMode', 'append', 'DelayTime', delayTime);
            end

            % Optional: free loaded data early
            clear cameraData pm imind cm frame;
        end

        close(fig);
    end

end
%%
tic;
skip_frames = []; % unique([correlation_data.droppedFrames{:}]); % these are frames to avoid
allCameras= struct();
allCameras.u = cell(150 - length(skip_frames), 5);
allCameras.v = cell(150 - length(skip_frames),5); % init empty cell to hold the camera information

d = dir(savePath); d= d(~ismember({d.name}, {'.', '..'}));
totalLoops = d;
tic;
numOfWindows = 5;
windowCenterCameras= struct(); % per loop need to reinit.
windowCenterCameras.winx1 = cell(1,5);
windowCenterCameras.winx2  =  cell(1,5);

for loopNo = 1: length(totalLoops) % might want to change the num of loops to just ensemble running
    for a = 1:length(cameraList) % loop through each camera
        fprintf("\t On %s loop %s.\n", cameraList(a),  totalLoops(loopNo).name);
        dropped_frames = [];
        counter = 0;
        saveFolder = [savePath, totalLoops(loopNo).name];
        % checking instantaneous frames - pixel displacement and GIF of peak_mag
        d = dir(fullfile(saveFolder, '\uncalibrated_piv\'));
        numImages = d(~ismember({d.name}, {'.', '..'}));          % remove . and ..

        base_dir = dir(fullfile(saveFolder, '//uncalibrated_piv//', numImages(1).name, ...
            ['//Cam' num2str(a)], '//instantaneous//', '*.mat')); % has all 150 images.
        
        for b = 1:length(base_dir)-1 % for num of images.
            if ismember(b, skip_frames) % check if the frame is in the skip list
                continue; % skip processing for this frame
            else
                counter = counter+1;
            end
            cameraData = load(fullfile([base_dir(b).folder, '\', base_dir(b).name])); % loads the stuff for SPECIFIC CAMERA !
            velocity_dpixel = {single(cameraData.piv_result(numOfWindows).uy), single(cameraData.piv_result(numOfWindows).ux)} ; % pixel displacement.
            window_centers = {single(cameraData.piv_result(numOfWindows).win_ctrs_y), single(cameraData.piv_result(numOfWindows).win_ctrs_x)}; % window centers in pixels
            % need to make a velocity grid the same same size as the px grid !

            % each vector is the displacement FROM that center.
            %
            [wincenter_x1 ,  wincenter_x2] = meshgrid(window_centers{2}, window_centers{1});

            px1_grid = wincenter_x1 + velocity_dpixel{1}; % this is uy
            px2_grid = wincenter_x2 + velocity_dpixel{2}; % this is ux


            cx_velgrid_x1 = polyvaln(calib.fit_cx1{a}.p, [px1_grid(:), px2_grid(:)]);
            cx_velgrid_x2 = polyvaln(calib.fit_cx2{a}.p, [px1_grid(:), px2_grid(:)]);

            cx_velgrid_x1 = reshape(cx_velgrid_x1, size(velocity_dpixel{1}));
            cx_velgrid_x2 = reshape(cx_velgrid_x2, size(velocity_dpixel{1}));

            % now need to convert the pixel centers in pix to mm
            cx_windowcenterx1 = polyvaln(calib.fit_cx1{a}.p, [wincenter_x1(:), wincenter_x2(:)]);
            cx_windowcenterx2 = polyvaln(calib.fit_cx2{a}.p, [wincenter_x1(:), wincenter_x2(:)]);

            cx_windowcenterx1 = reshape(cx_windowcenterx1, size(velocity_dpixel{1}));
            cx_windowcenterx2 = reshape(cx_windowcenterx2, size(velocity_dpixel{1}));

            x1_ms = ( cx_velgrid_x1 - cx_windowcenterx1) * 10^-3 /dt ; % make it into mm.
            x2_ms = (cx_velgrid_x2 - cx_windowcenterx2 )* 10^-3 /dt ;

            allCameras.u{counter, a} = x2_ms(2:end -1, 2:end-1); % truncate the edges on all 4 sides.
            allCameras.v{counter, a} = x1_ms(2:end -1, 2:end-1);

        end % end of 150 images / camera loop
        fprintf("------Finished processing %d images for camera %d\n", str2num(numImages.name)-length(skip_frames), a )
        windowCenterCameras.winx1{a} = cx_windowcenterx1(2:end -1, 2:end-1); % might need to end up truncating more ?
        windowCenterCameras.winx2{a} = cx_windowcenterx2(2:end -1, 2:end-1);

    end % end for cameras
    S = whos('allCameras');
    bytesUsed = S.bytes/10^9;
    fprintf("\tThe total cell is %.2f GB\n", bytesUsed); % chonk bois

    % Create a MAT-file object with write access and version 7.3 for large data
    filename = fullfile([saveFolder,'\processedvelocityfields.mat']);
    save(filename, 'allCameras', '-v7.3');
    fprintf("Finished saving for loop %d. \n", loopNo);
    clear allCameras; % save space.
end % end of loops

endTime = toc;
fprintf("Total time to process all the loops was %.2f min.\n", endTime/60);

%%
figure();
t = tiledlayout(length(cameraList), 1, 'TileSpacing', 'compact'); % nicer than subplot
title(t, "Velocity - U m/s");

for i = 1:length(cameraList)
    nexttile;
    U = allCameras.u(1,:);       % first frame cell? (see note below)
    Ui = U{i};                 % actual 2D matrix

    imagesc(Ui);               % 2D field, e.g. 574x663
    axis image;                % same as axis equal + square pixels
    colormap(redblue(10));
    clim([1.1*min(Ui(:)) 1.1*max(Ui(:))]);             % or caxis([-30 0]);
    colorbar;

    xlabel("pix");
    ylabel("pix");
    title(sprintf('Camera %d', i));
end

% clear allCameras;

%%

avgVelocityField = averageVelocityFields(savePath, totalLoops); 

function meanCameras = averageVelocityFields(savePath, totalLoops)
tic;
% Load FIRST dataset to init sums (150x5 cells)
firstPath = fullfile(savePath, totalLoops(1).name, "processedvelocityfields.mat");
load(firstPath);
sumU = allCameras.u;  % Direct init, no zeros needed
sumV = allCameras.v;
numDatasets = 1;  % Or numFrames = 150 if fixed
clear allCameras processedvelocityfield;
fprintf('\tInitialized with first loop (%d datasets).\n', numDatasets);

% Accumulate remaining
for k = 1:length(totalLoops)
    nextPath = fullfile(savePath, totalLoops(k).name, "processedvelocityfields.mat");
    load(nextPath);
    sumU = cellfun(@plus, sumU, allCameras.u, 'UniformOutput', false);
    sumV = cellfun(@plus, sumV, allCameras.v, 'UniformOutput', false);
    numDatasets = numDatasets + 1;
    clear allCameras processedvelocityfield;
    fprintf('\tAdded loop %d (%d total datasets)\n', k, numDatasets);
end

% Average across datasets: 150x5 -> divide each cell
sumU = cellfun(@(s) s / numDatasets, sumU, 'UniformOutput', false);
sumV = cellfun(@(s) s / numDatasets, sumV, 'UniformOutput', false);

% Final: average down 150 frames per camera -> 1x5
% sumU, sumV are 150x5 cells of 291x1325 arrays

colsU = num2cell(sumU, 1);    % 1x5, each cell is 150x1 column
colsV = num2cell(sumV, 1);

meanU = cellfun(@(col) mean(cat(3, col{:}), 3), colsU, 'UniformOutput', false);
meanV = cellfun(@(col) mean(cat(3, col{:}), 3), colsV, 'UniformOutput', false);
% meanU, meanV are 1x5 cells of 291x1325
meanCameras.u = meanU;
meanCameras.v = meanV;
endTime = toc;
fprintf('Averaged %d datasets to 1x5 fields in %.1f sec.\n', numDatasets, endTime);
S = whos('meanCameras');
bytesUsed = S.bytes/10^9;
fprintf("\tThe total cell is %.2f GB\n", bytesUsed); % chonk bois
avgFramesNo = sprintf('averagedvelfields_uv_%d', length(totalLoops)); 

% Create a MAT-file object with write access and version 7.3 for large data
filename = fullfile([savePath,avgFramesNo,'.mat']);
save(filename, 'meanCameras', '-v7.3');
fprintf("Finished saving %s", avgFramesNo);
end

%%
figure('Position', [100, 100, 1200, 800]);  % Wide figure for 5 subplots
meanU = avgVelocityField.u(1, :); 
for cam = 1:5
    subplot(2, 3, cam);  % 2 rows optional for colorbar
    
    % Plot raw meanU{cam} as image/quiver (adjust based on shape)
    if numel(meanU{cam}) == numel(meanU{cam}(:))  % 2D field
        imagesc(meanU{cam});  % Or pcolor/heatmap
        colormap(redblue(10))
        colorbar();clim([0 30])
        title(sprintf('Camera %d: Mean U', cam));
        xlabel('X'); ylabel('Y');
    else
        quiver(meanU{cam});  % Vector field
        title(sprintf('Camera %d: Mean U Vectors', cam));
    end
end

sgtitle('Raw Mean U Velocity Fields Across 5 Cameras', 'FontSize', 14);
% caxis([min(cellfun(@min, meanU)), max(cellfun(@max, meanU))]);  % Shared color limits

%% Prepare inputs
velocityU = avgVelocityField.u(1, :);  % Cell array of U fields
velocityV = avgVelocityField.v(1, :);  % Cell array of V fields
masks = {};         % Empty if no masks, or cell array of logical masks

% Run with Tukey window (Python default)
[worldX, worldY, U_tukey, ~] = merge_cameras_python_style_median(...
    windowCenterCameras, velocityU, velocityV, masks, 'tukey', 0.5);

% Or try Hann window
[~, ~, U_hann, ~] = merge_cameras_python_style_median(...
    windowCenterCameras, velocityU, velocityV, masks, 'hann', []);
%%
% Compare
figure(1);
hold on
subplot(2,1,1); imagesc(worldX(1,:), flipud(worldY(:, 1)), U_tukey); title('Tukey'); colorbar;  
subplot(2,1,2); imagesc(worldX(1,:), flipud(worldY(:, 1)),U_hann); title('Hann'); colorbar; 
colormap(parula); 
sgtitle("median dx dy"); 


%% Prepare inputs
velocityU = avgVelocityField.u(1, :);  % Cell array of U fields
velocityV = avgVelocityField.v(1, :);  % Cell array of V fields
masks = {};         % Empty if no masks, or cell array of logical masks

% Run with Tukey window (Python default)
[worldX, worldY, U_tukey, ~] = merge_cameras_python_style_mean(...
    windowCenterCameras, velocityU, velocityV, masks, 'tukey', 0.5);

% Or try Hann window
[~, ~, U_hann, ~] = merge_cameras_python_style_mean(...
    windowCenterCameras, velocityU, velocityV, masks, 'hann', []);
%%
% Compare
figure(2);
hold on
subplot(2,1,1); imagesc(worldX(1,:), flipud(worldY(:, 1)), U_tukey); title('Tukey'); colorbar;  
subplot(2,1,2); imagesc(worldX(1,:), flipud(worldY(:, 1)),U_hann); title('Hann'); colorbar; 
colormap(parula); 
sgtitle("mean dx dy"); 

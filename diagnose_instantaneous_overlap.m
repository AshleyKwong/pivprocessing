%% diagnose_instantaneous_overlap.m
% Diagnostic script to visualise the overlap region between cameras 4 and 5
% for a single loop, across a selection of instantaneous frames and the mean field.
%
% Plots:
%   - Instantaneous U fields from cam4 and cam5 side by side with overlap boundaries marked
%   - Mean U field from cam4 and cam5 with overlap boundaries marked
%
% Required inputs:
%   savePath    - Root folder containing all loop subfolders. Must contain at least
%                 one loop subfolder with a valid windowCenterCameras_mm.mat file.
%                 Example: 'C:\Users\ak1u24\Downloads\PG_fixedcal\'
%
%   loopName    - Name of the loop folder to load instantaneous fields from.
%                 Must contain processedvelocityfields.mat in the root of the folder.
%                 Example: 'loop=06'
%
%   framesToPlot - Vector of frame indices to plot from the chosen loop.
%                  Frames must exist within processedvelocityfields.mat.
%                  Example: [1, 5, 10, 20, 50]
%
% Expected folder structure:
%   savePath/
%   ├── loop=00/
%   │     ├── windowCenterCameras_mm.mat   ← physical grid per camera (mm)
%   │     └── processedvelocityfields.mat  ← allCameras.u, allCameras.v {nFrames x nCams}
%   ├── loop=01/
%   │     └── ...
%   └── ...
%
% Notes:
%   - windowCenterCameras_mm.mat is loaded from the FIRST loop folder in which it is found
%   - processedvelocityfields.mat is loaded from loopName specifically
%   - The mean field in the BONUS plot is computed from ALL frames in loopName only,
%     not across all loops — for a full mean use averagedvelfields_uv_*.mat instead
%   - Cameras are assumed to have 5 entries in windowCenterCameras_mm; only 4 and 5 are plotted
%   - YDir is set to normal so wall (y=0) appears at the bottom of each plot
clear; 

%% USER OPTIONS
clc; close all;
savePath    = 'C:\Users\ak1u24\Downloads\PG_fixedcal\';
loopName    = 'loop=06';   % pick any single looploopName  = 'loop = 0';
framesToPlot = [1, 5, 10, 20, 50];   % pick a spread of frames

%% LOAD
physFile = '';
d = dir(savePath); d = d([d.isdir]);
for k = 1:length(d)
    c = fullfile(savePath, d(k).name, 'windowCenterCameras_mm.mat');
    if isfile(c); physFile = c; break; end
end
load(physFile, 'windowCenterCameras_mm');
if ~exist('data','var')
    data = load(fullfile(savePath, loopName, 'processedvelocityfields.mat'), 'allCameras');
end 
% Grid vectors
x1_4 = windowCenterCameras_mm.x1_mm{4};
x2_4 = windowCenterCameras_mm.x2_mm{4};
x1_5 = windowCenterCameras_mm.x1_mm{5};
x2_5 = windowCenterCameras_mm.x2_mm{5};

x_vec4 = x1_4(1,:);
x_vec5 = x1_5(1,:);
y_vec4 = x2_4(:,1);
y_vec5 = x2_5(:,1);

% Ensure y ascending for imagesc orientation
flip4 = y_vec4(1) > y_vec4(end);
flip5 = y_vec5(1) > y_vec5(end);
if flip4; y_vec4 = flipud(y_vec4); end
if flip5; y_vec5 = flipud(y_vec5); end

% Overlap boundaries
x_overlap_min = max(min(x_vec4), min(x_vec5));
x_overlap_max = min(max(x_vec4), max(x_vec5));
fprintf('Cam4/5 overlap: x = [%.2f, %.2f] mm\n', x_overlap_min, x_overlap_max);

% Shared colour limits — compute from mean field
u_all = [];
for fr = framesToPlot
    u_all = [u_all; data.allCameras.u{fr,4}(:); data.allCameras.u{fr,5}(:)];
end
u_all_clean = u_all(~isnan(u_all));
clims = [prctile(u_all_clean, 2), prctile(u_all_clean, 98)];
%% PLOT: One figure per frame
for fr = framesToPlot
    u4 = data.allCameras.u{fr, 4};
    u5 = data.allCameras.u{fr, 5};
    if flip4; u4 = flipud(u4); end
    if flip5; u5 = flipud(u5); end

    figure('Name', sprintf('Frame %d', fr), 'Position', [50 50 1600 700]);

    %--- Cam4 ---
    subplot(2,1,1);
    imagesc(x_vec4, y_vec4, u4);
    clim(clims); colorbar; colormap(turbo);
    hold on;
    % Mark overlap boundaries as vertical lines
    xline(x_overlap_min, 'w--', 'LineWidth', 2, 'Label', 'overlap start');
    xline(x_overlap_max, 'w-',  'LineWidth', 2, 'Label', 'overlap end');
    % Shade overlap region
    patch([x_overlap_min x_overlap_max x_overlap_max x_overlap_min], ...
          [y_vec4(1) y_vec4(1) y_vec4(end) y_vec4(end)], ...
          'white', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    title(sprintf('Cam4 — Frame %d  |  U [m/s]', fr));
    xlabel('x [mm]'); ylabel('y [mm]'); set(gca, 'YDir', 'normal'); 

    %--- Cam5 ---
    subplot(2,1,2);
    imagesc(x_vec5, y_vec5, u5);
    clim(clims); colorbar; colormap(turbo);
    hold on;
    xline(x_overlap_min, 'w--', 'LineWidth', 2, 'Label', 'overlap start');
    xline(x_overlap_max, 'w-',  'LineWidth', 2, 'Label', 'overlap end');
    patch([x_overlap_min x_overlap_max x_overlap_max x_overlap_min], ...
          [y_vec5(1) y_vec5(1) y_vec5(end) y_vec5(end)], ...
          'white', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    title(sprintf('Cam5 — Frame %d  |  U [m/s]', fr));
    xlabel('x [mm]'); ylabel('y [mm]'); set(gca, 'YDir', 'normal'); 

    sgtitle(sprintf('Instantaneous frame %d — overlap region shaded', fr));
end

%% BONUS: Also plot the mean field for each camera with overlap marked
meanU4 = mean(cat(3, data.allCameras.u{:,4}), 3, 'omitnan');
meanU5 = mean(cat(3, data.allCameras.u{:,5}), 3, 'omitnan');
if flip4; meanU4 = flipud(meanU4); end
if flip5; meanU5 = flipud(meanU5); end

figure('Name', 'Mean U per camera with overlap', 'Position', [50 50 1600 700]);
subplot(2,1,1);
imagesc(x_vec4, y_vec4, meanU4);
clim(clims); colorbar; colormap(turbo); hold on;
xline(x_overlap_min, 'w--', 'LineWidth', 2);
xline(x_overlap_max, 'w-',  'LineWidth', 2);
title('Cam4 — Mean U [m/s]');
xlabel('x [mm]'); ylabel('y [mm]'); set(gca, 'YDir', 'normal'); 

subplot(2,1,2);
imagesc(x_vec5, y_vec5, meanU5);
clim(clims); colorbar; colormap(turbo); hold on;
xline(x_overlap_min, 'w--', 'LineWidth', 2);
xline(x_overlap_max, 'w-',  'LineWidth', 2);
title('Cam5 — Mean U [m/s]');
xlabel('x [mm]'); ylabel('y [mm]'); set(gca, 'YDir', 'normal'); 
sgtitle('Per-camera mean U — overlap boundaries marked');
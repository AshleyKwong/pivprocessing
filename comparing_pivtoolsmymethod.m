% =========================================================================
% compare_pivtools_vs_mymethod.m
% Ashley Kwong
%
% Reads uncalibrated_piv from PIVtools output, applies your calibration
% pipeline (identical to PIV_stitching.m lines 1-233), saves:
%   - processedvelocityfields.mat  (per loop)
%   - windowCenterCameras_mm.mat   (per loop)
%   - averagedvelfields_uv_<N>.mat (global mean, your method)
% Then loads PIVtools calibrated_piv/Merged fields and compares.
% =========================================================================

clear; clc; close all;

set(groot, 'defaultAxesFontName',  'Cambria Math');
set(groot, 'defaultAxesFontSize',  20);
set(groot, 'defaultTextFontName',  'Cambria Math');
set(groot, 'defaultTextFontSize',  12);
set(groot, 'defaultTextColor',     'k');

% -------------------------------------------------------------------------
%% USER SETTINGS  (edit these)
% -------------------------------------------------------------------------
savePath   = 'E:\ProcessedPIV_case2fullpipe';   % root folder
cameraList = ["Cam1","Cam2","Cam3","Cam4","Cam5"];

SKIP_LOOP_PROCESSING = false;   % true  → skip per-loop calibration (already done)
FORCE_RECOMPUTE_AVG  = false;   % true  → recompute mean even if file exists

% Calibration
addpath(genpath('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\calib_tools\'));
load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\calib_tools\poly_2d2c\Step_01_Calibration\y235_aoan11_aoafn11\calib_AK.mat');

% Raw image frame sizes (per camera) – needed for dewarping grid
frame = { zeros(3528,5312,'single'), ...
    zeros(3472,5312,'single'), ...
    zeros(3512,5312,'single'), ...
    zeros(3536,5312,'single'), ...
    zeros(3520,5312,'single') };

dt             = 49.6e-6;   % [s] inter-frame time
final_pass_index = 3;       % PIVtools multi-pass index to use
skip_frames    = [];        % frame indices to skip (leave empty if none)

windowsize         = [16, 16];
overlap_percentage = 20;

% -------------------------------------------------------------------------
%% DEWARPING  (needed to build the world grid – identical to PIV_stitching)
% -------------------------------------------------------------------------
cx_img = est_dewarp_size_2d(frame, calib);
[~, cx_img, ~] = dewarp_image_2d(frame, calib, cx_img);

x1_min = cellfun(@(c) min(c.x1(1,:)),  cx_img, 'UniformOutput', true);
x1_max = cellfun(@(c) max(c.x1(1,:)),  cx_img, 'UniformOutput', true);
x2_min = cellfun(@(c) min(c.x2(:,1)),  cx_img, 'UniformOutput', true);
x2_max = cellfun(@(c) max(c.x2(:,1)),  cx_img, 'UniformOutput', true);

magFactor = zeros(2,5);
for i = 1:5
    magFactor(1,i) = calib.var{i}.mmppix(1);
    magFactor(2,i) = calib.var{i}.mmppix(2);
end
mean_del_x1 = mean(magFactor(1,:));
mean_del_x2 = mean(magFactor(2,:));

dx1 = mean_del_x1 * windowsize(2) * (1 - overlap_percentage/100);
dx2 = mean_del_x1 * windowsize(1) * (1 - overlap_percentage/100);

[worldx1, worldx2] = meshgrid(min(x1_min):dx1:max(x1_max), ...
    min(x2_max):-dx2:max(x2_min));
worldx1 = single(worldx1);
worldx2 = single(worldx2);

worldGrid.x1  = worldx1;
worldGrid.x2  = worldx2;
worldGrid.dx1 = dx1;
worldGrid.dx2 = dx2;

% -------------------------------------------------------------------------
%% DISCOVER LOOP FOLDERS
% -------------------------------------------------------------------------
d = dir(savePath);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));

loopPattern = '^loop\s*=\s*\d+$';
validLoops  = false(size(d));
for i = 1:length(d)
    validLoops(i) = ~isempty(regexpi(d(i).name, loopPattern));
end
totalLoops = d(validLoops);

if isempty(totalLoops)
    error('No "loop=XX" folders found in %s', savePath);
else
    fprintf('Found %d loop folders:\n', length(totalLoops));
    for i = 1:length(totalLoops)
        fprintf('  %s\n', totalLoops(i).name);
    end
end

% -------------------------------------------------------------------------
%% PER-LOOP PROCESSING  (calibrate uncalibrated_piv → processedvelocityfields)
% -------------------------------------------------------------------------
if SKIP_LOOP_PROCESSING
    fprintf('\n⚠  SKIP_LOOP_PROCESSING = true. Skipping per-loop calibration.\n');
else
    tic;
    for loopNo = 1:length(totalLoops)
        fprintf('\n==== Processing %s  (%d/%d) ====\n', ...
            totalLoops(loopNo).name, loopNo, length(totalLoops));

        saveFolder = fullfile(savePath, totalLoops(loopNo).name);
        physFile   = fullfile(saveFolder, 'windowCenterCameras_mm.mat');
        velFile    = fullfile(saveFolder, 'processedvelocityfields.mat');

        % --- Physical window centres (skip if already saved) -------------
        if isfile(physFile)
            load(physFile, 'windowCenterCameras_mm');
            fprintf('  ✓ Loaded existing windowCenterCameras_mm\n');
            skipCalibration = true;
        else
            skipCalibration = false;
            windowCenterCameras_mm.x1_mm       = cell(1,5);
            windowCenterCameras_mm.x2_mm       = cell(1,5);
            windowCenterCameras_mm.description = ...
                'Physical window centre locations in mm for each camera';
            windowCenterCameras_mm.cameras     = cameraList;
        end

        % --- Determine frame count from Cam1 of first image subfolder ----
        imgDirs     = dir(fullfile(saveFolder, 'uncalibrated_piv'));
        imgDirs     = imgDirs(~ismember({imgDirs.name},{'.','..'}));
        sample_dir  = dir(fullfile(saveFolder, 'uncalibrated_piv', ...
            imgDirs(1).name, 'Cam1', 'instantaneous', '*.mat'));
        % Pre-allocate across all cameras
        nTotalFrames = length(sample_dir) - 1;
        nValidFrames = nTotalFrames - length(skip_frames);

        allCameras.u = cell(nValidFrames, 5);
        allCameras.v = cell(nValidFrames, 5);

        for a = 1:length(cameraList)
            fprintf('  Camera %s\n', cameraList(a));

            base_dir = dir(fullfile(saveFolder, 'uncalibrated_piv', ...
                imgDirs(1).name, ['Cam' num2str(a)], ...
                'instantaneous', '*.mat'));

            [camData, windowCenterCameras_mm] = calibrateAndConvertVelocities( ...
                base_dir, calib, windowCenterCameras_mm, a, ...
                final_pass_index, skip_frames, dt, skipCalibration);

            % Slot this camera's frames into the full allCameras struct
            allCameras.u(:, a) = camData.u;
            allCameras.v(:, a) = camData.v;
        end

        % --- Save -----------------------------------------------------------------
        if ~skipCalibration
            save(physFile, 'windowCenterCameras_mm', '-v7.3');
            fprintf('  ✓ Saved windowCenterCameras_mm\n');
        end

        S = whos('allCameras');
        fprintf('  Memory: %.2f GB\n', S.bytes/1e9);
        save(velFile, 'allCameras', '-v7.3');
        fprintf('  ✓ Saved processedvelocityfields.mat\n');
        clear allCameras;


        % --- Save ---------------------------------------------------------
        if ~skipCalibration
            save(physFile, 'windowCenterCameras_mm', '-v7.3');
            fprintf('  ✓ Saved windowCenterCameras_mm\n');
        end

        % S = whos('allCameras');
        % fprintf('  Memory: %.2f GB\n', S.bytes/1e9);
        % save(velFile, 'allCameras', '-v7.3');
        % fprintf('  ✓ Saved processedvelocityfields.mat\n');
        % clear allCameras;
    end
    fprintf('\n✓ Total loop time: %.2f min\n', toc/60);
end

% -------------------------------------------------------------------------
%% AVERAGE  (your method)
% -------------------------------------------------------------------------
avgFileName = sprintf('averagedvelfields_uv_%d.mat', length(totalLoops)*150);
avgFilePath  = fullfile(savePath, avgFileName);

if isfile(avgFilePath) && ~FORCE_RECOMPUTE_AVG
    fprintf('\n✓ Loading existing mean: %s\n', avgFileName);
    avgVelocityField = load(avgFilePath, 'meanCameras').meanCameras;
else
    fprintf('\nComputing mean velocity fields...\n');
    avgVelocityField = averageVelocityFields(savePath, totalLoops);
end

% -------------------------------------------------------------------------
%% LOAD PIVTOOLS CALIBRATED/MERGED FIELDS
% -------------------------------------------------------------------------
fprintf('\nLoading PIVtools calibrated_piv/Merged fields...\n');

% Collect all loops that have a Merged folder
pivtools_inst_all = {};   % instantaneous: {frame, loopNo} cell of structs
pivtools_mean_U   = [];
pivtools_mean_V   = [];
pivtools_mean_count = 0;

for loopNo = 1:length(totalLoops)
    mergedDir = fullfile(savePath, totalLoops(loopNo).name, ...
        'calibrated_piv', '150', 'Merged', 'instantaneous/');

    if ~isfolder(mergedDir)
        fprintf('  ⚠  No Merged folder in %s – skipping\n', ...
            totalLoops(loopNo).name);
        continue;
    end

    mergedFiles = dir(fullfile(mergedDir, '*.mat'));
    if isempty(mergedFiles)
        fprintf('  ⚠  Merged folder empty in %s – skipping\n', ...
            totalLoops(loopNo).name);
        continue;
    end

    loopU = [];
    loopV = [];

    for f = 1:length(mergedFiles)
        pivData = load(fullfile(mergedFiles(f).folder, mergedFiles(f).name));

        % --- Flexible field detection ------------------------------------
        % PIVtools may store fields under different variable names.
        % Common candidates: piv_result, U, V, u, v, Ux, Uy, vel
        if isfield(pivData, 'piv_result')
            % Multi-pass struct (same as uncalibrated)
            U_frame = single(pivData.piv_result(end).ux);
            V_frame = single(pivData.piv_result(end).uy);
        elseif isfield(pivData, 'U') && isfield(pivData, 'V')
            U_frame = single(pivData.U);
            V_frame = single(pivData.V);
        elseif isfield(pivData, 'u') && isfield(pivData, 'v')
            U_frame = single(pivData.u);
            V_frame = single(pivData.v);
        elseif isfield(pivData, 'Ux') && isfield(pivData, 'Uy')
            U_frame = single(pivData.Ux);
            V_frame = single(pivData.Uy);
        else
            fprintf('  ⚠  Frame %d loop %s: unrecognised field names. Fields: %s\n', ...
                f, totalLoops(loopNo).name, strjoin(fieldnames(pivData)', ', '));
            continue;
        end

        % Store instantaneous frame
        pivtools_inst_all{end+1} = struct('U', U_frame, 'V', V_frame, ...
            'loop', loopNo, 'frame', f); %#ok<SAGROW>

        % Accumulate for mean
        if isempty(loopU)
            loopU = double(U_frame);
            loopV = double(V_frame);
        else
            loopU = loopU + double(U_frame);
            loopV = loopV + double(V_frame);
        end
    end

    % Running mean accumulation across loops
    if ~isempty(loopU)
        loopU = loopU / length(mergedFiles);
        loopV = loopV / length(mergedFiles);
        if isempty(pivtools_mean_U)
            pivtools_mean_U = loopU;
            pivtools_mean_V = loopV;
        else
            pivtools_mean_U = pivtools_mean_U + loopU;
            pivtools_mean_V = pivtools_mean_V + loopV;
        end
        pivtools_mean_count = pivtools_mean_count + 1;
        fprintf('  ✓ Loop %s: loaded %d merged frames\n', ...
            totalLoops(loopNo).name, length(mergedFiles));
    end
end

if pivtools_mean_count > 0
    pivtools_mean_U = pivtools_mean_U / pivtools_mean_count;
    pivtools_mean_V = pivtools_mean_V / pivtools_mean_count;
    fprintf('✓ PIVtools mean computed from %d loops, %d total frames\n', ...
        pivtools_mean_count, length(pivtools_inst_all));
else
    warning('No PIVtools Merged data found. Check folder structure.');
end

% -------------------------------------------------------------------------
%% COMPARE: MEAN FIELDS  (your method vs PIVtools)
% -------------------------------------------------------------------------
fprintf('\n=== Comparing Mean Velocity Fields ===\n');


velocityU = avgVelocityField.u(1, :);
velocityV = avgVelocityField.v(1, :);
masks = {};
% addpath for the functions
addpath(genpath('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV')); 
% Load physical locations for merging (from first loop folder)
load(fullfile(savePath, totalLoops(1).name, 'windowCenterCameras_mm.mat'), 'windowCenterCameras_mm');

[worldX, worldY, U_hann, ~] = merge_cameras_python_style_median(...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'hann', []);

figure('Position',[100 100 1600 900]);

% -- Your method: camera-level means (before merging) --------------------
subplot(3,1,1);
imagesc(worldX(1,:), worldY(:,1), U_hann);
colormap(gca, jet); colorbar; clim([0 30]);
title(sprintf('My Method')); 
axis image; set(gca,'YDir','normal');

% -- PIVtools merged mean -------------------------------------------------
subplot(3,1,2);
if ~isempty(pivtools_mean_U)
    imagesc(pivData.coordinates(final_pass_index).x(1,:),pivData.coordinates(final_pass_index).y (:, 1), pivtools_mean_U);
    colormap(gca, jet); colorbar; clim([0 30]);
    title('PIVtools \langle U \rangle  (calibrated\_piv/Merged)');
    axis image; set(gca,'YDir','normal');
    xlabel('X (px or mm)'); ylabel('Y (px or mm)');
else
    text(0.5, 0.5, 'No PIVtools data loaded', ...
        'HorizontalAlignment','center','Units','normalized');
    axis off;
end

% -- Difference (only if grids match) -------------------------------------
% Attempt a simple side-by-side difference plot using cam-level mean U.
% A proper grid-matched diff requires merging first (see below).
subplot(3,1,3);
myMeanU = U_hann;   % example: Cam1 only
ptMeanU_sub  = pivtools_mean_U;

if ~isempty(ptMeanU_sub) && isequal(size(myMeanU), size(ptMeanU_sub))
    diff_U = myMeanU - single(ptMeanU_sub);
    imagesc(diff_U);
    colormap(gca, redblue(256)); colorbar;
    title('Difference: MyMethod\_Cam1 − PIVtools (same grid only)');
    axis image; set(gca,'YDir','normal');
else
    text(0.5, 0.5, sprintf( ...
        'Grid mismatch: MyMethod Cam1 [%dx%d]  vs  PIVtools [%dx%d]\n(merge cameras first for a valid diff)', ...
        size(myMeanU,1), size(myMeanU,2), ...
        size(ptMeanU_sub,1),  size(ptMeanU_sub,2)), ...
        'HorizontalAlignment','center','Units','normalized');
    axis off;
end

sgtitle('Mean U: Your Calibration vs PIVtools calibrated\_piv/Merged', 'FontSize',14);
% =========================================================================
%% GRID-MATCHED COMPARISON  (interpolate both onto common grid, mean fields)
% =========================================================================
fprintf('\n=== Grid-matched comparison ===\n');

% --- 1. Your method: merge cameras onto world grid -----------------------
%   (worldX, worldY already exist from merge_cameras_python_style_median above)
%   U_my and V_my are your merged mean fields on (worldX, worldY)
[worldX_my, worldY_my, U_my, V_my] = merge_cameras_python_style_mean( ...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'hann', []);

% --- 2. PIVtools: extract its native coordinate grid ---------------------
%   pivData.coordinates should still be in workspace from the load loop.
%   Adjust field names below if PIVtools uses a different layout.
pt_x = pivData.coordinates(final_pass_index).x;   % [Ny × Nx] mm
pt_y = pivData.coordinates(final_pass_index).y;   % [Ny × Nx] mm
U_pt = single(pivtools_mean_U);                    % [Ny × Nx] m/s
V_pt = single(pivtools_mean_V);

fprintf('  My method grid  : [%d × %d],  dx ≈ %.3f mm\n', ...
    size(worldX_my,1), size(worldX_my,2), median(diff(worldX_my(1,:))));
fprintf('  PIVtools grid   : [%d × %d],  dx ≈ %.3f mm\n', ...
    size(pt_x,1), size(pt_x,2), median(diff(pt_x(1,:))));

% --- 3. Build a common target grid (coarser spacing → no upsampling) -----
dx_my = median(diff(worldX_my(1,:)));
dy_my = median(diff(worldY_my(:,1)));
dx_pt = median(diff(pt_x(1,:)));
dy_pt = median(diff(pt_y(:,1)));

dx_common = max(abs(dx_my), abs(dx_pt));   % coarser x-spacing
dy_common = max(abs(dy_my), abs(dy_pt));   % coarser y-spacing

% Overlap extent only
x_lo = max(min(worldX_my(:)), min(pt_x(:)));
x_hi = min(max(worldX_my(:)), max(pt_x(:)));
y_lo = max(min(worldY_my(:)), min(pt_y(:)));
y_hi = min(max(worldY_my(:)), max(pt_y(:)));

if x_lo >= x_hi || y_lo >= y_hi
    error('No spatial overlap between your method and PIVtools grids. Check coordinate units / axes.');
end

[Xc, Yc] = meshgrid(x_lo : dx_common : x_hi, ...
                     y_lo : dy_common : y_hi);
fprintf('  Common grid     : [%d × %d],  dx = %.3f mm,  dy = %.3f mm\n', ...
    size(Xc,1), size(Xc,2), dx_common, dy_common);

% --- 4. Valid-data masks (NaN where a field has no data) -----------------
%   Your merged field: NaN where Hann weight was zero
valid_my = ~isnan(U_my);
valid_pt = ~isnan(U_pt);

% --- 5. Interpolate both onto Xc/Yc (linear, no extrapolation = NaN) ----
F_my_U = scatteredInterpolant(double(worldX_my(valid_my)), double(worldY_my(valid_my)), ...
    double(U_my(valid_my)), 'linear', 'none');
F_my_V = scatteredInterpolant(double(worldX_my(valid_my)), double(worldY_my(valid_my)), ...
    double(V_my(valid_my)), 'linear', 'none');

F_pt_U = scatteredInterpolant(double(pt_x(valid_pt)), double(pt_y(valid_pt)), ...
    double(U_pt(valid_pt)), 'linear', 'none');
F_pt_V = scatteredInterpolant(double(pt_x(valid_pt)), double(pt_y(valid_pt)), ...
    double(V_pt(valid_pt)), 'linear', 'none');

U_my_c  = single(F_my_U(double(Xc), double(Yc)));
V_my_c  = single(F_my_V(double(Xc), double(Yc)));
U_pt_c  = single(F_pt_U(double(Xc), double(Yc)));
V_pt_c  = single(F_pt_V(double(Xc), double(Yc)));

% --- 6. Overlap mask: only points valid in BOTH fields -------------------
overlap = ~isnan(U_my_c) & ~isnan(U_pt_c);
fprintf('  Overlap coverage: %.1f %% of common grid\n', ...
    100 * sum(overlap(:)) / numel(overlap));

U_diff = nan(size(Xc), 'single');
V_diff = nan(size(Xc), 'single');
U_diff(overlap) = U_my_c(overlap) - U_pt_c(overlap);
V_diff(overlap) = V_my_c(overlap) - V_pt_c(overlap);

% --- 7. Summary statistics (overlap region only) -------------------------
fprintf('\n  --- U (streamwise) difference ---\n');
fprintf('  Mean  : %+.4f m/s\n',  mean(U_diff(overlap), 'omitnan'));
fprintf('  Std   : %.4f m/s\n',   std( U_diff(overlap), 'omitnan'));
fprintf('  RMSE  : %.4f m/s\n',   rms( U_diff(overlap)));
fprintf('  |max| : %.4f m/s\n',   max( abs(U_diff(overlap))));

fprintf('\n  --- V (wall-normal) difference ---\n');
fprintf('  Mean  : %+.4f m/s\n',  mean(V_diff(overlap), 'omitnan'));
fprintf('  Std   : %.4f m/s\n',   std( V_diff(overlap), 'omitnan'));
fprintf('  RMSE  : %.4f m/s\n',   rms( V_diff(overlap)));
fprintf('  |max| : %.4f m/s\n',   max( abs(V_diff(overlap))));

% --- 8. Plots ------------------------------------------------------------
clim_vel  = [0 30];
clim_diff = [-3  3];   % ← adjust to your expected discrepancy range

figure('Position', [50 50 1600 1100]);

% Row 1 – U
subplot(3,3,1);
imagesc(Xc(1,:), Yc(:,1), U_my_c);
set(gca,'YDir','normal'); axis image; colormap(gca,jet); colorbar; clim(clim_vel);
title('U  –  My method (on common grid)'); xlabel('X (mm)'); ylabel('Y (mm)');

subplot(3,3,2);
imagesc(Xc(1,:), Yc(:,1), U_pt_c);
set(gca,'YDir','normal'); axis image; colormap(gca,jet); colorbar; clim(clim_vel);
title('U  –  PIVtools (on common grid)'); xlabel('X (mm)');

subplot(3,3,3);
imagesc(Xc(1,:), Yc(:,1), U_diff);
set(gca,'YDir','normal'); axis image; colormap(gca,redblue(256)); colorbar; clim(clim_diff);
title('\Delta U = My − PIVtools [m/s]'); xlabel('X (mm)');

% Row 2 – V
subplot(3,3,4);
imagesc(Xc(1,:), Yc(:,1), V_my_c);
set(gca,'YDir','normal'); axis image; colormap(gca,jet); colorbar; clim([-5 5]);
title('V  –  My method'); xlabel('X (mm)'); ylabel('Y (mm)');

subplot(3,3,5);
imagesc(Xc(1,:), Yc(:,1), V_pt_c);
set(gca,'YDir','normal'); axis image; colormap(gca,jet); colorbar; clim([-5 5]);
title('V  –  PIVtools'); xlabel('X (mm)');

subplot(3,3,6);
imagesc(Xc(1,:), Yc(:,1), V_diff);
set(gca,'YDir','normal'); axis image; colormap(gca,redblue(256)); colorbar; clim(clim_diff);
title('\Delta V = My − PIVtools [m/s]'); xlabel('X (mm)');

% Row 3 – diagnostics
subplot(3,3,7);
histogram(U_diff(overlap), 80, 'FaceColor', [0.2 0.5 0.8]);
xlabel('\Delta U (m/s)'); ylabel('Count');
title(sprintf('\\Delta U  (\\mu=%.3f, \\sigma=%.3f m/s)', ...
    mean(U_diff(overlap),'omitnan'), std(U_diff(overlap),'omitnan')));
xline(0,'k--','LineWidth',1.5);

subplot(3,3,8);
histogram(V_diff(overlap), 80, 'FaceColor', [0.8 0.3 0.2]);
xlabel('\Delta V (m/s)'); ylabel('Count');
title(sprintf('\\Delta V  (\\mu=%.3f, \\sigma=%.3f m/s)', ...
    mean(V_diff(overlap),'omitnan'), std(V_diff(overlap),'omitnan')));
xline(0,'k--','LineWidth',1.5);

subplot(3,3,9);
imagesc(Xc(1,:), Yc(:,1), double(overlap));
set(gca,'YDir','normal'); axis image; colormap(gca, gray(2)); colorbar;
title('Overlap mask (white = valid in both)'); xlabel('X (mm)'); ylabel('Y (mm)');

sgtitle('Mean field comparison: My Method vs PIVtools  (common grid, linear interp)', ...
    'FontSize', 14);

% -------------------------------------------------------------------------
% =========================================================================
%% GRID-MATCHED COMPARISON  (U component only)
% =========================================================================
fprintf('\n=== Grid-matched comparison ===\n');

% --- 1. Your method: merge cameras onto world grid -----------------------
[worldX_my, worldY_my, U_my, ~] = merge_cameras_python_style_mean( ...
    windowCenterCameras_mm, velocityU, velocityV, masks, 'hann', []);

% Clear variables no longer needed
clear velocityU velocityV avgVelocityField;
clear worldx1 worldx2 worldGrid cx_img;
clear magFactor mean_del_x1 mean_del_x2;
clear x1_min x1_max x2_min x2_max;
clear U_hann V_hann;   % from median merge earlier in script

% --- 2. PIVtools: extract its native coordinate grid ---------------------
pt_x = pivData.coordinates(final_pass_index).x;   % [Ny × Nx] mm
pt_y = pivData.coordinates(final_pass_index).y;   % [Ny × Nx] mm
U_pt = single(pivtools_mean_U);

% Clear PIVtools accumulator variables
clear pivtools_mean_V pivtools_inst_all;
clear loopU loopV;

fprintf('  My method grid  : [%d × %d],  dx ≈ %.3f mm\n', ...
    size(worldX_my,1), size(worldX_my,2), median(diff(worldX_my(1,:))));
fprintf('  PIVtools grid   : [%d × %d],  dx ≈ %.3f mm\n', ...
    size(pt_x,1), size(pt_x,2), median(diff(pt_x(1,:))));

% --- 3. Build common target grid (coarser spacing → no upsampling) -------
dx_my = median(diff(worldX_my(1,:)));
dy_my = median(diff(worldY_my(:,1)));
dx_pt = median(diff(pt_x(1,:)));
dy_pt = median(diff(pt_y(:,1)));

dx_common = max(abs(dx_my), abs(dx_pt));
dy_common = max(abs(dy_my), abs(dy_pt));

x_lo = max(min(worldX_my(:)), min(pt_x(:)));
x_hi = min(max(worldX_my(:)), max(pt_x(:)));
y_lo = max(min(worldY_my(:)), min(pt_y(:)));
y_hi = min(max(worldY_my(:)), max(pt_y(:)));

if x_lo >= x_hi || y_lo >= y_hi
    error('No spatial overlap between grids. Check coordinate units / axes.');
end

[Xc, Yc] = meshgrid(x_lo : dx_common : x_hi, ...
                     y_lo : dy_common : y_hi);
fprintf('  Common grid     : [%d × %d],  dx = %.3f mm,  dy = %.3f mm\n', ...
    size(Xc,1), size(Xc,2), dx_common, dy_common);

% Clear grid construction intermediates
clear dx_my dy_my dx_pt dy_pt dx_common dy_common;
clear x_lo x_hi y_lo y_hi;

% --- 4. Valid-data masks --------------------------------------------------
valid_my = ~isnan(U_my);
valid_pt = ~isnan(U_pt);

% --- 5. Interpolate U onto common grid (linear, no extrapolation) --------
F_my_U = scatteredInterpolant(double(worldX_my(valid_my)), double(worldY_my(valid_my)), ...
    double(U_my(valid_my)), 'linear', 'none');
U_my_c = single(F_my_U(double(Xc), double(Yc)));
clear F_my_U worldX_my worldY_my U_my valid_my;

F_pt_U = scatteredInterpolant(double(pt_x(valid_pt)), double(pt_y(valid_pt)), ...
    double(U_pt(valid_pt)), 'linear', 'none');
U_pt_c = single(F_pt_U(double(Xc), double(Yc)));
clear F_pt_U pt_x pt_y U_pt valid_pt;
clear pivtools_mean_U;

% --- 6. Overlap mask and difference field --------------------------------
overlap = ~isnan(U_my_c) & ~isnan(U_pt_c);
fprintf('  Overlap coverage: %.1f %% of common grid\n', ...
    100 * sum(overlap(:)) / numel(overlap));

U_diff          = nan(size(Xc), 'single');
U_diff(overlap) = U_my_c(overlap) - U_pt_c(overlap);

% --- 7. Summary statistics -----------------------------------------------
fprintf('\n  --- U (streamwise) difference ---\n');
fprintf('  Mean  : %+.4f m/s\n', mean(U_diff(overlap), 'omitnan'));
fprintf('  Std   : %.4f m/s\n',  std( U_diff(overlap), 'omitnan'));
fprintf('  RMSE  : %.4f m/s\n',  rms( U_diff(overlap)));
fprintf('  |max| : %.4f m/s\n',  max( abs(U_diff(overlap))));

% --- 8. Plots (U only) ---------------------------------------------------
clim_vel  = [0 30];
clim_diff = [-3  3];   % adjust to expected discrepancy range

figure('Position', [50 50 1600 900]);

subplot(2,3,1);
imagesc(Xc(1,:), Yc(:,1), U_my_c);
set(gca,'YDir','normal'); axis image; colormap(gca,jet); colorbar; clim(clim_vel);
title('U  –  My method (common grid)'); xlabel('X (mm)'); ylabel('Y (mm)');

subplot(2,3,2);
imagesc(Xc(1,:), Yc(:,1), U_pt_c);
set(gca,'YDir','normal'); axis image; colormap(gca,jet); colorbar; clim(clim_vel);
title('U  –  PIVtools (common grid)'); xlabel('X (mm)');

subplot(2,3,3);
imagesc(Xc(1,:), Yc(:,1), U_diff);
set(gca,'YDir','normal'); axis image; colormap(gca,redblue(256)); colorbar; clim(clim_diff);
title('\Delta U = My − PIVtools [m/s]'); xlabel('X (mm)');

subplot(2,3,4:5);
histogram(U_diff(overlap), 80, 'FaceColor', [0.2 0.5 0.8]);
xlabel('\Delta U (m/s)'); ylabel('Count');
title(sprintf('\\Delta U distribution  (\\mu = %.3f m/s,  \\sigma = %.3f m/s)', ...
    mean(U_diff(overlap),'omitnan'), std(U_diff(overlap),'omitnan')));
xline(0, 'k--', 'LineWidth', 1.5);

subplot(2,3,6);
imagesc(Xc(1,:), Yc(:,1), double(overlap));
set(gca,'YDir','normal'); axis image; colormap(gca, gray(2)); colorbar;
title('Overlap mask'); xlabel('X (mm)'); ylabel('Y (mm)');

sgtitle('Mean U: My Method vs PIVtools  (common grid, linear interp)', 'FontSize', 14);

% Final cleanup
clear U_diff overlap Xc Yc;

%% COMPARE: INSTANTANEOUS FIELDS  (first available frame)
% -------------------------------------------------------------------------
if ~isempty(pivtools_inst_all)
    fprintf('\n=== Comparing Instantaneous Fields (frame 1, loop 1) ===\n');

    % Load your method frame 1 from loop 1
    myVelFile = fullfile(savePath, totalLoops(1).name, 'processedvelocityfields.mat');
    myData    = load(myVelFile, 'allCameras');

    figure('Position',[100 100 1400 700]);
    for cam = 1:5
        subplot(2, 5, cam);
        imagesc(myData.allCameras.u{1, cam});
        colormap(gca, jet); colorbar; clim([-5 35]);
        title(sprintf('My Cam%d  f=1', cam));
        axis image; set(gca,'YDir','normal');
    end

    % PIVtools – first loaded frame
    ptFrame = pivtools_inst_all{1};
    subplot(2, 5, 6:10);
    imagesc(ptFrame.U);
    colormap(gca, jet); colorbar; clim([-5 35]);
    title(sprintf('PIVtools Merged  loop=%d  f=%d', ptFrame.loop, ptFrame.frame));
    axis image; set(gca,'YDir','normal');
    xlabel('X'); ylabel('Y');

    sgtitle('Instantaneous U: Your Method (per camera) vs PIVtools Merged','FontSize',14);
    clear myData;
end
%%
% =========================================================================
% Example:
% Extract profiles with all integral thicknesses
pt_x = pivData.coordinates(final_pass_index).x;   % [Ny × Nx] mm
pt_y = pivData.coordinates(final_pass_index).y;   % [Ny × Nx] mm
U_pt = single(ptMeanU_sub);                    % [Ny × Nx] m/s
figure(); imagesc(pt_x(1,:), pt_y(:,1), U_pt); axis image; 
yline(2, "r")
set(gca,'YDir','normal');




%% LOCAL FUNCTION: averageVelocityFields
% =========================================================================
function meanCameras = averageVelocityFields(savePath, totalLoops)
tic;
fprintf('\n=== Averaging Velocity Fields ===\n');

firstPath = fullfile(savePath, totalLoops(1).name, 'processedvelocityfields.mat');
firstData = load(firstPath, 'allCameras');
[nFrames, nCameras] = size(firstData.allCameras.u);
sumU = firstData.allCameras.u;
sumV = firstData.allCameras.v;
clear firstData;
fprintf('  Initialised from loop 1/%d\n', length(totalLoops));

for k = 2:length(totalLoops)
    nextPath = fullfile(savePath, totalLoops(k).name, 'processedvelocityfields.mat');
    nextData = load(nextPath, 'allCameras');
    for fr = 1:nFrames
        for cam = 1:nCameras
            sumU{fr,cam} = sumU{fr,cam} + nextData.allCameras.u{fr,cam};
            sumV{fr,cam} = sumV{fr,cam} + nextData.allCameras.v{fr,cam};
        end
    end
    clear nextData;
    fprintf('  Added loop %d/%d\n', k, length(totalLoops));
end

nL = length(totalLoops);
meanU = cell(1, nCameras);
meanV = cell(1, nCameras);
for cam = 1:nCameras
    stackU = cat(3, sumU{:,cam}) / nL;
    stackV = cat(3, sumV{:,cam}) / nL;
    meanU{cam} = mean(stackU, 3);
    meanV{cam} = mean(stackV, 3);
end

meanCameras.u = meanU;
meanCameras.v = meanV;

S = whos('meanCameras');
fprintf('✓ Averaged %d loops × %d frames  (%.2f GB) in %.1f s\n', ...
    nL, nFrames, S.bytes/1e9, toc);

tag      = sprintf('averagedvelfields_uv_%d', nL*150);
filename = fullfile(savePath, [tag '.mat']);
save(filename, 'meanCameras', '-v7.3');
fprintf('✓ Saved %s\n', tag);
end

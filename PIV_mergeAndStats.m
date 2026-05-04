% Ashley Kwong

% PIV_mergeInstantaneous.m
% Loads processed (unmerged) velocity fields per loop, merges all 5 cameras
% for each instantaneous frame, saves merged frames as a cell array.
% Also computes and saves the merged mean velocity field at the end.

clear; clc; close all;
addpath(genpath('/iridisfs/scratch/ak1u24/calib_tools'));

%% ======== USER OPTIONS ================================================
savePath   = '/iridisfs/scratch/ak1u24/ProcessedPIV_fullpipeline';
calibFile  = '/iridisfs/scratch/ak1u24/calib_allcamerascropped_floorpixelspecify_20260303_211711.mat';
cameraList = ["Cam1","Cam2","Cam3","Cam4","Cam5"];
dt         = 49.6e-6;
TEST_LOOP  = [];   % set to [] to process all loops
% ======================================================================

%% OUTPUT DIRECTORY
tstamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(savePath, sprintf('merge_instantaneous_%s', tstamp));
mkdir(outDir);
fprintf('Output directory: %s\n', outDir);

%% LOAD CALIBRATION
load(calibFile, 'calib');
fprintf('✓ Loaded calibration\n');

%% DISCOVER LOOP FOLDERS
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
    error('No loop folders found in %s', savePath);
end
fprintf('Found %d loop folders.\n', length(totalLoops));

% Restrict to test loop if specified
if ~isempty(TEST_LOOP)
    totalLoops = totalLoops(TEST_LOOP);
    fprintf('TEST_LOOP = %d: processing only %s\n', TEST_LOOP, totalLoops(1).name);
end

%% LOAD WINDOW CENTRES (from first available loop folder)
fprintf('\n=== Loading window centres ===\n');
physicalLocationFile = '';
for k = 1:length(totalLoops)
    candidate = fullfile(savePath, totalLoops(k).name, 'windowCenterCameras_mm.mat');
    if isfile(candidate)
        physicalLocationFile = candidate;
        break;
    end
end
if isempty(physicalLocationFile)
    error('No windowCenterCameras_mm.mat found in any loop folder.');
end
load(physicalLocationFile, 'windowCenterCameras_mm');
fprintf('✓ Loaded window centres from: %s\n', physicalLocationFile);
for a = 1:length(cameraList)
    x1 = windowCenterCameras_mm.x1_mm{a};
    x2 = windowCenterCameras_mm.x2_mm{a};
    fprintf('  Cam%d: x1=[%.2f, %.2f] mm | x2=[%.2f, %.2f] mm\n', a, ...
        min(x1(:)), max(x1(:)), min(x2(:)), max(x2(:)));
end

masks = {}; % empty = no masking

%% PROCESS EACH LOOP
for loopNo = 1:length(totalLoops)
    loopName   = totalLoops(loopNo).name;
    loopFolder = fullfile(savePath, loopName);
    fprintf('\n--- Processing %s (%d/%d) ---\n', loopName, loopNo, length(totalLoops));

    % Load processed (unmerged) velocity fields
    velFile = fullfile(loopFolder, 'processedvelocityfields.mat');
    if ~isfile(velFile)
        warning('processedvelocityfields.mat not found in %s — skipping.', loopName);
        continue;
    end
    data     = load(velFile, 'allCameras');
    nFrames  = size(data.allCameras.u, 1);
    nCams    = size(data.allCameras.u, 2);
    fprintf('  Loaded: %d frames x %d cameras\n', nFrames, nCams);

    % Preallocate output cell arrays
    mergedFrames.U = cell(nFrames, 1);  % merged instantaneous U per frame
    mergedFrames.V = cell(nFrames, 1);  % merged instantaneous V per frame
    worldX_stored  = [];
    worldY_stored  = [];

    % Merge each frame
    for fr = 1:nFrames
        % Extract this frame's data across all cameras into {1x5} cells
        frameU = data.allCameras.u(fr, :);  % {1x5}
        frameV = data.allCameras.v(fr, :);  % {1x5}

        [worldX, worldY, U_merged, V_merged] = merge_cameras_python_style_mean(...
            windowCenterCameras_mm, frameU, frameV, masks, 'hann', []);

        mergedFrames.U{fr} = U_merged;
        mergedFrames.V{fr} = V_merged;

        % Store grid from first frame (same for all)
        if fr == 1
            worldX_stored = worldX;
            worldY_stored = worldY;
            fprintf('  Grid: [%d x %d], X:[%.1f, %.1f] mm, Y:[%.1f, %.1f] mm\n', ...
                size(U_merged,1), size(U_merged,2), ...
                min(worldX(:)), max(worldX(:)), ...
                min(worldY(:)), max(worldY(:)));
        end

        if mod(fr, 25) == 0
            fprintf('  Frame %d/%d done\n', fr, nFrames);
        end
    end
    clear data;

    % NaN diagnostic on last frame as spot check
    nan_pct = 100 * sum(isnan(mergedFrames.U{end}(:))) / numel(mergedFrames.U{end});
    fprintf('  NaNs in last merged frame: %.1f%%\n', nan_pct);

    % Save merged instantaneous fields for this loop
    mergedFile = fullfile(loopFolder, sprintf('mergedvelocityfields_%s.mat', tstamp));
    worldX = worldX_stored;
    worldY = worldY_stored;
    save(mergedFile, 'mergedFrames', 'worldX', 'worldY', '-v7.3');
    fprintf('✓ Saved merged frames → %s\n', mergedFile);
    clear mergedFrames;

end % loopNo

%% MERGED MEAN VELOCITY (from pre-averaged per-camera fields)
fprintf('\n=== Computing merged mean velocity ===\n');
avgFileName = sprintf('averagedvelfields_uv_%d.mat', length(d(validLoops))*150);
avgFilePath  = fullfile(savePath, avgFileName);

if ~isfile(avgFilePath)
    warning('Averaged file not found: %s\nSkipping mean merge.', avgFilePath);
else
    avgVelocityField = load(avgFilePath, 'meanCameras').meanCameras;
    fprintf('✓ Loaded averaged fields from: %s\n', avgFileName);

    [worldX, worldY, U_hann_mean, V_hann_mean] = merge_cameras_python_style_mean(...
        windowCenterCameras_mm, avgVelocityField.u, avgVelocityField.v, masks, 'hann', []);

    fprintf('✓ Merged mean field: [%d x %d]\n', size(U_hann_mean,1), size(U_hann_mean,2));
    fprintf('  NaNs: %.1f%% | U range: [%.3f, %.3f] m/s\n', ...
        100*sum(isnan(U_hann_mean(:)))/numel(U_hann_mean), ...
        min(U_hann_mean(:),[],'omitnan'), max(U_hann_mean(:),[],'omitnan'));

    mergedMeanFile = fullfile(outDir, sprintf('merged_meanUV_%dloops_%s.mat', ...
        length(d(validLoops)), tstamp));
    save(mergedMeanFile, 'worldX', 'worldY', 'U_hann_mean', 'V_hann_mean', '-v7.3');
    fprintf('✓ Saved merged mean UV → %s\n', mergedMeanFile);
end

fprintf('\n=== Done: %s ===\n', tstamp);
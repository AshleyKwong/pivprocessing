%% Compute per-camera mean PIV velocity fields (PIVtools calibrated data)
% Runs over all discovered loop folders and saves mean U/V per camera.
% Intended to be run on Iridis 6 via SLURM.
clear; clc; close all;

%%
savePath   = '/iridisfs/scratch/ak1u24/ProcessedPIV_case2fullpipe';  % ← UPDATE if different
cameraList = ["Cam1","Cam2","Cam3","Cam4","Cam5"];

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
    error('No "loop=XX" folders found in %s', savePath);
else
    fprintf('Found %d loop folders:\n', length(totalLoops));
    for i = 1:length(totalLoops)
        fprintf('  %s\n', totalLoops(i).name);
    end
end

mean_U_per_cam = cell(1, numel(cameraList));
mean_V_per_cam = cell(1, numel(cameraList));
coords_per_cam = cell(1, numel(cameraList));

%% LOAD PIVTOOLS CALIBRATED/MERGED FIELDS
for cameraNo = 1:numel(cameraList)
    pivtools_mean_U     = [];
    pivtools_mean_V     = [];
    pivtools_mean_count = 0;

    for loopNo = 1:length(totalLoops)   % ← was hardcoded to 3; now runs all loops
        loopU = [];
        loopV = [];
        fprintf('\n%s\n', fullfile(savePath, totalLoops(loopNo).name, ...
            'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous'));
        camDir = fullfile(savePath, totalLoops(loopNo).name, ...
            'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous/');
        if ~isfolder(camDir)
            fprintf('  ⚠  No %s folder in %s – skipping\n', ...
                cameraList(cameraNo), totalLoops(loopNo).name);
            continue;
        end
        mergedFiles = dir(fullfile(camDir, '*.mat'));
        mergedFiles = mergedFiles(~strcmp({mergedFiles.name}, 'coordinates.mat'));
        if isempty(mergedFiles)
            fprintf('  ⚠  %s folder empty in %s – skipping\n', ...
                cameraList(cameraNo), totalLoops(loopNo).name);
            continue;
        end
        coordinates_file = fullfile(savePath, totalLoops(loopNo).name, ...
            'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous/coordinates.mat');
        if ~isfile(coordinates_file)
            fprintf('   ⚠  No coordinates.mat for %s in %s\n', ...
                cameraList(cameraNo), totalLoops(loopNo).name);
            continue;
        else
            tmp = load(coordinates_file).coordinates;
            if iscell(tmp)
                coords_per_cam{cameraNo} = tmp{end};
            else
                coords_per_cam{cameraNo} = tmp;
            end
        end

        for f = 1:length(mergedFiles)
            pivData = load(fullfile(mergedFiles(f).folder, mergedFiles(f).name));
            if isfield(pivData, 'piv_result')
                if iscell(pivData.piv_result)
                    U_frame = single(pivData.piv_result{end}.ux);
                    V_frame = single(pivData.piv_result{end}.uy);
                else
                    U_frame = single(pivData.piv_result(end).ux);
                    V_frame = single(pivData.piv_result(end).uy);
                end
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
            if isempty(loopU)
                loopU = double(U_frame);
                loopV = double(V_frame);
            else
                loopU = loopU + double(U_frame);
                loopV = loopV + double(V_frame);
            end
        end

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
    end % loopNo

    if pivtools_mean_count > 0
        mean_U_per_cam{cameraNo} = pivtools_mean_U / pivtools_mean_count;
        mean_V_per_cam{cameraNo} = pivtools_mean_V / pivtools_mean_count;
        fprintf('✓ PIVtools mean computed for %s, %d total frames\n', ...
            cameraList(cameraNo), length(mergedFiles));
    else
        warning('No PIVtools Merged data found. Check folder structure.');
    end
end % cameraNo

%% SAVE MEAN VELOCITY FIELDS
% -v7.3 uses HDF5 format — required for variables > 2 GB
outputFile = fullfile(savePath, 'pivtools_mean_UV_allcams.mat');
save(outputFile, 'mean_U_per_cam', 'mean_V_per_cam', 'coords_per_cam', 'cameraList', '-v7.3');
fprintf('\n✓ Saved mean U/V fields to:\n  %s\n', outputFile);

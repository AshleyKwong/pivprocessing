% PIV_applyCorrections_batch.m
% Runner script called by SLURM via: matlab -batch "PIV_applyCorrections_batch"
% Applies floor shift (coordinates.mat) and velocity rotation (frame .mat files)
% for ALL loops and ALL cameras in the save path.
%
% !! USER: fill in p_C_all from your offline PIV_detectFloor analysis !!

%% ====================================================================
%  USER SETTINGS  — edit these before submitting the job
% =====================================================================
savePath   = 'E:\ProcessedPIV_case2fullpipe'; %'/iridisfs/scratch/ak1u24/case1_fullpipelinetest';
cameraList = ["Cam1", "Cam2", "Cam3", "Cam4", "Cam5"];
subPath    = fullfile('calibrated_piv', '150');   % path segment inside each loop folder

% p_C = [slope, intercept] per camera — from your offline PIV_detectFloor run
p_C_all = {
    [0.0033,-12.3862],   % Cam1
    [0.0028,-6.2082],   % Cam2
    [0.0028,-14.5398],   % Cam3
    [0.0028,-21.7908],   % Cam4
    [0.0041,-3.5662],   % Cam5
};
% =====================================================================

%% Discover all loop folders
d           = dir(savePath);
d           = d([d.isdir]);
d           = d(~ismember({d.name}, {'.','..'}));
loopPattern = '^loop\s*=\s*\d+$';
validLoops  = arrayfun(@(x) ~isempty(regexpi(x.name, loopPattern)), d);
loopFolders = d(validLoops);

if isempty(loopFolders)
    error('PIV_applyCorrections_batch:noLoops', ...
          'No loop folders found in: %s', savePath);
end

fprintf('Found %d loop folders.\n', numel(loopFolders));

%% Apply corrections: loop over every loop × camera
for loopNo = 1%:numel(loopFolders) % running test on loop = 01

    loopName = loopFolders(loopNo).name;
    fprintf('\n========== %s ==========\n', loopName);

    for camNo = 1:numel(cameraList)

        instDir   = fullfile(savePath, loopName, subPath, ...
                             cameraList(camNo), 'instantaneous');
        coordFile = fullfile(instDir, 'coordinates.mat');

        fprintf('\n  Camera: %s\n', cameraList(camNo));

        % --- 1. Shift coordinates.mat ---
        if ~isfile(coordFile)
            warning('Missing coordinates.mat, skipping: %s', coordFile);
        else
            PIV_applyFloorShift(coordFile, p_C_all{camNo});
        end

        % --- 2. Rotate velocity vectors in all frame .mat files ---
        if ~isfolder(instDir)
            warning('Missing instantaneous folder, skipping: %s', instDir);
        else
            PIV_applyVelocityRotation(instDir, p_C_all{camNo});
        end

    end % camNo
end % loopNo

fprintf('\n[DONE] All floor shifts and velocity rotations applied.\n');

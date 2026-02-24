% PIV_applyCorrections_batch.m
% Runner script called by SLURM via: matlab -batch "PIV_applyCorrections_batch"
% Applies floor shift (coordinates.mat) and velocity rotation (frame .mat files)
% for ALL loops and ALL cameras in the save path.
%
% !! USER: fill in p_C_all from your offline PIV_detectFloor analysis !!

%% ====================================================================
%  USER SETTINGS  — edit these before submitting the job
% =====================================================================
savePath   = '/iridisfs/scratch/ak1u24/ProcessedPIV_fullpipeline'; %'/iridisfs/scratch/ak1u24/case1_fullpipelinetest';
cameraList = ["Cam1", "Cam2", "Cam3", "Cam4", "Cam5"];
subPath    = fullfile('calibrated_piv', '150');   % path segment inside each loop folder

% p_C = [slope, intercept] per camera — from your offline PIV_detectFloor run
p_C_all = {
    [0.0041,-12.6361],   % Cam1
    [0.0025,-6.1333],   % Cam2
    [0.0028,-14.5699],   % Cam3
    [0.0028,-21.7745],   % Cam4
    [0.0041,-3.5699],   % Cam5
};
% =====================================================================

%% Start parallel pool using all cores SLURM gave us
nCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nCores) || nCores < 1
    nCores = feature('numcores');   % fallback if env var not set
end
fprintf('Starting parpool with %d workers...\n', nCores);
if isempty(gcp('nocreate'))
    parpool('local', nCores);
end

%% Discover loop folders
d           = dir(savePath);
d           = d([d.isdir]);
d           = d(~ismember({d.name},{'.','..'}));
loopPattern = '^loop\s*=\s*\d+$';
validLoops  = arrayfun(@(x) ~isempty(regexpi(x.name, loopPattern)), d);
loopFolders = d(validLoops);

if isempty(loopFolders)
    error('No loop folders found in: %s', savePath);
end
fprintf('Found %d loop folders.\n', numel(loopFolders));

%% Apply corrections
% Cameras are serial here — parallelism is inside PIV_applyVelocityRotation
% (parfor over files). Nested parfor is not supported in MATLAB.
for loopNo = 1%:numel(loopFolders)
    fprintf("\n ASHLEY TESTING LOOP =01 BATCH !!! NEED TO FIX FOR REAL PROCESSING. \n"); 
    loopName = loopFolders(loopNo).name;
    fprintf('\n========== %s ==========\n', loopName);

    for camNo = 1:numel(cameraList)

        instDir   = fullfile(savePath, loopName, subPath, ...
                             cameraList(camNo), 'instantaneous');
        coordFile = fullfile(instDir, 'coordinates.mat');
        fprintf('\n  Camera: %s\n', cameraList(camNo));

        if ~isfile(coordFile)
            warning('Missing coordinates.mat, skipping: %s', coordFile);
        else
            PIV_applyFloorShift(coordFile, p_C_all{camNo});
        end

        if ~isfolder(instDir)
            warning('Missing instantaneous folder, skipping: %s', instDir);
        else
            PIV_applyVelocityRotation(instDir, p_C_all{camNo});
        end

    end
end

fprintf('\n[DONE] All corrections applied.\n');
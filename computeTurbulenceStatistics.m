function turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField, varargin)
% COMPUTETURBULENCESTATISTICS Calculate velocity fluctuations and turbulence statistics
%
% Syntax:
%   turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField)
%   turbStats = computeTurbulenceStatistics(..., 'Name', Value)
%
% Inputs:
%   savePath          - Base path where loop folders are located
%   totalLoops        - Structure array of loop directories (from dir())
%   avgVelocityField  - Structure with fields .u and .v (1x5 cell arrays of mean velocities)
%
% Optional Name-Value Pairs:
%   'SaveIndividualFluctuations' - Save all frames for PDF analysis (default: true)
%   'ForceRecompute'             - Recompute even if files exist (default: false)
%   'Verbose'                    - Print progress messages (default: true)
%   'ProcessFramesInChunks'      - Process frames in chunks to reduce peak RAM (default: false)
%   'ChunkSize'                  - Number of frames per chunk (default: 30)
%   'VelocityThreshold'          - [min, max] valid velocity in m/s (default: [0, 50])
%   'FluctuationThreshold'       - Max reasonable fluctuation in m/s (default: 10)
%   'MinValidFraction'           - Minimum fraction of valid frames per point (default: 0.7)
%
% Outputs:
%   turbStats - Structure containing:
%       .u_variance{cam}          - Variance field <u'^2> in (m/s)^2
%       .v_variance{cam}          - Variance field <v'^2> in (m/s)^2
%       .u_rms{cam}               - RMS field sqrt(<u'^2>) in m/s
%       .v_rms{cam}               - RMS field sqrt(<v'^2>) in m/s
%       .u_rms_mean_scalar(cam)   - Spatial-average RMS in m/s
%       .v_rms_mean_scalar(cam)   - Spatial-average RMS in m/s
%       .u_variance_mean_scalar   - Spatial-average variance in (m/s)^2
%
% Notes:
%   - Fluctuations: u' = u - U_mean (m/s)
%   - Variance: <u'^2> (m/s)^2
%   - RMS: sqrt(<u'^2>) (m/s)
%   - Turbulence Intensity: <u'^2> / U_inf^2 (dimensionless) - compute separately
%
% File Structure Created:
%   savePath/
%   ├── loop=X/
%   │   └── vel_fluctuations/
%   │       └── fluctuations_all_frames.mat (~2.58 GB per loop)
%   ├── loop_averaged_variance_statistics.mat
%   └── turbulence_statistics.mat
%
% Example:
%   turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField);
%   
%   turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField, ...
%       'VelocityThreshold', [0, 50], ...
%       'FluctuationThreshold', 10, ...
%       'MinValidFraction', 0.7);
%
% Author: Ashley Kwong
% Date: 02/16/2026

%% Parse inputs
p = inputParser;
addRequired(p, 'savePath', @ischar);
addRequired(p, 'totalLoops', @isstruct);
addRequired(p, 'avgVelocityField', @isstruct);
addParameter(p, 'SaveIndividualFluctuations', true, @islogical);
addParameter(p, 'ForceRecompute', false, @islogical);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'ProcessFramesInChunks', false, @islogical);
addParameter(p, 'ChunkSize', 30, @isnumeric);
addParameter(p, 'VelocityThreshold', [0, 50], @isnumeric);
addParameter(p, 'FluctuationThreshold', 10, @isnumeric);
addParameter(p, 'MinValidFraction', 0.7, @isnumeric);

parse(p, savePath, totalLoops, avgVelocityField, varargin{:});

SAVE_INDIVIDUAL_FLUCTUATIONS = p.Results.SaveIndividualFluctuations;
FORCE_RECOMPUTE = p.Results.ForceRecompute;
VERBOSE = p.Results.Verbose;
PROCESS_IN_CHUNKS = p.Results.ProcessFramesInChunks;
CHUNK_SIZE = p.Results.ChunkSize;
VEL_MIN = p.Results.VelocityThreshold(1);
VEL_MAX = p.Results.VelocityThreshold(2);
FLUCT_MAX = p.Results.FluctuationThreshold;
MIN_VALID_FRAC = p.Results.MinValidFraction;

%% Initialize
if VERBOSE
    fprintf('\n=== Computing Turbulence Fluctuations ===\n');
    fprintf('Estimated storage: ~%.1f GB total (~2.58 GB per loop)\n', 2.58 * length(totalLoops));
    fprintf('Quality control thresholds:\n');
    fprintf('  Velocity range: [%.1f, %.1f] m/s\n', VEL_MIN, VEL_MAX);
    fprintf('  Max fluctuation: %.1f m/s\n', FLUCT_MAX);
    fprintf('  Min valid fraction: %.1f%%\n', MIN_VALID_FRAC * 100);
    tic;
end

% Get mean velocity for each camera
meanU = avgVelocityField.u;
meanV = avgVelocityField.v;
nLoops = length(totalLoops);
nCameras = length(meanU);

% Create validity masks for mean fields
meanMask = cell(1, nCameras);
for cam = 1:nCameras
    % Mask where mean velocity is reasonable
    validU = ~isnan(meanU{cam}) & meanU{cam} >= VEL_MIN & meanU{cam} <= VEL_MAX;
    validV = ~isnan(meanV{cam}) & meanV{cam} >= VEL_MIN & meanV{cam} <= VEL_MAX;
    meanMask{cam} = validU & validV;
    
    if VERBOSE
        valid_pct = 100 * sum(meanMask{cam}(:)) / numel(meanMask{cam});
        fprintf('Camera %d mean field: %.1f%% valid\n', cam, valid_pct);
    end
end

% Preallocate structure for loop-averaged statistics
loopStats = struct();
loopStats.u_variance_mean = cell(nLoops, nCameras);    % Variance <u'^2> in (m/s)^2
loopStats.v_variance_mean = cell(nLoops, nCameras);    % Variance <v'^2> in (m/s)^2
loopStats.valid_fraction = cell(nLoops, nCameras);     % Track data quality

%% Process each loop
for loopNo = 1:nLoops
    if VERBOSE
        fprintf('\n Processing fluctuations for %s (%d/%d)\n', ...
            totalLoops(loopNo).name, loopNo, nLoops);
    end
    
    saveFolder = fullfile(savePath, totalLoops(loopNo).name);
    fluctFolder = fullfile(saveFolder, 'vel_fluctuations');
    if ~exist(fluctFolder, 'dir')
        mkdir(fluctFolder);
    end
    
    fluctFile = fullfile(fluctFolder, 'fluctuations_all_frames.mat');
    
    if isfile(fluctFile) && ~FORCE_RECOMPUTE
        if VERBOSE
            fprintf('  ⚠ Fluctuations already exist, loading...\n');
        end
        load(fluctFile, 'fluctuations');
        [nFrames, ~] = size(fluctuations.u_prime);
    else
        % Load instantaneous velocity fields
        velocityFile = fullfile(saveFolder, 'processedvelocityfields.mat');
        if VERBOSE
            fprintf('  Loading velocity fields... ');
        end
        load(velocityFile, 'allCameras');
        if VERBOSE
            fprintf('✓\n');
        end
        
        [nFrames, ~] = size(allCameras.u);
        
        % Preallocate fluctuation storage
        fluctuations = struct();
        fluctuations.u_prime = cell(nFrames, nCameras);
        fluctuations.v_prime = cell(nFrames, nCameras);
        fluctuations.valid_mask = cell(nFrames, nCameras);
        fluctuations.loop_name = totalLoops(loopNo).name;
        fluctuations.n_frames = nFrames;
        fluctuations.n_cameras = nCameras;
        fluctuations.description = 'Velocity fluctuations: u'' = u - U_mean, v'' = v - V_mean (m/s)';
        
        % Process each camera
        for cam = 1:nCameras
            if VERBOSE
                fprintf('  Camera %d: Computing %d fluctuation fields with QC... ', cam, nFrames);
            end
            
            rejected_count = 0;
            
            % Process frames with quality control
            for frame = 1:nFrames
                u_inst = allCameras.u{frame, cam};
                v_inst = allCameras.v{frame, cam};
                
                % Quality control mask for this frame
                instMask = ~isnan(u_inst) & ~isnan(v_inst) & ...
                          u_inst >= VEL_MIN & u_inst <= VEL_MAX & ...
                          v_inst >= VEL_MIN & v_inst <= VEL_MAX;
                
                % Calculate fluctuations: u' = u - U_mean (m/s)
                u_prime = single(u_inst - meanU{cam});
                v_prime = single(v_inst - meanV{cam});
                
                % Additional check: reject unreasonably large fluctuations
                fluctMask = abs(u_prime) <= FLUCT_MAX & abs(v_prime) <= FLUCT_MAX;
                
                % Combined validity mask
                validMask = meanMask{cam} & instMask & fluctMask;
                
                % Apply mask: set invalid regions to NaN
                u_prime(~validMask) = NaN;
                v_prime(~validMask) = NaN;
                
                rejected_count = rejected_count + sum(~validMask(:));
                
                fluctuations.u_prime{frame, cam} = u_prime;
                fluctuations.v_prime{frame, cam} = v_prime;
                fluctuations.valid_mask{frame, cam} = validMask;
            end
            
            if VERBOSE
                total_points = nFrames * numel(u_inst);
                reject_pct = 100 * rejected_count / total_points;
                fprintf('✓ (%.1f%% rejected)\n', reject_pct);
            end
        end
        
        % Save fluctuations
        if SAVE_INDIVIDUAL_FLUCTUATIONS
            if VERBOSE
                S = whos('fluctuations');
                fprintf('  Saving to vel_fluctuations/ (%.2f GB)... ', S.bytes/1e9);
            end
            save(fluctFile, 'fluctuations', '-v7.3');
            if VERBOSE
                fprintf('✓\n');
            end
        end
        
        clear allCameras;
    end
    
    % Calculate loop-averaged variance with proper masking
    if VERBOSE
        fprintf('  Computing loop-averaged variance statistics...\n');
    end
    
    for cam = 1:nCameras
        u_variance_sum = zeros(size(fluctuations.u_prime{1, cam}), 'single');
        v_variance_sum = zeros(size(fluctuations.v_prime{1, cam}), 'single');
        valid_count = zeros(size(fluctuations.u_prime{1, cam}), 'single');
        
        for frame = 1:nFrames
            u_sq = fluctuations.u_prime{frame, cam}.^2;  % u'^2 in (m/s)^2
            v_sq = fluctuations.v_prime{frame, cam}.^2;  % v'^2 in (m/s)^2
            
            valid = ~isnan(u_sq) & ~isnan(v_sq);
            
            % Accumulate squared fluctuations (variance)
            u_variance_sum(valid) = u_variance_sum(valid) + u_sq(valid);
            v_variance_sum(valid) = v_variance_sum(valid) + v_sq(valid);
            valid_count = valid_count + single(valid);
        end
        
        % Calculate valid data fraction
        valid_frac = valid_count / nFrames;
        
        % Only compute variance where we have sufficient valid data
        sufficient_data = valid_frac >= MIN_VALID_FRAC;
        
        u_variance = NaN(size(u_variance_sum), 'single');
        v_variance = NaN(size(v_variance_sum), 'single');
        
        % Mean variance: <u'^2> in (m/s)^2
        u_variance(sufficient_data) = u_variance_sum(sufficient_data) ./ valid_count(sufficient_data);
        v_variance(sufficient_data) = v_variance_sum(sufficient_data) ./ valid_count(sufficient_data);
        
        % Store variance (NOT RMS!)
        loopStats.u_variance_mean{loopNo, cam} = u_variance;  % <u'^2> in (m/s)^2
        loopStats.v_variance_mean{loopNo, cam} = v_variance;  % <v'^2> in (m/s)^2
        loopStats.valid_fraction{loopNo, cam} = valid_frac;
        
        if VERBOSE
            good_data_pct = 100 * sum(sufficient_data(:)) / numel(sufficient_data);
            fprintf('    Camera %d: %.1f%% with sufficient data\n', cam, good_data_pct);
        end
    end
    
    clear fluctuations u_variance_sum v_variance_sum valid_count;
    
    if loopNo < nLoops
        pause(0.1);
    end
end

% Save loop-averaged statistics
loopStatsFile = fullfile(savePath, 'loop_averaged_variance_statistics.mat');
save(loopStatsFile, 'loopStats', '-v7.3');
if VERBOSE
    fprintf('\n✓ Saved loop-averaged variance statistics\n');
end

%% Calculate overall turbulence statistics
if VERBOSE
    fprintf('\n=== Computing RMS and Turbulence Statistics ===\n');
end

turbStats = struct();
turbStats.description = 'Turbulence statistics: variance <u''^2>, RMS sqrt(<u''^2>), with quality control';
turbStats.n_loops = nLoops;
turbStats.n_cameras = nCameras;
turbStats.n_frames = nFrames;
turbStats.qc_params.vel_range = [VEL_MIN, VEL_MAX];
turbStats.qc_params.max_fluctuation = FLUCT_MAX;
turbStats.qc_params.min_valid_fraction = MIN_VALID_FRAC;

for cam = 1:nCameras
    if VERBOSE
        fprintf(' Camera %d: ', cam);
    end
    
    [nRows, nCols] = size(loopStats.u_variance_mean{1, cam});
    
    % Stack variance from all loops
    u_variance_all_loops = zeros(nRows, nCols, nLoops, 'single');
    v_variance_all_loops = zeros(nRows, nCols, nLoops, 'single');
    
    for loopNo = 1:nLoops
        u_variance_all_loops(:, :, loopNo) = loopStats.u_variance_mean{loopNo, cam};
        v_variance_all_loops(:, :, loopNo) = loopStats.v_variance_mean{loopNo, cam};
    end
    
    % Mean variance across all loops: <u'^2> in (m/s)^2
    u_variance_mean = mean(u_variance_all_loops, 3, 'omitnan');
    v_variance_mean = mean(v_variance_all_loops, 3, 'omitnan');
    
    % Standard deviation of variance across loops
    u_variance_std = std(u_variance_all_loops, 0, 3, 'omitnan');
    v_variance_std = std(v_variance_all_loops, 0, 3, 'omitnan');
    
    % RMS (standard deviation of fluctuations): sqrt(<u'^2>) in m/s
    u_rms = sqrt(u_variance_mean);  % m/s
    v_rms = sqrt(v_variance_mean);  % m/s
    
    % Store variance fields
    turbStats.u_variance{cam} = u_variance_mean;      % <u'^2> in (m/s)^2
    turbStats.v_variance{cam} = v_variance_mean;      % <v'^2> in (m/s)^2
    turbStats.u_variance_std{cam} = u_variance_std;   % Std of variance across loops
    turbStats.v_variance_std{cam} = v_variance_std;
    
    % Store RMS fields
    turbStats.u_rms{cam} = u_rms;                     % sqrt(<u'^2>) in m/s
    turbStats.v_rms{cam} = v_rms;                     % sqrt(<v'^2>) in m/s
    
    % Store mean mask for plotting
    turbStats.valid_mask{cam} = meanMask{cam};
    
    % Calculate spatial-average statistics
    turbStats.u_variance_mean_scalar(cam) = mean(u_variance_mean(:), 'omitnan');  % Mean variance in (m/s)^2
    turbStats.v_variance_mean_scalar(cam) = mean(v_variance_mean(:), 'omitnan');
    
    turbStats.u_rms_mean_scalar(cam) = mean(u_rms(:), 'omitnan');                 % Mean RMS in m/s
    turbStats.v_rms_mean_scalar(cam) = mean(v_rms(:), 'omitnan');
    turbStats.u_rms_std_scalar(cam) = std(u_rms(:), 'omitnan');                   % Spatial std of RMS
    turbStats.v_rms_std_scalar(cam) = std(v_rms(:), 'omitnan');
    
    % Data quality metrics
    valid_percent = 100 * sum(~isnan(u_rms(:))) / numel(u_rms);
    rms_max = max(u_rms(:), [], 'omitnan');
    variance_max = max(u_variance_mean(:), [], 'omitnan');
    
    if VERBOSE
        fprintf('u_rms = %.3f ± %.3f m/s, max = %.2f m/s (%.1f%% valid) ✓\n', ...
            turbStats.u_rms_mean_scalar(cam), turbStats.u_rms_std_scalar(cam), ...
            rms_max, valid_percent);
    end
    
    clear u_variance_all_loops v_variance_all_loops;
end

% Save turbulence statistics
turbStatsFile = fullfile(savePath, 'turbulence_statistics.mat');
save(turbStatsFile, 'turbStats', '-v7.3');

if VERBOSE
    endTime = toc;
    fprintf('\n✓ Fluctuation analysis completed in %.2f min\n', endTime/60);
    
    fprintf('\n=== Summary ===\n');
    fprintf('Fluctuation files: %d × ~2.58 GB in vel_fluctuations/ folders\n', nLoops);
    fprintf('Total fluctuation storage: ~%.1f GB\n', 2.58 * nLoops);
    fprintf('Statistics files: 2 files in root directory\n');
    fprintf('Peak RAM usage: ~2.6 GB per loop\n');
    fprintf('\nFile structure:\n');
    fprintf('  %s\n', savePath);
    fprintf('  ├── loop=0/vel_fluctuations/fluctuations_all_frames.mat (2.58 GB)\n');
    fprintf('  ├── loop=1/vel_fluctuations/fluctuations_all_frames.mat (2.58 GB)\n');
    fprintf('  ├── ...\n');
    fprintf('  ├── loop_averaged_variance_statistics.mat\n');
    fprintf('  └── turbulence_statistics.mat\n\n');
    
    fprintf('Stored statistics:\n');
    fprintf('  - u_variance{cam}: Variance <u''^2> field in (m/s)^2\n');
    fprintf('  - u_rms{cam}: RMS sqrt(<u''^2>) field in m/s\n');
    fprintf('  - u_rms_mean_scalar: Spatial-average RMS in m/s\n');
    fprintf('  - u_variance_mean_scalar: Spatial-average variance in (m/s)^2\n\n');
    fprintf('To compute turbulence intensity: TI = u_variance / U_inf^2\n');
end

end

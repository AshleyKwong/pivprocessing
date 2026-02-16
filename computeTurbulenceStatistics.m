function turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField, varargin)
% COMPUTETURBULENCESTATISTICS Calculate velocity fluctuations and turbulence intensity
%
% Syntax:
%   turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField)
%   turbStats = computeTurbulenceStatistics(..., 'Name', Value)

% Example
% turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField, ...
%     'VelocityThreshold', [0, 50], ...      % Valid velocity range in m/s
%     'FluctuationThreshold', 10, ...        % Max reasonable fluctuation
%     'MinValidFraction', 0.7);              % Need 70% valid frames per point
% 
% % For stricter quality control
% turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField, ...
%     'VelocityThreshold', [5, 40], ...      % Tighter range
%     'FluctuationThreshold', 5, ...         % More conservative
%     'MinValidFraction', 0.8);              % Need 80% valid data



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
%
% Outputs:
%   turbStats - Structure containing:
%       .TI_u{cam}           - Turbulence intensity field for u component (m/s)
%       .TI_v{cam}           - Turbulence intensity field for v component (m/s)
%       .TI_u_mean(cam)      - Spatial-average TI_u (scalar)
%       .TI_v_mean(cam)      - Spatial-average TI_v (scalar)
%       .u_rms_mean_field    - Mean RMS field across all loops
%       .u_rms_std           - Standard deviation of RMS across loops
%
% File Structure Created:
%   savePath/
%   ├── loop=X/
%   │   └── vel_fluctuations/
%   │       └── fluctuations_all_frames.mat
%   ├── loop_averaged_rms_statistics.mat
%   └── turbulence_intensity_statistics.mat
%
% Example:
%   turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField);
%   turbStats = computeTurbulenceStatistics(savePath, totalLoops, avgVelocityField, ...
%                                           'SaveIndividualFluctuations', true, ...
%                                           'ForceRecompute', false);
%
% Author: Ashley Kwong
% Date: 02/15/2026

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
addParameter(p, 'VelocityThreshold', [0, 50], @isnumeric); % [min, max] valid velocity in m/s
addParameter(p, 'FluctuationThreshold', 10, @isnumeric); % Max reasonable fluctuation in m/s
addParameter(p, 'MinValidFraction', 0.7, @isnumeric); % Minimum fraction of valid frames per point

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
% nLoops = length(totalLoops);
nLoops = 1; 
fprintf("Testing on loop 1 !!\n"); 
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
loopStats.u_rms_mean = cell(nLoops, nCameras);
loopStats.v_rms_mean = cell(nLoops, nCameras);
loopStats.valid_fraction = cell(nLoops, nCameras); % Track data quality

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
        fluctuations.valid_mask = cell(nFrames, nCameras); % Store validity per frame
        fluctuations.loop_name = totalLoops(loopNo).name;
        fluctuations.n_frames = nFrames;
        fluctuations.n_cameras = nCameras;
        
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
                
                % Calculate fluctuations
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
    
    % Calculate loop-averaged RMS with proper masking
    if VERBOSE
        fprintf('  Computing loop-averaged RMS statistics...\n');
    end
    
    for cam = 1:nCameras
        u_rms_sum = zeros(size(fluctuations.u_prime{1, cam}), 'single');
        v_rms_sum = zeros(size(fluctuations.v_prime{1, cam}), 'single');
        valid_count = zeros(size(fluctuations.u_prime{1, cam}), 'single');
 
        for frame = 1:nFrames
            u_sq = fluctuations.u_prime{frame, cam}.^2;
            v_sq = fluctuations.v_prime{frame, cam}.^2;
            
            valid = ~isnan(u_sq) & ~isnan(v_sq);
            
            u_rms_sum(valid) = u_rms_sum(valid) + u_sq(valid);
            v_rms_sum(valid) = v_rms_sum(valid) + v_sq(valid);
            valid_count = valid_count + single(valid);
        end
        
        % Calculate valid data fraction
        valid_frac = valid_count / nFrames;
        
        % Only compute RMS where we have sufficient valid data
        sufficient_data = valid_frac >= MIN_VALID_FRAC;
        
        u_rms = NaN(size(u_rms_sum), 'single');
        v_rms = NaN(size(v_rms_sum), 'single');
        
        u_rms(sufficient_data) = u_rms_sum(sufficient_data) ./ valid_count(sufficient_data);
        v_rms(sufficient_data) = v_rms_sum(sufficient_data) ./ valid_count(sufficient_data);
        
        loopStats.u_rms_mean{loopNo, cam} = u_rms;
        loopStats.v_rms_mean{loopNo, cam} = v_rms;
        loopStats.valid_fraction{loopNo, cam} = valid_frac;
        
        if VERBOSE
            good_data_pct = 100 * sum(sufficient_data(:)) / numel(sufficient_data);
            fprintf('    Camera %d: %.1f%% with sufficient data\n', cam, good_data_pct);
        end
    end
    
    clear fluctuations u_rms_sum v_rms_sum valid_count;
    
    if loopNo < nLoops
        pause(0.1);
    end
end

% Save loop-averaged statistics
loopStatsFile = fullfile(savePath, 'loop_averaged_rms_statistics.mat');
save(loopStatsFile, 'loopStats', '-v7.3');
if VERBOSE
    fprintf('\n✓ Saved loop-averaged RMS statistics\n');
end

%% Calculate overall turbulence intensity
if VERBOSE
    fprintf('\n=== Computing Turbulence Intensity ===\n');
end

turbStats = struct();
turbStats.description = 'Turbulence intensity statistics with quality control';
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
    
    [nRows, nCols] = size(loopStats.u_rms_mean{1, cam});
    
    u_rms_all_loops = zeros(nRows, nCols, nLoops, 'single');
    v_rms_all_loops = zeros(nRows, nCols, nLoops, 'single');
    
    for loopNo = 1:nLoops
        u_rms_all_loops(:, :, loopNo) = loopStats.u_rms_mean{loopNo, cam};
        v_rms_all_loops(:, :, loopNo) = loopStats.v_rms_mean{loopNo, cam};
    end
    
    u_rms_mean = mean(u_rms_all_loops, 3, 'omitnan');
    v_rms_mean = mean(v_rms_all_loops, 3, 'omitnan');
    
    u_rms_std = std(u_rms_all_loops, 0, 3, 'omitnan');
    v_rms_std = std(v_rms_all_loops, 0, 3, 'omitnan');
    
    turbStats.TI_u{cam} = sqrt(u_rms_mean);
    turbStats.TI_v{cam} = sqrt(v_rms_mean);
    
    turbStats.u_rms_std{cam} = u_rms_std;
    turbStats.v_rms_std{cam} = v_rms_std;
    turbStats.u_rms_mean_field{cam} = u_rms_mean;
    turbStats.v_rms_mean_field{cam} = v_rms_mean;
    
    % Store mean mask for plotting
    turbStats.valid_mask{cam} = meanMask{cam};
    
    turbStats.TI_u_mean(cam) = mean(turbStats.TI_u{cam}(:), 'omitnan');
    turbStats.TI_v_mean(cam) = mean(turbStats.TI_v{cam}(:), 'omitnan');
    turbStats.TI_u_std(cam) = std(turbStats.TI_u{cam}(:), 'omitnan');
    turbStats.TI_v_std(cam) = std(turbStats.TI_v{cam}(:), 'omitnan');
    
    valid_percent = 100 * sum(~isnan(turbStats.TI_u{cam}(:))) / numel(turbStats.TI_u{cam});
    ti_max = max(turbStats.TI_u{cam}(:), [], 'omitnan');
    
    if VERBOSE
        fprintf('TI_u = %.3f ± %.3f m/s, max = %.2f m/s (%.1f%% valid) ✓\n', ...
            turbStats.TI_u_mean(cam), turbStats.TI_u_std(cam), ti_max, valid_percent);
    end
    
    clear u_rms_all_loops v_rms_all_loops;
end

turbStatsFile = fullfile(savePath, 'turbulence_intensity_statistics.mat');
save(turbStatsFile, 'turbStats', '-v7.3');

if VERBOSE
    endTime = toc;
    fprintf('\n✓ Fluctuation analysis completed in %.2f min\n', endTime/60);
end

end

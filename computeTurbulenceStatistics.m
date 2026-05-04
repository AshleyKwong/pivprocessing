function fluctStats = computeTurbulenceStatistics(savePath, totalLoops, meanFile, varargin)
% COMPUTETURBULENCESTATISTICS
% Computes turbulence statistics from merged instantaneous velocity fields.
% Operates on mergedvelocityfields.mat (world-grid) rather than per-camera fields.
%
% Inputs:
%   savePath    - Base path where loop folders are located
%   totalLoops  - Structure array of loop directories (from dir())
%   meanFile    - Full path to merged mean UV file (contains U_hann_mean, V_hann_mean)
%                 e.g. '.../merge_stats_20260308_154918/merged_meanUV_14loops_20260308_154918.mat'
%
% Optional Name-Value Pairs:
%   'ForceRecompute'       - Recompute even if output exists (default: false)
%   'Verbose'              - Print progress (default: true)
%   'VelocityThreshold'    - [min, max] valid velocity in m/s (default: [-10, 35])
%   'FluctuationThreshold' - Max abs fluctuation in m/s (default: 20)
%   'MinValidFraction'     - Min fraction of valid frames per point (default: 0.7)
%
% Outputs:
%   fluctStats.u_prime        - {N_frames x 1} cell, instantaneous u' (m/s)
%   fluctStats.v_prime        - {N_frames x 1} cell, instantaneous v' (m/s)
%   fluctStats.u_variance     - [Ny x Nx] time-averaged <u'^2> (m/s)^2
%   fluctStats.v_variance     - [Ny x Nx] time-averaged <v'^2> (m/s)^2
%   fluctStats.uv_stress      - [Ny x Nx] time-averaged <u'v'> (m/s)^2
%   fluctStats.u_rms          - [Ny x Nx] sqrt(<u'^2>) m/s
%   fluctStats.v_rms          - [Ny x Nx] sqrt(<v'^2>) m/s
%   fluctStats.valid_fraction - [Ny x Nx] fraction of valid frames per point
%   fluctStats.worldX         - world grid x (mm)
%   fluctStats.worldY         - world grid y (mm)
%
% Author: Ashley Kwong
% Date:   08/03/2026

%% Parse inputs
p = inputParser;
addRequired(p, 'savePath',   @ischar);
addRequired(p, 'totalLoops', @isstruct);
addRequired(p, 'meanFile',   @ischar);   % ← added
addParameter(p, 'ForceRecompute',       false,    @islogical);
addParameter(p, 'Verbose',              true,     @islogical);
addParameter(p, 'VelocityThreshold',    [-10,35], @isnumeric);
addParameter(p, 'FluctuationThreshold', 20,       @isnumeric);
addParameter(p, 'MinValidFraction',     0.7,      @isnumeric);
parse(p, savePath, totalLoops, meanFile, varargin{:});


FORCE_RECOMPUTE = p.Results.ForceRecompute;
VERBOSE         = p.Results.Verbose;
VEL_MIN         = p.Results.VelocityThreshold(1);
VEL_MAX         = p.Results.VelocityThreshold(2);
FLUCT_MAX       = p.Results.FluctuationThreshold;
MIN_VALID_FRAC  = p.Results.MinValidFraction;

%% Check if output already exists
outFile = fullfile(savePath, 'turbulence_statistics.mat');
if isfile(outFile) && ~FORCE_RECOMPUTE
    if VERBOSE
        fprintf('✓ turbulence_statistics.mat already exists — loading.\n');
        fprintf('  Set ForceRecompute=true to recompute.\n');
    end
    loaded = load(outFile, 'fluctStats');
    fluctStats = loaded.fluctStats;
    return;
end

%% Load merged mean velocity field
if ~isfile(meanFile)
    error('Mean file not found: %s', meanFile);
end
avg    = load(meanFile, 'U_hann_mean', 'V_hann_mean');
U_mean = avg.U_hann_mean;
V_mean = avg.V_hann_mean;
clear avg;
fprintf('✓ Loaded merged mean from: %s\n', meanFile);


if VERBOSE
    fprintf('\n=== Computing Turbulence Statistics (merged field) ===\n');
    fprintf('QC thresholds: vel=[%.1f, %.1f] m/s | fluct<=%.1f m/s | min_valid=%.0f%%\n', ...
        VEL_MIN, VEL_MAX, FLUCT_MAX, MIN_VALID_FRAC*100);
end

%% First pass — collect all frames across all loops, accumulate statistics
% We do one pass to get variance/stress, storing u_prime cells as we go

nLoops      = length(totalLoops);
allU_prime  = {};   % will grow to {N_total x 1}
allV_prime  = {};
frameCount  = 0;

% Preallocate accumulators (size determined from first frame)
u_var_sum   = [];
v_var_sum   = [];
uv_sum      = [];
valid_count = [];

tic;
for loopNo = 1:nLoops
    loopName   = totalLoops(loopNo).name;

    % Find mergedvelocityfields file in this loop folder (timestamp-agnostic)
    mergedFileSearch = dir(fullfile(savePath, loopName, 'mergedvelocityfields_*.mat'));
    if isempty(mergedFileSearch)
        warning('No mergedvelocityfields_*.mat found in %s — skipping.', loopName);
        continue;
    end
    if length(mergedFileSearch) > 1
        warning('Multiple mergedvelocityfields files found in %s — using most recent.', loopName);
        [~, idx] = max([mergedFileSearch.datenum]);
        mergedFileSearch = mergedFileSearch(idx);
    end
    mergedFile = fullfile(mergedFileSearch.folder, mergedFileSearch.name);
    
    data   = load(mergedFile, 'mergedFrames', 'worldX', 'worldY');
    nFrames = numel(data.mergedFrames.U);

    % Grab grid from first loop
    if loopNo == 1
        worldX      = data.worldX;
        worldY      = data.worldY;
        [Ny, Nx]    = size(data.mergedFrames.U{1});
        u_var_sum   = zeros(Ny, Nx, 'single');
        v_var_sum   = zeros(Ny, Nx, 'single');
        uv_sum      = zeros(Ny, Nx, 'single');
        valid_count = zeros(Ny, Nx, 'single');
        if VERBOSE
            fprintf('Grid: [%d x %d], worldX:[%.1f, %.1f] mm, worldY:[%.1f, %.1f] mm\n', ...
                Ny, Nx, min(worldX(:)), max(worldX(:)), min(worldY(:)), max(worldY(:)));
        end
    end

    if VERBOSE
        fprintf('  Loop %d/%d (%s): %d frames\n', loopNo, nLoops, loopName, nFrames);
    end

    for fr = 1:nFrames
        frameCount = frameCount + 1;
        U_inst = data.mergedFrames.U{fr};
        V_inst = data.mergedFrames.V{fr};

        % QC mask on raw velocity
        qcMask = ~isnan(U_inst) & ~isnan(V_inst) & ...
                 U_inst >= VEL_MIN & U_inst <= VEL_MAX & ...
                 V_inst >= VEL_MIN & V_inst <= VEL_MAX;

        % Fluctuations
        u_prime = single(U_inst - U_mean);
        v_prime = single(V_inst - V_mean);

        % QC mask on fluctuation magnitude
        fluctMask = abs(u_prime) <= FLUCT_MAX & abs(v_prime) <= FLUCT_MAX;
        validMask = qcMask & fluctMask;

        % Apply mask
        u_prime(~validMask) = NaN;
        v_prime(~validMask) = NaN;

        % Store instantaneous fluctuation fields
        allU_prime{frameCount, 1} = u_prime; %#ok<AGROW>
        allV_prime{frameCount, 1} = v_prime; %#ok<AGROW>

        % Accumulate statistics
        valid = ~isnan(u_prime) & ~isnan(v_prime);
        u_var_sum(valid)   = u_var_sum(valid)   + u_prime(valid).^2;
        v_var_sum(valid)   = v_var_sum(valid)   + v_prime(valid).^2;
        uv_sum(valid)      = uv_sum(valid)       + u_prime(valid) .* v_prime(valid);
        valid_count(valid) = valid_count(valid)  + 1;
    end

    clear data;

    if VERBOSE
        fprintf('    Cumulative frames: %d | RAM: %.2f GB\n', ...
            frameCount, getfield(whos('allU_prime'), 'bytes')/1e9);
    end
end

%% Compute time-averaged statistics
if VERBOSE
    fprintf('\n=== Finalising statistics over %d frames ===\n', frameCount);
end

valid_fraction = valid_count / frameCount;
sufficient     = valid_fraction >= MIN_VALID_FRAC;

u_variance = NaN(Ny, Nx, 'single');
v_variance = NaN(Ny, Nx, 'single');
uv_stress  = NaN(Ny, Nx, 'single');

u_variance(sufficient) = u_var_sum(sufficient) ./ valid_count(sufficient);
v_variance(sufficient) = v_var_sum(sufficient) ./ valid_count(sufficient);
uv_stress(sufficient)  = uv_sum(sufficient)    ./ valid_count(sufficient);

u_rms = sqrt(u_variance);
v_rms = sqrt(v_variance);

if VERBOSE
    fprintf('  u_rms: mean=%.4f m/s | max=%.4f m/s\n', ...
        mean(u_rms(:),'omitnan'), max(u_rms(:),[],'omitnan'));
    fprintf('  v_rms: mean=%.4f m/s | max=%.4f m/s\n', ...
        mean(v_rms(:),'omitnan'), max(v_rms(:),[],'omitnan'));
    fprintf('  Valid fraction: %.1f%% of points above threshold\n', ...
        100*sum(sufficient(:))/numel(sufficient));
end

%% Pack output
fluctStats.u_prime        = allU_prime;
fluctStats.v_prime        = allV_prime;
fluctStats.u_variance     = u_variance;
fluctStats.v_variance     = v_variance;
fluctStats.uv_stress      = uv_stress;
fluctStats.u_rms          = u_rms;
fluctStats.v_rms          = v_rms;
fluctStats.valid_fraction = valid_fraction;
fluctStats.worldX         = worldX;
fluctStats.worldY         = worldY;
fluctStats.n_frames       = frameCount;
fluctStats.n_loops        = nLoops;
fluctStats.qc_params.vel_range        = [VEL_MIN, VEL_MAX];
fluctStats.qc_params.max_fluctuation  = FLUCT_MAX;
fluctStats.qc_params.min_valid_frac   = MIN_VALID_FRAC;
fluctStats.description = 'Turbulence statistics on merged world-grid field. u_prime = u_merged - U_mean_merged.';

%% Save
save(outFile, 'fluctStats', '-v7.3');
if VERBOSE
    S = whos('fluctStats');
    fprintf('✓ Saved turbulence_statistics.mat (%.2f GB) in %.1f min\n', ...
        S.bytes/1e9, toc/60);
end

end

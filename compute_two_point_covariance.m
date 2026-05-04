% Ashley Kwong
% compute_two_point_covariance.m
%
% Computes the two-point covariance field C(x,y) = <u'(x_ref,y_ref) * u'(x,y)>
% per camera domain, then merges the resulting covariance fields using the
% Hanning-weighted merge. The normalised correlation R is computed via
% twopointcorr.m using pre-computed merged variance/std fields.
%
% Pipeline:
%   1. For each (x_ref, y_ref) pair, identify the owning camera
%      (closest domain centre).
%   2. Loop over all loops and frames, accumulating the covariance
%      sum per camera using an online accumulator (memory efficient —
%      only one snapshot in memory at a time).
%   3. After all snapshots, divide by N to get C_avg per camera.
%   4. Merge the 5 per-camera covariance fields with Hanning blending.
%   5. Call twopointcorr to normalise and plot R.

clear; clc; close all;
addpath(genpath('/iridisfs/scratch/ak1u24/calib_tools'));

%% ======== USER OPTIONS ================================================
savePath   = '/iridisfs/scratch/ak1u24/case1_fullmask';
masks      = {};   % empty = no masking
TEST_LOOP  = 1;   % set to [] to process all loops, integer to test single

% Pre-computed merged std fields (Ny x Nx) on the merged world grid
% used as denominator in twopointcorr normalisation
mergedStdFile = '/iridisfs/scratch/ak1u24/case1_fullmask/merge_instantaneousavg_20260309_213859/merged_turbrms_hann_14loops_20260309_213859.mat';% Expected variables inside: turbStats_merged.u_std, turbStats_merged.v_std
% and worldX_merged, worldY_merged

% Fluctuation file name within each loop folder
fluctFileName = 'fluctuations_all_frames.mat';   % struct: fluctuations.u_prime, .v_prime

% Reference point pairs [x_ref_mm, y_ref_mm] — one row per pair
% Choose points well within a single camera interior, away from overlaps
refPairs = [
    100.0,  5.0;   % pair 1: (x_ref, y_ref) in mm
    100.0, 15.0;   % pair 2
];

% Correlation type: 'uu', 'vv', 'uv', or 'vu'
%   'uu' : <u'(ref) * u'(x,y)>
%   'vv' : <v'(ref) * v'(x,y)>
%   'uv' : <u'(ref) * v'(x,y)>
%   'vu' : <v'(ref) * u'(x,y)>
corrType = 'uu';


% ======================================================================

%% TIMESTAMP AND OUTPUT DIRECTORY
tstamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(savePath, sprintf('two_point_covariance_%s', tstamp));
mkdir(outDir);
fprintf('Output directory: %s\n', outDir);

%% DISCOVER LOOP FOLDERS (same pattern as PIV_mergeAndStats)
d = dir(savePath);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));
loopPattern  = '^loop\s*=\s*\d+$';
validLoops   = false(size(d));
for i = 1:length(d)
    validLoops(i) = ~isempty(regexpi(d(i).name, loopPattern));
end
totalLoops = d(validLoops);
if isempty(totalLoops)
    error('No loop folders found in %s', savePath);
end
fprintf('Found %d loop folders.\n', length(totalLoops));

if ~isempty(TEST_LOOP)
    totalLoops = totalLoops(TEST_LOOP);
    fprintf('TEST_LOOP = %d: processing only %s\n', TEST_LOOP, totalLoops(1).name);
end

%% LOAD WINDOW CENTRES (from first available loop — same grid for all loops)
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
nCams = length(windowCenterCameras_mm.x1_mm);
fprintf('Loaded window centres: %d cameras\n', nCams);

% Print camera domain extents
for c = 1:nCams
    x1 = windowCenterCameras_mm.x1_mm{c};
    x2 = windowCenterCameras_mm.x2_mm{c};
    fprintf('  Cam%d: x1=[%.2f, %.2f] mm | x2=[%.2f, %.2f] mm\n', c, ...
        min(x1(:)), max(x1(:)), min(x2(:)), max(x2(:)));
end

%% IDENTIFY OWNING CAMERA FOR EACH REFERENCE PAIR
nPairs   = size(refPairs, 1);
ownerCam = zeros(nPairs, 1);   % which camera owns each ref pair
iy_refs  = zeros(nPairs, 1);   % row index within owning camera grid
ix_refs  = zeros(nPairs, 1);   % col index within owning camera grid

for p = 1:nPairs
    xr = refPairs(p, 1);
    yr = refPairs(p, 2);

    % Compute distance from (xr,yr) to the domain centre of each camera
    dist_to_centre = inf(1, nCams);
    for c = 1:nCams
        x1c = windowCenterCameras_mm.x1_mm{c};
        x2c = windowCenterCameras_mm.x2_mm{c};
        cx  = mean(x1c(:));
        cy  = mean(x2c(:));
        dist_to_centre(c) = sqrt((xr - cx)^2 + (yr - cy)^2);
    end

    % Pick camera whose centre is closest — consistent owning camera
    [~, camIdx] = min(dist_to_centre);
    ownerCam(p) = camIdx;

    % Find nearest grid point within that camera
    x1c = windowCenterCameras_mm.x1_mm{camIdx};   % (Ny x Nx) or vector
    x2c = windowCenterCameras_mm.x2_mm{camIdx};
    x_vec = x1c(1, :);   % streamwise vector
    y_vec = x2c(:, 1);   % wall-normal vector

    [~, ix_refs(p)] = min(abs(x_vec - xr));
    [~, iy_refs(p)] = min(abs(y_vec - yr));

    fprintf('Pair %d: (x_ref=%.2f, y_ref=%.2f) mm → Cam%d, grid index (iy=%d, ix=%d), nearest=(%.3f, %.3f) mm\n', ...
        p, xr, yr, camIdx, iy_refs(p), ix_refs(p), ...
        y_vec(iy_refs(p)), x_vec(ix_refs(p)));
end

%% PREALLOCATE ONLINE ACCUMULATORS
% C_accum{p, cam} = running sum of u'(xref,yref) * u'(x,y) for pair p, camera cam
% Shape of each accumulator matches the camera's own grid
C_accum  = cell(nPairs, nCams);
N_accum  = zeros(nPairs, 1);   % total frame count per pair (same for all cams)

for p = 1:nPairs
    for c = 1:nCams
        x1c = windowCenterCameras_mm.x1_mm{c};
        [Ny_c, Nx_c] = size(x1c);
        C_accum{p, c} = zeros(Ny_c, Nx_c, 'double');
    end
end

%% MAIN ACCUMULATION LOOP
fprintf('\n=== Accumulating covariance ===\n');

for loopNo = 1:length(totalLoops)
    loopName   = totalLoops(loopNo).name;
    loopFolder = fullfile(savePath, loopName);
    fluctFile  = fullfile(loopFolder, 'vel_fluctuations', fluctFileName);

    if ~isfile(fluctFile)
        warning('Fluctuation file not found in %s — skipping.', loopName);
        continue;
    end

    data     = load(fluctFile, 'fluctuations');
    nFrames  = size(data.fluctuations.u_prime, 1);
    fprintf('\n--- %s: %d frames ---\n', loopName, nFrames);

    for fr = 1:nFrames
        for p = 1:nPairs
            oc  = ownerCam(p);
            iyr = iy_refs(p);
            ixr = ix_refs(p);

            % Extract reference scalar — always from the FIRST corrType of corrType
            switch corrType(1)
                case 'u';  u_ref_field = data.fluctuations.u_prime{fr, oc};
                case 'v';  u_ref_field = data.fluctuations.v_prime{fr, oc};
            end
            u_scalar = u_ref_field(iyr, ixr);

            if isnan(u_scalar); continue; end

            % Accumulate against SECOND corrType of corrType
            for c = 1:nCams
                switch corrType(2)
                    case 'u';  u_field = data.fluctuations.u_prime{fr, c};
                    case 'v';  u_field = data.fluctuations.v_prime{fr, c};
                end
                C_accum{p, c} = C_accum{p, c} + u_scalar .* u_field;
            end

            N_accum(p) = N_accum(p) + 1;
        end

        if mod(fr, 50) == 0
            fprintf('  Frame %d/%d\n', fr, nFrames);
        end
    end

    clear data;
end

fprintf('\n=== Accumulation complete ===\n');
for p = 1:nPairs
    fprintf('  Pair %d: N = %d snapshots\n', p, N_accum(p));
end

%% AVERAGE AND MERGE PER-CAMERA COVARIANCE FIELDS
fprintf('\n=== Averaging and merging covariance fields ===\n');

% Load merged std fields for normalisation
stdData = load(mergedStdFile, 'U_rms', 'V_rms', 'worldX', 'worldY');
switch corrType(1)
    case 'u'; std_ref_field = stdData.U_rms;
    case 'v'; std_ref_field = stdData.V_rms;
end
switch corrType(2)
    case 'u'; std_field = stdData.U_rms;
    case 'v'; std_field = stdData.V_rms;
end
worldX_merged = stdData.worldX;
worldY_merged = stdData.worldY;

for p = 1:nPairs
    fprintf('\nPair %d: (x_ref=%.2f, y_ref=%.2f) mm\n', p, refPairs(p,1), refPairs(p,2));

    % Average each per-camera covariance field
    C_avg_cams = cell(1, nCams);
    for c = 1:nCams
        C_avg_cams{c} = C_accum{p, c} ./ N_accum(p);
    end

    % Merge per-camera covariance fields onto world grid using Hanning blend
    % Reuse merge function: pass C_avg as "U", zeros as "V" (unused)
    V_dummy = cell(1, nCams);
    for c = 1:nCams
        V_dummy{c} = zeros(size(C_avg_cams{c}));
    end

    [worldX_out, worldY_out, C_merged, ~] = merge_cameras_python_style_mean( ...
        windowCenterCameras_mm, C_avg_cams, V_dummy, masks, 'hann', []);

    fprintf('  Merged covariance field: [%d x %d]\n', size(C_merged,1), size(C_merged,2));
    fprintf('  NaNs: %.1f%%\n', 100*sum(isnan(C_merged(:)))/numel(C_merged));

    % Save merged covariance field
    xr  = refPairs(p, 1);
    yr  = refPairs(p, 2);
    outFile = fullfile(outDir, sprintf('covariance_%s_xref%.1f_yref%.1f_%s.mat', ...
        corrType, xr, yr, tstamp));
    save(outFile, 'C_merged', 'worldX_out', 'worldY_out', ...
        'xr', 'yr', 'N_accum', '-v7.3');
    fprintf('  ✓ Saved merged covariance → %s\n', outFile);
    
    assert(isequal(size(C_merged), size(std_field)), ...
    'Grid mismatch: C_merged and std_field have different sizes. Check merge grids.');

    % --- Normalise and plot via twopointcorr ---
    % std_ref: scalar std at the reference point on the merged grid
    [~, iy_ref_merged] = min(abs(worldY_merged(:,1)   - yr));
    [~, ix_ref_merged] = min(abs(worldX_merged(1,:)   - xr));


    % std_ref: scalar std of the REFERENCE component at (x_ref, y_ref)
    std_ref_scalar = std_ref_field(iy_ref_merged, ix_ref_merged);

    % Normalise: C(x,y) / ( std_ref_scalar * std_field(x,y) )
    R = C_merged ./ (std_ref_scalar .* std_field);
    R(abs(R) > 1.5) = NaN;
   
    % Save R field
    outFileR = fullfile(outDir, sprintf('R_%s_xref%.1f_yref%.1f_%s.mat', ...
        corrType, xr, yr, tstamp));
    save(outFileR, 'R', 'worldX_merged', 'worldY_merged', ...
        'xr', 'yr', '-v7.3');
    fprintf('  ✓ Saved R → %s\n', outFileR);
end

fprintf('\n=== Done: %s ===\n', tstamp);

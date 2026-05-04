% =========================================================================
% compute_quad_spectra_skewkurt.m
%
% For a set of user-specified xref locations, accumulates turbulence
% statistics across all loop folders in a single pass.
%
% At each xref, extracts a column of pixels from each camera that covers
% that location, applies Hann-weighted blending in the overlap region,
% and accumulates:
%   - u'rms, v'rms
%   - Reynolds shear stress -<u'v'>
%   - Skewness  S_u = <u'^3> / <u'^2>^(3/2)
%   - Kurtosis  K_u = <u'^4> / <u'^2>^2
%   - Quadrant fractions Q1-Q4
%   - Quadrant RS contributions per hole size H (Lu & Willmarth 1973)
%
% NaN protection: NaNs zero-filled before accumulation.
% No plotting — saves results only.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc;

%% ===================== USER INPUTS =====================

caseDir       = '/iridisfs/scratch/ak1u24/case1_fullmask/';
outDir        = '/iridisfs/scratch/ak1u24/case1_fullmask/turbulence_statistics_PIV/';
caseLabel     = 'case2';
fluctFileName = 'fluctuations_all_frames.mat';

% Merged RMS field — used as fixed reference for Lu & Willmarth hole
mergedRmsFile = ['/iridisfs/scratch/ak1u24/case1_fullmask/' ...
    'merge_instantaneousavg_20260413_103619/' ...
    'merged_turbrms_hann_14loops_20260413_103619.mat'];

% xref locations to analyse (mm)
xref_targets = [350.0, 550.0, 711.0, 1000.0];

% Hole size parameter — Lu & Willmarth (1973)
H_vals = [0, 1.0];
nH     = numel(H_vals);

%% ===================== SETUP =====================

if ~exist(outDir, 'dir'), mkdir(outDir); end

nXref = numel(xref_targets);

% --- Find window centres ---
d = dir(caseDir);
d = d([d.isdir]);
physicalLocationFile = '';
for k = 1:numel(d)
    candidate = fullfile(caseDir, d(k).name, 'windowCenterCameras_mm.mat');
    if isfile(candidate)
        physicalLocationFile = candidate;
        break;
    end
end
assert(~isempty(physicalLocationFile), 'windowCenterCameras_mm.mat not found');
load(physicalLocationFile, 'windowCenterCameras_mm');
nCams = numel(windowCenterCameras_mm.x1_mm);
fprintf('Cameras: %d\n', nCams);

% --- Find all loop folders ---
loopDirs = dir(fullfile(caseDir, 'loop*'));
loopDirs = loopDirs([loopDirs.isdir]);
nLoops   = numel(loopDirs);
assert(nLoops > 0, 'No loop folders found in %s', caseDir);
fprintf('Found %d loop folders\n', nLoops);
fprintf('xref targets: ');
fprintf('%.1f  ', xref_targets);
fprintf('\n\n');

%% ===================== PRE-COMPUTE COLUMN WEIGHTS =====================

fprintf('Computing column weights for each xref...\n');

xref_cam_info = struct();

for iX = 1:nXref
    xr = xref_targets(iX);
    contributing = [];

    for c = 1:nCams
        x_cam   = windowCenterCameras_mm.x1_mm{c}(1,:);
        x_start = min(x_cam);
        x_end   = max(x_cam);

        if xr < x_start || xr > x_end
            continue;
        end

        [~, colIdx] = min(abs(x_cam - xr));
        t      = (x_cam(colIdx) - x_start) / (x_end - x_start);
        w_hann = 0.5 * (1 - cos(pi * t));
        contributing(end+1, :) = [c, colIdx, w_hann]; %#ok<AGROW>
    end

    if isempty(contributing)
        warning('xref = %.1f mm not covered by any camera', xr);
        xref_cam_info(iX).xr           = xr;
        xref_cam_info(iX).contributing = [];
        continue;
    end

    contributing(:,3) = contributing(:,3) / sum(contributing(:,3));
    xref_cam_info(iX).xr           = xr;
    xref_cam_info(iX).contributing = contributing;

    fprintf('  xref = %7.1f mm | %d camera(s): ', xr, size(contributing,1));
    for row = 1:size(contributing,1)
        fprintf('cam%d (col%d, w=%.3f)  ', ...
            contributing(row,1), contributing(row,2), contributing(row,3));
    end
    fprintf('\n');
end

%% ===================== GET y AXIS FROM FIRST CAMERA =====================

fluctFile0 = '';
for iLoop = 1:nLoops
    f = fullfile(caseDir, loopDirs(iLoop).name, ...
        'vel_fluctuations', fluctFileName);
    if isfile(f), fluctFile0 = f; break; end
end
assert(~isempty(fluctFile0), 'No fluctuation file found');

tmp    = load(fluctFile0, 'fluctuations');
y_axis = windowCenterCameras_mm.x2_mm{1}(:,1);
clear tmp;

fprintf('\ny_axis: %d points, range [%.1f, %.1f] mm\n\n', ...
    numel(y_axis), min(y_axis), max(y_axis));

%% ===================== LOAD RMS REFERENCE FOR LU & WILLMARTH =====================

fprintf('Loading RMS reference for Lu & Willmarth hole criterion...\n');

Mrms       = load(mergedRmsFile, 'U_rms', 'V_rms', 'worldX', 'worldY');
urms_world = double(Mrms.U_rms);
vrms_world = double(Mrms.V_rms);
worldX_rms = double(Mrms.worldX);
worldY_rms = double(Mrms.worldY);

ref_LW = cell(nXref, 1);

for iX = 1:nXref
    contrib = xref_cam_info(iX).contributing;
    if isempty(contrib), continue; end

    xr          = xref_targets(iX);
    x_world_vec = worldX_rms(1,:);
    [~, colW]   = min(abs(x_world_vec - xr));

    urms_col    = urms_world(:, colW);
    vrms_col    = vrms_world(:, colW);
    y_world_col = worldY_rms(:, colW);

    [~, iDom]  = max(contrib(:,3));
    domCam     = contrib(iDom, 1);
    y_cam_col  = windowCenterCameras_mm.x2_mm{domCam}(:,1);

    urms_interp = interp1(y_world_col, urms_col, y_cam_col, 'linear', 'extrap');
    vrms_interp = interp1(y_world_col, vrms_col, y_cam_col, 'linear', 'extrap');

    urms_interp(~isfinite(urms_interp)) = eps;
    vrms_interp(~isfinite(vrms_interp)) = eps;

    ref_LW{iX} = max(urms_interp .* vrms_interp, eps);

    fprintf('  xref = %.1f mm | ref range: [%.4f, %.4f] m^2/s^2\n', ...
        xr, min(ref_LW{iX}), max(ref_LW{iX}));
end

%% ===================== PREALLOCATE GLOBAL ACCUMULATORS =====================

acc_u2    = cell(nXref, 1);
acc_u3    = cell(nXref, 1);
acc_u4    = cell(nXref, 1);
acc_v2    = cell(nXref, 1);
acc_uv    = cell(nXref, 1);
acc_uv2   = cell(nXref, 1);
acc_Q1    = cell(nXref, 1);
acc_Q2    = cell(nXref, 1);
acc_Q3    = cell(nXref, 1);
acc_Q4    = cell(nXref, 1);
acc_RS_Q1 = cell(nXref, 1);
acc_RS_Q2 = cell(nXref, 1);
acc_RS_Q3 = cell(nXref, 1);
acc_RS_Q4 = cell(nXref, 1);

for iX = 1:nXref
    contrib = xref_cam_info(iX).contributing;
    if isempty(contrib), continue; end

    [~, iDom] = max(contrib(:,3));
    domCam    = contrib(iDom, 1);
    Ny_c      = size(windowCenterCameras_mm.x1_mm{domCam}, 1);

    acc_u2{iX}    = zeros(Ny_c, 1);
    acc_u3{iX}    = zeros(Ny_c, 1);
    acc_u4{iX}    = zeros(Ny_c, 1);
    acc_v2{iX}    = zeros(Ny_c, 1);
    acc_uv{iX}    = zeros(Ny_c, 1);
    acc_uv2{iX}   = zeros(Ny_c, 1);
    acc_Q1{iX}    = zeros(Ny_c, 1);
    acc_Q2{iX}    = zeros(Ny_c, 1);
    acc_Q3{iX}    = zeros(Ny_c, 1);
    acc_Q4{iX}    = zeros(Ny_c, 1);
    acc_RS_Q1{iX} = zeros(Ny_c, nH);
    acc_RS_Q2{iX} = zeros(Ny_c, nH);
    acc_RS_Q3{iX} = zeros(Ny_c, nH);
    acc_RS_Q4{iX} = zeros(Ny_c, nH);
end

totalFrames = 0;

%% ===================== MAIN ACCUMULATION LOOP =====================

tTotal = tic;

for iLoop = 1:nLoops

    loopName  = loopDirs(iLoop).name;
    fluctFile = fullfile(caseDir, loopName, 'vel_fluctuations', fluctFileName);

    fprintf('[%02d/%02d] %s\n', iLoop, nLoops, loopName);

    if ~isfile(fluctFile)
        warning('  Fluctuation file not found — skipping');
        continue;
    end

    % --- Load and pre-extract ---
    data    = load(fluctFile, 'fluctuations');
    nFrames = size(data.fluctuations.u_prime, 1);
    fprintf('  %d frames\n', nFrames);

    u_fields = cell(1, nCams);
    v_fields = cell(1, nCams);

    for c = 1:nCams
        [Ny_c, Nx_c] = size(data.fluctuations.u_prime{1,c});
        u_arr = zeros(Ny_c, Nx_c, nFrames, 'double');
        v_arr = zeros(Ny_c, Nx_c, nFrames, 'double');
        for fr = 1:nFrames
            u_tmp = double(data.fluctuations.u_prime{fr,c});
            v_tmp = double(data.fluctuations.v_prime{fr,c});
            u_tmp(~isfinite(u_tmp)) = 0;
            v_tmp(~isfinite(v_tmp)) = 0;
            u_arr(:,:,fr) = u_tmp;
            v_arr(:,:,fr) = v_tmp;
        end
        u_fields{c} = u_arr;
        v_fields{c} = v_arr;
    end
    clear data u_arr v_arr;

    % --- Accumulate per xref ---
    for iX = 1:nXref

        contrib = xref_cam_info(iX).contributing;
        if isempty(contrib), continue; end

        % Fixed Lu & Willmarth reference
        ref = ref_LW{iX};

        for fr = 1:nFrames

            % Weighted column blend
            u_col = zeros(size(acc_u2{iX}));
            v_col = zeros(size(acc_u2{iX}));

            for row = 1:size(contrib, 1)
                c      = contrib(row, 1);
                colIdx = contrib(row, 2);
                w      = contrib(row, 3);
                u_col  = u_col + w * u_fields{c}(:, colIdx, fr);
                v_col  = v_col + w * v_fields{c}(:, colIdx, fr);
            end

            % Moments
            acc_u2{iX}  = acc_u2{iX}  + u_col.^2;
            acc_u3{iX}  = acc_u3{iX}  + u_col.^3;
            acc_u4{iX}  = acc_u4{iX}  + u_col.^4;
            acc_v2{iX}  = acc_v2{iX}  + v_col.^2;
            acc_uv{iX}  = acc_uv{iX}  + u_col.*v_col;
            acc_uv2{iX} = acc_uv2{iX} + (u_col.*v_col).^2;

            % Quadrant event fractions
            acc_Q1{iX} = acc_Q1{iX} + double(u_col > 0 & v_col > 0);
            acc_Q2{iX} = acc_Q2{iX} + double(u_col < 0 & v_col > 0);
            acc_Q3{iX} = acc_Q3{iX} + double(u_col < 0 & v_col < 0);
            acc_Q4{iX} = acc_Q4{iX} + double(u_col > 0 & v_col < 0);

            % Lu & Willmarth RS contributions
            uv_col = u_col .* v_col;

            for iH = 1:nH
                H        = H_vals(iH);
                holeMask = abs(uv_col) > H * ref;

                acc_RS_Q1{iX}(:,iH) = acc_RS_Q1{iX}(:,iH) + ...
                    uv_col .* double(u_col > 0 & v_col > 0 & holeMask);
                acc_RS_Q2{iX}(:,iH) = acc_RS_Q2{iX}(:,iH) + ...
                    uv_col .* double(u_col < 0 & v_col > 0 & holeMask);
                acc_RS_Q3{iX}(:,iH) = acc_RS_Q3{iX}(:,iH) + ...
                    uv_col .* double(u_col < 0 & v_col < 0 & holeMask);
                acc_RS_Q4{iX}(:,iH) = acc_RS_Q4{iX}(:,iH) + ...
                    uv_col .* double(u_col > 0 & v_col < 0 & holeMask);
            end

        end % frame loop

    end % xref loop

    totalFrames = totalFrames + nFrames;
    clear u_fields v_fields;

    fprintf('  Running total: %d frames | elapsed: %.1f min\n', ...
        totalFrames, toc(tTotal)/60);

end % loop folder

fprintf('\nTotal frames accumulated: %d\n\n', totalFrames);

%% ===================== NORMALISE AND DERIVE =====================

results = struct();

for iX = 1:nXref

    contrib = xref_cam_info(iX).contributing;
    if isempty(contrib)
        results(iX).xr    = xref_targets(iX);
        results(iX).valid = false;
        continue;
    end

    u2  = acc_u2{iX}  / totalFrames;
    u3  = acc_u3{iX}  / totalFrames;
    u4  = acc_u4{iX}  / totalFrames;
    v2  = acc_v2{iX}  / totalFrames;
    uv  = acc_uv{iX}  / totalFrames;
    uv2 = acc_uv2{iX} / totalFrames;
    Q1  = acc_Q1{iX}  / totalFrames;
    Q2  = acc_Q2{iX}  / totalFrames;
    Q3  = acc_Q3{iX}  / totalFrames;
    Q4  = acc_Q4{iX}  / totalFrames;

    urms    = sqrt(u2);
    vrms    = sqrt(v2);
    safe_u2 = max(u2, eps);
    skew    = u3 ./ safe_u2.^(3/2);
    kurt    = u4 ./ safe_u2.^2;
    RS      = -uv;

    results(iX).xr          = xref_targets(iX);
    results(iX).valid       = true;
    results(iX).y_mm        = windowCenterCameras_mm.x2_mm{...
                                  contrib(find(contrib(:,3)==max(contrib(:,3)),1),1)}(:,1);
    results(iX).urms        = urms;
    results(iX).vrms        = vrms;
    results(iX).skew        = skew;
    results(iX).kurt        = kurt;
    results(iX).RS          = RS;
    results(iX).Q1          = Q1;
    results(iX).Q2          = Q2;
    results(iX).Q3          = Q3;
    results(iX).Q4          = Q4;
    results(iX).uv2         = uv2;
    results(iX).RS_Q1       = acc_RS_Q1{iX} / totalFrames;
    results(iX).RS_Q2       = acc_RS_Q2{iX} / totalFrames;
    results(iX).RS_Q3       = acc_RS_Q3{iX} / totalFrames;
    results(iX).RS_Q4       = acc_RS_Q4{iX} / totalFrames;
    results(iX).H_vals      = H_vals;
    results(iX).ref_LW      = ref_LW{iX};

    fprintf('xref = %.1f mm — done\n', xref_targets(iX));

end

%% ===================== SAVE =====================

saveFile = fullfile(outDir, sprintf('turbstats_%s.mat', caseLabel));
save(saveFile, 'results', 'xref_targets', 'caseLabel', ...
    'totalFrames', 'H_vals', '-v7.3');

fprintf('\nSaved -> %s\n', saveFile);
fprintf('Total elapsed: %.1f min\n', toc(tTotal)/60);
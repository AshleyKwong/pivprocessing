% compute_fluctuations_SW.m
% Ashley Kwong
%
% Computes velocity fluctuations u' = u - U_mean, v' = v - V_mean
% for all positions (Pos_1 to Pos_4), saving ONE combined file at the
% baseDir level: baseDir/vel_fluctuations/fluctuations_all_frames.mat
%
% Output cell array is {nFrames x nPos} with ALL columns populated —
% directly compatible with compute_two_point_covariance_v5.m where
% each Pos_n maps to one camera (nCams = nPos).
%
% File structure:
%   baseDir/
%   ├── Pos_1/  Mean_Flow_Feilds.mat, coords.mat, 0001.mat ...
%   ├── Pos_2/  ...
%   ├── Pos_3/  ...
%   ├── Pos_4/  ...
%   └── vel_fluctuations/
%       └── fluctuations_all_frames.mat  ← single combined output

clear; clc;

%% ── USER SETTINGS ────────────────────────────────────────────────────────
baseDir         = '/iridisfs/scratch/ak1u24/Tom_aoan08_400mm/-8';
nPos            = 4; 
FORCE_RECOMPUTE = false;   % set true to overwrite existing file

% Files to ignore when scanning for instantaneous frames
IGNORE_FILES = {'Mean_Flow_Feilds.mat', 'coords.mat'};
%% ─────────────────────────────────────────────────────────────────────────

fprintf('\n=== compute_fluctuations_SW.m ===\n');
fprintf('Base dir: %s\n', baseDir);
fprintf('Positions: %d\n\n', nPos);

%% Output location — single combined file at baseDir level
fluctDir  = fullfile(baseDir, 'vel_fluctuations');
fluctFile = fullfile(fluctDir, 'fluctuations_all_frames.mat');

if ~exist(fluctDir, 'dir')
    mkdir(fluctDir);
end

if isfile(fluctFile) && ~FORCE_RECOMPUTE
    fprintf('Combined fluctuations file already exists:\n  %s\n', fluctFile);
    fprintf('Set FORCE_RECOMPUTE = true to overwrite.\n');
    return;
end

%% STEP 1: Load mean fields and discover instantaneous files per position
meanU      = cell(1, nPos);
meanV      = cell(1, nPos);
framePaths = cell(1, nPos);

for p = 1:nPos
    posDir = fullfile(baseDir, sprintf('Pos_%d', p));

    % --- Load mean fields ------------------------------------------------
    meanFile = fullfile(posDir, 'Mean_Flow_Feilds.mat');
    mf = load(meanFile, 'U_mean', 'V_mean');
    meanU{p} = double(mf.U_mean);
    meanV{p} = double(mf.V_mean);
    fprintf('Pos_%d: mean field loaded  [%d x %d]\n', p, size(mf.U_mean,1), size(mf.U_mean,2));

    % --- Discover instantaneous .mat files (exclude non-frame files) -----
    allMats = dir(fullfile(posDir, '*.mat'));
    keep    = true(size(allMats));
    for k = 1:numel(allMats)
        if ismember(allMats(k).name, IGNORE_FILES)
            keep(k) = false;
        end
    end
    allMats = allMats(keep);

    % Sort numerically by filename
    [~, sortIdx]  = sort({allMats.name});
    allMats       = allMats(sortIdx);
    framePaths{p} = fullfile(posDir, {allMats.name});
    fprintf('Pos_%d: %d instantaneous frames found\n', p, numel(framePaths{p}));
end

% All positions must have the same number of frames
nFrames_per_pos = cellfun(@numel, framePaths);
if numel(unique(nFrames_per_pos)) > 1
    warning('Frame counts differ across positions: %s', mat2str(nFrames_per_pos));
end
nFrames = nFrames_per_pos(1);
fprintf('\nFrames per position: %d\n', nFrames);

%% STEP 2: Preallocate combined {nFrames x nPos} cell arrays
fprintf('\nPreallocating {%d x %d} cell arrays...\n', nFrames, nPos);
u_prime_all = cell(nFrames, nPos);
v_prime_all = cell(nFrames, nPos);

%% STEP 3: Loop over positions, fill all columns
for p = 1:nPos
    fprintf('\n--- Pos_%d (%d/%d) ---\n', p, p, nPos);
    nF           = numel(framePaths{p});
    U_mean_p     = meanU{p};
    V_mean_p     = meanV{p};
    report_every = max(1, round(nF / 10));

    for fr = 1:nF
        inst   = load(framePaths{p}{fr}, 'U', 'V');
        u_pr   = single(double(inst.U) - U_mean_p);
        v_pr   = single(double(inst.V) - V_mean_p);

        u_prime_all{fr, p} = u_pr;
        v_prime_all{fr, p} = v_pr;

        if mod(fr, report_every) == 0
            fprintf('  Frame %d / %d (%.0f%%)\n', fr, nF, 100*fr/nF);
        end
    end
    fprintf('  Pos_%d complete.\n', p);
end

%% STEP 4: Build and save combined fluctuations struct
fluctuations.u_prime     = u_prime_all;
fluctuations.v_prime     = v_prime_all;
fluctuations.n_frames    = nFrames;
fluctuations.n_cameras   = nPos;
fluctuations.description = ['Combined velocity fluctuations for all positions. ' ...
    'u'' = u - U_mean, v'' = v - V_mean (m/s). ' ...
    'Cell array {nFrames x nPos}: col 1=Pos_1, col 2=Pos_2, col 3=Pos_3, col 4=Pos_4.'];

S = whos('fluctuations');
fprintf('\nSaving combined fluctuations (%.2f GB):\n  %s\n', S.bytes/1e9, fluctFile);
save(fluctFile, 'fluctuations', '-v7.3');
fprintf('\n=== Done ===\n');
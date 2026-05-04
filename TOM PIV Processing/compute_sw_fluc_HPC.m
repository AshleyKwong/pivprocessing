% compute_fluctuations_SW.m
% Ashley Kwong
%
% Computes velocity fluctuations u' = u - U_mean, v' = v - V_mean
% for each position independently from 3D instantaneous fields.
%
% Input per Pos_n:
%   Mean_Flow_Feilds.mat        — U_mean, V_mean, Y (as before)
%   Instantance_Flow_Feilds.mat — U, V [Ny x Nx x nFrames]
%
% Output per Pos_n:
%   Pos_n/loop=00/vel_fluctuations/fluctuations_all_frames.mat
%   {nFrames x 1} cell array — single-camera self-contained domain

clear; clc;

%% ── USER SETTINGS ────────────────────────────────────────────────────────
baseDir         = '/iridisfs/scratch/ak1u24/TOMDATAn08';
nPos            = 4;
FORCE_RECOMPUTE = false;
%% ─────────────────────────────────────────────────────────────────────────

fprintf('\n=== compute_fluctuations_SW.m ===\n');
fprintf('Base dir: %s\n', baseDir);
fprintf('Positions: %d\n\n', nPos);

%% STEP 1: Load mean fields per position
meanU     = cell(1, nPos);
meanV     = cell(1, nPos);
belowWall = cell(1, nPos);

for p = 1:nPos
    posDir   = fullfile(baseDir, sprintf('Pos_%d', p));
    meanFile = fullfile(posDir, 'mean_fields_blsweep.mat');

    mf = load(meanFile, 'U_hann_mean', 'V_hann_mean', 'worldY');
    meanU{p}     = double(mf.U_hann_mean);
    meanV{p}     = double(mf.V_hann_mean);
    belowWall{p} = mf.worldY <= 0;   % already in mm, mask wall and below

    fprintf('Pos_%d: mean field loaded  [%d x %d]\n', p, ...
        size(mf.U_hann_mean,1), size(mf.U_hann_mean,2));
end

%% STEP 2: Compute and save fluctuations per position
for p = 1:nPos
    posDir    = fullfile(baseDir, sprintf('Pos_%d', p));
    loopDir   = fullfile(posDir, 'loop=00');
    fluctDir  = fullfile(loopDir, 'vel_fluctuations');
    fluctFile = fullfile(fluctDir, 'fluctuations_all_frames.mat');

    if ~exist(loopDir,  'dir'), mkdir(loopDir);  end
    if ~exist(fluctDir, 'dir'), mkdir(fluctDir); end

    if isfile(fluctFile) && ~FORCE_RECOMPUTE
        fprintf('\nPos_%d: already exists — skipping (set FORCE_RECOMPUTE=true to overwrite)\n', p);
        continue;
    end

    fprintf('\n--- Pos_%d (%d/%d) ---\n', p, p, nPos);

    % --- Load 3D instantaneous fields ------------------------------------
    instFile = fullfile(posDir, 'Instantance_Flow_Feilds.mat');
    fprintf('  Loading %s ...\n', instFile);
    inst   = load(instFile, 'U', 'V');
    U_inst = double(inst.U);   % [Ny x Nx x nFrames]
    V_inst = double(inst.V);
    clear inst;

    nFrames = size(U_inst, 3);
    fprintf('  Loaded: [%d x %d x %d]  (%.2f GB)\n', ...
        size(U_inst,1), size(U_inst,2), nFrames, ...
        (numel(U_inst) + numel(V_inst)) * 8 / 1e9);

    U_mean_p     = meanU{p};
    V_mean_p     = meanV{p};
    below_wall_p = belowWall{p};
    report_every = max(1, round(nFrames / 10));

    % --- Preallocate ---------------------------------------------------
    u_prime_all = cell(nFrames, 1);
    v_prime_all = cell(nFrames, 1);

    % --- Compute fluctuations frame by frame ---------------------------
    for fr = 1:nFrames
        u_pr = single(U_inst(:,:,fr) - U_mean_p);
        v_pr = single(V_inst(:,:,fr) - V_mean_p);

        % NaN-mask below-wall points
        u_pr(below_wall_p) = NaN;
        v_pr(below_wall_p) = NaN;

        u_prime_all{fr} = u_pr;
        v_prime_all{fr} = v_pr;

        if mod(fr, report_every) == 0
            fprintf('  Frame %d / %d (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
        end
    end

    clear U_inst V_inst;

    % --- Build and save struct -----------------------------------------
    fluctuations.u_prime     = u_prime_all;
    fluctuations.v_prime     = v_prime_all;
    fluctuations.n_frames    = nFrames;
    fluctuations.n_cameras   = 1;
    fluctuations.loop_name   = sprintf('Pos_%d', p);
    fluctuations.description = sprintf(['Pos_%d velocity fluctuations. ' ...
        'u'' = u - U_mean, v'' = v - V_mean (m/s). ' ...
        '{nFrames x 1} — single-camera self-contained domain.'], p);

    S = whos('fluctuations');
    fprintf('  Saving (%.2f GB) -> %s\n', S.bytes/1e9, fluctFile);
    save(fluctFile, 'fluctuations', '-v7.3');
    fprintf('  Pos_%d done.\n', p);

    clear fluctuations u_prime_all v_prime_all;
end

fprintf('\n=== All positions complete ===\n');
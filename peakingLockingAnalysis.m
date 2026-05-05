% Ashley Kwong
% peakLockingAnalysis.m
% Peak locking diagnostic — two stages:
%   Stage 1: Raw pixel displacement PDF (floating range, 0.1px bins)
%            Spikes at integers = pixel peak locking
%   Stage 2: Sub-pixel distribution following Hearst & Ganapathisubramani
%            (2015) convention — mod folded to [-0.5, +0.5], stored as
%            both counts and probability
% ux and uy treated separately throughout.
% Cameras processed in parallel (parfor).
% No figures generated (HPC safe). Saves peakLocking.mat to savePath.
% Intended to run on Iridis6 via SLURM.

clear; clc;

%% ======== USER OPTIONS ================================================
savePath      = '/iridisfs/scratch/ak1u24/case2_smallerwindows/';
cameraList    = ["Cam1","Cam2","Cam3","Cam4","Cam5"];
binWidth_raw  = 0.1;    % Stage 1 bin width [px] — floating range
nBins_subpx   = 64;     % Stage 2 bins over [-0.5, +0.5]
rawRange_ux   = [-30, 5];   % generous range for streamwise (negative flow)
rawRange_uy   = [-10, 10];  % generous range for wall-normal
% ======================================================================
% Note: raw range is intentionally wide — empty bins at edges are
% trimmed during local plotting. Adjust if your displacement exceeds this.

%% DISCOVER LOOP FOLDERS (identical pattern to PIV_calibrateAndMeanUV.m)
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
end
fprintf('Found %d loop folders\n', length(totalLoops));

%% STAGE 1 BIN EDGES — fixed wide range, 0.1px bins
edges_raw_ux  = rawRange_ux(1) : binWidth_raw : rawRange_ux(2);
edges_raw_uy  = rawRange_uy(1) : binWidth_raw : rawRange_uy(2);
ctrs_raw_ux   = 0.5 * (edges_raw_ux(1:end-1) + edges_raw_ux(2:end));
ctrs_raw_uy   = 0.5 * (edges_raw_uy(1:end-1) + edges_raw_uy(2:end));

%% STAGE 2 BIN EDGES — [-0.5, +0.5], Hearst convention
edges_subpx   = linspace(-0.5, 0.5, nBins_subpx + 1);
ctrs_subpx    = 0.5 * (edges_subpx(1:end-1) + edges_subpx(2:end));

%% START PARALLEL POOL — read worker count from SLURM, fall back to detected cores
nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nWorkers) || nWorkers < 1
    nWorkers = feature('numcores');
end
if isempty(gcp('nocreate'))
    parpool('threads', nWorkers);
end
pool = gcp('nocreate');
fprintf('Parallel pool ready with %d workers.\n', pool.NumWorkers);

%% PRE-ALLOCATE CELL ARRAYS FOR parfor RESULTS
nCams = length(cameraList);

% Stage 1 — raw pixel displacement
raw_counts_ux_all = cell(1, nCams);
raw_counts_uy_all = cell(1, nCams);
raw_pdf_ux_all    = cell(1, nCams);
raw_pdf_uy_all    = cell(1, nCams);

% Stage 2 — sub-pixel [-0.5, +0.5]
subpx_counts_ux_all = cell(1, nCams);
subpx_counts_uy_all = cell(1, nCams);
subpx_prob_ux_all   = cell(1, nCams);
subpx_prob_uy_all   = cell(1, nCams);

nVec_ux_all = cell(1, nCams);
nVec_uy_all = cell(1, nCams);

%% MAIN LOOP — parfor over cameras
parfor a = 1:nCams

    camName = cameraList(a);
    fprintf('\n=== %s ===\n', camName);

    % Accumulators — raw displacements and sub-pixel parts
    raw_ux_all   = [];
    raw_uy_all   = [];
    subpx_ux_all = [];
    subpx_uy_all = [];

    % --- Loop over all loop folders ------------------------------------
    for loopNo = 1:length(totalLoops)

        loopFolder = fullfile(savePath, totalLoops(loopNo).name);

        base_dir = dir(fullfile(loopFolder, 'uncalibrated_piv', ...
            '150', char(camName), 'instantaneous', '*.mat'));
        base_dir = base_dir(~strcmp({base_dir.name}, 'coordinates.mat'));

        if isempty(base_dir)
            warning('No frames found for %s in %s — skipping', ...
                camName, totalLoops(loopNo).name);
            continue;
        end

        % --- Loop over all frames -------------------------------------
        for b = 1:length(base_dir)

            frameData = load(fullfile(base_dir(b).folder, base_dir(b).name), ...
                'piv_result');
            ux = double(frameData.piv_result(end).ux);   % streamwise px
            uy = double(frameData.piv_result(end).uy);   % wall-normal px
            frameData = [];   % free memory — clear not allowed in parfor

            % Strip NaNs
            ux = ux(~isnan(ux(:)));
            uy = uy(~isnan(uy(:)));
            ux = ux(:);
            uy = uy(:);

            % Stage 1 — raw displacements
            raw_ux_all = [raw_ux_all; ux]; %#ok<AGROW>
            raw_uy_all = [raw_uy_all; uy]; %#ok<AGROW>

            % Stage 2 — sub-pixel part folded to [-0.5, +0.5]
            % mod(u,1) gives [0,1); subtracting round() folds to [-0.5,+0.5]
            frac_ux  = mod(ux, 1);
            frac_uy  = mod(uy, 1);
            subpx_ux = frac_ux - round(frac_ux);
            subpx_uy = frac_uy - round(frac_uy);

            subpx_ux_all = [subpx_ux_all; subpx_ux]; %#ok<AGROW>
            subpx_uy_all = [subpx_uy_all; subpx_uy]; %#ok<AGROW>

        end % frames

        fprintf('  %s — Loop %d/%d done (%d frames)\n', ...
            camName, loopNo, length(totalLoops), length(base_dir));

    end % loops

    fprintf('  %s — Total vectors: %d (ux),  %d (uy)\n', ...
        camName, length(raw_ux_all), length(raw_uy_all));

    %% STAGE 1 HISTOGRAMS — raw pixel displacement
    raw_counts_ux = histcounts(raw_ux_all, edges_raw_ux, 'Normalization', 'count');
    raw_counts_uy = histcounts(raw_uy_all, edges_raw_uy, 'Normalization', 'count');
    raw_pdf_ux    = histcounts(raw_ux_all, edges_raw_ux, 'Normalization', 'pdf');
    raw_pdf_uy    = histcounts(raw_uy_all, edges_raw_uy, 'Normalization', 'pdf');

    %% STAGE 2 HISTOGRAMS — sub-pixel [-0.5, +0.5], Hearst convention
    subpx_counts_ux = histcounts(subpx_ux_all, edges_subpx, 'Normalization', 'count');
    subpx_counts_uy = histcounts(subpx_uy_all, edges_subpx, 'Normalization', 'count');
    subpx_prob_ux   = histcounts(subpx_ux_all, edges_subpx, 'Normalization', 'probability');
    subpx_prob_uy   = histcounts(subpx_uy_all, edges_subpx, 'Normalization', 'probability');

    %% STORE IN CELL ARRAYS
    raw_counts_ux_all{a} = raw_counts_ux;
    raw_counts_uy_all{a} = raw_counts_uy;
    raw_pdf_ux_all{a}    = raw_pdf_ux;
    raw_pdf_uy_all{a}    = raw_pdf_uy;

    subpx_counts_ux_all{a} = subpx_counts_ux;
    subpx_counts_uy_all{a} = subpx_counts_uy;
    subpx_prob_ux_all{a}   = subpx_prob_ux;
    subpx_prob_uy_all{a}   = subpx_prob_uy;

    nVec_ux_all{a} = length(raw_ux_all);
    nVec_uy_all{a} = length(raw_uy_all);

end % parfor cameras

%% ASSEMBLE OUTPUT STRUCT
peakLocking = struct();
for a = 1:nCams
    camName = cameraList(a);

    % Stage 1
    peakLocking.(camName).raw_counts_ux  = raw_counts_ux_all{a};
    peakLocking.(camName).raw_counts_uy  = raw_counts_uy_all{a};
    peakLocking.(camName).raw_pdf_ux     = raw_pdf_ux_all{a};
    peakLocking.(camName).raw_pdf_uy     = raw_pdf_uy_all{a};
    peakLocking.(camName).edges_raw_ux   = edges_raw_ux;
    peakLocking.(camName).edges_raw_uy   = edges_raw_uy;
    peakLocking.(camName).ctrs_raw_ux    = ctrs_raw_ux;
    peakLocking.(camName).ctrs_raw_uy    = ctrs_raw_uy;

    % Stage 2
    peakLocking.(camName).subpx_counts_ux = subpx_counts_ux_all{a};
    peakLocking.(camName).subpx_counts_uy = subpx_counts_uy_all{a};
    peakLocking.(camName).subpx_prob_ux   = subpx_prob_ux_all{a};
    peakLocking.(camName).subpx_prob_uy   = subpx_prob_uy_all{a};
    peakLocking.(camName).edges_subpx     = edges_subpx;
    peakLocking.(camName).ctrs_subpx      = ctrs_subpx;

    % Metadata
    peakLocking.(camName).n_vectors_ux = nVec_ux_all{a};
    peakLocking.(camName).n_vectors_uy = nVec_uy_all{a};
end

%% SAVE
outPath = fullfile(savePath, 'peakLocking.mat');
save(outPath, 'peakLocking', '-v7.3');
fprintf('\n✓ Saved peakLocking.mat → %s\n', outPath);
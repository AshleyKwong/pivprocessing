% Ashley Kwong
% peakLockingAnalysis.m
% Subpixel peak locking diagnostic via the modulus method.
% Pools fractional pixel displacements across all loops x 150 frames
% per camera, builds histogram, fits b(eps) = A*sin(2*pi*eps + phi).
% Saves peakLocking.mat to savePath — no figures generated (HPC safe).
% Intended to run on Iridis6 via SLURM.

clear; clc;

%% ======== USER OPTIONS ================================================
savePath   = '/iridisfs/scratch/ak1u24/case2_smallerwindows/';
cameraList = ["Cam1","Cam2","Cam3","Cam4","Cam5"];
nBins      = 64;          % histogram bins over [0, 1)
% ======================================================================

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

%% BIN EDGES
binEdges = linspace(0, 1, nBins + 1);   % nBins+1 edges → nBins bins

%% INITIALISE OUTPUT STRUCT
peakLocking = struct();

%% MAIN LOOP — per camera
for a = 1:length(cameraList)

    camName = cameraList(a);
    fprintf('\n=== %s ===\n', camName);

    % Accumulators for fractional parts across all loops x frames
    frac_ux_all = [];   % streamwise (will be large vector)
    frac_uy_all = [];   % wall-normal

    % --- Loop over all loop folders ------------------------------------
    for loopNo = 1:length(totalLoops)

        loopFolder = fullfile(savePath, totalLoops(loopNo).name);

        % File discovery — identical to calibration script
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

            frameData  = load(fullfile(base_dir(b).folder, base_dir(b).name), ...
                'piv_result');
            ux = double(frameData.piv_result(end).ux);   % streamwise px displacement
            uy = double(frameData.piv_result(end).uy);   % wall-normal px displacement
            clear frameData;

            % Fractional part — mod handles negative ux correctly in
            % MATLAB (always returns value in [0,1))
            frac_ux = mod(ux(:), 1);
            frac_uy = mod(uy(:), 1);

            % Remove NaNs (masked / invalid vectors)
            frac_ux = frac_ux(~isnan(frac_ux));
            frac_uy = frac_uy(~isnan(frac_uy));

            frac_ux_all = [frac_ux_all; frac_ux]; %#ok<AGROW>
            frac_uy_all = [frac_uy_all; frac_uy]; %#ok<AGROW>

        end % frames

        fprintf('  Loop %d/%d done (%d frames)\n', ...
            loopNo, length(totalLoops), length(base_dir));

    end % loops

    fprintf('  Total vectors pooled: %d (ux),  %d (uy)\n', ...
        length(frac_ux_all), length(frac_uy_all));

    %% HISTOGRAMS
    counts_u = histcounts(frac_ux_all, binEdges, 'Normalization', 'probability');
    counts_v = histcounts(frac_uy_all, binEdges, 'Normalization', 'probability');

    % Bin centres for fitting
    binCentres = 0.5 * (binEdges(1:end-1) + binEdges(2:end));   % [1 x nBins]

    %% FIT SINUSOID:  b(eps) = A * sin(2*pi*eps + phi)
    % Fit to (histogram - expected uniform level)
    % Under no peak locking, counts_u should be flat at 1/nBins
    uniformLevel = 1 / nBins;

    bias_u = counts_u - uniformLevel;   % deviation from flat
    bias_v = counts_v - uniformLevel;

    fitResult_u = fitSinusoid(binCentres, bias_u);
    fitResult_v = fitSinusoid(binCentres, bias_v);

    fprintf('  ux → A = %.4f px,  phi = %.4f rad\n', fitResult_u(1), fitResult_u(2));
    fprintf('  uy → A = %.4f px,  phi = %.4f rad\n', fitResult_v(1), fitResult_v(2));

    %% STORE IN OUTPUT STRUCT
    peakLocking.(camName).hist_u    = counts_u;
    peakLocking.(camName).hist_v    = counts_v;
    peakLocking.(camName).bias_u    = bias_u;
    peakLocking.(camName).bias_v    = bias_v;
    peakLocking.(camName).bin_edges = binEdges;
    peakLocking.(camName).bin_centres = binCentres;
    peakLocking.(camName).fit_u     = fitResult_u;   % [A, phi]
    peakLocking.(camName).fit_v     = fitResult_v;   % [A, phi]
    peakLocking.(camName).n_vectors_u = length(frac_ux_all);
    peakLocking.(camName).n_vectors_v = length(frac_uy_all);

end % cameras

%% SAVE
outPath = fullfile(savePath, 'peakLocking.mat');
save(outPath, 'peakLocking', '-v7.3');
fprintf('\n✓ Saved peakLocking.mat → %s\n', outPath);

%% ======== LOCAL FUNCTION ==============================================
function coeffs = fitSinusoid(x, y)
% Fits y = A * sin(2*pi*x + phi) using nonlinear least squares.
% Returns coeffs = [A, phi].
% x: bin centres [0,1), y: bias (histogram - uniform level)

    % Initial guess: amplitude from std, phase = 0
    A0   = 2 * std(y);
    phi0 = 0;

    modelFun = @(p, x) p(1) .* sin(2 * pi * x + p(2));

    opts = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8, ...
        'MaxFunEvals', 1e4, 'MaxIter', 1e4);

    try
        coeffs = lsqcurvefit(modelFun, [A0, phi0], x, y, ...
            [-1, -2*pi], [1, 2*pi], opts);
    catch
        warning('lsqcurvefit failed for this camera — returning zeros');
        coeffs = [0, 0];
    end

end
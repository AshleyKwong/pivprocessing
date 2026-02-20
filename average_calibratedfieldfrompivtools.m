% =========================================================================
%  PIV_AverageVelField.m
%
%  For each camera in a given loop, averages all piv_result mat files
%  (150 per camera) over the final PIV pass to produce mean Ux and Uy
%  fields. Saves result as:
%    <loopPath>/calibrated_piv/<nImages>/average_velfield_loop.mat
%
%  Saved variable: avg_velfield  (cell array, 1 x nCams)
%    avg_velfield{camIdx}.Ux  [nY x nX] mean streamwise velocity
%    avg_velfield{camIdx}.Uy  [nY x nX] mean wall-normal velocity
%    avg_velfield{camIdx}.camName  string
%
%  Ashley Kwong — University of Southampton, February 2026
% =========================================================================

clc; clear; close all;

%% -------------------------------------------------------------------------
%  User Inputs
% -------------------------------------------------------------------------
savedPath  = input('Enter savedPath (e.g. E:\\ProcessedPIV_case2fullpipe): ', 's');
nImages    = input('Enter number of images (e.g. 150): ');
nPasses    = input('Enter number of PIV passes (final pass will be used): ');
nCams      = input('Enter number of cameras: ');

methodChar = '';
while ~ismember(methodChar, {'i','e'})
    methodChar = lower(strtrim(input('Method: instantaneous (i) or ensemble (e)? ', 's')));
end
methodStr = struct('i','instantaneous','e','ensemble');
methodStr = methodStr.(methodChar);

fprintf('\nSettings:\n');
fprintf('  Path    : %s\n', savedPath);
fprintf('  Images  : %d\n', nImages);
fprintf('  Pass    : %d\n', nPasses);
fprintf('  Cameras : %d\n', nCams);
fprintf('  Method  : %s\n\n', methodStr);

%% -------------------------------------------------------------------------
%  Find loop folder
% -------------------------------------------------------------------------
loopDir = dir(fullfile(savedPath, 'loop=*'));
if isempty(loopDir)
    error('No loop= folder found in: %s', savedPath);
end
if numel(loopDir) > 1
    warning('Multiple loop= folders found — using first: %s', loopDir(1).name);
end
loopPath = fullfile(savedPath, loopDir(1).name);
fprintf('Loop folder: %s\n\n', loopPath);

% Output save directory: <loopPath>/calibrated_piv/<nImages>/
saveDir = fullfile(loopPath, 'calibrated_piv', num2str(nImages));
if ~exist(saveDir, 'dir')
    error('Save directory does not exist: %s', saveDir);
end

%% -------------------------------------------------------------------------
%  Loop over cameras — accumulate and average
% -------------------------------------------------------------------------
avg_velfield = cell(1, nCams);

for camIdx = 1:nCams

    camName  = sprintf('Cam%d', camIdx);
    camPath  = fullfile(loopPath, 'calibrated_piv', num2str(nImages), methodStr, camName);

    fprintf('========== %s ==========\n', camName);
    fprintf('  Looking in: %s\n', camPath);

    if ~exist(camPath, 'dir')
        warning('Camera folder not found: %s\n  Skipping.', camPath);
        continue;
    end

    % Find all piv_result mat files in this camera folder
    matFiles = dir(fullfile(camPath, '*.mat'));
    % Exclude coordinates.mat — we only want piv_result files
    isCoord  = arrayfun(@(f) strcmpi(f.name, 'coordinates.mat'), matFiles);
    matFiles = matFiles(~isCoord);

    nFiles = numel(matFiles);
    if nFiles == 0
        warning('  No .mat files found in %s. Skipping.', camPath);
        continue;
    end
    fprintf('  Found %d mat files.\n', nFiles);

    % ---- Load first file to get array dimensions -------------------------
    tmp    = load(fullfile(camPath, matFiles(1).name));
    piv    = tmp.piv_result{nPasses};
    [nY, nX] = size(piv.ux);

    sumUx  = zeros(nY, nX, 'double');
    sumUy  = zeros(nY, nX, 'double');
    nValid = zeros(nY, nX, 'double');   % track non-NaN counts per pixel

    % ---- Accumulate ------------------------------------------------------
    fprintf('  Averaging %d snapshots ', nFiles);
    for fIdx = 1:nFiles
        if mod(fIdx, 25) == 0
            fprintf('.');   % progress dots every 25 files
        end

        raw = load(fullfile(camPath, matFiles(fIdx).name));

        if ~isfield(raw, 'piv_result')
            warning('\n  piv_result not found in %s — skipping file.', matFiles(fIdx).name);
            continue;
        end

        piv_f = raw.piv_result{nPasses};

        ux_f  = double(piv_f.ux);
        uy_f  = double(piv_f.uy);

        % NaN-safe accumulation
        validMask = ~isnan(ux_f) & ~isnan(uy_f);
        ux_f(~validMask) = 0;
        uy_f(~validMask) = 0;

        sumUx  = sumUx  + ux_f;
        sumUy  = sumUy  + uy_f;
        nValid = nValid + double(validMask);
    end
    fprintf(' done.\n');

    % ---- Divide by valid count (NaN-safe mean) ---------------------------
    nValid(nValid == 0) = NaN;   % avoid division by zero
    meanUx = sumUx ./ nValid;
    meanUy = sumUy ./ nValid;

    avg_velfield{camIdx}.Ux      = single(meanUx);
    avg_velfield{camIdx}.Uy      = single(meanUy);
    avg_velfield{camIdx}.camName = camName;
    avg_velfield{camIdx}.nFiles  = nFiles;
    avg_velfield{camIdx}.nPass   = nPasses;

    fprintf('  Mean field size: [%d x %d],  NaN fraction: %.2f%%\n', ...
        nY, nX, 100*mean(isnan(meanUx(:))));
end

%% -------------------------------------------------------------------------
%  Save
% -------------------------------------------------------------------------
savePath = fullfile(saveDir, 'average_velfield_loop.mat');
fprintf('\nSaving to: %s\n', savePath);
save(savePath, 'avg_velfield', '-v7.3');
fprintf('Saved successfully.\n');

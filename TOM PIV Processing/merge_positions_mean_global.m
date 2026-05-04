%% merge_positions_mean.m
% Loads mean velocity and variance fields from multi-position PIV data,
% merges them onto a unified world grid, and saves outputs with variable
% names that bl_sweep.m expects directly.
%
% Folder structure expected:
%   baseDir\
%     Pos_1\  Mean_Flow_Feilds.mat  (X, Y, U_mean, V_mean, uu_mean, vv_mean)
%     Pos_2\  ...
%% merge_positions_mean.m
% Loads mean velocity and variance fields from multi-position PIV data,
% merges them onto a unified world grid, and saves outputs with variable
% names that bl_sweep.m expects directly.
%
% Folder structure expected:
%   baseDir\
%     Pos_1\  Mean_Flow_Feilds.mat  (X, Y, U_mean, V_mean, uu_mean, vv_mean)
%     Pos_2\  ...
%     Pos_3\  ...
%     Pos_4\  ...
%
% Output variables (saved to outputFile):
%   worldX, worldY      — unified grid [mm]
%   U_hann_mean         — merged streamwise mean velocity
%   V_hann_mean         — merged wall-normal mean velocity
%   U_rms               — sqrt(uu_merged), streamwise RMS
%   V_rms               — sqrt(vv_merged), wall-normal RMS

%% ---- USER SETTINGS -------------------------------------------------------
baseDir     = 'G:\SW 500mm Mean Flow Fields\';
nPos        = 4;               % number of positions
window_type = 'hann';          % blending window: 'tukey','hann','cosine','linear','flat'
taper_param = [];              % tukey alpha (ignored for hann)
timestamp   = datestr(now, 'yyyymmdd_HHMMSS');
outputFile  = fullfile(baseDir, sprintf('globalmerged_mean_field_%s.mat', timestamp));
%% -------------------------------------------------------------------------

%% STEP 1: Load per-position data
windowCenters.x1_mm = cell(1, nPos);
windowCenters.x2_mm = cell(1, nPos);
velocityU  = cell(1, nPos);
velocityV  = cell(1, nPos);
varianceUU = cell(1, nPos);   % uu_mean — streamwise velocity variance
varianceVV = cell(1, nPos);   % vv_mean — wall-normal velocity variance

for p = 1:nPos
    posDir   = fullfile(baseDir, sprintf('Pos_%d', p));
    meanFile = fullfile(posDir, 'Mean_Flow_Feilds.mat');   % typo in filename is intentional

    mf = load(meanFile, 'X', 'Y', 'U_mean', 'V_mean', 'uu_mean', 'vv_mean');

    % X, Y in metres — convert to mm
    windowCenters.x1_mm{p} = mf.X * 1e3;
    windowCenters.x2_mm{p} = mf.Y * 1e3;

    fprintf('Pos_%d  X:[%.2f, %.2f] mm   Y:[%.2f, %.2f] mm\n', p, ...
        min(mf.X(:))*1e3, max(mf.X(:))*1e3, ...
        min(mf.Y(:))*1e3, max(mf.Y(:))*1e3);

    velocityU{p}  = mf.U_mean;
    velocityV{p}  = mf.V_mean;
    varianceUU{p} = mf.uu_mean;
    varianceVV{p} = mf.vv_mean;

    fprintf('  Field size: %d x %d\n', size(mf.U_mean, 1), size(mf.U_mean, 2));
end

%% STEP 2: Empty masks
masks = {};

%% STEP 3: Merge mean velocity fields (grid built automatically)
fprintf('\n--- Merging U_mean / V_mean ---\n');
[worldX, worldY, U_hann_mean, V_hann_mean] = merge_cameras_python_style_mean( ...
    windowCenters, velocityU, velocityV, masks, window_type, taper_param);

%% STEP 4: Merge variance fields on the SAME grid
% Pass worldGrid so the merge function skips grid construction and uses
% exactly the same points — guarantees identical array sizes.
worldGrid.x1 = worldX;
worldGrid.x2 = worldY;

fprintf('\n--- Merging uu_mean / vv_mean (variance) ---\n');
[~, ~, uu_merged, vv_merged] = merge_cameras_python_style_mean( ...
    windowCenters, varianceUU, varianceVV, masks, window_type, taper_param, worldGrid);

% RMS = sqrt(variance); clamp any small negatives from interpolation to 0
U_rms = sqrt(max(uu_merged, 0));
V_rms = sqrt(max(vv_merged, 0));

%% STEP 4b: NaN-mask all points at or below the wall (y <= 0)
% Zero-filling below the wall biases covariance averages — NaN is safer.
% Downstream nan_mask logic in the covariance code excludes these points
% from N_pix accumulation, so they never dilute the correlation average.
below_wall = worldY <= 0;
U_hann_mean(below_wall) = NaN;
V_hann_mean(below_wall) = NaN;
U_rms(below_wall)       = NaN;
V_rms(below_wall)       = NaN;
fprintf('Below-wall points NaN-masked: %d / %d (%.1f%%)\n', ...
    sum(below_wall(:)), numel(below_wall), 100*mean(below_wall(:)));

%% STEP 5: Save — variable names match bl_sweep.m inputs exactly
fprintf('\nSaving merged fields to:\n  %s\n', outputFile);
save(outputFile, 'worldX', 'worldY', ...
    'U_hann_mean', 'V_hann_mean', ...
    'U_rms',       'V_rms',       ...
    '-v7.3');
fprintf('Done.\n');


%% STEP 6: Quick diagnostic plot (2 Figures)
figure('Name', 'Merged Fields', 'Position', [100 100 1200 800]);

subplot(2,1,1);
pcolor(worldX, worldY, U_hann_mean); shading interp; colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); title('U_{mean}  [m/s]');
axis image; clim([ 0 30]); 
xline([5866 6480 6921 7320 8000 8573], "k", LineWidth=1.5) 

subplot(2,1,2);
pcolor(worldX, worldY, V_hann_mean); shading interp; colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); title('V_{mean}  [m/s]');
axis image; clim([ -2 5]); 
axis equal tight;
sgtitle(sprintf('Merged %d-position fields  |  %s', nPos, timestamp));


figure(); 
subplot(2,1,1);
pcolor(worldX, worldY, U_rms); shading interp; colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); title('U_{rms} = \surd\overline{u^2}  [m/s]');
axis image; clim([ 0 3]); 

subplot(2,1,2);
pcolor(worldX, worldY, V_rms); shading interp; colorbar;
xlabel('X (mm)'); ylabel('Y (mm)'); title('V_{rms} = \surd\overline{v^2}  [m/s]');
axis image; clim([ 0 2]); 

sgtitle(sprintf('Merged %d-position fields  |  %s', nPos, timestamp));
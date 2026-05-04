% wall_detection_diagnostic.m
% Ashley Kwong
%
% Standalone local script to:
%   1. Detect wall position per column from U_hann_mean == 0 mask boundary
%   2. Smooth y_wall_raw with a median filter (robust, no camera assumptions)
%   3. Compute column-by-column correction — no clamping, camera 1 included
%   4. Save smoothed correction vector for upload to IRIDIS
%
% Files needed locally (download from IRIDIS):
%   - merged_meanUV_*.mat
%   - windowCenterCameras_mm.mat  (from any loop folder)
%
% NOTE: nothing in this script modifies any HPC data.

clear; clc; close all;

%% ===== USER INPUTS =====================================================

mergedMeanFile    = 'C:\Users\ak1u24\Downloads\merged_meanUV_1loops_20260409_233423.mat';
windowCentresFile = 'C:\Users\ak1u24\Downloads\windowCenterCameras_mm.mat';

y_ref_nominal    = 2.36;    % mm — nominal wall-normal reference height for overlay

% Median filter half-width in mm — controls smoothing of floor profile.
% Wider = smoother, less sensitive to per-column noise.
% Start with 50 mm; increase if the smoothed curve still looks noisy.
smoothing_halfwidth_mm = 50;

%% ===== LOAD =============================================================

fprintf('Loading mean field...\n');
d = load(mergedMeanFile);
if ~isfield(d, 'U_hann_mean')
    error('Expected U_hann_mean. Found: %s', strjoin(fieldnames(d), ', '));
end
U_mean = double(d.U_hann_mean);
x_axis = d.worldX(1, :);    % 1 x Nx
y_axis = d.worldY(:, 1);    % Ny x 1
Nx     = length(x_axis);
dx     = mean(diff(x_axis));

fprintf('Loading window centres...\n');
wc    = load(windowCentresFile, 'windowCenterCameras_mm');
wc    = wc.windowCenterCameras_mm;
nCams = length(wc.x1_mm);

fprintf('Grid: [%d x %d] | X: [%.1f, %.1f] mm | Y: [%.1f, %.1f] mm\n', ...
        size(U_mean,1), size(U_mean,2), ...
        min(x_axis), max(x_axis), min(y_axis), max(y_axis));
fprintf('Cameras: %d | dx = %.3f mm\n\n', nCams, dx);

%% ===== WALL DETECTION: topmost U==0 row per column =====================
% Where U_hann_mean == 0 the region is masked — this is the floor.
% Take the highest y value among zero rows: last zero row before flow starts.
% No y >= 0 clamping — camera 1 may have a genuine negative floor offset.

y_wall_raw = NaN(1, Nx);
for ix = 1:Nx
    col       = U_mean(:, ix);
    zero_rows = find(col == 0);
    if isempty(zero_rows); continue; end
    y_wall_raw(ix) = max(y_axis(zero_rows));
end

fprintf('Raw detection range:  [%.4f, %.4f] mm\n', ...
        min(y_wall_raw,[], 'omitnan'), max(y_wall_raw,[], 'omitnan'));

%% ===== SMOOTH y_wall_raw ===============================================
% Median filter in physical units — window width = 2 * halfwidth / dx columns.
% Applied only over the valid (non-NaN) range; NaN columns are filled by
% nearest-neighbour extrapolation so the filter has no edge gaps.

win_cols = 2 * round(smoothing_halfwidth_mm / dx) + 1;   % always odd
fprintf('Smoothing window: %d columns (%.1f mm)\n', win_cols, win_cols * dx);

% Fill NaNs with nearest valid value before filtering (avoids edge dropout)
y_filled = y_wall_raw;
nan_idx  = find(isnan(y_filled));
ok_idx   = find(~isnan(y_filled));
if ~isempty(nan_idx) && ~isempty(ok_idx)
    y_filled(nan_idx) = interp1(ok_idx, y_filled(ok_idx), nan_idx, 'nearest', 'extrap');
end

y_wall_smooth = medfilt1(y_filled, win_cols);

fprintf('Smoothed range:       [%.4f, %.4f] mm\n\n', ...
        min(y_wall_smooth), max(y_wall_smooth));

%% ===== PER-CAMERA SUMMARY (informational only) =========================
% Print the median smoothed offset within each camera's x-domain.
% This is for reference — the actual correction is column-by-column.

cam_colours = lines(nCams);
fprintf('%-8s  %-14s  %-14s  %-16s\n', ...
        'Camera', 'x_min (mm)', 'x_max (mm)', 'median offset (mm)');
fprintf('%s\n', repmat('-', 1, 56));
x_cam_min = zeros(1, nCams);
x_cam_max = zeros(1, nCams);
for c = 1:nCams
    x_cam_min(c) = min(wc.x1_mm{c}(:));
    x_cam_max(c) = max(wc.x1_mm{c}(:));
    in_cam = x_axis >= x_cam_min(c) & x_axis <= x_cam_max(c);
    fprintf('  Cam%-4d  %-14.2f  %-14.2f  %.4f\n', ...
            c, x_cam_min(c), x_cam_max(c), median(y_wall_smooth(in_cam), 'omitnan'));
end
fprintf('\nCorrection is applied column-by-column from y_wall_smooth.\n');
fprintf('This will be subtracted from x2_mm{cam} in the covariance script.\n\n');

%% ===== SAVE CORRECTION VECTOR ==========================================

save('wall_correction_smooth.mat', 'x_axis', 'y_wall_raw', 'y_wall_smooth', ...
     'x_cam_min', 'x_cam_max', 'nCams', 'smoothing_halfwidth_mm');
fprintf('Saved wall_correction_smooth.mat\n');
fprintf('Upload this file to IRIDIS alongside the covariance script.\n\n');

%% ===== FIGURE 1: Raw detection + smoothed correction ===================

figure('Name','Wall Detection & Smoothed Correction','Position',[50 50 1300 540]);

subplot(2,1,1);
plot(x_axis, y_wall_raw,    'b.', 'MarkerSize', 3, 'DisplayName', 'Raw: topmost U=0'); hold on;
plot(x_axis, y_wall_smooth, 'r-', 'LineWidth', 2, 'DisplayName', ...
     sprintf('Smoothed (%.0f mm window)', win_cols*dx));
yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', 'y=0 datum');
xline(0, 'k:',  'LineWidth', 1.2, 'HandleVisibility', 'off');
for c = 1:nCams
    xline(x_cam_min(c), '--', 'Color', cam_colours(c,:), 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Cam%d boundary', c));
end
xlabel('x (mm)'); ylabel('y_{wall} (mm)');
title('Wall detection: raw mask boundary + smoothed correction');
legend Location best; grid on;

subplot(2,1,2);
plot(x_axis, y_wall_smooth, 'k-', 'LineWidth', 2); hold on;
yline(0, 'b--', 'LineWidth', 1.2);
xline(0, 'k:',  'LineWidth', 1.2);
for c = 1:nCams
    xline(x_cam_min(c), '--', 'Color', cam_colours(c,:), 'LineWidth', 1.2);
end
xlabel('x (mm)'); ylabel('Correction subtracted from x2\_mm (mm)');
title('Column-by-column smoothed correction — negative = floor below datum (Cam1)');
grid on;

sgtitle('Wall correction diagnostic — verify before uploading to IRIDIS');

%% ===== FIGURE 2: Mean U field with wall overlay ========================

figure('Name','Mean U + wall overlay','Position',[50 650 1400 420]);
imagesc(x_axis, y_axis, U_mean); set(gca, 'YDir', 'normal');
colormap(turbo); colorbar;
valid_vals = U_mean(U_mean > 0);
if ~isempty(valid_vals); clim([0, prctile(valid_vals, 95)]); end
hold on;
plot(x_axis, y_wall_raw,    'w.', 'MarkerSize', 2,   'DisplayName', 'Raw boundary');
plot(x_axis, y_wall_smooth, 'r-', 'LineWidth', 2,    'DisplayName', 'Smoothed correction');
yline(y_ref_nominal, 'g--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('y_{ref} nominal = %.1f mm', y_ref_nominal));
for c = 1:nCams
    xline(x_cam_min(c), '--', 'Color', cam_colours(c,:), 'LineWidth', 1.2, ...
          'DisplayName', sprintf('Cam%d', c));
end
xlabel('x (mm)'); ylabel('y (mm)');
title('U_{hann,mean} with smoothed wall correction overlay');
legend Location northeast; grid on;

%% ===== FIGURE 3: Corrected y_ref across FOV ============================
% figure('Name','Corrected y_ref','Position',[50 1130 1100 360]);
figure(3); 
x_q    = x_axis;
corr_q = y_wall_smooth;
plot(x_q, repmat(y_ref_nominal, size(x_q)), 'b--', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Nominal y_{ref} = %.1f mm', y_ref_nominal)); hold on;
plot(x_q, y_ref_nominal + corr_q, 'r-', 'LineWidth', 2, ...
     'DisplayName', 'Corrected y_{ref} (nominal + floor rise)');
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
for c = 1:nCams
    xline(x_cam_min(c), '--', 'Color', cam_colours(c,:), 'LineWidth', 1.2, ...
          'HandleVisibility', 'off');
end
xlabel('x (mm)'); ylabel('y_{ref} (mm)');
title('Effect of smoothed column-by-column wall correction on y_{ref}');
legend Location best; grid on;

fprintf('Max correction magnitude: %.4f mm\n', max(abs(y_wall_smooth)));
fprintf('Done. Check figures then upload wall_correction_smooth.mat to IRIDIS.\n');

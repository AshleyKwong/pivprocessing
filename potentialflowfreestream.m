% load in the .mat file : load("C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\BGroup Multielement Vortex Panel Code\panelmethod_velocity_chord_dimensional.mat")
% this will generate the following struct: potential_flowsoln with fields
% u, x_m and y_m where x and y are in meters. 

%% Compare Hanning U_mean with Potential Flow Solution
% Load potential flow solution (if not already in workspace)
if ~exist('potential_flowsoln', 'var')
    load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\BGroup Multielement Vortex Panel Code\panelmethod_velocity_chord_dimensional_CASE1.mat');
end 
% Extract potential flow data
pf_u = potential_flowsoln.u;      % Potential flow U velocity
pf_x = potential_flowsoln.x_m;    % X coordinates in meters
pf_y = potential_flowsoln.y_m;    % Y coordinates in meters

% Extract Hanning mean solution (from your merged data) - need to run piv
% stitching first. 
mean_field = load("D:\merge_instantaneousavg_20260413_103619\merged_meanUV_14loops_20260413_103619.mat"); 
hanning_u = mean_field.U_hann_mean;          % From your camera merging-->
hanning_x = mean_field.worldX;               % In mm from merging
hanning_y = mean_field.worldY;               % In mm from merging

% Convert worldX, worldY from mm to meters to match potential flow
hanning_x_m = hanning_x * 1e-3;
hanning_y_m = hanning_y * 1e-3;

%% Interpolate Potential Flow onto Hanning Grid
fprintf('\n=== Interpolating Potential Flow Solution ===\n');

% Check if potential flow is meshgrid or vectors
if isvector(pf_x) && isvector(pf_y)
    % If vectors, create meshgrid
    [pf_X, pf_Y] = meshgrid(pf_x, pf_y);
else
    % Already meshgrid
    pf_X = pf_x;
    pf_Y = pf_y;
end

% Interpolate potential flow U onto Hanning grid
pf_u_interp = interp2(pf_X, pf_Y, pf_u, hanning_x_m, hanning_y_m, 'linear', NaN);

fprintf('✓ Interpolated potential flow onto %d × %d Hanning grid\n', ...
    size(hanning_u, 1), size(hanning_u, 2));

%% Compute Difference and Identify Freestream
% Calculate absolute difference
u_diff = (hanning_u - pf_u_interp);

% Calculate relative error (normalized by potential flow)
u_relative_error = abs(u_diff) ./ abs(pf_u_interp) * 100;  % Always positive;  % Percentage

% Define freestream criterion (adjust threshold as needed)
freestream_threshold = 1;  % 5% relative error
freestream_mask = u_relative_error < freestream_threshold;

fprintf('\n=== Freestream Region Analysis ===\n');
fprintf('Threshold: %.1f%% relative error\n', freestream_threshold);
fprintf('Freestream points: %d / %d (%.1f%%)\n', ...
    sum(freestream_mask(:), 'omitnan'), numel(freestream_mask), ...
    100*sum(freestream_mask(:), 'omitnan')/numel(freestream_mask));

% Get Y-coordinate range where freestream is detected
[row_indices, ~] = find(freestream_mask);
if ~isempty(row_indices)
    y_freestream_min = min(hanning_y_m(row_indices, 1));
    y_freestream_max = max(hanning_y_m(row_indices, 1));
    fprintf('Freestream Y-range: %.4f to %.4f m\n', y_freestream_min, y_freestream_max);
end

%% Visualization: 4-Panel Comparison
close all; 
figure('Position', [100, 100, 1400, 1000], 'Visible', 'on');

% Panel 1: Potential Flow U
subplot(2,2,1);
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), pf_u_interp);
axis image; colorbar;
colormap(parula); clim([0 30])
title('Potential Flow U (interpolated) [m/s]');
xlabel('X [m]'); ylabel('Y [m]'); clim([0 30]);
set(gca, 'YDir', 'normal');

% Panel 2: Hanning Mean U
subplot(2,2,2);
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), hanning_u);
axis image; colorbar;
colormap(parula);
clim([0 30]); 
title('Hanning Mean U (PIV) [m/s]');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'YDir', 'normal');

% Panel 3: Absolute Difference
subplot(2,2,3);
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), u_diff);
axis image; colorbar; clim([-5 5])
colormap(hot);
title('|U_{PIV} - U_{potential}| [m/s]');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'YDir', 'normal');

% Panel 4: Freestream Region Overlay
subplot(2,2,4);
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), u_relative_error);
hold on;
% Overlay freestream contour
contour(hanning_x_m, hanning_y_m, u_relative_error, [freestream_threshold freestream_threshold], ...
    'LineColor', 'r', 'LineWidth', 2);
axis image; colorbar;
colormap(jet);
clim([0 20]);  % 0-20% error range
title(sprintf('Relative Error (< %.0f%% = Freestream)', freestream_threshold));
xlabel('X [m]'); ylabel('Y [m]');
legend('Freestream Boundary', 'Location', 'best');
set(gca, 'YDir', 'normal');

sgtitle('Hanning Mean vs Potential Flow Comparison', 'FontSize', 16);
%%
% Panel 3: Absolute Difference
figure(1); 
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), u_diff);
axis image; colorbar;
colormap(hot);
title('|U_{PIV} - U_{potential}| [m/s]');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'YDir', 'normal');
% Panel 4: Freestream Region Overlay
figure(2)
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), u_relative_error);
hold on;
% Overlay freestream contour
contour(hanning_x_m, hanning_y_m, u_relative_error, [1 1], ...
    'LineColor', 'r', 'LineWidth', 2);
axis image; colorbar;
colormap(parula);
clim([0 5]);  % 0-5% error range
title(sprintf('Relative Error (< %.0f%% = Freestream)', freestream_threshold));
xlabel('X [m]'); ylabel('Y [m]');
legend('Freestream Boundary', 'Location', 'best');
set(gca, 'YDir', 'normal');

%% Extract Freestream Velocity Statistics
% Mask NaN values
valid_mask = ~isnan(hanning_u) & ~isnan(pf_u_interp) & freestream_mask;

if any(valid_mask(:))
    u_freestream_piv = hanning_u(valid_mask);
    u_freestream_potential = pf_u_interp(valid_mask);
    
    fprintf('\n=== Freestream Velocity Statistics ===\n');
    fprintf('PIV Freestream U:        %.3f ± %.3f m/s\n', ...
        mean(u_freestream_piv), std(u_freestream_piv));
    fprintf('Potential Freestream U:  %.3f ± %.3f m/s\n', ...
        mean(u_freestream_potential), std(u_freestream_potential));
    fprintf('Mean difference:         %.3f m/s (%.2f%%)\n', ...
        mean(u_freestream_piv - u_freestream_potential), ...
        100*mean(abs(u_freestream_piv - u_freestream_potential))./mean(u_freestream_potential));
end

%% Optional: Line Profile Comparison at Specific X
% Choose X location to compare vertical profiles
x_profile = 1;  % meters, adjust as needed
[~, x_idx] = min(abs(hanning_x_m(1,:) - x_profile));

figure('Position', [100, 100, 800, 600], 'Visible', 'on');
plot(pf_u_interp(:, x_idx), hanning_y_m(:, x_idx), 'b-', 'LineWidth', 2, 'DisplayName', 'Potential Flow');
hold on;
plot(hanning_u(:, x_idx), hanning_y_m(:, x_idx), 'r--', 'LineWidth', 2, 'DisplayName', 'PIV (Hanning)');
xlabel('U [m/s]');
ylabel('Y [m]');
title(sprintf('Velocity Profile at X = %.3f m', x_profile));
legend('Location', 'best');
grid on;
%% ── Extract Freestream Boundary Line ─────────────────────────────────────
fprintf('\n=== Extracting Freestream Boundary Line ===\n');

Nx = size(hanning_u, 2);
freestream_line = nan(Nx, 3);   % [x_m, y_m, U_mean]
N_consec = 3;

for col = 1:Nx

    x_col   = hanning_x_m(1, col);
    y_col   = hanning_y_m(:, col);
    err_col = u_relative_error(:, col);
    u_col   = hanning_u(:, col);

    below   = err_col < freestream_threshold & ~isnan(err_col);

    % Find first run of N_consec consecutive true values
    first_fs_idx = [];
    all_idxs = find(below);
    for k = 1 : length(all_idxs) - (N_consec - 1)
        if all_idxs(k + N_consec - 1) - all_idxs(k) == N_consec - 1
            first_fs_idx = all_idxs(k);
            break
        end
    end

    if isempty(first_fs_idx) || first_fs_idx == 1
        continue
    end

    % Interpolate precisely at the threshold crossing
    e1 = err_col(first_fs_idx - 1);
    e2 = err_col(first_fs_idx);
    y1 = y_col(first_fs_idx - 1);
    y2 = y_col(first_fs_idx);

    if ~isnan(e1) && ~isnan(e2) && e1 ~= e2
        y_threshold    = interp1([e1, e2], [y1, y2], freestream_threshold);
        u_at_threshold = interp1([y1, y2], [u_col(first_fs_idx-1), u_col(first_fs_idx)], y_threshold);
    else
        y_threshold    = y2;
        u_at_threshold = u_col(first_fs_idx);
    end

    freestream_line(col, 1) = x_col;
    freestream_line(col, 2) = y_threshold;
    freestream_line(col, 3) = u_at_threshold;
end

% Remove columns with no valid crossing
valid_rows      = ~isnan(freestream_line(:, 1));
freestream_line = freestream_line(valid_rows, :);

fprintf('  Extracted %d valid freestream boundary points\n', size(freestream_line, 1));
fprintf('  x range:     [%.4f, %.4f] m\n', min(freestream_line(:,1)), max(freestream_line(:,1)));
fprintf('  y range:     [%.4f, %.4f] m\n', min(freestream_line(:,2)), max(freestream_line(:,2)));
fprintf('  U_inf range: [%.3f, %.3f] m/s\n', min(freestream_line(:,3)), max(freestream_line(:,3)));


%% ── Quick sanity plot ────────────────────────────────────────────────────
figure('Position', [100, 100, 1200, 500]);
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), u_relative_error);
axis xy equal tight;
colormap(parula); clim([0, 10]); colorbar;
hold on;
plot(freestream_line(:,1), freestream_line(:,2), 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Freestream edge (%.0f%% threshold)', freestream_threshold));
xlabel('X [m]'); ylabel('Y [m]');
title('Relative Error Field with Extracted Freestream Boundary');
legend('Location', 'best'); grid off;

% %% ── Save ─────────────────────────────────────────────────────────────────
% save('freestream_line.mat', 'freestream_line');
% fprintf('\n✓ Saved → freestream_line.mat  [Nx3: x_m, y_m, U_inf]\n');
% 
% %% Save Results
% comparison_results = struct();
% comparison_results.hanning_u = hanning_u;
% comparison_results.potential_u_interp = pf_u_interp;
% comparison_results.difference = u_diff;
% comparison_results.relative_error = u_relative_error;
% comparison_results.freestream_mask = freestream_mask;
% comparison_results.x_m = hanning_x_m;
% comparison_results.y_m = hanning_y_m;
% comparison_results.threshold = freestream_threshold;
% 
% save(fullfile(savePath, 'hanning_vs_potential_comparison.mat'), 'comparison_results', '-v7.3');
% fprintf('\n✓ Saved comparison results to: hanning_vs_potential_comparison.mat\n');

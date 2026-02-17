% load in the .mat file : load("C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\BGroup Multielement Vortex Panel Code\panelmethod_velocity_chord_dimensional.mat")
% this will generate the following struct: potential_flowsoln with fields
% u, x_m and y_m where x and y are in meters. 

%% Compare Hanning U_mean with Potential Flow Solution
% Load potential flow solution (if not already in workspace)
if ~exist('potential_flowsoln', 'var')
    load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\BGroup Multielement Vortex Panel Code\panelmethod_velocity_chord_dimensional.mat');
end 
% Extract potential flow data
pf_u = potential_flowsoln.u;      % Potential flow U velocity
pf_x = potential_flowsoln.x_m;    % X coordinates in meters
pf_y = potential_flowsoln.y_m;    % Y coordinates in meters

% Extract Hanning mean solution (from your merged data) - need to run piv
% stitching first. 
hanning_u = U_hann_mean;          % From your camera merging-->
hanning_x = worldX;               % In mm from merging
hanning_y = worldY;               % In mm from merging

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
u_diff = abs(hanning_u - pf_u_interp);

% Calculate relative error (normalized by potential flow)
u_relative_error = u_diff ./ abs(pf_u_interp) * 100;  % Percentage

% Define freestream criterion (adjust threshold as needed)
freestream_threshold = 10;  % 5% relative error
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
title('Hanning Mean U (PIV) [m/s]');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'YDir', 'normal');

% Panel 3: Absolute Difference
subplot(2,2,3);
imagesc(hanning_x_m(1,:), hanning_y_m(:,1), u_diff);
axis image; colorbar;
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
contour(hanning_x_m, hanning_y_m, u_relative_error, [freestream_threshold freestream_threshold], ...
    'LineColor', 'r', 'LineWidth', 2);
axis image; colorbar;
colormap(jet);
clim([0 20]);  % 0-20% error range
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
x_profile = 0.5;  % meters, adjust as needed
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

%% Save Results
comparison_results = struct();
comparison_results.hanning_u = hanning_u;
comparison_results.potential_u_interp = pf_u_interp;
comparison_results.difference = u_diff;
comparison_results.relative_error = u_relative_error;
comparison_results.freestream_mask = freestream_mask;
comparison_results.x_m = hanning_x_m;
comparison_results.y_m = hanning_y_m;
comparison_results.threshold = freestream_threshold;

save(fullfile(savePath, 'hanning_vs_potential_comparison.mat'), 'comparison_results', '-v7.3');
fprintf('\n✓ Saved comparison results to: hanning_vs_potential_comparison.mat\n');

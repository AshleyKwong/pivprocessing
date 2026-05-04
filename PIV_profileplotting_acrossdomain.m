% Ashley Kwong
% PIV Boundary Layer Profile Comparison — Multiple x-locations, no HWA
clear; clc; close all;
set(groot, 'defaultAxesFontName', 'Cambria Math');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontName',  'Cambria Math');
set(groot, 'defaultTextFontSize',  16);

%% 1. File paths & user inputs
mean_field     = 'G:\Y235_AOAN04_AOAFN06_SmallerWindows\merge_instantaneousavg_20260429_083207\merged_meanUV_14loops_20260429_083207.mat';
Cf_ofi         = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\OFI RESULTS\OFI_Cf_results_20260426_165929.mat';
blSweep        = load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\16x96init\blSweep_20260430_084136.mat').blSweep;

case_name_query = 'Case 2';
x_locations     = [54, 550, 1000];   % <-- set your x-locations (mm)
y0_wall         = 0;
nu_air          = 1.5e-5;
n_gray          = 6;   % number of near-wall points to gray out

%% 2. Load fields
U_mean     = load(mean_field).U_hann_mean;
worldX_mean = load(mean_field).worldX;
worldY_mean = load(mean_field).worldY;

mag_factor = (max(worldX_mean(1,:)) - min(worldX_mean(1,:))) / size(worldX_mean, 2);

%% 3. Cf matching
Cf_data_raw = load(Cf_ofi);
case_idx    = find(strcmp({Cf_data_raw.results.caseName}, case_name_query), 1);
if isempty(case_idx)
    error('Case "%s" not found. Available: %s', ...
        case_name_query, strjoin({Cf_data_raw.results.caseName}, ', '));
end
Cf_data = Cf_data_raw.results(case_idx);
fprintf('Using Cf data for: %s (index %d)\n', Cf_data.caseName, case_idx);

%% 4. Coordinate vectors
x_vec = worldX_mean(1, :);
y_vec = worldY_mean(:, 1);
fprintf('X range: [%.2f, %.2f] mm (%d pts)\n', min(x_vec), max(x_vec), length(x_vec));
fprintf('Y range: [%.2f, %.2f] mm (%d pts)\n', min(y_vec), max(y_vec), length(y_vec));

%% 5. Extract profiles
profileData = struct();
for i = 1:length(x_locations)
    x_target = x_locations(i);

    % --- Cf / utau at this x -----------------------------------------
    x_global_query = x_target + 7200;   % local PIV mm → global tunnel mm
    [~, matchIdx]  = min(abs(Cf_data.x_centers - x_global_query));
    Cf_half        = Cf_data.Cf(matchIdx) / 2;
    utau_piv       = sqrt(Cf_half * Cf_data.U_inf_local(matchIdx)^2);

    % --- Column averaging --------------------------------------------
    cols = find(x_vec >= x_target - mag_factor & x_vec <= x_target + mag_factor);
    if isempty(cols)
        fprintf('⚠ No data near x = %.2f mm — skipping.\n', x_target);
        continue
    end
    x_actual = mean(x_vec(cols));

    U_profile = mean(U_mean(:, cols), 2, 'omitnan');
    cum_z     = double(y_vec - y0_wall);
    Um        = double(U_profile);

    % --- NaN removal & sort ------------------------------------------
    valid_mask     = ~isnan(Um);
    Um             = Um(valid_mask);
    cum_z          = cum_z(valid_mask);
    [cum_z, sidx]  = sort(cum_z, 'ascend');
    Um             = Um(sidx);

    % --- Wall detection ----------------------------------------------
    U_inf       = max(Um);
    wall_thresh = 0.0005 * U_inf;
    valid_above = find(Um > wall_thresh, 1, 'first');
    if isempty(valid_above)
        warning('No valid above-wall points at x = %.2f mm — skipping.', x_target);
        continue
    end
    y_offset  = cum_z(valid_above);
    cum_z_corr = cum_z - y_offset + 0.2;

    % --- Trim first 3 (reflections) then store -----------------------
    trim_idx       = valid_above + 3;
    cum_z_plot     = cum_z_corr(trim_idx:end);
    Um_plot        = Um(trim_idx:end);

    % --- delta99 -----------------------------------------------------
    delta99 = interp1(blSweep.x_mm, blSweep.delta99_hybrid_mm, x_target, 'linear', 'extrap');

    % --- Store -------------------------------------------------------
    profileData(i).x_target    = x_target;
    profileData(i).x_actual    = x_actual;
    profileData(i).y_corrected = cum_z_plot;
    profileData(i).U           = Um_plot;
    profileData(i).utau_piv    = utau_piv;
    profileData(i).U_inf       = U_inf;
    profileData(i).y_offset    = y_offset;
    profileData(i).delta99     = delta99;

    fprintf('x = %.1f mm | wall offset = %.3f mm | delta99 = %.3f mm | utau = %.4f m/s\n', ...
        x_actual, y_offset, delta99, utau_piv);
end

%% 6. Plot — physical and inner-scaled overlaid
n_profiles = length(profileData);
gray_ramp  = linspace(0.75, 0.1, n_profiles);   % light → dark
colors     = repmat(gray_ramp', 1, 3);            % RGB: equal channels = gray
discard_vectors = [0.98 0.65 0.60];   % salmon pink

figure('Position', [100 100 1300 560]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile; hold on; box on; grid on;
ax2 = nexttile; hold on; box on; grid on;

for i = 1:length(profileData)
    if isempty(profileData(i).U), continue; end

    y_mm      = profileData(i).y_corrected;          % mm
    U_ms      = profileData(i).U;                    % m/s
    utau_piv  = profileData(i).utau_piv;
    yplus     = (y_mm / 1000) * utau_piv / nu_air;
    uplus     = U_ms / utau_piv;
    lbl       = sprintf('x = %.0f mm', profileData(i).x_actual);

    % Split into gray (first n_gray points) and coloured (rest)
    n         = length(y_mm);
    idx_gray  = 1 : min(n_gray, n);
    idx_color = min(n_gray, n) + 1 : n;

    % Physical
    axes(ax1);
    plot(y_mm(idx_gray),  U_ms(idx_gray),  'o', 'Color', discard_vectors, ...
        'MarkerFaceColor', discard_vectors, 'MarkerSize', 4, 'HandleVisibility', 'off');
    if ~isempty(idx_color)
        plot(y_mm(idx_color), U_ms(idx_color), 'o', 'Color', colors(i,:), ...
            'MarkerFaceColor', colors(i,:), 'MarkerSize', 4, 'DisplayName', lbl);
    end

    % Inner-scaled
    axes(ax2);
    plot(yplus(idx_gray),  uplus(idx_gray),  'o', 'Color', discard_vectors, ...
        'MarkerFaceColor', discard_vectors, 'MarkerSize', 4, 'HandleVisibility', 'off');
    if ~isempty(idx_color)
        plot(yplus(idx_color), uplus(idx_color), 'o', 'Color', colors(i,:), ...
            'MarkerFaceColor', colors(i,:), 'MarkerSize', 4, 'DisplayName', lbl);
    end
end

% Reference curves on inner-scaled panel (using last valid profile's utau)
last = find(~cellfun(@isempty, {profileData.U}), 1, 'last');
utau_ref   = profileData(last).utau_piv;
delta99_ref = profileData(last).delta99 / 1000;
yp_max     = delta99_ref * utau_ref / nu_air;

kappa = 0.39; B = 4.3;
yp_visc = linspace(1, 5, 50);
yp_log  = logspace(log10(30), log10(yp_max), 200);
up_log  = (1/kappa) .* log(yp_log) + B;

axes(ax2);
plot(yp_visc, yp_visc,  'k-',  'LineWidth', 1.5, 'DisplayName', 'U^+ = y^+');
plot(yp_log,  up_log,   'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Log law (\\kappa=%.2f, B=%.1f)', kappa, B));
set(ax2, 'XScale', 'log');
xlim(ax2, [1, yp_max * 1.3]);
xlabel(ax2, '$y^+$',  'Interpreter', 'latex');
ylabel(ax2, '$U^+$',  'Interpreter', 'latex');

title(ax2, 'Inner-scaled');
legend(ax2, 'Location', 'northwest');

axes(ax1);
set(ax1, 'XScale', 'log');
xlabel(ax1, '$y$ (mm)', 'Interpreter', 'latex');
ylabel(ax1, '$U$ (m/s)', 'Interpreter', 'latex');
title(ax1, 'Physical profile');
legend(ax1, 'Location', 'best');

sgtitle(sprintf('%s — PIV profiles, %d x-locations', case_name_query, length(x_locations)), ...
    'FontWeight', 'bold');
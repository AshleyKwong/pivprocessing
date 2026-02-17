%% Compare HWA and PIV - Match Based on Boundary Layer Thickness
% Load HWA data
if ~exist('caseData_run2', 'var')
    load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\HWA\TSFP_Abstract\Case 6\xd=2\case6_xd2_bumpprofile_wH_fixedofi.mat');
    % loads in as caseData_run2
end

set(groot, 'defaultFigureVisible', 'remove');
set(groot, 'defaultFigureVisible', 'on');

%% USER INPUTS
x_nominal = 1050;  % Nominal HWA measurement location [mm] - USER DEFINED
threshold = 0.06;  % Matching tolerance (e.g., 0.05 = 5%, 0.08 = 8%)

%% Load HWA data
HWA_y_delta = caseData_run2.delta_BL;  % in m
HWA_y = caseData_run2.y_cumulative_corrected; % in m
HWA_U = caseData_run2.meanvel_corrected; 
HWA_Uinf = caseData_run2.U_infty_new; 

fprintf('=== HWA Reference Data ===\n');
fprintf('  U_inf = %.3f m/s\n', HWA_Uinf);
fprintf('  δ_99 = %.2f mm\n', HWA_y_delta * 1000);
fprintf('  Location: x = %.0f mm (nominal)\n\n', x_nominal);

%% Extract PIV profiles along a range of x-locations
x_search_range = [1000:5:1200];  % Search from 900 to 1400 mm in 10mm steps

fprintf('Extracting %d PIV profiles...\n', length(x_search_range));
profiles = extractBoundaryLayerProfiles(x_search_range, U_hann_mean, worldX, worldY, ...
    'WallPosition', 0, ...
    'Verbose', false);

% Extract incoming freestream (for normalization)
incoming_Uinf = extractBoundaryLayerProfiles([100], U_hann_mean, worldX, worldY, ...
    'WallPosition', 0, ...
    'Verbose', false);

fprintf('✓ Extracted profiles\n');
fprintf('  Incoming U_inf (PIV) = %.3f m/s\n\n', incoming_Uinf.U_max);

%% Find x-locations and extract delta_99 values
x_positions = [];
delta_99_array = [];

for i = 1:length(profiles)
    if ~isempty(profiles(i).U) && ~isnan(profiles(i).delta_99)
        x_positions(end+1) = profiles(i).x_actual;
        delta_99_array(end+1) = profiles(i).delta_99;  % in mm
    end
end

fprintf('Valid profiles: %d / %d\n', length(x_positions), length(profiles));

%% Find matching location based on delta_99
HWA_delta_mm = HWA_y_delta * 1000;  % Convert to mm

% Find PIV profiles within threshold
delta_error = abs(delta_99_array - HWA_delta_mm) / HWA_delta_mm;
within_threshold = delta_error < threshold;

if any(within_threshold)
    x_matching = x_positions(within_threshold);
    delta_matching = delta_99_array(within_threshold);
    
    % Use the closest match
    [~, best_idx] = min(delta_error);
    x_best = x_positions(best_idx);
    delta_best = delta_99_array(best_idx);
    error_best = delta_error(best_idx) * 100;
    
    fprintf('\n=== Matching Results (δ₉₉ criterion) ===\n');
    fprintf('Target δ₉₉: %.2f mm (HWA)\n', HWA_delta_mm);
    fprintf('Number of matches within %.0f%%: %d\n', threshold*100, sum(within_threshold));
    fprintf('X-range of matches: [%.1f, %.1f] mm\n', min(x_matching), max(x_matching));
    fprintf('\nBest match:\n');
    fprintf('  x = %.1f mm\n', x_best);
    fprintf('  δ₉₉ = %.2f mm (error: %.1f%%)\n', delta_best, error_best);
    
    x_suggested = x_best;
else
    % No exact match - find closest
    [min_error, closest_idx] = min(delta_error);
    x_suggested = x_positions(closest_idx);
    delta_closest = delta_99_array(closest_idx);
    
    fprintf('\n⚠ No profiles within %.0f%% of HWA δ₉₉\n', threshold*100);
    fprintf('Closest match:\n');
    fprintf('  x = %.1f mm\n', x_suggested);
    fprintf('  δ₉₉ = %.2f mm (error: %.1f%%)\n', delta_closest, min_error*100);
end

%% Extract the matching profile
% Find the profile closest to x_suggested
[~, profile_idx] = min(abs([profiles.x_actual] - x_suggested));
profile_matched = profiles(profile_idx);

fprintf('\n=== Profile Comparison ===\n');
fprintf('HWA:\n');
fprintf('  x = %.0f mm (nominal)\n', x_nominal);
fprintf('  δ₉₉ = %.2f mm\n', HWA_delta_mm);
fprintf('  U_inf = %.3f m/s\n', HWA_Uinf);

fprintf('\nPIV (matched):\n');
fprintf('  x = %.1f mm\n', profile_matched.x_actual);
fprintf('  δ₉₉ = %.2f mm\n', profile_matched.delta_99);
fprintf('  U_max = %.3f m/s\n', profile_matched.U_max);

fprintf('\nDifferences:\n');
fprintf('  Δx = %.1f mm\n', profile_matched.x_actual - x_nominal);
fprintf('  Δδ₉₉ = %.2f mm (%.1f%%)\n', profile_matched.delta_99 - HWA_delta_mm, ...
    100*(profile_matched.delta_99 - HWA_delta_mm)/HWA_delta_mm);
fprintf('  ΔU_inf = %.3f m/s (%.1f%%)\n', profile_matched.U_max - HWA_Uinf, ...
    100*(profile_matched.U_max - HWA_Uinf)/HWA_Uinf);

%% Visualization
close all;
figure('Position', [100, 100, 1800, 1000]);

% FIRST ROW: Boundary Layer Growth, Outer Scaling, Outer Scaling Difference

% Subplot 1: Boundary layer growth
subplot(2,3,1);
plot(x_positions, delta_99_array, 'b-', 'LineWidth', 2);
hold on;
yline(HWA_delta_mm, 'r-', 'LineWidth', 2, 'DisplayName', 'HWA δ₉₉');
yline(HWA_delta_mm * (1+threshold), 'r--', 'LineWidth', 1, ...
    'DisplayName', sprintf('±%.0f%%', threshold*100));
yline(HWA_delta_mm * (1-threshold), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
if any(within_threshold)
    plot(x_matching, delta_matching, 'go', 'MarkerSize', 10, 'LineWidth', 2, ...
        'DisplayName', 'Matches');
end
plot(x_suggested, profile_matched.delta_99, 'ks', 'MarkerSize', 15, 'LineWidth', 3, ...
    'DisplayName', sprintf('Best: x=%.0f mm', x_suggested));
xlabel('x [mm]');
ylabel('δ₉₉ [mm]');
title('Boundary Layer Growth');
legend('Location', 'best', 'FontSize', 9);
grid on;

% Subplot 2: Outer scaling
subplot(2,3,2);
% PIV profile
y_delta_PIV = profile_matched.y ./ profile_matched.delta_99;
Uinfty_PIV = (0.99*max(profile_matched.U));
U_outer_PIV = profile_matched.U ./ Uinfty_PIV;

% HWA profile
HWA_y_delta_normalized = HWA_y ./ HWA_y_delta;
HWA_U_normalized = HWA_U ./ HWA_Uinf;

hold on;
plot(y_delta_PIV, U_outer_PIV, "k-o", 'MarkerFaceColor', "k", 'MarkerEdgeColor', "k", ...
    'LineWidth', 1.5, 'MarkerSize', 4, ...
    'DisplayName', sprintf('PIV (x=%.0f mm)', profile_matched.x_actual));
plot(HWA_y_delta_normalized, HWA_U_normalized, "r-o", 'MarkerFaceColor', "r", ...
    'MarkerEdgeColor', "r", 'LineWidth', 1.5, 'MarkerSize', 4, ...
    'DisplayName', sprintf('HWA (x=%.0f mm)', x_nominal));
xlabel('$y/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$U/U_{\infty}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Outer Scaling');
set(gca, 'XScale', 'log');
legend('Location', 'southeast', 'FontSize', 9);
grid on;
xlim([0.001, 10]);
ylim([0, 1.2]);

% Subplot 3: OUTER SCALING DIFFERENCE
subplot(2,3,3);
hold on;

% Interpolate HWA data onto PIV y/δ coordinates
HWA_U_outer_interp = interp1(HWA_y_delta_normalized, HWA_U_normalized, y_delta_PIV, 'linear', NaN);

% Calculate difference: PIV - HWA
U_outer_diff = U_outer_PIV - HWA_U_outer_interp;

% Plot difference vs y/δ
valid_outer_diff = ~isnan(U_outer_diff);
plot(y_delta_PIV(valid_outer_diff), U_outer_diff(valid_outer_diff), 'b-o', ...
    'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'PIV - HWA');

% Add zero reference line
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero difference');

% Shade ±threshold% error band
y_delta_band = logspace(-3, 1, 100);
outer_band_upper = threshold * ones(size(y_delta_band));
outer_band_lower = -outer_band_upper;
fill([y_delta_band, fliplr(y_delta_band)], [outer_band_upper, fliplr(outer_band_lower)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('±%.0f%% band', threshold*100));

xlabel('$y/\delta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Delta(U/U_{\infty})$', 'Interpreter', 'latex', 'FontSize', 14);
title('Difference in Outer Scaling');
legend('Location', 'best', 'FontSize', 9);
grid on;
set(gca, 'XScale', 'log');
xlim([0.001, 10]);
ylim([-0.15, 0.15]);

% Calculate statistics for outer layer
outer_region = y_delta_PIV > 0.1 & y_delta_PIV < 1.0 & valid_outer_diff;
if any(outer_region)
    rms_outer = sqrt(mean(U_outer_diff(outer_region).^2));
    bias_outer = mean(U_outer_diff(outer_region));
    text(0.05, 0.25, sprintf('Outer region (0.1 < y/\\delta < 1.0):\nRMS = %.4f\nBias = %.4f', ...
        rms_outer, bias_outer), 'Units', 'normalized', 'FontSize', 9, ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');
end

% SECOND ROW: Inner Scaling, Inner Scaling Difference, Velocity Field

% Get inner scaling variables
nidaq_cond = load("F:\LAB7 COMPUTER\AK PIV NIDAQ\AK_PIV_NIDaq\U20_y235_aoan11_aoafn11.mat");
Uinf = sqrt(nidaq_cond.qu*2*caseData_run2.rho_air);

ofi_casedata = load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\PLOTTING\case6cf.mat');
Cf_ofi = interp1(ofi_casedata.case2cfdata.("sorted global x mm"), ofi_casedata.case2cfdata.("Cf Avg window 16"), 8250);

q_1 = check.qu;
P1 = check.P0 / 1013;
T1 = check.T_plc;
q_2 = check2.qu;
P2 = check2.P0 / 1013;
T2 = check2.T_plc;

[rho_air_1, visc_air_1] = fluid_prop(T1, P1);
[rho_air_2, visc_air_2] = fluid_prop(T2, P2);

U1 = sqrt(2*q_1/rho_air_1);
U2 = sqrt(2*q_2/rho_air_2);
Uinf_pitot = mean([U1 U2]);

Cf_ofi_corrected = Cf_ofi*(Uinf_pitot/ Uinfty_PIV)^2;

utau_musker = caseData_run2.utau_musker;
Cf_musker = 2*utau_musker^2 / (HWA_Uinf^2);

fprintf("\tThe error bw Cf_musker and Cf_ofi = %.2f %% \n", (abs(Cf_musker - Cf_ofi_corrected)/ Cf_ofi_corrected)*100 )

utau_piv = sqrt(Cf_musker*Uinfty_PIV^2*0.5);
kin_vis = caseData_run2.nu_air;

% Calculate inner scaling coordinates
y_inner_m = (profile_matched.y./1000);
y_inner_piv = (y_inner_m) .* utau_piv ./ kin_vis;
U_inner_PIV = profile_matched.U ./ utau_piv;

HWA_y_inner_normalized = (HWA_y .* utau_musker) ./ kin_vis;
HWA_U_inner_normalized = HWA_U ./ utau_musker;

% Subplot 4: Inner scaling
subplot(2,3,4);
hold on;

% Plot both profiles
plot(y_inner_piv, U_inner_PIV, "k-o", 'MarkerFaceColor', "k", 'MarkerEdgeColor', "k", ...
    'LineWidth', 1.5, 'MarkerSize', 4, ...
    'DisplayName', sprintf('PIV (x=%.0f mm)', profile_matched.x_actual));
plot(HWA_y_inner_normalized, HWA_U_inner_normalized, "r-o", 'MarkerFaceColor', "r", ...
    'MarkerEdgeColor', "r", 'LineWidth', 1.5, 'MarkerSize', 4, ...
    'DisplayName', sprintf('HWA (x=%.0f mm)', x_nominal));

% Add log-law and viscous sublayer references
y_plus_ref = logspace(0, 4, 100);
U_plus_viscous = y_plus_ref;
U_plus_log = (1/0.41) * log(y_plus_ref) + 5.0;

plot(y_plus_ref, U_plus_viscous, 'b--', 'LineWidth', 1, 'DisplayName', 'u^+ = y^+');
plot(y_plus_ref, U_plus_log, 'g--', 'LineWidth', 1, 'DisplayName', 'Log law');

xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$U^+$', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('Inner Scaling (u_\\tau = %.3f m/s)', utau_musker));
legend('Location', 'southeast', 'FontSize', 9);
grid on;
set(gca, 'XScale', 'log');
xlim([1, 15000]);
ylim([0, 35]);

% Subplot 5: INNER SCALING DIFFERENCE
subplot(2,3,5);
hold on;

% Interpolate HWA data onto PIV y+ coordinates
HWA_U_interp = interp1(HWA_y_inner_normalized, HWA_U_inner_normalized, y_inner_piv, 'linear', NaN);

% Calculate difference: PIV - HWA
U_diff = U_inner_PIV - HWA_U_interp;

% Plot difference vs y+
valid_diff = ~isnan(U_diff);
plot(y_inner_piv(valid_diff), U_diff(valid_diff), 'b-o', 'LineWidth', 2, 'MarkerSize', 4, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'PIV - HWA');

% Add zero reference line
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero difference');

% Shade ±threshold% error band
y_plus_band = logspace(0, 4, 100);
U_plus_band_upper = threshold * ((1/0.41) * log(y_plus_band) + 5.0);
U_plus_band_lower = -U_plus_band_upper;
fill([y_plus_band, fliplr(y_plus_band)], [U_plus_band_upper, fliplr(U_plus_band_lower)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('±%.0f%% band', threshold*100));

xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Delta U^+ = U^+_{PIV} - U^+_{HWA}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Difference in Inner Scaling');
legend('Location', 'best', 'FontSize', 9);
ylim([-6, 6]);
grid on;
set(gca, 'XScale', 'log');
xlim([1, 10000]);

% Calculate RMS difference and bias in overlap region
overlap_region = y_inner_piv > 30 & y_inner_piv < 300 & valid_diff;
if any(overlap_region)
    rms_diff = sqrt(mean(U_diff(overlap_region).^2));
    bias_diff = mean(U_diff(overlap_region));
    text(0.05, 0.25, sprintf('Overlap (30 < y^+ < 300):\nRMS = %.2f\nBias = %.2f', ...
        rms_diff, bias_diff), 'Units', 'normalized', 'FontSize', 9, ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');
end

% Subplot 6: PIV Velocity field with profile location
subplot(2,3,6);
imagesc(worldX(1,:), worldY(:,1), U_hann_mean);
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
colormap(gca, jet);
hold on;

% Mark the matched PIV profile location
xline(x_suggested, 'k--', 'LineWidth', 2, ...
    'DisplayName', sprintf('PIV x=%.0f mm', x_suggested));
yline(profile_matched.delta_99, 'k--', 'LineWidth', 2, ...
    'DisplayName', sprintf('δ₉₉=%.3f mm', profile_matched.delta_99));

% Mark ONLY the nominal HWA location (user-defined)
% xline(x_nominal, 'k--', 'LineWidth', 3, ...
%     'DisplayName', sprintf('HWA x=%.0f mm', x_nominal));
ax = gca; 
ax.FontSize = 12; 
xlabel('X [mm]');
ylabel('Y [mm]');
title('PIV Velocity Field with Profile Locations');
legend('Location', 'northeast', 'FontSize', 9);
clim([0, 25]);

sgtitle(sprintf('PIV-HWA Comparison: δ₉₉ Matching (PIV x=%.0f mm, HWA x=%.0f mm)', ...
    x_suggested, x_nominal), 'FontSize', 16, 'FontWeight', 'bold');
%% 
figure(); 

%% Print summary statistics
fprintf('\n=== Difference Statistics ===\n');

% OUTER SCALING ANALYSIS
fprintf('Outer Scaling (0.1 < y/δ < 1.0):\n');
if any(outer_region)
    fprintf('  Mean difference: %.4f (bias)\n', bias_outer);
    fprintf('  RMS difference: %.4f\n', rms_outer);
    fprintf('  Relative error: %.1f%%\n', 100*rms_outer);
    
    % Find last point within threshold % in outer scaling
    within_threshold_outer = abs(U_outer_diff) <= threshold & valid_outer_diff;
    if any(within_threshold_outer)
        last_idx_outer = find(within_threshold_outer, 1, 'last');
        y_last_outer_normalized = y_delta_PIV(last_idx_outer);
        y_last_outer_mm = profile_matched.y(last_idx_outer);
        fprintf('  Last point within ±%.0f%%: y/δ = %.4f (y = %.2f mm)\n', ...
            threshold*100, y_last_outer_normalized, y_last_outer_mm);
    else
        fprintf('  ⚠ No points within ±%.0f%%\n', threshold*100);
    end
end

% INNER SCALING ANALYSIS
fprintf('\nInner Scaling (30 < y+ < 300):\n');
if any(overlap_region)
    fprintf('  Mean difference: %.2f U+ (bias)\n', bias_diff);
    fprintf('  RMS difference: %.2f U+\n', rms_diff);
    fprintf('  Relative error: %.1f%% (based on U+ ~ 15)\n', 100*rms_diff/15);
    
    % Find last point within threshold% in inner scaling
    U_plus_reference = (1/0.41) * log(y_inner_piv) + 5.0;
    error_threshold_inner = threshold * U_plus_reference;
    
    within_threshold_inner = abs(U_diff) <= error_threshold_inner & valid_diff;
    if any(within_threshold_inner)
        last_idx_inner = find(within_threshold_inner, 1, 'last');
        y_last_inner_plus = y_inner_piv(last_idx_inner);
        y_last_inner_mm = profile_matched.y(last_idx_inner);
        fprintf('  Last point within ±%.0f%%: y+ = %.1f (y = %.2f mm)\n', ...
            threshold*100, y_last_inner_plus, y_last_inner_mm);
    else
        fprintf('  ⚠ No points within ±%.0f%%\n', threshold*100);
    end
end

%% Summary table
fprintf('\n=== Summary Table ===\n');
fprintf('%-20s %12s %12s %12s\n', 'Parameter', 'HWA', 'PIV', 'Difference');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('%-20s %12.1f %12.1f %12.1f\n', 'x-location [mm]', x_nominal, profile_matched.x_actual, profile_matched.x_actual - x_nominal);
fprintf('%-20s %12.2f %12.2f %12.2f\n', 'δ₉₉ [mm]', HWA_delta_mm, profile_matched.delta_99, profile_matched.delta_99 - HWA_delta_mm);
fprintf('%-20s %12.3f %12.3f %12.3f\n', 'U_inf [m/s]', HWA_Uinf, profile_matched.U_max, profile_matched.U_max - HWA_Uinf);
fprintf('%-20s %12s %12.2f %12s\n', 'Error in δ₉₉ [%]', '-', 100*delta_error(profile_idx), '-');

% Add agreement quality section
fprintf('\n=== Agreement Quality (±%.0f%% threshold) ===\n', threshold*100);
if exist('y_last_outer_normalized', 'var')
    fprintf('Outer scaling: Agreement within ±%.0f%% from y/δ = %.4f (%.2f mm)\n', ...
        threshold*100, y_last_outer_normalized, y_last_outer_mm);
else
    fprintf('Outer scaling: No agreement within ±%.0f%%\n', threshold*100);
end

if exist('y_last_inner_plus', 'var')
    fprintf('Inner scaling: Agreement within ±%.0f%% from y+ = %.1f (%.2f mm)\n', ...
        threshold*100, y_last_inner_plus, y_last_inner_mm);
else
    fprintf('Inner scaling: No agreement within ±%.0f%%\n', threshold*100);
end

%% Save matching data
matching_data = struct();
matching_data.method = 'delta_99_matching';
matching_data.x_HWA_nominal = x_nominal;
matching_data.x_PIV = profile_matched.x_actual;
matching_data.delta_HWA = HWA_delta_mm;
matching_data.delta_PIV = profile_matched.delta_99;
matching_data.Uinf_HWA = HWA_Uinf;
matching_data.Umax_PIV = profile_matched.U_max;
matching_data.error_delta_percent = 100*delta_error(profile_idx);
matching_data.threshold = threshold;
matching_data.profile_PIV = profile_matched;
matching_data.x_search_range = x_search_range;
matching_data.all_x_positions = x_positions;
matching_data.all_delta_99 = delta_99_array;

% Add agreement information
if exist('y_last_outer_normalized', 'var')
    matching_data.outer_agreement_y_delta = y_last_outer_normalized;
    matching_data.outer_agreement_y_mm = y_last_outer_mm;
end
if exist('y_last_inner_plus', 'var')
    matching_data.inner_agreement_y_plus = y_last_inner_plus;
    matching_data.inner_agreement_y_mm = y_last_inner_mm;
end

save('PIV_HWA_delta99_matching.mat', 'matching_data');
fprintf('\n✓ Saved matching data to PIV_HWA_delta99_matching.mat\n');

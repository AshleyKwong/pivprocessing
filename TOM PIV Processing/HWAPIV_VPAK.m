%% compare_VP_AK_profiles_HWA.m
% VP (grayscale circles) vs AK (red diamonds), PIV + HWA comparison
clear; close all; clc;
set(groot, 'defaultAxesFontName',  'Cambria Math');
set(groot, 'defaultAxesFontSize',  14);
set(groot, 'defaultTextFontName',  'Cambria Math');
set(groot, 'defaultTextFontSize',  14);

%% ========================================================================
%% USER SETTINGS
%% ========================================================================

%% --- VP PIV file paths --------------------------------------------------
VP_mergedFile  = 'G:\SW 500mm Mean Flow Fields\Mean_Flow_Feilds_AOA_-8.mat';
VP_Cf_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\Cf_minus8.mat';
VP_Re_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\Re_minus8.mat';
VP_PIVsum_file = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\SW_PIV_Summary_-8.mat';

%% --- VP HWA file --------------------------------------------------------
VP_HWA_file    = 'C:\Users\ak1u24\OneDrive - University of Southampton\Thomas Preskett Old Work\PRF_Data\Turbulence_Profiles\SW_Vel_20_AOA_-8.csv';

%% --- AK PIV file paths --------------------------------------------------
AK_mean_field  = 'G:\Y235_AOAN04_AOAFN06_SmallerWindows\merge_instantaneousavg_20260429_083207\merged_meanUV_14loops_20260429_083207.mat';
AK_Cf_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\OFI RESULTS\OFI_Cf_results_20260426_165929.mat';
AK_blSweep     = load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\16x96init\blSweep_20260430_084136.mat').blSweep;
AK_case_name   = 'Case 2';

%% --- AK HWA file --------------------------------------------------------
AK_HWA_file    = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\HWA\TSFP_Abstract\Case 2\xd=2\case2_xd2_HWA_20260429_115708.mat';
AK_HWA_mode    = 'poly';   % field suffix in caseData

%% --- x-locations (matched by index) ------------------------------------
VP_x_extract_mm = [ 9000];   % VP global x (mm)
AK_x_extract_mm = [ 1000];   % AK local PIV x (mm)

%% --- Shared parameters --------------------------------------------------
VP_x_tol_mm  = 5;
AK_nu_air    = 1.5e-5;
VP_nu_air    = 1.51e-5;
VP_Uinf0     = 19.2085724;
n_gray       = 6;        % AK near-wall salmon points
marker_skip  = 1;
font_size    = 13;
kappa        = 0.39;
B            = 4.3;

%% --- Colour definitions -------------------------------------------------
n_x = numel(VP_x_extract_mm);

% VP PIV: dark → light gray ramp
vp_dark  = [0.15 0.15 0.15];
vp_light = [0.70 0.70 0.70];
vp_dark      =  [0.15 0.15 0.15];   % near-black

VP_colors = vp_dark; % [linspace(vp_dark(1),  vp_light(1),  n_x)', ...
             % linspace(vp_dark(2),  vp_light(2),  n_x)', ...
             % linspace(vp_dark(3),  vp_light(3),  n_x)'];
VP_HWA_color = [0.55 0.55 0.55];   % mid-gray, clearly distinct from PIV

% VP HWA: single light gray
VP_HWA_color = [0.75 0.75 0.75];

% AK PIV: dark → bright red ramp
% red_dark   = [0.6, 0.0, 0.0];
red_dark     = [0.6 0.00 0.00];   % deep crimson
red_bright = [1.0, 0.4, 0.4];
AK_colors  = red_dark ; %[linspace(red_dark(1), red_bright(1), n_x)', ...
              % linspace(red_dark(2), red_bright(2), n_x)', ...
              % linspace(red_dark(3), red_bright(3), n_x)'];

% AK HWA: single light red
% AK_HWA_color = [1.0 0.7 0.7];
AK_HWA_color = [1.00 0.55 0.55];   % soft coral, distinct from PIV red

% Near-wall discard colour (AK PIV only)
salmon_raw   = [1.0, 0.4, 0.4];
alpha_equiv  = 0.4;   % 0 = white, 1 = full colour — adjust to taste
salmon       = salmon_raw * alpha_equiv + [1 1 1] * (1 - alpha_equiv);
%% ========================================================================
%% LOAD VP PIV DATA
%% ========================================================================
fprintf('=== Loading VP PIV data ===\n');
Cf_data    = load(VP_Cf_file);
Re_data    = load(VP_Re_file);
SW_PIVData = load(VP_PIVsum_file);
vp_piv     = load(VP_mergedFile, 'X', 'Y', 'U_mean');

Cf_vec1 = Cf_data.Cf_vec1;
Re_vec1 = Re_data.Re_vec1;

Virgilio_Rex       = Re_vec1(:, 2);
VP_x_OFI_m        = Virgilio_Rex .* VP_nu_air ./ VP_Uinf0;
VP_Cf_raw          = mean(Cf_vec1, 2);
Virgilio_x_PIV     = SW_PIVData.X_m(:);
Virgilio_localUinf = SW_PIVData.U99 .* 0.99;

[~, matched_idx] = min(abs(VP_x_OFI_m(:) - Virgilio_x_PIV(:)'), [], 2);
VP_Uinf_local    = Virgilio_localUinf(matched_idx);
if size(VP_Cf_raw,1) ~= size(VP_Uinf_local,1)
    VP_Uinf_local = VP_Uinf_local';
end

VP_x_PIV_mm_all = Virgilio_x_PIV .* 1000;
VP_xf_all       = max(VP_x_PIV_mm_all);
VP_Uout_mask    = VP_x_PIV_mm_all >= (VP_xf_all - 4);
Virg_Uinfoutlet = mean(Virgilio_localUinf(VP_Uout_mask), 'omitnan');

cf_threshold  = 0.0004 * 2;
VP_Cf_outlet  = VP_Cf_raw .* 2 .* (VP_Uinf0 ./ Virg_Uinfoutlet)^2;
valid_vp      = VP_Cf_outlet >= cf_threshold & ~isnan(VP_Cf_outlet);
VP_x_OFI_m   = VP_x_OFI_m(valid_vp);
VP_Cf_raw     = VP_Cf_raw(valid_vp);
VP_Uinf_local = VP_Uinf_local(valid_vp);
[VP_x_OFI_m, vp_sort] = sort(VP_x_OFI_m, 'ascend');
VP_Cf_raw     = VP_Cf_raw(vp_sort);
VP_Uinf_local = VP_Uinf_local(vp_sort);

VP_x_BL_m   = SW_PIVData.X_m(:);
VP_delta99_m = SW_PIVData.delta_m(:);
VP_Ue_m      = SW_PIVData.U99(:);
[VP_x_BL_m, bl_sort] = sort(VP_x_BL_m, 'ascend');
VP_delta99_m = VP_delta99_m(bl_sort);
VP_Ue_m      = VP_Ue_m(bl_sort);

vp_x_mm = vp_piv.X(1,:) .* 1e3;
vp_y_mm = vp_piv.Y(:,1) .* 1e3;

%% ========================================================================
%% LOAD VP HWA DATA (CSV)
%% ========================================================================
fprintf('=== Loading VP HWA data ===\n');
VP_HWA_raw  = readtable(VP_HWA_file);
% Columns: index | y(m) | U(m/s) | uu | delta | U_99 | utau | nu | Cf | theta | deltastar | PI | Re_tau
VP_HWA_y_m  = VP_HWA_raw{:, 2};    % y (m)
VP_HWA_U_ms = VP_HWA_raw{:, 3};    % U (m/s)
VP_HWA_utau = VP_HWA_raw{1, 7};    % single utau value (first row)
VP_HWA_nu   = VP_HWA_raw{1, 8};    % nu (m2/s)

% Sort ascending in y
[VP_HWA_y_m, hwa_sort] = sort(VP_HWA_y_m, 'ascend');
VP_HWA_U_ms = VP_HWA_U_ms(hwa_sort);

VP_HWA_y_mm  = VP_HWA_y_m  .* 1e3;                        % m → mm
VP_HWA_yplus = VP_HWA_y_m  .* VP_HWA_utau ./ VP_HWA_nu;  % y+
VP_HWA_uplus = VP_HWA_U_ms ./ VP_HWA_utau;                % U+

fprintf('  VP HWA: %d points | utau=%.4f m/s | y range=[%.3f, %.1f] mm\n', ...
    numel(VP_HWA_y_mm), VP_HWA_utau, min(VP_HWA_y_mm), max(VP_HWA_y_mm));

%% ========================================================================
%% LOAD AK PIV DATA
%% ========================================================================
fprintf('=== Loading AK PIV data ===\n');
AK_U_mean   = load(AK_mean_field).U_hann_mean;
AK_worldX   = load(AK_mean_field).worldX;
AK_worldY   = load(AK_mean_field).worldY;
AK_mag      = (max(AK_worldX(1,:)) - min(AK_worldX(1,:))) / size(AK_worldX,2);
AK_x_vec    = AK_worldX(1,:);
AK_y_vec    = AK_worldY(:,1);

AK_Cf_raw   = load(AK_Cf_file);
ak_case_idx = find(strcmp({AK_Cf_raw.results.caseName}, AK_case_name), 1);
if isempty(ak_case_idx)
    error('AK case "%s" not found.', AK_case_name);
end
AK_Cf_data  = AK_Cf_raw.results(ak_case_idx);

%% ========================================================================
%% LOAD AK HWA DATA (.mat)
%% ========================================================================
fprintf('=== Loading AK HWA data ===\n');
AK_HWA_raw   = load(AK_HWA_file);
AK_HWA_y_mm  = AK_HWA_raw.caseData.(sprintf('y_corrected_%s',       AK_HWA_mode)) .* 1e3;
AK_HWA_U_ms  = AK_HWA_raw.caseData.(sprintf('meanvel_%s_corrected', AK_HWA_mode));
AK_HWA_utau  = AK_HWA_raw.caseData.utau_ofi;  % adjust field name if needed
AK_HWA_nu    = AK_HWA_raw.caseData.nu_air;

% Sort ascending
[AK_HWA_y_mm, ak_hwa_sort] = sort(AK_HWA_y_mm, 'ascend');
AK_HWA_U_ms  = AK_HWA_U_ms(ak_hwa_sort);
AK_HWA_yplus = (AK_HWA_y_mm ./ 1e3) .* AK_HWA_utau ./ AK_HWA_nu;
AK_HWA_uplus = AK_HWA_U_ms ./ AK_HWA_utau;

fprintf('  AK HWA: %d points | utau=%.4f m/s | y range=[%.3f, %.1f] mm\n', ...
    numel(AK_HWA_y_mm), AK_HWA_utau, min(AK_HWA_y_mm), max(AK_HWA_y_mm));

%% ========================================================================
%% EXTRACT VP PIV PROFILES
%% ========================================================================
fprintf('=== Extracting VP PIV profiles ===\n');
VP_prof(n_x) = struct();

for k = 1:n_x
    x_target_mm = VP_x_extract_mm(k);
    x_target_m  = x_target_mm / 1000;

    [dist_x, col_idx] = min(abs(vp_x_mm - x_target_mm));
    if dist_x > VP_x_tol_mm
        warning('VP station %d: %.1f mm from target', k, dist_x);
    end

    U_col = vp_piv.U_mean(:, col_idx);
    valid = ~isnan(U_col) & vp_y_mm > 0;
    y_mm  = vp_y_mm(valid);
    U_ms  = U_col(valid);
    [y_mm, ys] = sort(y_mm, 'ascend');
    U_ms       = U_ms(ys);

    Cf_raw_here  = interp1(VP_x_OFI_m, VP_Cf_raw,     x_target_m, 'linear', NaN);
    Uinf_loc_k   = interp1(VP_x_OFI_m, VP_Uinf_local, x_target_m, 'linear', NaN);
    u_tau        = sqrt(Cf_raw_here) * VP_Uinf0 * (Uinf_loc_k / VP_Uinf0);
    delta99_here = interp1(VP_x_BL_m,  VP_delta99_m,  x_target_m, 'linear', NaN);

    y_m   = y_mm / 1000;
    yplus = y_m  .* u_tau ./ VP_nu_air;
    uplus = U_ms ./ u_tau;

    VP_prof(k).x_mm    = x_target_mm;
    VP_prof(k).y_mm    = y_mm;
    VP_prof(k).U_ms    = U_ms;
    VP_prof(k).yplus   = yplus;
    VP_prof(k).uplus   = uplus;
    VP_prof(k).u_tau   = u_tau;
    VP_prof(k).delta99 = delta99_here * 1e3;
    VP_prof(k).Re_tau  = u_tau * delta99_here / VP_nu_air;
    VP_prof(k).color   = VP_colors(k,:);
    VP_prof(k).label   = sprintf('VP PIV  $x=%d$ mm,  $Re_\\tau=%.0f$', ...
                             x_target_mm, VP_prof(k).Re_tau);

    fprintf('  VP x=%4d mm | u_tau=%.4f | delta99=%.2f mm | Re_tau=%.0f\n', ...
        x_target_mm, u_tau, delta99_here*1e3, VP_prof(k).Re_tau);
end

%% ========================================================================
%% EXTRACT AK PIV PROFILES
%% ========================================================================
fprintf('=== Extracting AK PIV profiles ===\n');
AK_prof(n_x) = struct();

for i = 1:n_x
    x_target = AK_x_extract_mm(i);

    x_global  = x_target + 7200;
    [~, midx] = min(abs(AK_Cf_data.x_centers - x_global));
    Cf_half   = AK_Cf_data.Cf(midx) / 2;
    utau_piv  = sqrt(Cf_half * AK_Cf_data.U_inf_local(midx)^2);

    cols = find(AK_x_vec >= x_target - AK_mag & AK_x_vec <= x_target + AK_mag);
    if isempty(cols)
        warning('AK: no data near x=%.1f mm', x_target); continue
    end

    U_profile = mean(AK_U_mean(:, cols), 2, 'omitnan');
    cum_z     = double(AK_y_vec);
    Um        = double(U_profile);

    valid_mask    = ~isnan(Um);
    Um            = Um(valid_mask);
    cum_z         = cum_z(valid_mask);
    [cum_z, sidx] = sort(cum_z, 'ascend');
    Um            = Um(sidx);

    U_inf       = max(Um);
    wall_thresh = 0.0005 * U_inf;
    first_valid = find(Um > wall_thresh, 1, 'first');
    if isempty(first_valid), warning('AK: no valid wall at x=%.1f', x_target); continue; end

    y_offset   = cum_z(first_valid);
    cum_z_corr = cum_z - y_offset + 0.2;
    trim_idx   = first_valid + 3;

    y_mm_plot = cum_z_corr(trim_idx:end);
    U_ms_plot = Um(trim_idx:end);
    delta99   = interp1(AK_blSweep.x_mm, AK_blSweep.delta99_hybrid_mm, x_target, 'linear', 'extrap');
    yplus     = (y_mm_plot / 1000) .* utau_piv ./ AK_nu_air;
    uplus     = U_ms_plot ./ utau_piv;

    AK_prof(i).x_mm    = x_target;
    AK_prof(i).x_actual = mean(AK_x_vec(cols));
    AK_prof(i).y_mm    = y_mm_plot;
    AK_prof(i).U_ms    = U_ms_plot;
    AK_prof(i).yplus   = yplus;
    AK_prof(i).uplus   = uplus;
    AK_prof(i).utau    = utau_piv;
    AK_prof(i).delta99 = delta99;
    AK_prof(i).color   = AK_colors(i,:);
    AK_prof(i).label   = sprintf('AK PIV  $x=%d$ mm', x_target);

    fprintf('  AK x=%4d mm | u_tau=%.4f | delta99=%.3f mm\n', x_target, utau_piv, delta99);
end

%% ========================================================================
%% REFERENCE LOG-LAW CURVES  (use last AK profile for y+ range)
%% ========================================================================
last_ak    = find(~cellfun(@isempty, {AK_prof.U_ms}), 1, 'last');
utau_ref   = AK_prof(last_ak).utau;
delta99_ref = AK_prof(last_ak).delta99 / 1000;
yp_max     = delta99_ref * utau_ref / AK_nu_air;

yp_visc = linspace(1, 5, 50);
yp_log  = logspace(log10(30), log10(yp_max), 300);
up_log  = (1/kappa) .* log(yp_log) + B;

%% ========================================================================
%% HELPER: add one PIV profile pair to axes (circles=VP, diamonds=AK)
%% ========================================================================
function add_piv_pair(ax_phys, ax_inner, vp, ak, n_gray, salmon, marker_skip)
    idx_mk = 1:marker_skip:numel(vp.y_mm);
    n      = numel(ak.y_mm);
    ig     = 1:min(n_gray, n);
    ic     = min(n_gray,n)+1:n;

    % VP PIV — circles
    plot(ax_phys,  vp.y_mm,          vp.U_ms,  '-',  'Color', vp.color, 'LineWidth', 1.4, 'DisplayName', vp.label);
    plot(ax_phys,  vp.y_mm(idx_mk),  vp.U_ms(idx_mk), 'o', ...
        'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(ax_inner, vp.yplus,         vp.uplus, '-',  'Color', vp.color, 'LineWidth', 1.4, 'DisplayName', vp.label);
    plot(ax_inner, vp.yplus(idx_mk), vp.uplus(idx_mk), 'o', ...
        'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 5, 'HandleVisibility', 'off');

    % AK PIV — diamonds, salmon near-wall
    plot(ax_phys,  ak.y_mm(ig),  ak.U_ms(ig),  'd', 'Color', salmon, ...
        'MarkerFaceColor', salmon, 'MarkerSize', 5, 'HandleVisibility', 'off');
    if ~isempty(ic)
        plot(ax_phys,  ak.y_mm(ic),  ak.U_ms(ic),  '-d', 'Color', ak.color, ...
            'MarkerFaceColor', ak.color, 'MarkerSize', 5, 'LineWidth', 1.4, 'DisplayName', ak.label);
        plot(ax_inner, ak.yplus(ig),  ak.uplus(ig),  'd', 'Color', salmon, ...
            'MarkerFaceColor', salmon, 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(ax_inner, ak.yplus(ic), ak.uplus(ic), '-d', 'Color', ak.color, ...
            'MarkerFaceColor', ak.color, 'MarkerSize', 5, 'LineWidth', 1.4, 'DisplayName', ak.label);
    end
end

%% ========================================================================
%% HELPER: format a panel pair
%% ========================================================================
function fmt_panels(ax_phys, ax_inner, yp_max, font_size)
    axes_list = [ax_phys, ax_inner];
    for ax = axes_list
        set(ax, 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
        grid(ax, 'on'); box(ax, 'on');
    end
    xlabel(ax_phys,  '$y$ (mm)', 'Interpreter', 'latex');
    ylabel(ax_phys,  '$U$ (m/s)', 'Interpreter', 'latex');
    title(ax_phys,   'Physical',  'Interpreter', 'latex');
    xlabel(ax_inner, '$y^+$',    'Interpreter', 'latex');
    ylabel(ax_inner, '$U^+$',    'Interpreter', 'latex');
    title(ax_inner,  'Inner-scaled', 'Interpreter', 'latex');
    xlim(ax_inner, [1, yp_max * 1.4]);
    legend(ax_phys,  'Interpreter', 'latex', 'Location', 'best');
    legend(ax_inner, 'Interpreter', 'latex', 'Location', 'northwest');
end

%% ========================================================================
%% FIGURE 1 — Mega overlay: all x-locations, PIV only
%% ========================================================================
fig1 = figure('Name', 'Fig1: All PIV profiles — VP vs AK', 'Position', [50 50 1300 560]);
tl1  = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1_phys  = nexttile(tl1); hold on;
ax1_inner = nexttile(tl1); hold on;

for k = 1:n_x
    add_piv_pair(ax1_phys, ax1_inner, VP_prof(k), AK_prof(k), n_gray, salmon, marker_skip);
end

% Log-law refs on inner panel
plot(ax1_inner, yp_visc, yp_visc, 'k-',  'LineWidth', 1.5, 'DisplayName', '$U^+=y^+$');
plot(ax1_inner, yp_log,  up_log,  'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('$\\kappa=%.2f,\\,B=%.1f$', kappa, B));
fmt_panels(ax1_phys, ax1_inner, yp_max, font_size);
sgtitle('VP vs AK — all PIV x-locations', 'FontWeight', 'bold');

%% ========================================================================
%% FIGURE 2 — HWA comparison: 2 panels, 4 lines each
%% ========================================================================
% Use darkest VP and AK colors for the PIV lines on this figure
vp_col_hwa = VP_colors(end,:);   % darkest gray
ak_col_hwa = AK_colors(end,:);   % darkest red

fig2 = figure('Name', 'Fig2: HWA vs PIV — VP and AK', 'Position', [80 80 1300 560]);
tl2  = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax2_phys  = nexttile(tl2); hold on;
ax2_inner = nexttile(tl2); hold on;

% VP PIV (circle, dark gray, solid) — use closest x-location to HWA
[~, vp_hwa_k] = min(abs(VP_x_extract_mm - VP_HWA_raw{1,5}*1e3));  % delta col ~ x proxy; just use last
vp_hwa_k = n_x;  % default: use last (typically closest to HWA x)
idx_mk = 1:marker_skip:numel(VP_prof(vp_hwa_k).y_mm);
plot(ax2_phys,  VP_prof(vp_hwa_k).y_mm,          VP_prof(vp_hwa_k).U_ms,  '-',  'Color', vp_col_hwa, 'LineWidth', 1.5, 'DisplayName', 'VP PIV');
plot(ax2_phys,  VP_prof(vp_hwa_k).y_mm(idx_mk),  VP_prof(vp_hwa_k).U_ms(idx_mk), 'o', ...
    'Color', vp_col_hwa, 'MarkerFaceColor', vp_col_hwa, 'MarkerSize', 5, 'HandleVisibility', 'off');
plot(ax2_inner, VP_prof(vp_hwa_k).yplus,          VP_prof(vp_hwa_k).uplus, '-',  'Color', vp_col_hwa, 'LineWidth', 1.5, 'DisplayName', 'VP PIV');
plot(ax2_inner, VP_prof(vp_hwa_k).yplus(idx_mk),  VP_prof(vp_hwa_k).uplus(idx_mk), 'o', ...
    'Color', vp_col_hwa, 'MarkerFaceColor', vp_col_hwa, 'MarkerSize', 5, 'HandleVisibility', 'off');

% VP HWA (circle, light gray, dashed)
idx_mk_vph = 1:marker_skip:numel(VP_HWA_y_mm);
plot(ax2_phys,  VP_HWA_y_mm,           VP_HWA_U_ms,  '--', 'Color', VP_HWA_color, 'LineWidth', 1.5, 'DisplayName', 'VP HWA');
plot(ax2_phys,  VP_HWA_y_mm(idx_mk_vph), VP_HWA_U_ms(idx_mk_vph), 'o', ...
    'Color', VP_HWA_color, 'MarkerFaceColor', VP_HWA_color, 'MarkerSize', 5, 'HandleVisibility', 'off');
plot(ax2_inner, VP_HWA_yplus,          VP_HWA_uplus, '--', 'Color', VP_HWA_color, 'LineWidth', 1.5, 'DisplayName', 'VP HWA');
plot(ax2_inner, VP_HWA_yplus(idx_mk_vph), VP_HWA_uplus(idx_mk_vph), 'o', ...
    'Color', VP_HWA_color, 'MarkerFaceColor', VP_HWA_color, 'MarkerSize', 5, 'HandleVisibility', 'off');

% AK PIV (diamond, dark red, solid) — use closest AK x to AK HWA x
ak_hwa_i = n_x;  % default last; adjust if needed
n_ak = numel(AK_prof(ak_hwa_i).y_mm);
ig   = 1:min(n_gray, n_ak);
ic   = min(n_gray, n_ak)+1:n_ak;
idx_mk_ak = 1:marker_skip:numel(AK_prof(ak_hwa_i).y_mm);
plot(ax2_phys,  AK_prof(ak_hwa_i).y_mm(ig),  AK_prof(ak_hwa_i).U_ms(ig),  'd', ...
    'Color', salmon, 'MarkerFaceColor', salmon, 'MarkerSize', 5, 'HandleVisibility', 'off');
if ~isempty(ic)
    plot(ax2_phys,  AK_prof(ak_hwa_i).y_mm(ic),  AK_prof(ak_hwa_i).U_ms(ic),  '-d', ...
        'Color', ak_col_hwa, 'MarkerFaceColor', ak_col_hwa, 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'AK PIV');
    plot(ax2_inner, AK_prof(ak_hwa_i).yplus(ig),  AK_prof(ak_hwa_i).uplus(ig),  'd', ...
        'Color', salmon, 'MarkerFaceColor', salmon, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(ax2_inner, AK_prof(ak_hwa_i).yplus(ic), AK_prof(ak_hwa_i).uplus(ic), '-d', ...
        'Color', ak_col_hwa, 'MarkerFaceColor', ak_col_hwa, 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'AK PIV');
end

% AK HWA (diamond, light red, dashed)
idx_mk_akh = 1:marker_skip:numel(AK_HWA_y_mm);
plot(ax2_phys,  AK_HWA_y_mm,             AK_HWA_U_ms,  '--', 'Color', AK_HWA_color, 'LineWidth', 1.5, 'DisplayName', 'AK HWA');
plot(ax2_phys,  AK_HWA_y_mm(idx_mk_akh), AK_HWA_U_ms(idx_mk_akh), 'd', ...
    'Color', AK_HWA_color, 'MarkerFaceColor', AK_HWA_color, 'MarkerSize', 5, 'HandleVisibility', 'off');
plot(ax2_inner, AK_HWA_yplus,            AK_HWA_uplus, '--', 'Color', AK_HWA_color, 'LineWidth', 1.5, 'DisplayName', 'AK HWA');
plot(ax2_inner, AK_HWA_yplus(idx_mk_akh), AK_HWA_uplus(idx_mk_akh), 'd', ...
    'Color', AK_HWA_color, 'MarkerFaceColor', AK_HWA_color, 'MarkerSize', 5, 'HandleVisibility', 'off');

% Log-law refs
plot(ax2_inner, yp_visc, yp_visc, 'k-',  'LineWidth', 1.5, 'DisplayName', '$U^+=y^+$');
plot(ax2_inner, yp_log,  up_log,  'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('$\\kappa=%.2f,\\,B=%.1f$', kappa, B));

fmt_panels(ax2_phys, ax2_inner, yp_max, font_size);
sgtitle('HWA vs PIV — VP (gray) and AK (red)', 'FontWeight', 'bold');

%% ========================================================================
%% FIGURES 3..N+2 — Per x-location: 2 panels, VP+AK PIV only
%% ========================================================================
for k = 1:n_x
    fk = figure('Name', sprintf('Fig%d: VP x=%d / AK x=%d mm', k+2, ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), ...
        'Position', [100+k*25, 120+k*20, 1200, 520]);
    tlk       = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_k_phys = nexttile(tlk); hold on;
    ax_k_in   = nexttile(tlk); hold on;

    add_piv_pair(ax_k_phys, ax_k_in, VP_prof(k), AK_prof(k), n_gray, salmon, marker_skip);

    plot(ax_k_in, yp_visc, yp_visc, 'k-',  'LineWidth', 1.5, 'DisplayName', '$U^+=y^+$');
    plot(ax_k_in, yp_log,  up_log,  'k--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('$\\kappa=%.2f,\\,B=%.1f$', kappa, B));
    fmt_panels(ax_k_phys, ax_k_in, yp_max, font_size);
    sgtitle(sprintf('VP $x=%d$ mm  |  AK $x=%d$ mm', ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), ...
        'Interpreter', 'latex', 'FontWeight', 'bold');
end

%% ========================================================================
%% FIGURE N+3 — Grid: 2 rows x n_x cols (top=inner, bottom=physical)
%% ========================================================================
fig_grid = figure('Name', sprintf('Fig%d: Grid', n_x+3), ...
    'Position', [60 60 370*n_x 880]);
tl_grid = tiledlayout(2, n_x, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_top = gobjects(1, n_x);
ax_bot = gobjects(1, n_x);
for k = 1:n_x
    ax_top(k) = nexttile(tl_grid, k);       hold on;
end
for k = 1:n_x
    ax_bot(k) = nexttile(tl_grid, n_x + k); hold on;
end

for k = 1:n_x
    add_piv_pair(ax_bot(k), ax_top(k), VP_prof(k), AK_prof(k), n_gray, salmon, marker_skip);

    plot(ax_top(k), yp_visc, yp_visc, 'k-',  'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(ax_top(k), yp_log,  up_log,  'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    set(ax_top(k), 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size-1);
    set(ax_bot(k), 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size-1);
    grid(ax_top(k), 'on'); box(ax_top(k), 'on');
    grid(ax_bot(k), 'on'); box(ax_bot(k), 'on');
    xlim(ax_top(k), [1, yp_max * 1.3]);

    title(ax_top(k), sprintf('$x_{VP}=%d$, $x_{AK}=%d$ mm', ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), 'Interpreter', 'latex', 'FontSize', font_size-2);
    xlabel(ax_top(k), '$y^+$',    'Interpreter', 'latex');
    xlabel(ax_bot(k), '$y$ (mm)', 'Interpreter', 'latex');

    if k == 1
        ylabel(ax_top(k), '$U^+$',    'Interpreter', 'latex', 'FontSize', font_size);
        ylabel(ax_bot(k), '$U$ (m/s)', 'Interpreter', 'latex', 'FontSize', font_size);
    end
    legend(ax_top(k), 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', font_size-3);
    legend(ax_bot(k), 'Interpreter', 'latex', 'Location', 'best',      'FontSize', font_size-3);
end

sgtitle('Inner-scaled (top) | Physical (bottom) — VP \circ vs AK \diamond', 'FontWeight', 'bold');
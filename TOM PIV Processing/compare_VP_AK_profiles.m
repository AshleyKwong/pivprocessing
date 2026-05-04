%% compare_VP_AK_profiles.m
% Overlays Virgilio (VP, grayscale) and AK (red scale) BL profiles
% Produces 3 figures:
%   Fig 1 — mega overlay: all x-locations, 2 panels (physical | inner)
%   Fig 2..N+1 — per x-location: 2 panels (physical | inner), 2 lines each
%   Fig N+2 — grid: 2 rows x n_x cols (top=inner, bottom=physical)
clear; close all; clc;

%% ========================================================================
%% USER SETTINGS
%% ========================================================================

%% --- VP (Virgilio) file paths -------------------------------------------
VP_mergedFile  = 'G:\SW 500mm Mean Flow Fields\Mean_Flow_Feilds_AOA_-8.mat';
VP_Cf_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\Cf_minus8.mat';
VP_Re_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\Re_minus8.mat';
VP_PIVsum_file = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\SW_PIV_Summary_-8.mat';

%% --- AK file paths ------------------------------------------------------
AK_mean_field  = 'G:\Y235_AOAN04_AOAFN06_SmallerWindows\merge_instantaneousavg_20260429_083207\merged_meanUV_14loops_20260429_083207.mat';
AK_Cf_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\OFI RESULTS\OFI_Cf_results_20260426_165929.mat';
AK_blSweep     = load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\16x96init\blSweep_20260430_084136.mat').blSweep;
AK_case_name   = 'Case 2';

%% --- Shared parameters --------------------------------------------------
% x-locations must be the same length and correspond index-by-index
VP_x_extract_mm = [ 6850 ];   % VP global x (mm) 5853, 6850, 7755, 9030 
AK_x_extract_mm = [ 550];   % AK local PIV x (mm)  54,  550, 1000, 1000

VP_x_tol_mm   = 5;
AK_nu_air     = 1.5e-5;
VP_nu_air     = 1.51e-5;
VP_Uinf0      = 19.2085724;

n_gray        = 6;       % AK near-wall points to salmon-out
marker_skip   = 1;
font_size     = 13;

%% --- Colour maps --------------------------------------------------------
n_x    = numel(VP_x_extract_mm);   % same as AK

% VP: dark→light gray
vp_gray_dark  = [0.15 0.15 0.15];
vp_gray_light = [0.70 0.70 0.70];
VP_colors = [linspace(vp_gray_dark(1), vp_gray_light(1), n_x)', ...
             linspace(vp_gray_dark(2), vp_gray_light(2), n_x)', ...
             linspace(vp_gray_dark(3), vp_gray_light(3), n_x)'];

% AK: dark→bright red
red_dark   = [0.6, 0.0, 0.0];
red_bright = [1.0, 0.4, 0.4];
AK_colors  = [linspace(red_dark(1), red_bright(1), n_x)', ...
              linspace(red_dark(2), red_bright(2), n_x)', ...
              linspace(red_dark(3), red_bright(3), n_x)'];

salmon = [0.98 0.65 0.60];   % AK near-wall marker colour

%% ========================================================================
%% LOAD VP DATA
%% ========================================================================
fprintf('=== Loading VP data ===\n');
Cf_data    = load(VP_Cf_file);
Re_data    = load(VP_Re_file);
SW_PIVData = load(VP_PIVsum_file);
vp_piv     = load(VP_mergedFile, 'X', 'Y', 'U_mean');

Cf_vec1 = Cf_data.Cf_vec1;
Re_vec1 = Re_data.Re_vec1;

Virgilio_Rex   = Re_vec1(:, 2);
VP_x_OFI_m    = Virgilio_Rex .* VP_nu_air ./ VP_Uinf0;
VP_Cf_raw      = mean(Cf_vec1, 2);
Virgilio_x_PIV     = SW_PIVData.X_m(:);
Virgilio_localUinf = SW_PIVData.U99 .* 0.99;

[~, matched_idx] = min(abs(VP_x_OFI_m(:) - Virgilio_x_PIV(:)'), [], 2);
VP_Uinf_local    = Virgilio_localUinf(matched_idx);
if size(VP_Cf_raw, 1) ~= size(VP_Uinf_local, 1)
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
%% LOAD AK DATA
%% ========================================================================
fprintf('=== Loading AK data ===\n');
AK_U_mean    = load(AK_mean_field).U_hann_mean;
AK_worldX    = load(AK_mean_field).worldX;
AK_worldY    = load(AK_mean_field).worldY;
AK_mag       = (max(AK_worldX(1,:)) - min(AK_worldX(1,:))) / size(AK_worldX, 2);
AK_x_vec     = AK_worldX(1,:);
AK_y_vec     = AK_worldY(:,1);

AK_Cf_raw    = load(AK_Cf_file);
ak_case_idx  = find(strcmp({AK_Cf_raw.results.caseName}, AK_case_name), 1);
if isempty(ak_case_idx)
    error('AK case "%s" not found.', AK_case_name);
end
AK_Cf_data   = AK_Cf_raw.results(ak_case_idx);

%% ========================================================================
%% EXTRACT PROFILES — VP
%% ========================================================================
fprintf('=== Extracting VP profiles ===\n');
VP_prof(n_x) = struct();

for k = 1:n_x
    x_target_mm = VP_x_extract_mm(k);
    x_target_m  = x_target_mm / 1000;

    [dist_x, col_idx] = min(abs(vp_x_mm - x_target_mm));
    if dist_x > VP_x_tol_mm
        warning('VP station %d: %.1f mm away from target', k, dist_x);
    end

    U_col = vp_piv.U_mean(:, col_idx);
    y_col = vp_y_mm;
    valid = ~isnan(U_col) & y_col > 0;
    y_mm  = y_col(valid);
    U_ms  = U_col(valid);
    [y_mm, ys] = sort(y_mm, 'ascend');
    U_ms       = U_ms(ys);

    Cf_raw_here  = interp1(VP_x_OFI_m, VP_Cf_raw,     x_target_m, 'linear', NaN);
    Uinf_loc_k   = interp1(VP_x_OFI_m, VP_Uinf_local, x_target_m, 'linear', NaN);
    u_tau        = sqrt(Cf_raw_here) * VP_Uinf0 * (Uinf_loc_k / VP_Uinf0);
    delta99_here = interp1(VP_x_BL_m, VP_delta99_m, x_target_m, 'linear', NaN);
    Ue_here      = interp1(VP_x_BL_m, VP_Ue_m,      x_target_m, 'linear', NaN);

    y_m    = y_mm / 1000;
    yplus  = y_m * u_tau / VP_nu_air;
    uplus  = U_ms / u_tau;
    VP_prof(k).Ue        = Ue_here;           % ADD this line (m/s)
    VP_prof(k).x_mm      = x_target_mm;
    VP_prof(k).y_mm      = y_mm;
    VP_prof(k).U_ms      = U_ms;
    VP_prof(k).yplus     = yplus;
    VP_prof(k).uplus     = uplus;
    VP_prof(k).u_tau     = u_tau;
    VP_prof(k).delta99   = delta99_here * 1e3;   % mm
    VP_prof(k).Re_tau    = u_tau * delta99_here / VP_nu_air;
    VP_prof(k).color     = VP_colors(k,:);
    VP_prof(k).label     = sprintf('VP  $x=%d$ mm', ...
                               x_target_mm);


    fprintf('  VP x=%4d mm | u_tau=%.4f | delta99=%.2f mm | Re_tau=%.0f\n', ...
        x_target_mm, u_tau, delta99_here*1e3, VP_prof(k).Re_tau);
end

%% ========================================================================
%% EXTRACT PROFILES — AK
%% ========================================================================
fprintf('=== Extracting AK profiles ===\n');
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

    U_profile  = mean(AK_U_mean(:, cols), 2, 'omitnan');
    cum_z      = double(AK_y_vec - 0);
    Um         = double(U_profile);

    valid_mask = ~isnan(Um);
    Um         = Um(valid_mask);
    cum_z      = cum_z(valid_mask);
    [cum_z, sidx] = sort(cum_z, 'ascend');
    Um         = Um(sidx);

    U_inf       = max(Um);
    wall_thresh = 0.0005 * U_inf;
    first_valid = find(Um > wall_thresh, 1, 'first');
    if isempty(first_valid), warning('AK: no valid wall point at x=%.1f', x_target); continue; end

    y_offset   = cum_z(first_valid);
    cum_z_corr = cum_z - y_offset + 0.2;
    trim_idx   = first_valid + 3;

    y_mm_plot  = cum_z_corr(trim_idx:end);
    U_ms_plot  = Um(trim_idx:end);

    delta99    = interp1(AK_blSweep.x_mm, AK_blSweep.delta99_hybrid_mm, x_target, 'linear', 'extrap');
    yplus      = (y_mm_plot / 1000) * utau_piv / AK_nu_air;
    uplus      = U_ms_plot / utau_piv;

    AK_prof(i).x_mm      = x_target;
    AK_prof(i).Ue = interp1(AK_blSweep.x_mm, AK_blSweep.U99, x_target, 'linear', 'extrap');
    AK_prof(i).x_actual  = mean(AK_x_vec(cols));
    AK_prof(i).y_mm      = y_mm_plot;
    AK_prof(i).U_ms      = U_ms_plot;
    AK_prof(i).yplus     = yplus;
    AK_prof(i).uplus     = uplus;
    AK_prof(i).utau      = utau_piv;
    AK_prof(i).delta99   = delta99;
    AK_prof(i).color     = AK_colors(i,:);
    AK_prof(i).label     = sprintf('AK  $x=%d$ mm', x_target);

    fprintf('  AK x=%4d mm | u_tau=%.4f | delta99=%.3f mm\n', ...
        x_target, utau_piv, delta99);
end

%% ========================================================================
%% HELPER — shared plot function
% %% ========================================================================
% function plot_pair(ax_phys, ax_inner, vp, ak, n_gray, salmon, marker_skip)
%     % VP — physical
%     axes(ax_phys); hold on;
%     idx_mk = 1:marker_skip:numel(vp.y_mm);
%     plot(ax_phys, vp.y_mm, vp.U_ms, '-', 'Color', vp.color, ...
%         'LineWidth', 1.5, 'DisplayName', vp.label);
%     plot(ax_phys, vp.y_mm(idx_mk), vp.U_ms(idx_mk), 'o', ...
%         'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 4, ...
%         'HandleVisibility', 'off');
% 
%     % AK — physical (salmon near-wall, red rest)
%     n = numel(ak.y_mm);
%     ig = 1:min(n_gray, n);
%     ic = min(n_gray,n)+1:n;
%     plot(ax_phys, ak.y_mm(ig), ak.U_ms(ig), 'o', 'Color', salmon, ...
%         'MarkerFaceColor', salmon, 'MarkerSize', 4, 'HandleVisibility', 'off');
%     if ~isempty(ic)
%         plot(ax_phys, ak.y_mm(ic), ak.U_ms(ic), '-o', 'Color', ak.color, ...
%             'MarkerFaceColor', ak.color, 'MarkerSize', 4, ...
%             'LineWidth', 1.5, 'DisplayName', ak.label);
%     end
%     set(ax_phys, 'XScale', 'log');
%     xlabel(ax_phys, '$y$ (mm)', 'Interpreter', 'latex');
%     ylabel(ax_phys, '$U$ (m/s)', 'Interpreter', 'latex');
%     grid(ax_phys, 'on'); box(ax_phys, 'on');
%     legend(ax_phys, 'Interpreter', 'latex', 'Location', 'best');
% 
%     % VP — inner
%     axes(ax_inner); hold on;
%     plot(ax_inner, vp.yplus, vp.uplus, '-', 'Color', vp.color, ...
%         'LineWidth', 1.5, 'DisplayName', vp.label);
%     plot(ax_inner, vp.yplus(idx_mk), vp.uplus(idx_mk), 'o', ...
%         'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 4, ...
%         'HandleVisibility', 'off');
% 
%     % AK — inner
%     if ~isempty(ic)
%         plot(ax_inner, ak.yplus(ig),  ak.uplus(ig),  'o', 'Color', salmon, ...
%             'MarkerFaceColor', salmon, 'MarkerSize', 4, 'HandleVisibility', 'off');
%         plot(ax_inner, ak.yplus(ic), ak.uplus(ic), '-o', 'Color', ak.color, ...
%             'MarkerFaceColor', ak.color, 'MarkerSize', 4, ...
%             'LineWidth', 1.5, 'DisplayName', ak.label);
%     end
%     set(ax_inner, 'XScale', 'log');
%     xlabel(ax_inner, '$y^+$',  'Interpreter', 'latex');
%     ylabel(ax_inner, '$U^+$',  'Interpreter', 'latex');
%     grid(ax_inner, 'on'); box(ax_inner, 'on');
%     legend(ax_inner, 'Interpreter', 'latex', 'Location', 'northwest');
% end
%% ========================================================================
%% HELPER — shared plot function  (UPDATED: optional ax_delta)
%% ========================================================================
function plot_pair(ax_phys, ax_inner, vp, ak, n_gray, salmon, marker_skip, ax_delta)

    if nargin < 8, ax_delta = []; end

    idx_mk = 1:marker_skip:numel(vp.y_mm);
    n  = numel(ak.y_mm);
    ig = 1:min(n_gray, n);
    ic = min(n_gray, n)+1:n;

    % ── VP & AK — physical ───────────────────────────────────────────────
    if ~isempty(ax_phys)
        axes(ax_phys); hold on;
        plot(ax_phys, vp.y_mm, vp.U_ms, '-', 'Color', vp.color, ...
            'LineWidth', 1.5, 'DisplayName', vp.label);
        plot(ax_phys, vp.y_mm(idx_mk), vp.U_ms(idx_mk), 'o', ...
            'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 4, ...
            'HandleVisibility', 'off');
        plot(ax_phys, ak.y_mm(ig), ak.U_ms(ig), 'o', 'Color', salmon, ...
            'MarkerFaceColor', salmon, 'MarkerSize', 4, 'HandleVisibility', 'off');
        if ~isempty(ic)
            plot(ax_phys, ak.y_mm(ic), ak.U_ms(ic), '-o', 'Color', ak.color, ...
                'MarkerFaceColor', ak.color, 'MarkerSize', 4, ...
                'LineWidth', 1.5, 'DisplayName', ak.label);
        end
        set(ax_phys, 'XScale', 'log');
        xlabel(ax_phys, '$y$ (mm)',  'Interpreter', 'latex');
        ylabel(ax_phys, '$U$ (m/s)', 'Interpreter', 'latex');
        grid(ax_phys, 'on'); box(ax_phys, 'on');
        legend(ax_phys, 'Interpreter', 'latex', 'Location', 'best');
    end

    % ── VP & AK — inner ──────────────────────────────────────────────────
    if ~isempty(ax_inner)
        axes(ax_inner); hold on;
        plot(ax_inner, vp.yplus, vp.uplus, '-', 'Color', vp.color, ...
            'LineWidth', 1.5, 'DisplayName', vp.label);
        plot(ax_inner, vp.yplus(idx_mk), vp.uplus(idx_mk), 'o', ...
            'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 4, ...
            'HandleVisibility', 'off');
        plot(ax_inner, ak.yplus(ig), ak.uplus(ig), 'o', 'Color', salmon, ...
            'MarkerFaceColor', salmon, 'MarkerSize', 4, 'HandleVisibility', 'off');
        if ~isempty(ic)
            plot(ax_inner, ak.yplus(ic), ak.uplus(ic), '-o', 'Color', ak.color, ...
                'MarkerFaceColor', ak.color, 'MarkerSize', 4, ...
                'LineWidth', 1.5, 'DisplayName', ak.label);
        end
        set(ax_inner, 'XScale', 'log');
        xlabel(ax_inner, '$y^+$', 'Interpreter', 'latex');
        ylabel(ax_inner, '$U^+$', 'Interpreter', 'latex');
        grid(ax_inner, 'on'); box(ax_inner, 'on');
        legend(ax_inner, 'Interpreter', 'latex', 'Location', 'northwest');
    end

    % ── Delta-scaled (optional) ──────────────────────────────────────────
    if ~isempty(ax_delta)
        axes(ax_delta); hold on;

        vp_eta   = (vp.y_mm / 1000) / (vp.delta99 / 1000);
        vp_Unorm = vp.U_ms / vp.Ue;
        plot(ax_delta, vp_eta, vp_Unorm, '-', 'Color', vp.color, ...
            'LineWidth', 1.5, 'DisplayName', vp.label);
        plot(ax_delta, vp_eta(idx_mk), vp_Unorm(idx_mk), 'o', ...
            'Color', vp.color, 'MarkerFaceColor', vp.color, 'MarkerSize', 4, ...
            'HandleVisibility', 'off');

        ak_eta   = (ak.y_mm / 1000) / (ak.delta99 / 1000);
        ak_Unorm = ak.U_ms / ak.Ue;
        plot(ax_delta, ak_eta(ig), ak_Unorm(ig), 'o', 'Color', salmon, ...
            'MarkerFaceColor', salmon, 'MarkerSize', 4, 'HandleVisibility', 'off');
        if ~isempty(ic)
            plot(ax_delta, ak_eta(ic), ak_Unorm(ic), '-o', 'Color', ak.color, ...
                'MarkerFaceColor', ak.color, 'MarkerSize', 4, ...
                'LineWidth', 1.5, 'DisplayName', ak.label);
        end

        xlabel(ax_delta, '$y/\delta_{99}$',        'Interpreter', 'latex');
        ylabel(ax_delta, '$U/U_{\mathrm{e}}$',     'Interpreter', 'latex');
        grid(ax_delta, 'on'); box(ax_delta, 'on');
        legend(ax_delta, 'Interpreter', 'latex', 'Location', 'southeast');
    end
end

%% ========================================================================
%% FIGURE 1 — Mega overlay: all x-locations, 2 panels
%% ========================================================================
fig1 = figure('Name', 'All profiles — VP vs AK', 'Position', [50 50 1300 560]);
tl1  = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1_phys  = nexttile(tl1); hold on; box on; grid on;
ax1_inner = nexttile(tl1); hold on; box on; grid on;

for k = 1:n_x
    plot_pair(ax1_phys, ax1_inner, VP_prof(k), AK_prof(k), n_gray, salmon, marker_skip);
end

set(ax1_phys,  'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
set(ax1_inner, 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
title(ax1_phys,  'Physical',     'Interpreter', 'latex');
title(ax1_inner, 'Inner-scaled', 'Interpreter', 'latex');
sgtitle('VP vs AK — all x-locations', 'FontWeight', 'bold');

%% ========================================================================
%% FIGURES 2..N+1 — Per x-location, 2 panels
%% ========================================================================
for k = 1:n_x
    fig_k = figure('Name', sprintf('x VP=%d mm / AK=%d mm', ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), ...
        'Position', [100+k*30, 100+k*20, 1200, 520]);
    tl_k      = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_k_phys = nexttile(tl_k); hold on; box on; grid on;
    ax_k_in   = nexttile(tl_k); hold on; box on; grid on;

    plot_pair(ax_k_phys, ax_k_in, VP_prof(k), AK_prof(k), n_gray, salmon, marker_skip);

    set(ax_k_phys, 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
    set(ax_k_in,   'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
    title(ax_k_phys,  'Physical',     'Interpreter', 'latex');
    title(ax_k_in,    'Inner-scaled', 'Interpreter', 'latex');
    sgtitle(sprintf('VP $x=%d$ mm  |  AK $x=%d$ mm', ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), ...
        'Interpreter', 'latex', 'FontWeight', 'bold');
end

%% ========================================================================
%% FIGURE N+2 — Grid: 2 rows x n_x cols
%%   Row 1 = inner-scaled,  Row 2 = physical
%%   Each column = one x-location pair
%% ========================================================================
fig_grid = figure('Name', 'Grid — inner (top) / physical (bottom)', ...
    'Position', [80 80 380*n_x 900]);
tl_grid  = tiledlayout(2, n_x, 'TileSpacing', 'compact', 'Padding', 'compact');

% Pre-allocate axes in row-major order expected by tiledlayout
ax_inner_row = gobjects(1, n_x);
ax_phys_row  = gobjects(1, n_x);

for k = 1:n_x
    % Top row: inner-scaled  (tile index k)
    ax_inner_row(k) = nexttile(tl_grid, k);
    hold on; box on; grid on;
end
for k = 1:n_x
    % Bottom row: physical  (tile index n_x + k)
    ax_phys_row(k) = nexttile(tl_grid, n_x + k);
    hold on; box on; grid on;
end

for k = 1:n_x
    plot_pair(ax_phys_row(k), ax_inner_row(k), VP_prof(k), AK_prof(k), ...
        n_gray, salmon, marker_skip);

    set(ax_inner_row(k), 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size-1);
    set(ax_phys_row(k),  'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', font_size-1);

    title(ax_inner_row(k), sprintf('$x_{VP}=%d$, $x_{AK}=%d$ mm', ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), 'Interpreter', 'latex');
end

% Row labels via ylabel on leftmost panels only
ylabel(ax_inner_row(1), '$U^+$',    'Interpreter', 'latex', 'FontSize', font_size);
ylabel(ax_phys_row(1),  '$U$ (m/s)','Interpreter', 'latex', 'FontSize', font_size);

sgtitle('Inner-scaled (top) | Physical (bottom) — VP vs AK', 'FontWeight', 'bold');

%% ========================================================================
%% FIGURE N+3 — Grid: delta-scaled (y/δ₉₉ vs U/Uₑ), 1 row x n_x cols
%% ========================================================================
fig_delta = figure('Name', 'Grid — outer/delta-scaled', ...
    'Position', [80 80 380*n_x 480]);
tl_delta  = tiledlayout(1, n_x, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_delta_row = gobjects(1, n_x);
for k = 1:n_x
    ax_delta_row(k) = nexttile(tl_delta);
    hold on; box on; grid on;
end

for k = 1:n_x
    plot_pair([], [], VP_prof(k), AK_prof(k), n_gray, salmon, marker_skip, ax_delta_row(k));

    set(ax_delta_row(k), 'TickLabelInterpreter', 'latex', 'FontSize', font_size-1);
    title(ax_delta_row(k), sprintf('$x_{\\mathrm{VP}}=%d$, $x_{\\mathrm{AK}}=%d$ mm', ...
        VP_x_extract_mm(k), AK_x_extract_mm(k)), 'Interpreter', 'latex');
end
xline(0.05, 'k--', HandleVisibility= 'off'); 
xline(0.5, 'k--', HandleVisibility= 'off');
set(gca, 'Xscale', 'log')

ylabel(ax_delta_row(1), '$U/U_\mathrm{e}$', 'Interpreter', 'latex', 'FontSize', font_size);
sgtitle('Outer-scaled ($y/\delta$ vs $U/U_\infty$) -- VP vs AK', ...
    'Interpreter', 'latex', 'FontWeight', 'bold');
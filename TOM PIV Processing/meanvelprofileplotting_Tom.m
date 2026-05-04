%% plot_VP_velocity_profiles_standalone.m
clear; close all; clc;

%% ========================================================================
%% USER SETTINGS
%% ========================================================================
%% ---- File paths ---------------------------------------------------------
mergedFile  = 'G:\SW 500mm Mean Flow Fields\Mean_Flow_Feilds_AOA_-8.mat';Cf_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\Cf_minus8.mat';
Re_file     = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\Re_minus8.mat';
PIVsum_file = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\SW_PIV_Summary_-8.mat';
savePath    = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\meanprofiles\';

Virg_Uinf0    = 19.2085724; %19.7;

Virgilio_nuair = 1.51e-5;

x_extract_mm = [5853, 6850, 7755,  9000];
x_tol_mm     = 5;

font_size     = 12;
fig_width_in  = 8;
fig_height_in = 6;
timestamp     = datestr(now, 'yyyymmdd_HHMMSS');

gray_dark   = [0.3 0.3 0.3];
gray_light  = [0.7 0.7 0.7];
marker_type = 'o';
marker_skip = 1;      % must be > 0; controls marker density along profile

%% ========================================================================
%% STEP 1: Load raw VP data
%% ========================================================================
fprintf('Loading VP OFI data...\n');
Cf_data    = load(Cf_file);
Re_data    = load(Re_file);
SW_PIVData = load(PIVsum_file);

Cf_vec1 = Cf_data.Cf_vec1;
Re_vec1 = Re_data.Re_vec1;

%% ========================================================================
%% STEP 2: Reconstruct VP_x_OFI_m, VP_Cf_raw, VP_Uinf_local
%% ========================================================================
Virgilio_Rex   = Re_vec1(:, 2);
VP_x_OFI_m    = Virgilio_Rex .* Virgilio_nuair ./ Virg_Uinf0;
VP_Cf_raw      = mean(Cf_vec1, 2); % this is (utau/Uinf)^2 already... and so this is Cf/2 --> so this is consistent
Virgilio_x_PIV     = SW_PIVData.X_m(:);
Virgilio_localUinf = SW_PIVData.U99 .* 0.99;

[~, matched_idx] = min(abs(VP_x_OFI_m(:) - Virgilio_x_PIV(:)'), [], 2);
VP_Uinf_local    = Virgilio_localUinf(matched_idx);
if size(VP_Cf_raw, 1) ~= size(VP_Uinf_local, 1)
    VP_Uinf_local = VP_Uinf_local';
end

% Outlet reference (kept for reference/diagnostics only)
VP_x_PIV_mm_all = Virgilio_x_PIV .* 1000;
VP_xf_all       = max(VP_x_PIV_mm_all);
VP_Uout_mask    = VP_x_PIV_mm_all >= (VP_xf_all - 4);
Virg_Uinfoutlet = mean(Virgilio_localUinf(VP_Uout_mask), 'omitnan');
fprintf('VP: Uoutlet = %.4f m/s\n', Virg_Uinfoutlet);

% Threshold filter and sort — applied to raw Cf and local Uinf together
% u_tau will be computed as sqrt(VP_Cf_raw/2) * Virg_Uinf0 which is
% reference-invariant (the Uinf rescaling cancels exactly)
cf_threshold  = 0.0004 * 2;
VP_Cf_outlet  = VP_Cf_raw .* 2 .* (Virg_Uinf0 ./ Virg_Uinfoutlet)^2;  % for threshold only
valid_vp      = VP_Cf_outlet >= cf_threshold & ~isnan(VP_Cf_outlet);

VP_x_OFI_m   = VP_x_OFI_m(valid_vp);
VP_Cf_raw     = VP_Cf_raw(valid_vp);
VP_Uinf_local = VP_Uinf_local(valid_vp);

[VP_x_OFI_m, vp_sort] = sort(VP_x_OFI_m, 'ascend');
VP_Cf_raw              = VP_Cf_raw(vp_sort);
VP_Uinf_local          = VP_Uinf_local(vp_sort);

fprintf('VP OFI: %d valid points, x = [%.4f, %.4f] m\n', ...
    numel(VP_x_OFI_m), min(VP_x_OFI_m), max(VP_x_OFI_m));

%% ========================================================================
%% STEP 3: Prepare SW_PIVData BL vectors
%% ========================================================================
VP_x_BL_m   = SW_PIVData.X_m(:);
VP_delta99_m = SW_PIVData.delta_m(:);
VP_Ue_m      = SW_PIVData.U99(:);

[VP_x_BL_m, bl_sort] = sort(VP_x_BL_m, 'ascend');
VP_delta99_m          = VP_delta99_m(bl_sort);
VP_Ue_m               = VP_Ue_m(bl_sort);

%% ========================================================================
%% STEP 4: Load PIV field
%% ========================================================================
fprintf('Loading PIV field...\n');
piv = load(mergedFile, 'X', 'Y', 'U_mean', 'uu_mean');
% X, Y in metres, 554(y) x 9179(x) — convert to mm
x_vec_mm = piv.X(1, :) .* 1e3;    % [1 x 9179]
y_vec_mm = piv.Y(:, 1) .* 1e3;    % [554 x 1]

fprintf('Grid : %d y-points x %d x-points\n', numel(y_vec_mm), numel(x_vec_mm));
fprintf('X range : [%.1f, %.1f] mm\n', min(x_vec_mm), max(x_vec_mm));
fprintf('Y range : [%.3f, %.1f] mm\n', min(y_vec_mm), max(y_vec_mm));
%% ========================================================================
%% STEP 5: Extract profiles
%% ========================================================================
n_prof = numel(x_extract_mm);
profiles(n_prof) = struct();
prof_colors = [linspace(gray_dark(1), gray_light(1), n_prof)', ...
               linspace(gray_dark(2), gray_light(2), n_prof)', ...
               linspace(gray_dark(3), gray_light(3), n_prof)'];

fprintf('\n--- BL parameters at extraction stations ---\n');

for k = 1:n_prof

    x_target_mm = x_extract_mm(k);
    x_target_m  = x_target_mm / 1000;

    % Nearest PIV column
    [dist_x, col_idx] = min(abs(x_vec_mm - x_target_mm));
    if dist_x > x_tol_mm
        warning('Station %d: nearest column %.1f mm is %.1f mm away from target %.1f mm', ...
            k, x_vec_mm(col_idx), dist_x, x_target_mm);
    end

    % Wall-normal profile
    U_col = piv.U_mean(:, col_idx);    % [554 x 1]
    y_col = y_vec_mm;
    valid = ~isnan(U_col) & y_col > 0;
    y_mm  = y_col(valid);
    U_ms  = U_col(valid);
    [y_mm, y_sort] = sort(y_mm, 'ascend');
    U_ms           = U_ms(y_sort);

    % u_tau — reference-invariant form: sqrt(Cf_raw/2) * Uinf0
    % Derivation: Cf_local = Cf_raw*(Uinf0/Uinf_local)^2
    %             u_tau = sqrt(Cf_local/2)*Uinf_local
    %                   = sqrt(Cf_raw/2)*Uinf0   [Uinf_local cancels]
    Cf_raw_here = interp1(VP_x_OFI_m, VP_Cf_raw,     x_target_m, 'linear', NaN); % this is tech Cf/2
    Uinf_loc_k  = interp1(VP_x_OFI_m, VP_Uinf_local, x_target_m, 'linear', NaN);
    u_tau       = sqrt(Cf_raw_here) * Virg_Uinf0 * (Uinf_loc_k/Virg_Uinf0); % check this 

% 
    % delta99 and U_e
    delta99_here = interp1(VP_x_BL_m, VP_delta99_m, x_target_m, 'linear', NaN);
    Ue_here      = interp1(VP_x_BL_m, VP_Ue_m,      x_target_m, 'linear', NaN);

    fprintf('  x = %4d mm | Cf_raw = %.5f | U_inf_local = %.3f m/s | u_tau = %.4f m/s | delta99 = %.2f mm | Re_tau = %.0f\n', ...
        x_target_mm, Cf_raw_here, Uinf_loc_k, u_tau, delta99_here*1e3, ...
        u_tau * delta99_here / Virgilio_nuair);

    % Scaled coordinates
    y_m     = y_mm / 1000;
    y_plus  = y_m * u_tau / Virgilio_nuair;
    U_plus  = U_ms / u_tau;
    y_delta = y_m / delta99_here;
    U_Ue    = U_ms / Ue_here;

    profiles(k).x_target_mm = x_target_mm;
    profiles(k).x_actual_mm = x_vec_mm(col_idx);
    profiles(k).y_mm        = y_mm;
    profiles(k).U_ms        = U_ms;
    profiles(k).y_plus      = y_plus;
    profiles(k).U_plus      = U_plus;
    profiles(k).y_delta     = y_delta;
    profiles(k).U_Ue        = U_Ue;
    profiles(k).Cf_raw      = Cf_raw_here;
    profiles(k).Uinf_local  = Uinf_loc_k;
    profiles(k).u_tau       = u_tau;
    profiles(k).delta99_mm  = delta99_here * 1e3;
    profiles(k).Ue_ms       = Ue_here;
    profiles(k).Re_tau      = u_tau * delta99_here / Virgilio_nuair;
    profiles(k).color       = prof_colors(k, :);
    profiles(k).label       = sprintf('$x = %d$ mm,  $Re_\\tau = %.0f$', ...
                                  x_target_mm, u_tau * delta99_here / Virgilio_nuair);
end

%% ========================================================================
%% STEP 6: Reference lines
%% ========================================================================
yp_ref  = logspace(log10(1), log10(1e4), 300);
kappa   = 0.41;
B       = 5.0;
Up_log  = (1/kappa) .* log(yp_ref) + B;
Up_visc = yp_ref;

%% ========================================================================
%% STEP 7: Figure 1 — Inner scaling (U+ vs y+)
%% ========================================================================
figure('Name', 'VP Inner-Scaled Profiles', 'NumberTitle', 'off');
ax_in = axes; hold on;

% semilogx(ax_in, yp_ref(yp_ref <= 12), Up_visc(yp_ref <= 12), ...
%     'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
% semilogx(ax_in, yp_ref(yp_ref >= 30), Up_log(yp_ref >= 30), ...
%     'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

for k = 1:n_prof
    idx_mk = 1:marker_skip:numel(profiles(k).y_plus);
    % Line carries the legend entry
    semilogx(ax_in, profiles(k).y_plus, profiles(k).U_plus, '-', ...
        'Color',        profiles(k).color, ...
        'LineWidth',    1.5, ...
        'DisplayName',  profiles(k).label);
    % Markers — no legend entry to avoid duplication
    semilogx(ax_in, profiles(k).y_plus(idx_mk), profiles(k).U_plus(idx_mk), ...
        'LineStyle',        'none', ...
        'Marker',           marker_type, ...
        'MarkerSize',       5, ...
        'MarkerFaceColor',  profiles(k).color, ...
        'MarkerEdgeColor',  profiles(k).color, ...
        'HandleVisibility', 'off');
end

% text(ax_in, 300, (1/kappa)*log(300) + B + 1.8, ...
%     sprintf('$(1/%.2f)\\ln y^+ + %.1f$', kappa, B), ...
%     'Interpreter', 'latex', 'FontSize', font_size-1, 'Color', [0.4 0.4 0.4]);

xlabel(ax_in, '$y^+$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel(ax_in, '$U^+$', 'Interpreter', 'latex', 'FontSize', font_size);
set(ax_in, 'XScale', 'log');
ax_in.FontSize             = font_size;
ax_in.TickLabelInterpreter = 'latex';
grid(ax_in, 'on');
legend(ax_in, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', font_size);

% drawnow;
% set(gcf, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
% set(gcf, 'PaperUnits', 'inches', 'PaperSize', [fig_width_in fig_height_in], ...
%     'PaperPosition', [0 0 fig_width_in fig_height_in]);
% figName_in = fullfile(savePath, sprintf('VP_profiles_inner_%s.pdf', timestamp));
% exportgraphics(gcf, figName_in, 'ContentType', 'vector', 'BackgroundColor', 'white');
% fprintf('\nFigure saved: %s\n', figName_in);

%% ========================================================================
%% STEP 8: Figure 2 — Outer scaling (U/Ue vs y/delta99)
%% ========================================================================
figure('Name', 'VP Outer-Scaled Profiles', 'NumberTitle', 'off');
ax_out = axes; hold on;

for k = 1:n_prof
    idx_mk = 1:marker_skip:numel(profiles(k).y_delta);
    % Line carries the legend entry

    plot(ax_out, profiles(k).y_delta, profiles(k).U_Ue, '-', ...
        'Color',       profiles(k).color, ...
        'LineWidth',   1.5, ...
        'DisplayName', profiles(k).label);
    % Markers
    plot(ax_out, profiles(k).y_delta(idx_mk), profiles(k).U_Ue(idx_mk), ...
        'LineStyle',        'none', ...
        'Marker',           marker_type, ...
        'MarkerSize',       5, ...
        'MarkerFaceColor',  profiles(k).color, ...
        'MarkerEdgeColor',  profiles(k).color, ...
        'HandleVisibility', 'off');
end

xlabel(ax_out, '$y / \delta_{99}$', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel(ax_out, '$U / U_e$',         'Interpreter', 'latex', 'FontSize', font_size);
set(ax_out, 'XScale', 'log');
ax_out.FontSize             = font_size;
ax_out.TickLabelInterpreter = 'latex';
grid(ax_out, 'on');
legend(ax_out, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', font_size);
%% ========================================================================
%% STEP 10: Turbulence intensity profiles
%% ========================================================================

%% ---- Figure 3: Inner scaling (u_rms/u_tau vs y+) -----------------------
figure('Name', 'VP Turbulence Intensity — Inner Scaled', 'NumberTitle', 'off');
ax_ti_in = axes; hold on;

for k = 1:n_prof

    x_target_mm = x_extract_mm(k);
    [~, col_idx] = min(abs(x_vec_mm - x_target_mm));

    uu_col = piv.uu_mean(:, col_idx);    % [554 x 1], variance [m^2/s^2]
    y_col  = y_vec_mm;

    valid  = ~isnan(uu_col) & y_col > 0 & uu_col >= 0;
    y_mm   = y_col(valid);
    uu_ms2 = uu_col(valid);

    [y_mm, y_sort] = sort(y_mm, 'ascend');
    uu_ms2         = uu_ms2(y_sort);

    u_rms     = sqrt(uu_ms2);
    y_m       = y_mm / 1e3;
    y_plus    = y_m * profiles(k).u_tau / Virgilio_nuair;
    urms_plus = u_rms / profiles(k).u_tau;

    idx_mk = 1:marker_skip:numel(y_plus);

    semilogx(ax_ti_in, y_plus, urms_plus, '-', ...
        'Color',            profiles(k).color, ...
        'LineWidth',        1.5, ...
        'DisplayName',      profiles(k).label);
    semilogx(ax_ti_in, y_plus(idx_mk), urms_plus(idx_mk), ...
        'LineStyle',        'none', ...
        'Marker',           marker_type, ...
        'MarkerSize',       5, ...
        'MarkerFaceColor',  profiles(k).color, ...
        'MarkerEdgeColor',  profiles(k).color, ...
        'HandleVisibility', 'off');
end

xlabel(ax_ti_in, '$y^+$',                       'Interpreter', 'latex', 'FontSize', font_size);
ylabel(ax_ti_in, '$u_{\mathrm{rms}} / u_\tau$',  'Interpreter', 'latex', 'FontSize', font_size);
set(ax_ti_in, 'XScale', 'log');
ax_ti_in.FontSize             = font_size;
ax_ti_in.TickLabelInterpreter = 'latex';
grid(ax_ti_in, 'on');
legend(ax_ti_in, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', font_size);

% drawnow;
% set(gcf, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
% set(gcf, 'PaperUnits', 'inches', 'PaperSize', [fig_width_in fig_height_in], ...
%     'PaperPosition', [0 0 fig_width_in fig_height_in]);
% figName_ti_in = fullfile(savePath, sprintf('VP_TI_inner_%s.pdf', timestamp));
% exportgraphics(gcf, figName_ti_in, 'ContentType', 'vector', 'BackgroundColor', 'white');
% fprintf('Figure saved: %s\n', figName_ti_in);

%% ---- Figure 4: Outer scaling (u_rms/U_e vs y/delta99) ------------------
figure('Name', 'VP Turbulence Intensity — Outer Scaled', 'NumberTitle', 'off');
ax_ti_out = axes; hold on;

for k = 1:n_prof

    x_target_mm = x_extract_mm(k);
    [~, col_idx] = min(abs(x_vec_mm - x_target_mm));

    uu_col = piv.uu_mean(:, col_idx);
    y_col  = y_vec_mm;

    valid  = ~isnan(uu_col) & y_col > 0 & uu_col >= 0;
    y_mm   = y_col(valid);
    uu_ms2 = uu_col(valid);

    [y_mm, y_sort] = sort(y_mm, 'ascend');
    uu_ms2         = uu_ms2(y_sort);

    u_rms   = sqrt(uu_ms2);
    y_m     = y_mm / 1e3;
    y_delta = y_m / (profiles(k).delta99_mm / 1e3);
    urms_Ue = u_rms / profiles(k).Ue_ms;

    idx_mk = 1:marker_skip:numel(y_delta);

    plot(ax_ti_out, y_delta, urms_Ue, '-', ...
        'Color',            profiles(k).color, ...
        'LineWidth',        1.5, ...
        'DisplayName',      profiles(k).label);
    plot(ax_ti_out, y_delta(idx_mk), urms_Ue(idx_mk), ...
        'LineStyle',        'none', ...
        'Marker',           marker_type, ...
        'MarkerSize',       5, ...
        'MarkerFaceColor',  profiles(k).color, ...
        'MarkerEdgeColor',  profiles(k).color, ...
        'HandleVisibility', 'off');
end

xlabel(ax_ti_out, '$y / \delta_{99}$',          'Interpreter', 'latex', 'FontSize', font_size);
ylabel(ax_ti_out, '$u_{\mathrm{rms}} / U_e$',   'Interpreter', 'latex', 'FontSize', font_size);
set(ax_ti_out, 'XScale', 'log');
ax_ti_out.FontSize             = font_size;
ax_ti_out.TickLabelInterpreter = 'latex';
grid(ax_ti_out, 'on');
legend(ax_ti_out, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', font_size);

drawnow;
set(gcf, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [fig_width_in fig_height_in], ...
    'PaperPosition', [0 0 fig_width_in fig_height_in]);
figName_ti_out = fullfile(savePath, sprintf('VP_TI_outer_%s.pdf', timestamp));
exportgraphics(gcf, figName_ti_out, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Figure saved: %s\n', figName_ti_out);
drawnow;
set(gcf, 'Units', 'inches', 'Position', [1 1 fig_width_in fig_height_in]);
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [fig_width_in fig_height_in], ...
    'PaperPosition', [0 0 fig_width_in fig_height_in]);
figName_out = fullfile(savePath, sprintf('VP_profiles_outer_%s.pdf', timestamp));
exportgraphics(gcf, figName_out, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Figure saved: %s\n', figName_out);


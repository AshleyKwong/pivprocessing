% =========================================================================
% plot_turbstats.m
%
% Visualises single-point turbulence statistics computed by
% compute_quad_spectra_skewkurt.m:
%
%   Figure 1 — Quadrant fractions Q1-Q4 vs y/delta (all xref overlaid)
%   Figure 2 — Reynolds shear stress -<u'v'> vs y/delta
%   Figure 3 — u'rms, v'rms vs y/delta
%   Figure 4 — Premultiplied spectra kx*Phi_uu(lambda_x/delta, y/delta)
%              one subplot per xref
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

turbStatsFile = ['D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\turbulence_statistics_PIV\turbstats_case1.mat'];
blSweepFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case1_PIVresults\blSweep_20260419_193106.mat'];

% Normalise y axis by delta99
normalise_y = true;

% y/delta range to plot — set [] to use full range
% e.g. [0.01 1.2] to focus on boundary layer interior
ydelta_lim = [0.01 1.5];

% Spectral plot settings
lambda_lim  = [];     % x axis limits for lambda_x/delta — [] = auto
kxPuu_lim   = [];     % colour axis limits for spectrum  — [] = auto

%% ===================== LOAD =====================

S    = load(turbStatsFile, 'results', 'xref_targets');
res  = S.results;
xref_targets = S.xref_targets;
nXref = numel(res);

B      = load(blSweepFile, 'blSweep');
bl_x   = B.blSweep.x_mm;
bl_d99 = B.blSweep.delta99_hybrid_mm;
valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
bl_x_clean   = bl_x(valid_bl);
bl_d99_clean = bl_d99(valid_bl);

delta_vals = interp1(bl_x_clean, bl_d99_clean, xref_targets, 'linear', NaN);

fprintf('delta99 at each xref:\n');
for i = 1:nXref
    fprintf('  x = %.1f mm  ->  delta99 = %.2f mm\n', ...
        xref_targets(i), delta_vals(i));
end

%% ===================== COLOUR SCHEME =====================

% xref: light blue (upstream) -> dark navy (downstream)
cmapX = zeros(nXref, 3);
for i = 1:nXref
    t = (i-1) / max(nXref-1, 1);
    cmapX(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
end

% Quadrant colours — standard convention
cQ1 = [0.60 0.60 0.60];   % Q1 outward   — grey
cQ2 = [0.13 0.47 0.71];   % Q2 ejection  — blue
cQ3 = [0.80 0.80 0.80];   % Q3 inward    — light grey
cQ4 = [0.85 0.33 0.10];   % Q4 sweep     — red/orange

%% ===================== HELPER: GET y AXIS =====================

    function [y_plot, ylabel_str] = get_yaxis(y_mm, delta, normalise)
        if normalise && isfinite(delta)
            y_plot    = y_mm / delta;
            ylabel_str = 'y / \delta_{99}';
        else
            y_plot    = y_mm;
            ylabel_str = 'y  (mm)';
        end
    end

% %% ===================== FIGURE 1: ALL QUADRANT RS CONTRIBUTIONS =====================
% 
% figure('Color','w','Position',[100 80 700 600], ...
%     'Name','Quadrant RS contributions vs y/delta');
% ax = axes; hold on; box on;
% 
% for iX = 1:nXref
%     if ~res(iX).valid, continue; end
% 
%     [y_plot, yl_str] = get_yaxis(res(iX).y_mm, delta_vals(iX), normalise_y);
%     col = cmapX(iX,:);
% 
%     RS_tot  = res(iX).RS;
%     safe_RS = max(abs(RS_tot), eps);
%     iH0     = find(res(iX).H_vals == 0, 1);
% 
%     S_Q1 = -res(iX).RS_Q1(:, iH0) ./ safe_RS;
%     S_Q2 = -res(iX).RS_Q2(:, iH0) ./ safe_RS;
%     S_Q3 = -res(iX).RS_Q3(:, iH0) ./ safe_RS;
%     S_Q4 = -res(iX).RS_Q4(:, iH0) ./ safe_RS;
% 
%     S_sum = S_Q1 + S_Q2 + S_Q3 + S_Q4;
%     fprintf('xref = %.0f mm | S_Q sum range: [%.3f, %.3f]\n', ...
%         xref_targets(iX), min(S_sum,[],'omitnan'), max(S_sum,[],'omitnan'));
% 
%     % x = y/delta, y = statistic
%     plot(ax, y_plot, S_Q1, '^:', ...
%         'Color',           cQ1, ...
%         'MarkerFaceColor', cQ1, ...
%         'MarkerEdgeColor', col, ...
%         'MarkerSize',      4, ...
%         'LineWidth',       1.0, ...
%         'DisplayName', sprintf('Q1 (outward)  x_{ref}=%.0f mm', xref_targets(iX)));
% 
%     plot(ax, y_plot, S_Q2, 'o-', ...
%         'Color',           cQ2, ...
%         'MarkerFaceColor', cQ2, ...
%         'MarkerEdgeColor', col, ...
%         'MarkerSize',      5, ...
%         'LineWidth',       1.8, ...
%         'DisplayName', sprintf('Q2 (ejection)  x_{ref}=%.0f mm', xref_targets(iX)));
% 
%     plot(ax, y_plot, S_Q3, 'v:', ...
%         'Color',           cQ3, ...
%         'MarkerFaceColor', cQ3, ...
%         'MarkerEdgeColor', col, ...
%         'MarkerSize',      4, ...
%         'LineWidth',       1.0, ...
%         'DisplayName', sprintf('Q3 (inward)  x_{ref}=%.0f mm', xref_targets(iX)));
% 
%     plot(ax, y_plot, S_Q4, 's--', ...
%         'Color',           cQ4, ...
%         'MarkerFaceColor', cQ4, ...
%         'MarkerEdgeColor', col, ...
%         'MarkerSize',      5, ...
%         'LineWidth',       1.8, ...
%         'DisplayName', sprintf('Q4 (sweep)  x_{ref}=%.0f mm', xref_targets(iX)));
% end
% 
% yline(ax, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
% 
% set(ax, 'XScale', 'log');
% if ~isempty(ydelta_lim) && normalise_y, xlim(ax, ydelta_lim); end
% xlabel(ax, yl_str);
% ylabel(ax, 'S_{Q_i} = -\langle u''v'' \rangle_{Q_i} / -\langle u''v'' \rangle');
% title(ax, 'Quadrant RS contributions  [H = 0]');
% legend(ax, 'Location', 'best', 'FontSize', 7);
% grid(ax, 'on');

%% ===================== FIGURE 2: Q2 - Q4 IMBALANCE =====================
markerStyles = {'o', 's', '^', 'd'};   % circle, square, triangle, diamond
figure('Color','w','Position',[150 80 700 600], ...
    'Name','Ejection-sweep imbalance vs y/delta');
ax2 = axes; hold on; box on;

for iX = 1:nXref
    if ~res(iX).valid, continue; end

    [y_plot, yl_str] = get_yaxis(res(iX).y_mm, delta_vals(iX), normalise_y);
    col = cmapX(iX,:);

    RS_tot  = res(iX).RS;
    safe_RS = max(abs(RS_tot), eps);
    iH0     = find(res(iX).H_vals == 0, 1);

    S_Q2 = -res(iX).RS_Q2(:, iH0) ./ safe_RS;
    S_Q4 = -res(iX).RS_Q4(:, iH0) ./ safe_RS;
    dQ   = S_Q2 - S_Q4;

    plot(ax2, y_plot, dQ, 'o-', ...
        'Color',           col, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize',      5, ...
        'LineWidth',       1.8, ...
        'DisplayName', sprintf('x_{ref} = %.0f mm', xref_targets(iX)));
end

yline(ax2, 0, 'k--', 'LineWidth', 1.0, ...
    'Label', 'Q2 = Q4', 'LabelHorizontalAlignment', 'left');

set(ax2, 'XScale', 'log');
if ~isempty(ydelta_lim) && normalise_y, xlim(ax2, ydelta_lim); end
xlabel(ax2, yl_str);
ylabel(ax2, 'S_{Q2} - S_{Q4}');
title(ax2, 'Ejection-sweep RS imbalance  (S_{Q2} - S_{Q4})  [H = 0]');
legend(ax2, 'Location', 'best');
grid(ax2, 'on');

%% ===================== FIGURE 3: REYNOLDS SHEAR STRESS =====================
% xref colours and markers — one per xref
cmapX = [0.13 0.47 0.71;   % blue      — x_ref = 350
         0.47 0.67 0.19;   % green     — x_ref = 550
         0.85 0.33 0.10;   % orange    — x_ref = 711
         0.50 0.00 0.50];  % purple    — x_ref = 1000
figure('Color','w','Position',[200 80 700 600], ...
    'Name','Reynolds shear stress vs y/delta');
ax3 = axes; hold on; box on;

for iX = 1:nXref
    if ~res(iX).valid, continue; end

    [y_plot, yl_str] = get_yaxis(res(iX).y_mm, delta_vals(iX), normalise_y);
    col = cmapX(iX,:);

    plot(ax3, y_plot, res(iX).RS, '-', ...
        'Color',           col, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize',      5, ...
        'LineWidth',       1.8, ...
        'DisplayName', sprintf('x_{ref} = %.0f mm', xref_targets(iX)));
end

set(ax3, 'XScale', 'log');
if ~isempty(ydelta_lim) && normalise_y, xlim(ax3, ydelta_lim); end
xlabel(ax3, yl_str);
ylabel(ax3, '-\langle u''v'' \rangle  (m^2/s^2)');
title(ax3, 'Reynolds shear stress -\langle u''v'' \rangle');
legend(ax3, 'Location', 'best');
grid(ax3, 'on');

%% ===================== FIGURE 4: u'rms, v'rms =====================

figure('Color','w','Position',[250 80 700 600], ...
    'Name','RMS profiles vs y/delta');
ax4 = axes; hold on; box on;

for iX = 1:nXref
    if ~res(iX).valid, continue; end

    [y_plot, yl_str] = get_yaxis(res(iX).y_mm, delta_vals(iX), normalise_y);
    col = cmapX(iX,:);

    plot(ax4, y_plot, res(iX).urms, 'o-', ...
        'Color',           col, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize',      5, ...
        'LineWidth',       1.8, ...
        'DisplayName', sprintf('u''_{rms}  x_{ref}=%.0f mm', xref_targets(iX)));

    plot(ax4, y_plot, res(iX).vrms, 's--', ...
        'Color',           col, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize',      5, ...
        'LineWidth',       1.2, ...
        'DisplayName', sprintf('v''_{rms}  x_{ref}=%.0f mm', xref_targets(iX)));
end

set(ax4, 'XScale', 'log');
if ~isempty(ydelta_lim) && normalise_y, xlim(ax4, ydelta_lim); end
xlabel(ax4, yl_str);
ylabel(ax4, 'u''_{rms},  v''_{rms}  (m/s)');
title(ax4, 'Turbulence intensity profiles');
legend(ax4, 'Location', 'best', 'FontSize', 8);
grid(ax4, 'on');

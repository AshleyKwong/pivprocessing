% =========================================================================
% diagnose_ydelta_availability.m
%
% Three-panel subplot:
%   Panel 1 — dCp/dx vs x  (pressure gradient driver)
%   Panel 2 — delta99 vs x (boundary layer response)
%   Panel 3 — y/delta availability at each x_target
%
% Matching x_target markers across all three panels lets you visually
% identify any lag between the pressure gradient and the BL response.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

y_vals    = [2.0, 2.6, 4.6, 8.0, 14.0, 24.5, 89.1];   % mm
x_targets = [300, 500, 700, 900, 1100, 1300];           % mm

blSweepFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\' ...
    'MATLAB\Experimental Campaign 1\PIV\PIV_results\' ...
    'Case2_PIVresults\blSweep_20260415_125916.mat'];

pivResultsFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\' ...
    'MATLAB\Experimental Campaign 1\PIV\PIV_results\' ...
    'Case2_PIVresults\PIV_vs_taps_results.mat'];

%% ===================== LOAD DATA =====================

% --- delta99 ---
B      = load(blSweepFile, 'blSweep');
bl_x   = B.blSweep.x_mm;
bl_d99 = B.blSweep.delta99_hybrid_mm;

valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
bl_x_clean   = bl_x(valid_bl);
bl_d99_clean = bl_d99(valid_bl);

% --- dCp/dx ---
P         = load(pivResultsFile, 'piv_results');
pr        = P.piv_results;
dCp_x_mm  = pr.x_m * 1000;            % convert m -> mm
dCp_dx    = pr.dCp_dxd;

valid_cp     = isfinite(dCp_x_mm) & isfinite(dCp_dx);
dCp_x_clean  = dCp_x_mm(valid_cp);
dCp_dx_clean = dCp_dx(valid_cp);

% --- Interpolate delta99 at x_targets ---
nX = numel(x_targets);
delta_at_xtarget = interp1(bl_x_clean, bl_d99_clean, x_targets, ...
    'linear', NaN);

%% ===================== COMPUTE y/delta TABLE =====================

nY           = numel(y_vals);
ydelta_table = nan(nY, nX);
for iY = 1:nY
    for iX = 1:nX
        if isfinite(delta_at_xtarget(iX))
            ydelta_table(iY, iX) = y_vals(iY) / delta_at_xtarget(iX);
        end
    end
end

%% ===================== PRINT TABLE =====================

fprintf('\n=== y/delta availability ===\n\n');
fprintf('%-12s', 'y_ref (mm)');
for iX = 1:nX
    fprintf('  x=%6.0fmm', x_targets(iX));
end
fprintf('\n%s\n', repmat('-', 1, 12 + nX*14));
for iY = 1:nY
    fprintf('%-12.1f', y_vals(iY));
    for iX = 1:nX
        if isfinite(ydelta_table(iY, iX))
            fprintf('  %10.3f', ydelta_table(iY, iX));
        else
            fprintf('  %10s', 'NaN');
        end
    end
    fprintf('\n');
end

fprintf('\ndelta99 at each x_target:\n');
for iX = 1:nX
    fprintf('  x = %6.1f mm  ->  delta99 = %.2f mm\n', ...
        x_targets(iX), delta_at_xtarget(iX));
end

%% ===================== COLOUR SCHEME =====================

cmapX = zeros(nX, 3);
for i = 1:nX
    t = (i-1) / max(nX-1, 1);
    cmapX(i,:) = (1-t)*[0.95 0.70 0.30] + t*[0.50 0.05 0.05];
end

%% ===================== FIGURE: 3-panel subplot =====================

figure('Color','w', 'Position', [100 50 900 850], ...
    'Name', 'Pressure gradient, delta99 and y/delta availability');

%% --- Panel 1: dCp/dx ---
ax1 = subplot(3, 1, 1);
hold(ax1, 'on'); box(ax1, 'on');

plot(ax1, dCp_x_clean, dCp_dx_clean, 'k-', 'LineWidth', 1.5);
yline(ax1, 0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

for iX = 1:nX
    dcp_at_x = interp1(dCp_x_clean, dCp_dx_clean, x_targets(iX), ...
        'linear', NaN);
    if isfinite(dcp_at_x)
        plot(ax1, x_targets(iX), dcp_at_x, 'v', ...
            'MarkerFaceColor', cmapX(iX,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize',      8, ...
            'HandleVisibility', 'off');
    end
    xline(ax1, x_
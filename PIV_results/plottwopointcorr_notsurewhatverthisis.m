

% plot_R_contour_with_dCpdx.m
% Plots two R fields (top and bottom) with dCp/dx in the middle.
% R subplots in delta-normalised separation coords, linked to dCp/dx x-axis.

clear; clc; close all;

%% ======== USER OPTIONS ================================================
RFiles = {
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\R_uu_xref450.0_yref6.0_20260310_143454.mat';
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\R_uu_xref1050.0_yref26.0_20260310_143454.mat'

    };
%    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\R_uu_xref750.0_yref16.6_20260310_143454.mat';

subplotLabels = {'FPG', 'APG', 'ZPG'};
delta_mm      = [30, 83, 136];     % boundary layer thickness per case (mm)
climSym       = 0.8;
nLevels       = 5;
load('C:\Users\ak1u24\OneDrive - University of Southampton\Thomas Preskett Old Work\TOM MODEL SHI\Pressure Data - Upload\case_pressuredata.mat')

% --- dCp/dx USER INPUT ---
% x locations in mm (streamwise), same coordinate system as worldX_merged
x_dCpdx   = (case_pdata(2).xloc-720).*10;         % e.g. [600, 650, 700, 750, 800, 850, 900] in cm  <-- FILL IN
dCpdx_vals = case_pdata(2).mean_dpdx;        % e.g. [-0.02, -0.01, 0, 0.01, 0.02, 0.01, 0.0]  <-- FILL IN
% ======================================================================

%% LOAD R DATA (needed to get xr values for axis linking)
xr_vals = zeros(1, 2);
yr_vals = zeros(1, 2);
for n = 1:2
    tmp = load(RFiles{n}, 'xr', 'yr');
    xr_vals(n) = tmp.xr;
    yr_vals(n) = tmp.yr;
end

%% SET UP FIGURE AND SUBPLOTS
figure('Units', 'centimeters', 'Position', [2 2 20 26]);

ax1 = subplot(3, 1, 1);   % FPG R field
ax2 = subplot(3, 1, 2);   % dCp/dx
ax3 = subplot(3, 1, 3);   % APG R field


%% SUBPLOT 1 & 3 — R FIELDS
axR = {ax1, ax3};
for n = 1:2
    d   = load(RFiles{n}, 'R', 'worldX_merged', 'worldY_merged', 'xr', 'yr');
    R   = d.R;
    x   = d.worldX_merged(1, :);
    y   = d.worldY_merged(:, 1);
    xr  = d.xr;
    yr  = d.yr;
    del = delta_mm(n);

    % Separation coords normalised by delta
    dx = (x - xr) ./ del;
    dy = (y - yr)  ./ del;

    axes(axR{n}); %#ok<LAXES>
    contourf(dx, dy, R, nLevels, 'LineColor', 'none');
    colormap(axR{n}, jet(10));
    cb = colorbar;
    cb.Label.String = 'R_{uu}';
    cb.FontSize     = 10;
    clim([-climSym climSym]);
    hold on;
    scatter(0, 0, 80, 'k', 'filled', 'MarkerEdgeColor', 'w', ...
        'DisplayName', sprintf('(%.1f, %.1f) mm', xr, yr));

    xline(0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.5, 'DisplayName','off');
    yline(0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.5, 'DisplayName','off');
    ylabel('\Deltay/\delta', 'FontSize', 11);
    title(sprintf('%s  |  R_{uu}  |  x_{ref} = %.1f mm, y_{ref} = %.1f mm, \\delta = %d mm', ...
        subplotLabels{n}, xr, yr, del), 'FontSize', 12);
    legend('Location', 'northeast', 'FontSize', 9);
    grid off; box on;
    axis image;

    % --- Dual x-axis tick labels ---
    % Bottom ticks show Δx/δ (already the axis units)
    % Add a second x-axis on top showing raw mm via manual tick labels
    xtk     = get(axR{n}, 'XTick');              % current Δx/δ ticks
    xtk_mm  = xtk .* del + xr;                   % convert back to mm
    % Replace bottom tick labels with Δx/δ values (already correct)
    % % Add annotation-style top labels in mm — using a twin axis
    % axTop = axes('Position', get(axR{n}, 'Position'), ...
    %     'XAxisLocation', 'top', ...
    %     'YAxisLocation', 'right', ...
    %     'Color', 'none', ...
    %     'XLim', get(axR{n}, 'XLim'), ...
    %     'YLim', get(axR{n}, 'YLim'), ...
    %     'XTick', xtk, ...
    %     'XTickLabel', arrayfun(@(v) sprintf('%.0f', v), xtk_mm, 'UniformOutput', false), ...
    %     'YTick', [], ...
    %     'FontSize', 8);
    % xlabel(axTop, 'x (mm)', 'FontSize', 9);
    axes(axR{n});   % return focus to main axis
end
xlabel(axR{2}, '\Deltax/\delta', 'FontSize', 11);   % only bottom subplot gets x label

%% SUBPLOT 2 — dCp/dx
axes(ax2);
if ~isempty(x_dCpdx) && ~isempty(dCpdx_vals)
    plot(x_dCpdx, dCpdx_vals, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    hold on;
    yline(0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.5);
    % Mark the two reference x locations
    xline(xr_vals(1)-100, '--', 'Color', [0.2 0.5 0.9], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('x_{ref,FPG} = %.0f mm', xr_vals(1)));
    xline(xr_vals(2)-100, '--', 'Color', [0.9 0.3 0.2], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('x_{ref,APG} = %.0f mm', xr_vals(2)));
    legend('Location', 'best', 'FontSize', 9);
else
    text(0.5, 0.5, 'INSERT dCp/dx DATA', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');
end
ylabel('dC_p/dx', 'FontSize', 11);
xlabel('x (mm)', 'FontSize', 11);
title('Pressure Gradient', 'FontSize', 12);
grid on; box on;

%% LINK x-AXES (dCp/dx mm axis linked to R subplot mm axis via XLim)
% The R subplots use Δx/δ so we link by setting matching XLim on ax2
% based on the delta-space limits of subplot 1
xlim_dx  = get(ax1, 'XLim');                         % in Δx/δ units of FPG
xlim_mm1 = xlim_dx .* delta_mm(1) + xr_vals(1);      % convert to mm
set(ax2, 'XLim', xlim_mm1);

%% SAVE
[fPath, ~, ~] = fileparts(RFiles{1});
tstamp  = datestr(now, 'yyyymmdd_HHMMSS');
outName = fullfile(fPath, sprintf('R_uu_dCpdx_comparison_%s.png', tstamp));
saveas(gcf, outName);
fprintf('Saved → %s\n', outName);
 
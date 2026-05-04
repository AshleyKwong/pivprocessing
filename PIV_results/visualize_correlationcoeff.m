% plot_R_contour.m
% Plots two R fields as stacked subplots (2x1), axes in delta-normalised separation coordinates.

clear; clc; close all;

%% ======== USER OPTIONS ================================================
RFiles = {
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\R_uu_xref450.0_yref6.0_20260310_143454.mat';
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\R_uu_xref750.0_yref16.6_20260310_143454.mat';
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\R_uu_xref1050.0_yref26.0_20260310_143454.mat'
};
subplotLabels  = {'FPG', 'APG', 'ZPG'};
delta_mm       = [23, 86 , 136];     % boundary layer thickness per case (mm)
climSym        = 1.0;
nLevels        = 10;
% ======================================================================

%% PLOT
figure('Units', 'centimeters', 'Position', [2 2 18 20]);

for n = 1:length(RFiles)

    d   = load(RFiles{n}, 'R', 'worldX_merged', 'worldY_merged', 'xr', 'yr');
    R   = d.R;
    x   = d.worldX_merged(1, :);   % (1 x Nx) in mm
    y   = d.worldY_merged(:, 1);   % (Ny x 1) in mm
    xr  = d.xr;
    yr  = d.yr;
    del = delta_mm(n);

    % Separation coordinates normalised by delta
    dx = (x - xr) ./ del;         % (1 x Nx)  Δx/δ
    dy = (y - yr)  ./ del;         % (Ny x 1)  Δy/δ

    subplot(2, 1, n);
    contourf(dx, dy, R, nLevels, 'LineColor', 'none');
    colormap(jet(10));
    cb = colorbar;
    cb.Label.String = 'R_{uu}';
    cb.FontSize     = 11;
    clim([-climSym climSym]);

    hold on;
    % Reference point is always at (0,0) in separation coords
    scatter(0, 0, 80, 'k', 'filled', 'MarkerEdgeColor', 'w', ...
        'DisplayName', sprintf('(x_{ref}, y_{ref}) = (%.1f, %.1f) mm', xr, yr));
    xline(0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.5);
    yline(0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.5);

    xlabel('\Deltax/\delta', 'FontSize', 12);
    ylabel('\Deltay/\delta', 'FontSize', 12);
    title(sprintf('%s  |  R_{uu}  |  x_{ref} = %.1f mm, y_{ref} = %.1f mm, \\delta = %d mm', ...
        subplotLabels{n}, xr, yr, del), 'FontSize', 13);
    legend('Location', 'northeast', 'FontSize', 10);
    grid off; box on;
    axis image;

end

%% SAVE
[fPath, ~, ~] = fileparts(RFiles{1});
tstamp  = datestr(now, 'yyyymmdd_HHMMSS');
outName = fullfile(fPath, sprintf('R_uu_comparison_%s.png', tstamp));
saveas(gcf, outName);
fprintf('Saved → %s\n', outName);

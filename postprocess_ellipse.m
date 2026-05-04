% =========================================================================
% postprocess_ellipse.m
%
% Lightweight post-processing script. Loads bbox_results_*.mat files
% produced by analyse_Ruu_maps.m and computes ellipse parameters:
%
%   a     = (Lxu + Lxd) / 2          semi-major axis (mm)
%   b     = (Lyt + Lyb) / 2          semi-minor axis (mm)
%   cx    = (Lxd - Lxu) / 2          ellipse centre x offset from xref (mm)
%   cy    = (Lyt - Lyb) / 2          ellipse centre y offset from yref (mm)
%   theta = atand(cy / cx)            inclination angle (deg)
%                                     = angle of major axis from horizontal,
%                                       measured at the ellipse centre (x0,y0)
%                                       per Ganapathisubramani et al. (2005)
%
% theta is computed per rho level independently — no LS fit or tip
% coordinates required. Report theta at rho = 1/e as the canonical value.
%
% No mask operations — pure algebra on already-saved bbox values.
% Flexible to any number of rho levels.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc;

%% ===================== USER INPUTS ======================================

% Directory containing bbox_results_*.mat files
resultsDir = 'G:\two_point_covariance_20260426_120726\bbox_analysis';

% Output directory (can be same as resultsDir)
outDir = resultsDir;

%% ===================== SETUP ============================================

if ~exist(outDir, 'dir'), mkdir(outDir); end

matFiles = dir(fullfile(resultsDir, 'bbox_results_*.mat'));
if isempty(matFiles)
    error('No bbox_results_*.mat files found in:\n  %s', resultsDir);
end
nFiles = numel(matFiles);
fprintf('Found %d bbox_results file(s)\n\n', nFiles);

%% ===================== PROCESS EACH FILE ================================

for iF = 1:nFiles

    fname = fullfile(resultsDir, matFiles(iF).name);
    fprintf('[%d/%d] %s\n', iF, nFiles, matFiles(iF).name);

    D          = load(fname);
    rho_levels = D.rho_levels;
    nLevels    = numel(rho_levels);
    nPoints    = numel(D.xr_vec);

    fprintf('  rho levels : %s\n', mat2str(round(rho_levels, 3)));
    fprintf('  N positions: %d\n', nPoints);

    %% --- Semi-axes ---
    a_mat = D.a_ellipse_mat;   % semi-major from tip method
    b_mat = D.b_ellipse_mat;   % semi-minor from tip method
    %% --- Ellipse centre offset from reference point ---
    % cx > 0 means centre is downstream of xref (typical for Ruu)
    % cy > 0 means centre is above yref
    cx_mat = (D.Lxd_mat - D.Lxu_mat) / 2;   % [nPoints x nLevels]
    cy_mat = (D.Lyt_mat - D.Lyb_mat) / 2;
 
    %% --- Inclination angle per position per rho level ---
    % Load theta directly from saved .mat — computed via tip-based method
    % in analyse_Ruu_maps.m. Do NOT recompute from bbox quantities.
    theta_mat = D.theta_mat;   % [nPoints x nLevels]
    fprintf('  theta range (all levels): %.1f to %.1f deg\n', ...
        min(theta_mat(:), [], 'omitnan'), max(theta_mat(:), [], 'omitnan'));

    %% --- Append to .mat ---
    % save(fname, 'a_mat', 'b_mat', 'cx_mat', 'cy_mat', 'theta_mat', '-append');
    % fprintf('  Appended ellipse results -> %s\n', matFiles(iF).name);

    %% --- Independent variable for plotting ---
    if strcmp(D.corrMode, 'point')
        xplot      = D.yr_vec;
        xlabelStr  = 'y_{ref}  (mm)';
    else
        xplot      = D.xr_vec;
        xlabelStr  = 'x_{ref}  (mm)';
    end

    %% --- Save .csv per rho level ---
    for iL = 1:nLevels
        rhoTag  = sprintf('rho%02.0f', rho_levels(iL)*100);
        csvFile = fullfile(outDir, ...
            sprintf('ellipse_%s_%s.csv', D.folderLabel, rhoTag));

        T = table(...
            D.xr_vec,           D.yr_vec, ...
            D.Lxu_mat(:,iL),    D.Lxd_mat(:,iL), ...
            D.Lyt_mat(:,iL),    D.Lyb_mat(:,iL), ...
            D.Lx_mat(:,iL),     D.Ly_mat(:,iL), ...
            D.aspect_mat(:,iL), ...
            a_mat(:,iL),         b_mat(:,iL), ...
            cx_mat(:,iL),        cy_mat(:,iL), ...
            theta_mat(:,iL), ...
            'VariableNames', {'xr_mm','yr_mm', ...
                              'Lxu_mm','Lxd_mm','Lyt_mm','Lyb_mm', ...
                              'Lx_mm','Ly_mm','aspect', ...
                              'a_mm','b_mm','cx_mm','cy_mm','theta_deg'});

        writetable(T, csvFile);
        fprintf('  Saved .csv -> ellipse_%s_%s.csv\n', D.folderLabel, rhoTag);
    end

    %% --- Summary figure ---
    cmap = redblue(nLevels);

    figure('Color','w','Position',[50 50 1100 350], ...
        'Name', sprintf('Ellipse params — %s', D.folderLabel));

    % Panel 1: semi-major axis a
    ax1 = subplot(1,3,1); hold on; box on;
    for iL = 1:nLevels
        plot(ax1, xplot, a_mat(:,iL), 'o-', ...
            'Color', cmap(iL,:), 'MarkerFaceColor', cmap(iL,:), ...
            'MarkerSize', 5, 'LineWidth', 1.4);
    end
    xlabel(ax1, xlabelStr); ylabel(ax1, 'a  (mm)');
    title(ax1, 'Semi-major axis');
    legend(ax1, arrayfun(@(r) sprintf('\\rho=%.2f',r), rho_levels, ...
        'UniformOutput',false), 'Location','best','FontSize',7);
    grid(ax1,'on');

    % Panel 2: semi-minor axis b
    ax2 = subplot(1,3,2); hold on; box on;
    for iL = 1:nLevels
        plot(ax2, xplot, b_mat(:,iL), 'o-', ...
            'Color', cmap(iL,:), 'MarkerFaceColor', cmap(iL,:), ...
            'MarkerSize', 5, 'LineWidth', 1.4);
    end
    xlabel(ax2, xlabelStr); ylabel(ax2, 'b  (mm)');
    title(ax2, 'Semi-minor axis');
    grid(ax2,'on');

    % Panel 3: theta per rho level
    ax3 = subplot(1,3,3); hold on; box on;
    for iL = 1:nLevels
        plot(ax3, xplot, theta_mat(:,iL), 'o-', ...
            'Color', cmap(iL,:), 'MarkerFaceColor', cmap(iL,:), ...
            'MarkerSize', 5, 'LineWidth', 1.4);
    end
    xlabel(ax3, xlabelStr); ylabel(ax3, '\theta  (deg)');
    title(ax3, 'Inclination angle \theta  (per \rho level)');
    legend(ax3, arrayfun(@(r) sprintf('\\rho=%.2f',r), rho_levels, ...
        'UniformOutput',false), 'Location','best','FontSize',7);
    grid(ax3,'on');

    sgtitle(sprintf('Ellipse parameters — %s', strrep(D.folderLabel,'_',' ')), ...
        'FontSize',12);

    saveas(gcf, fullfile(outDir, ...
        sprintf('ellipse_summary_%s.png', D.folderLabel)));
    fprintf('  Saved figure -> ellipse_summary_%s.png\n\n', D.folderLabel);

end

fprintf('========================================\n');
fprintf('Done. Results in: %s\n', outDir);
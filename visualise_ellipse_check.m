% =========================================================================
% visualise_ellipse_check.m
%
% Visual confirmation of ellipse geometry.
%
% Ellipse parameters (a, b, theta) are derived from the saved bbox_results
% .mat files. Contour boundaries are recomputed on the fly from the raw
% R_uu_*.mat files so they are always accurate.
%
% USER INPUTS
%   bboxDir    — folder containing bbox_results_*.mat files
%   fileIdx    — which bbox file to load  (0 = list and exit)
%   rawFolders — cell array of raw R_uu_xref*/ folder paths, one per xref
%   xref_val   — xref (mm) identifying which rawFolder to use for this
%                bbox file; matched to the closest entry in rawFolders
%   gridFile   — grid.mat for coordinate conversion
%   posIdx     — position index within the bbox file (0 = all)
%   rho_show   — rho levels to show ([] = all)
%   dx_back, dx_fwd — search box, must match analyse_Ruu_maps.m
%
% Left panel : real contour outlines in redblue gradient per rho level
%              + greyscale ellipse overlays with a/b axes and theta arc,
%              centred on ellipse centre (x0, y0).
% Right panel: bar chart of a and b per rho level, theta annotated.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all ; 

%% ===================== USER INPUTS ======================================

% Folder containing bbox_results_*.mat files
bboxDir = 'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260424_231444\bbox_analysis';
% bboxDir = 'G:\two_point_covariance_20260426_120726\bbox_analysis';

% Which bbox file to load. Set to 0 to list available files and exit.
fileIdx = 6;

% Raw R_uu folder(s) — one entry per xref position processed.
% Add as many as needed; xref_val selects the right one.
rawFolders = {
    'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260424_231444\R_uu_xref711\'
    
};

% xref value (mm) for the bbox file being inspected.
% Must match the xref encoded in one of the rawFolders entries.
xref_val = 711;

% Grid file (same one used in analyse_Ruu_maps.m)
gridFile = 'G:\two_point_covariance_20260426_120726\grid.mat';

% Position index within the bbox file. Set to 0 to loop over all.
posIdx = 0;

% Rho levels to show. [] = show all levels in the file.
%   rho_show = [];                  % show all
%   rho_show = [1/exp(1), 0.5];    % two levels only
rho_show = [0.3 0.6];

% Search box — must match what was used in analyse_Ruu_maps.m
dx_back = 700;
dx_fwd  = 600;
% Correlation type — must match the R_uu/R_uv filenames in rawDir
corrType = 'uu';   % 'uu' or 'uv'

%% ===================== DISCOVER + SELECT BBOX FILE ======================

bboxFiles = dir(fullfile(bboxDir, 'bbox_results_*.mat'));
if isempty(bboxFiles)
    error('No bbox_results_*.mat files found in:\n  %s', bboxDir);
end

fprintf('Available bbox files:\n');
for f = 1:numel(bboxFiles)
    fprintf('  [%d] %s\n', f, bboxFiles(f).name);
end
fprintf('\n');

if fileIdx == 0
    fprintf('Set fileIdx to select a file and rerun.\n');
    return;
end
if fileIdx > numel(bboxFiles)
    error('fileIdx=%d but only %d file(s) found.', fileIdx, numel(bboxFiles));
end

%% ===================== LOAD BBOX FILE ===================================

fname = fullfile(bboxDir, bboxFiles(fileIdx).name);
fprintf('Loading bbox: %s\n', bboxFiles(fileIdx).name);
D = load(fname);

rho_levels = D.rho_levels;
nLevels    = numel(rho_levels);
nPositions = numel(D.xr_vec);

fprintf('Positions  : %d\n', nPositions);
fprintf('Rho levels : %s\n\n', mat2str(round(rho_levels, 3)));

%% ===================== SELECT RAW FOLDER ================================
% Match xref_val to the correct rawFolders entry by parsing the xref
% number from each folder name.

rawXrefs = nan(numel(rawFolders), 1);
for f = 1:numel(rawFolders)
    tok = regexp(rawFolders{f}, 'R_uu_xref([\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        rawXrefs(f) = str2double(tok{1});
    end
end

[minDiff, fMatch] = min(abs(rawXrefs - xref_val));
if minDiff > 1
    error(['xref_val=%.1f did not match any rawFolders entry within 1 mm.\n' ...
           'Available xrefs: %s'], xref_val, mat2str(rawXrefs));
end
rawDir = rawFolders{fMatch};
fprintf('Using raw folder: %s\n\n', rawDir);

%% ===================== LOAD GRID ========================================

G      = load(gridFile, 'worldX_merged', 'worldY_merged');
worldX = double(G.worldX_merged);
worldY = double(G.worldY_merged);

%% ===================== RHO SHOW FILTER ==================================

if isempty(rho_show)
    showMask = true(1, nLevels);
else
    showMask = ismember(round(rho_levels,6), round(rho_show(:)',6));
    if ~any(showMask)
        warning('rho_show matched nothing — showing all levels.');
        showMask = true(1, nLevels);
    end
end
showIdx = find(showMask);
nShow   = numel(showIdx);

% Contour colours — redblue gradient, one per shown level
cmap_contour = redblue(nShow);

% Ellipse colours — greyscale, dark (outermost) to light (innermost)
grey_vals = linspace(0.15, 0.60, nShow);
grey_ell  = repmat(grey_vals(:), 1, 3);

%% ===================== INDEPENDENT VARIABLE =============================

if strcmp(D.corrMode, 'point')
    indepVar   = D.yr_vec;
    indepLabel = 'y_{ref}';
else
    indepVar   = D.xr_vec;
    indepLabel = 'x_{ref}';
end

plotList = posIdx;
if posIdx == 0, plotList = 1:nPositions; end

fprintf('Plotting %d position(s)...\n\n', numel(plotList));

%% ===================== PLOT LOOP ========================================

t = linspace(0, 2*pi, 361);

for kk = 1:numel(plotList)

    k  = plotList(kk);
    xr = D.xr_vec(k);
    yr = D.yr_vec(k);

    % get_contour handles file loading, grid computation and mask
    % internally — just pass xr/yr as the target reference point.
    fprintf('  pos %d/%d: xr=%.1f yr=%.1f mm\n', k, nPositions, xr, yr);

    % Ellipse parameters computed per level via regionprops inside the
    % plot loop. Preallocate here for the right panel bar chart.
    a_vec     = nan(1, nLevels);
    b_vec     = nan(1, nLevels);
    cx_vec    = nan(1, nLevels);
    cy_vec    = nan(1, nLevels);
    theta_vec = nan(1, nLevels);

    % ================================================================
    % FIGURE
    % ================================================================
    figure('Color','w','Position',[50 50 1100 500], ...
        'Name', sprintf('Ellipse check | file %d | pos %d/%d | %s=%.2f mm', ...
            fileIdx, k, nPositions, indepLabel, indepVar(k)));

    % ================================================================
    % LEFT PANEL: contour outlines + ellipses
    % ================================================================
    ax1 = subplot(1,2,1);
    hold(ax1,'on'); box(ax1,'on'); axis(ax1,'equal');

    for si = nShow:-1:1   % outermost (lowest rho) drawn first
        iL    = showIdx(si);
        col_c = cmap_contour(si,:);
        col_e = grey_ell(si,:);

        % --- Get contour boundary and ellipse via shared helper ----------
        [bnd_dx, bnd_dy, ell] = get_contour_2(rawDir, corrType, xr, yr, ...
            worldX, worldY, dx_back, dx_fwd, rho_levels(iL));

        if isempty(bnd_dx) || ~ell.valid, continue; end

        a     = ell.a;
        b     = ell.b;
        theta = ell.theta;
        cx    = ell.cx;
        cy    = ell.cy;

        % Store for right panel
        a_vec(iL)     = a;
        b_vec(iL)     = b;
        cx_vec(iL)    = cx;
        cy_vec(iL)    = cy;
        theta_vec(iL) = theta;

        % --- Contour outline in colour (no fill) -------------------------
        plot(ax1, bnd_dx, bnd_dy, '-', ...
            'Color', col_c, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('\\rho = %.2f', rho_levels(iL)));

        % --- Greyscale ellipse overlay -----------------------------------
        %   x = cx + a*cos(t)*cos(theta) - b*sin(t)*sin(theta)
        %   y = cy + a*cos(t)*sin(theta) + b*sin(t)*cos(theta)
        xe = cx + a.*cos(t).*cosd(theta) - b.*sin(t).*sind(theta);
        ye = cy + a.*cos(t).*sind(theta) + b.*sin(t).*cosd(theta);

        plot(ax1, xe, ye, '-', 'Color', col_e, 'LineWidth', 2.2, ...
            'DisplayName', sprintf('\\rho = %.2f  ellipse', rho_levels(iL)));

        % a-axis through ellipse centre
        plot(ax1, [cx - a*cosd(theta),  cx + a*cosd(theta)], ...
                  [cy - a*sind(theta),  cy + a*sind(theta)], ...
            '-', 'Color', col_e, 'LineWidth', 1.6, 'HandleVisibility','off');

        % b-axis through ellipse centre (perpendicular)
        plot(ax1, [cx - b*cosd(theta+90),  cx + b*cosd(theta+90)], ...
                  [cy - b*sind(theta+90),  cy + b*sind(theta+90)], ...
            '-', 'Color', col_e, 'LineWidth', 1.6, 'HandleVisibility','off');

        % Axis labels on innermost shown level only
        if si == nShow
            text(ax1, cx + a*cosd(theta) + 3, cy + a*sind(theta), ...
                sprintf('a = %.0f mm', a), 'Color', col_e, ...
                'FontSize', 8, 'FontWeight','bold', ...
                'HorizontalAlignment','left', 'HandleVisibility','off');
            text(ax1, cx + b*cosd(theta+90), cy + b*sind(theta+90) + 3, ...
                sprintf('b = %.0f mm', b), 'Color', col_e, ...
                'FontSize', 8, 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'HandleVisibility','off');
        end

        % Theta arc at ellipse centre
        r_arc = a * 0.28;
        t_arc = linspace(0, deg2rad(theta), 60);
        plot(ax1, cx + r_arc.*cos(t_arc), cy + r_arc.*sin(t_arc), ...
            '-', 'Color', col_e, 'LineWidth', 1.4, 'HandleVisibility','off');
        text(ax1, cx + r_arc*1.25*cosd(theta/2), ...
                  cy + r_arc*1.25*sind(theta/2), ...
            sprintf('%.1f°', theta), ...
            'FontSize', 8, 'Color', col_e, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'HandleVisibility','off');

        % Ellipse centre marker
        % plot(ax1, cx, cy, '+', 'Color', col_e, ...
        %     'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility','off');

        % Dotted line: reference point → ellipse centre
        plot(ax1, [0, cx], [0, cy], ':', ...
            'Color', col_e, 'LineWidth', 0.9, 'HandleVisibility','off');

    end

    % Reference point
    plot(ax1, 0, 0, 'w+', 'MarkerSize', 12, 'LineWidth', 2, ...
        'DisplayName', 'x_{ref}, y_{ref}');

    xlabel(ax1, '\Deltax  (mm)');
    ylabel(ax1, '\Deltay  (mm)');
    title(ax1, sprintf('Contours + ellipses  |  %s = %.2f mm', ...
        indepLabel, indepVar(k)));
    legend(ax1, 'Location','best', 'FontSize', 7);
    set(ax1, 'YDir','normal');
    grid(ax1,'on');

    % ================================================================
    % RIGHT PANEL: bar chart of a, b and theta
    % ================================================================
    ax2 = subplot(1,2,2);
    hold(ax2,'on'); box(ax2,'on');

    bw = 0.35;
    for si = 1:nShow
        iL  = showIdx(si);
        col = cmap_contour(si,:);
        if isnan(a_vec(iL)), continue; end

        bar(ax2, si - bw/2, a_vec(iL), bw, ...
            'FaceColor', col, 'FaceAlpha', 0.90, 'EdgeColor','none');
        bar(ax2, si + bw/2, b_vec(iL), bw, ...
            'FaceColor', col, 'FaceAlpha', 0.40, 'EdgeColor', col, 'LineWidth', 1.2);

        text(ax2, si - bw/2, a_vec(iL) + 0.5, sprintf('%.0f', a_vec(iL)), ...
            'HorizontalAlignment','center', 'FontSize', 7, 'Color', col);
        text(ax2, si + bw/2, b_vec(iL) + 0.5, sprintf('%.0f', b_vec(iL)), ...
            'HorizontalAlignment','center', 'FontSize', 7, 'Color', col);

        if ~isnan(theta_vec(iL))
            text(ax2, si, 2, sprintf('\\theta=%.1f°', theta_vec(iL)), ...
                'HorizontalAlignment','center', 'FontSize', 7, ...
                'Color', col, 'FontWeight','bold');
        end
    end

    h_a = bar(ax2, NaN, NaN, 'FaceColor',[0.4 0.4 0.4], ...
        'FaceAlpha',0.9,'EdgeColor','none');
    h_b = bar(ax2, NaN, NaN, 'FaceColor',[0.4 0.4 0.4], ...
        'FaceAlpha',0.4,'EdgeColor',[0.4 0.4 0.4],'LineWidth',1.2);
    legend(ax2, [h_a, h_b], {'a  (semi-major)','b  (semi-minor)'}, ...
        'Location','best','FontSize',8);

    set(ax2, 'XTick', 1:nShow, ...
        'XTickLabel', arrayfun(@(r) sprintf('\\rho=%.2f',r), ...
            rho_levels(showIdx), 'UniformOutput',false), ...
        'XTickLabelRotation', 30);
    ylabel(ax2, 'Semi-axis  (mm)');
    title(ax2, 'a, b and \theta per \rho level');
    grid(ax2,'on');

    sgtitle(sprintf('Ellipse check  |  x_{ref}=%.1f mm,  y_{ref}=%.1f mm  |  pos %d/%d', ...
        xr, yr, k, nPositions), 'FontSize', 11);

    drawnow;

end

fprintf('Done.\n');
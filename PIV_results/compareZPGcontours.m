% =========================================================================
% plot_ellipse_compare_2x2.m
%
% Produces a nLevels x nYtargets figure comparing R_uu isocontour ellipse
% fits between two cases (AK and VP) across matched y_ref/delta positions.
%
% Layout:
%   Rows    -> rho levels             (e.g. [0.3, 0.6])
%   Columns -> y_ref / delta targets  (e.g. [0.1, 0.5])
%
% This organisation means all subplots in a row share the same rho level,
% so y-axes link cleanly across columns (same correlation magnitude, same
% contour scale). x-limits are shared globally.
%
% Within each panel, nXref streamwise stations are overlaid using a
% light-to-dark colour gradient (light = first station, dark = last).
%
%   AK  : red family,  solid ellipse,  semi-transparent solid raw contour
%   VP  : blue family, dashed ellipse, semi-transparent dashed raw contour
%
% Requires: get_contour.m  -> [bnd_dx, bnd_dy, ell]
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ========================= USER INPUTS ==================================

corrType   = 'uu';

% --- Font size (applied to all text in the figure) -----------------------
fs = 14;

% --- Contour levels and display weights -----------------------------------
rho_levels  = [0.3, 0.6];          % two levels -> two rows
rho_lw      = [1.5, 1.5];          % linewidth per rho level [low, high]
rho_alpha   = 0.30;                 % transparency of raw contour underlay

% --- Wall-normal targets (y_ref / delta) -> one column per entry ---------
ydelta_targets = [0.05, 0.5];       % change freely; nCols = numel(ydelta_targets)

% --- Correlation search window (mm, in raw pixel-space coordinates) ------
dx_back = 300;
dx_fwd  = 600;

% =========================================================================
% CASE AK
% =========================================================================
ak.label       = 'Case 2';
ak.blSweepFile = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';
ak.gridFile    = 'G:\two_point_covariance_20260426_120726\grid.mat';

% Representative streamwise stations — change nXref freely (3 or 5 etc.)
ak.xref_noms  = [54, 550, 1000];   % mm
ak.covFolders = { ...
    'G:\two_point_covariance_20260426_120726\R_uu_xref54',   ...
    'G:\two_point_covariance_20260426_120726\R_uu_xref550',  ...
    'G:\two_point_covariance_20260426_120726\R_uu_xref1000', ...
};

% Station labels for the legend (one per xref)
ak.stationLabels = {'Inlet', 'Crossover', 'Recovery'};

% =========================================================================
% CASE VP
% =========================================================================
vp.label       = 'VP2025';
vp.blSweepFiles = { ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215729_Pos1.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215751_Pos2.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215811_Pos3.mat', ...
};
vp.xref_noms  = [5866, 6921, 8000];
vp.covFolders = { ...
    'G:\SW 500mm Mean Flow Fields\Pos_1\two_point_covariance_20260424_225929\R_uu_xref5866', ...
    'G:\SW 500mm Mean Flow Fields\Pos_2\two_point_covariance_20260424_230729\R_uu_xref6921', ...
    'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260424_231444\R_uu_xref8000', ...
};
vp.stationLabels = {'Inlet', 'Crossover', 'Recovery'};

%% ========================= COLOUR SCHEME ================================
% Light -> dark gradient across nXref stations.
% First station = lightest, last station = darkest.
% Both cases built from the same t_vals so the gradient pacing is identical.

nXref  = numel(ak.xref_noms);       % must equal numel(vp.xref_noms)
% t_vals = linspace(0, 1, nXref);     % 0 = light, 1 = dark
t_vals = [0.0, 0.5, 1.0];   % light, mid, dark — widen the gap at the dark end
 
% AK: rose-to-crimson
% ak_light = [0.95, 0.72, 0.72];
ak_light = [0.88, 0.45, 0.45]; 

% ak_dark  = [0.50, 0.00, 0.05];
ak_dark  = [0.50, 0.00, 0.05];
ak_dark = [0.65, 0.00, 0.00];

% VP: sky-to-navy
% vp_light = [0.65, 0.82, 0.96];
vp_light = [0.30, 0.55, 0.85]; 

vp_dark  = [0.00, 0.05, 0.50];
vp_dark = [0.00, 0.00, 0.65];

ak_colors = arrayfun(@(t) (1-t)*ak_light + t*ak_dark, t_vals, ...
    'UniformOutput', false);
vp_colors = arrayfun(@(t) (1-t)*vp_light + t*vp_dark, t_vals, ...
    'UniformOutput', false);

%% ========================= LOAD AK GRID & delta99 =======================

G_ak      = load(ak.gridFile, 'worldX_merged', 'worldY_merged');
worldX_ak = double(G_ak.worldX_merged);
worldY_ak = double(G_ak.worldY_merged);

B_ak     = load(ak.blSweepFile, 'blSweep');
bl_x     = B_ak.blSweep.x_mm;
bl_d99   = B_ak.blSweep.delta99_hybrid_mm;
valid    = isfinite(bl_x) & isfinite(bl_d99);
bl_x_ak  = bl_x(valid);
bl_d_ak  = bl_d99(valid);
if any(diff(bl_x_ak) == 0)
    [bl_x_ak, ~, ic] = unique(bl_x_ak, 'sorted');
    bl_d_ak = accumarray(ic, bl_d_ak, [], @mean);
end
delta_ak = interp1(bl_x_ak, bl_d_ak, ak.xref_noms, 'linear', NaN);

fprintf('\nAK delta99 at each xref:\n');
for i = 1:nXref
    fprintf('  x = %5.0f mm  ->  delta99 = %.2f mm\n', ak.xref_noms(i), delta_ak(i));
end

%% ========================= LOAD VP GRIDS & delta99 ======================

delta_vp  = nan(1, nXref);
worldX_vp = cell(1, nXref);
worldY_vp = cell(1, nXref);

fprintf('\nVP delta99 at each xref:\n');
for iX = 1:nXref
    B      = load(vp.blSweepFiles{iX}, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;
    valid  = isfinite(bl_x) & isfinite(bl_d99);
    bl_xc  = bl_x(valid);
    bl_dc  = bl_d99(valid);
    if any(diff(bl_xc) == 0)
        [bl_xc, ~, ic] = unique(bl_xc, 'sorted');
        bl_dc = accumarray(ic, bl_dc, [], @mean);
    end
    delta_vp(iX) = interp1(bl_xc, bl_dc, vp.xref_noms(iX), 'linear', NaN);
    fprintf('  x = %5.0f mm  ->  delta99 = %.2f mm\n', vp.xref_noms(iX), delta_vp(iX));

    % Locate grid.mat: first try inside the covFolder, then one level up
    g1 = fullfile(vp.covFolders{iX}, 'grid.mat');
    g2 = fullfile(fileparts(vp.covFolders{iX}), 'grid.mat');
    if     isfile(g1); G = load(g1, 'worldX_merged', 'worldY_merged');
    elseif isfile(g2); G = load(g2, 'worldX_merged', 'worldY_merged');
    else;  error('grid.mat not found for VP xref %d.\n  %s\n  %s', ...
                 vp.xref_noms(iX), g1, g2);
    end
    worldX_vp{iX} = double(G.worldX_merged);
    worldY_vp{iX} = double(G.worldY_merged);
end

%% ========================= CACHE CONTOURS & ELLIPSES ====================

nYtargets = numel(ydelta_targets);
nLevels   = numel(rho_levels);

% Cell arrays: (iX, iY, iL)
bnd_ak = cell(nXref, nYtargets, nLevels);
bnd_vp = cell(nXref, nYtargets, nLevels);
ell_ak = cell(nXref, nYtargets, nLevels);
ell_vp = cell(nXref, nYtargets, nLevels);

fprintf('\nFetching contours and ellipses...\n');
for iX = 1:nXref
    for iY = 1:nYtargets
        for iL = 1:nLevels
            rho = rho_levels(iL);

            % -- AK --
            if isfinite(delta_ak(iX))
                yr = ydelta_targets(iY) * delta_ak(iX);
                [bdx, bdy, ell] = get_contour(ak.covFolders{iX}, corrType, ...
                    ak.xref_noms(iX), yr, worldX_ak, worldY_ak, ...
                    dx_back, dx_fwd, rho);
                if ~isempty(bdx)
                    bnd_ak{iX,iY,iL} = [bdx(:)/delta_ak(iX), bdy(:)/delta_ak(iX)];
                end
                ell_ak{iX,iY,iL} = ell;
            end

            % -- VP --
            if isfinite(delta_vp(iX))
                yr = ydelta_targets(iY) * delta_vp(iX);
                [bdx, bdy, ell] = get_contour(vp.covFolders{iX}, corrType, ...
                    vp.xref_noms(iX), yr, worldX_vp{iX}, worldY_vp{iX}, ...
                    dx_back, dx_fwd, rho);
                if ~isempty(bdx)
                    bnd_vp{iX,iY,iL} = [bdx(:)/delta_vp(iX), bdy(:)/delta_vp(iX)];
                end
                ell_vp{iX,iY,iL} = ell;
            end
        end
    end
    fprintf('  xref %d / %d done\n', iX, nXref);
end
fprintf('Done.\n\n');

%% ========================= AXIS LIMITS ==================================
% x: symmetric about zero, shared globally.
% y: symmetric about zero, fitted per ROW (= per rho level) with 20% margin.
%    All y/delta columns in the same row share scale -> clean y-axis linking.

y_row_lim = nan(nLevels, 2);
x_abs_max = 0;

% for iL = 1:nLevels
%     y_abs_max_row = 0;
%     for iX = 1:nXref
%         for iY = 1:nYtargets
%             for src = 1:2
%                 if src == 1; b = bnd_ak{iX,iY,iL}; else; b = bnd_vp{iX,iY,iL}; end
%                 if isempty(b), continue; end
%                 x_abs_max     = max(x_abs_max,     max(abs(b(:,1))));
%                 y_abs_max_row = max(y_abs_max_row, max(abs(b(:,2))));
%             end
%         end
%     end
%     if y_abs_max_row > 0
%         margin = y_abs_max_row * 0.20;
%         y_row_lim(iL,:) = [-(y_abs_max_row + margin), y_abs_max_row + margin];
% 
%     end
% end
for iL = 1:nLevels
    y_min_row =  inf;
    y_max_row = -inf;
    for iX = 1:nXref
        for iY = 1:nYtargets
            for src = 1:2
                if src == 1; b = bnd_ak{iX,iY,iL}; else; b = bnd_vp{iX,iY,iL}; end
                if isempty(b), continue; end
                x_abs_max = max(x_abs_max, max(abs(b(:,1))));
                y_min_row = min(y_min_row, min(b(:,2)));
                y_max_row = max(y_max_row, max(b(:,2)));
            end
        end
    end
    if isfinite(y_min_row)
        margin = (y_max_row - y_min_row) * 0.20;
        y_row_lim(iL,:) = [y_min_row - margin, y_max_row + margin];
    end
end
x_margin  = x_abs_max * 0.08;
x_lim_all = [-(x_abs_max + x_margin), x_abs_max + x_margin];

%% ========================= FIGURE =======================================
% Layout: nLevels rows x nYtargets columns
% Rows -> rho level;  Columns -> y_ref/delta
% y-axes link cleanly within each row (same rho = same correlation scale).

t_ell = linspace(0, 2*pi, 361);

fig = figure('Color', 'w', ...
    'Position', [60 60 400*nYtargets 340*nLevels], ...
    'Name', sprintf('R_{%s}  AK vs VP  |  ellipses + raw contours', corrType));

% ax indexed as (iL, iY): row = rho level, col = y/delta target
ax = gobjects(nLevels, nYtargets);
for iL = 1:nLevels
    for iY = 1:nYtargets
        ax(iL,iY) = subplot(nLevels, nYtargets, (iL-1)*nYtargets + iY);
        hold on; box on; grid on;
        set(ax(iL,iY), 'GridAlpha', 0.15, 'GridLineStyle', ':');
    end
end

% ---- Draw data ----------------------------------------------------------
for iX = 1:nXref
    col_ak = ak_colors{iX};
    col_vp = vp_colors{iX};

    for iL = 1:nLevels
        lw = rho_lw(iL);
        for iY = 1:nYtargets
            a_ax = ax(iL, iY);

            % -- Raw contour underlay (semi-transparent) ------------------
            b = bnd_ak{iX,iY,iL};
            if ~isempty(b)
                plot(a_ax, b(:,1), b(:,2), '-', ...
                    'Color', [col_ak, rho_alpha], ...
                    'LineWidth', lw * 0.7, ...
                    'HandleVisibility', 'off');
            end
            b = bnd_vp{iX,iY,iL};
            if ~isempty(b)
                plot(a_ax, b(:,1), b(:,2), '--', ...
                    'Color', [col_vp, rho_alpha], ...
                    'LineWidth', lw * 0.7, ...
                    'HandleVisibility', 'off');
            end

            % -- AK ellipse (solid, full opacity) ------------------------
            e = ell_ak{iX,iY,iL};
            if ~isempty(e) && e.valid && isfinite(delta_ak(iX))
                d   = delta_ak(iX);
                a_s = e.a/d;  b_s = e.b/d;
                cx  = e.cx/d; cy  = e.cy/d; th = e.theta;
                xe  = cx + a_s.*cos(t_ell).*cosd(th) - b_s.*sin(t_ell).*sind(th);
                ye  = cy + a_s.*cos(t_ell).*sind(th) + b_s.*sin(t_ell).*cosd(th);
                plot(a_ax, xe, ye, '-', 'Color', col_ak, 'LineWidth', lw, ...
                    'HandleVisibility', 'off');
            end

            % -- VP ellipse (dashed, full opacity) -----------------------
            e = ell_vp{iX,iY,iL};
            if ~isempty(e) && e.valid && isfinite(delta_vp(iX))
                d   = delta_vp(iX);
                a_s = e.a/d;  b_s = e.b/d;
                cx  = e.cx/d; cy  = e.cy/d; th = e.theta;
                xe  = cx + a_s.*cos(t_ell).*cosd(th) - b_s.*sin(t_ell).*sind(th);
                ye  = cy + a_s.*cos(t_ell).*sind(th) + b_s.*sin(t_ell).*cosd(th);
                plot(a_ax, xe, ye, '--', 'Color', col_vp, 'LineWidth', lw, ...
                    'HandleVisibility', 'off');
            end
        end
    end
end

% ---- Reference lines and cosmetics --------------------------------------
panel_labels = 'abcdefghij';
label_idx    = 0;

for iL = 1:nLevels
    for iY = 1:nYtargets
        a_ax = ax(iL, iY);
        xline(a_ax, 0, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, ...
            'LineStyle', ':', 'HandleVisibility', 'off');
        yline(a_ax, 0, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, ...
            'LineStyle', ':', 'HandleVisibility', 'off');

        % Panel label (a), (b), ... top-left corner in axes-normalised units
        label_idx = label_idx + 1;
        text(a_ax, 0.03, 0.97, sprintf('(%s)', panel_labels(label_idx)), ...
            'Units',               'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment',   'top', ...
            'Interpreter',         'latex', ...
            'FontSize',            fs, ...
            'FontWeight',          'bold');

        % Column header: y/delta target (top row only)
        if iL == 1
            title(a_ax, ...
                sprintf('$y/\\delta = %.2f$', ydelta_targets(iY)), ...
                'Interpreter', 'latex', 'FontSize', fs, 'FontWeight', 'bold');
        end

        % x-axis label (bottom row only)
        if iL == nLevels
            xlabel(a_ax, '$\Delta x \,/\, \delta$', ...
                'Interpreter', 'latex', 'FontSize', fs);
        else
            set(a_ax, 'XTickLabel', []);
        end

        % y-axis label: Delta y / delta on left column only
        % rho tag is placed outside the right edge of each row (see below)
        if iY == 1
            ylabel(a_ax, '$\Delta y \,/\, \delta$', ...
                'Interpreter', 'latex', 'FontSize', fs);
        else
            set(a_ax, 'YTickLabel', []);
        end

        set(a_ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    end
end

% ---- Apply limits -------------------------------------------------------
for iL = 1:nLevels
    row_ax = ax(iL, :);
    set(row_ax, 'XLim', x_lim_all);
    if all(isfinite(y_row_lim(iL,:)))
        set(row_ax, 'YLim', y_row_lim(iL,:));
    end
    linkaxes(row_ax, 'y');
end

% ---- rho level tags: rotated, outside left edge of each row ------------
drawnow;
for iL = 1:nLevels
    a_ref = ax(iL, 1);
    pos   = get(a_ref, 'Position');        % [left bottom width height]
    x_fig = pos(3) + 0.6;                 % left of first column
    y_fig = pos(2) + pos(4) * 0.5;        % vertically centred on row

    % Invisible axes pinned to figure-normalised coords
   ax_tag = axes(fig, 'Position', [x_fig, pos(2), 0.01, pos(4)], ...
        'Visible', 'off');
    text(ax_tag, 0.5, 0.5, ...
        sprintf('$\\rho_{uu} = %.1f$', rho_levels(iL)), ...
        'Units',               'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'middle', ...
        'Interpreter',         'latex', ...
        'FontSize',            fs, ...
        'FontWeight',          'bold', ...
        'Rotation',            90);
end

% ---- Legend (top-left panel) -------------------------------------------
a_leg = ax(1,1);

h_ak = plot(a_leg, NaN, NaN, '-',  'Color', ak_dark, 'LineWidth', 1.8, ...
    'DisplayName', ak.label);
h_vp = plot(a_leg, NaN, NaN, '--', 'Color', vp_dark, 'LineWidth', 1.8, ...
    'DisplayName', vp.label);
% 
% h_st = gobjects(nXref, 1);
% for iX = 1:nXref
%     col_mid = 0.5 * (ak_colors{iX} + vp_colors{iX});
%     h_st(iX) = plot(a_leg, NaN, NaN, 's', ...
%         'MarkerFaceColor', col_mid, ...
%         'MarkerEdgeColor', col_mid * 0.6, ...
%         'MarkerSize', 8, ...
%         'DisplayName', ak.stationLabels{iX});
% end
legend(a_leg, [h_ak; h_vp], ...
    'Interpreter', 'latex', ...
    'Location',    'best', ...
    'FontSize',    fs - 2, ...
    'Box',         'on');


% legend(a_leg, [h_ak; h_vp; h_st], ...
%     'Interpreter', 'latex', ...
%     'Location',    'best', ...
%     'FontSize',    fs - 2, ...
%     'Box',         'on');

% ---- Tighten layout -----------------------------------------------------
set(fig, 'Units', 'normalized');
tightfig_safe();

%% ========================= LOCAL FUNCTIONS ===============================

function tightfig_safe()
% Minimal tight-layout without requiring external toolboxes.
try
    set(gcf, 'Units', 'normalized');
    sp = findobj(gcf, 'Type', 'axes');
    for k = 1:numel(sp)
        sp(k).Position(1) = max(sp(k).Position(1), 0.06);
        sp(k).Position(2) = max(sp(k).Position(2), 0.08);
    end
catch
end
end
%% ========================= EXPORT =======================================

% --- Target dimensions (cm) — set to match journal column width ----------
fig_width_cm  = 17.0;   % full-width double column (adjust to journal spec)
fig_height_cm = 14.0;   % tune to taste; ~7 cm per row is a good starting point

% --- Convert to inches for MATLAB's PaperSize ----------------------------
cm2in = 1 / 2.54;
fig_width_in  = fig_width_cm  * cm2in;
fig_height_in = fig_height_cm * cm2in;

% --- Force figure to exact paper size so no scaling occurs ---------------
set(fig, 'Units',     'inches', ...
         'Position',  [1, 1, fig_width_in, fig_height_in]);
set(fig, 'PaperUnits',     'inches', ...
         'PaperSize',      [fig_width_in, fig_height_in], ...
         'PaperPosition',  [0, 0, fig_width_in, fig_height_in]);

% --- Timestamped filename -------------------------------------------------
timestamp  = datestr(now, 'yyyymmdd_HHMMSS');
save_dir   = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\figures\';
filename   = fullfile(save_dir, sprintf('Ruu_ellipse_compare_%s.pdf', timestamp));

% --- Export ---------------------------------------------------------------
exportgraphics(fig, filename, ...
    'ContentType',  'vector', ...   % vector PDF, not rasterised
    'BackgroundColor', 'white');

fprintf('Figure saved to:\n  %s\n', filename);

% =========================================================================
% ComparetTomCase2_isocontourgeometries.m
%
% Compares R_uu (or R_uv) isocontours between two cases across matched
% flow regimes. Produces three figures with identical subplot layouts:
%
%   Figure 1 — raw isocontours
%   Figure 2 — fitted ellipses
%   Figure 3 — contours + ellipses overlaid
%
% Layout: nTargets rows x nXref columns
%   AK:  red family   (solid lines)
%   SW:  gray family  (dashed lines)
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS ======================================

corrType      = 'uu';
rho_levels    = [0.3 0.6];
rho_linewidth = [2.0, 2.0];
rho_linestyle = {'-.', '-'};
dx_back       = 300;
dx_fwd        = 600;
fontSize      = 12;

ydelta_targets = [0.05, 0.5];

%% ===================== COLOUR SCHEME ====================================

nXref  = 3;

ak_dark  = [0.6 0 0.00];
ak_light = [1.0 0.4 0.4];
t = 1/2;
ak_mid = (1-t)*ak_dark + t*ak_light;

ak_colors = cell(1, nXref);
for k = 1:nXref
    ak_colors{k} = ak_mid;
end

sw_dark  = [0.3 0.3 0.3];
sw_colors = cell(1, nXref);
for k = 1:nXref
    sw_colors{k} = sw_dark;
end

%% ===================== CASE DEFINITIONS =================================

% --- Case AK ---
ak.blSweepFile = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\128x128init\blSweep_20260415_125916.mat';
ak.gridFile    = 'G:\two_point_covariance_20260426_120726\grid.mat';
ak.xref_noms   = [54, 550, 1000];
ak.covFolders  = { ...
    'G:\two_point_covariance_20260426_120726\R_uu_xref54',   ...
    'G:\two_point_covariance_20260426_120726\R_uu_xref550',  ...
    'G:\two_point_covariance_20260426_120726\R_uu_xref1000', ...
};
ak.regime_idx = [1 2 3];

% --- Case SW ---
sw.blSweepFiles = { ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215729_Pos1.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215751_Pos2.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215811_Pos3.mat', ...
};
sw.xref_noms  = [5853, 6850, 7755];
sw.covFolders = { ...
    'G:\SW 500mm Mean Flow Fields\Pos_1\two_point_covariance_20260428_214739\R_uu_xref5853', ...
    'G:\SW 500mm Mean Flow Fields\Pos_2\two_point_covariance_20260428_215235\R_uu_xref6850', ...
    'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260428_220246\R_uu_xref7755', ...
};
sw.regime_idx = [1, 2, 3];

regimeLabels = {'Inlet', 'Crossover', 'Recovery'};

%% ===================== LOAD AK GRID AND delta99 =========================

G_ak      = load(ak.gridFile, 'worldX_merged', 'worldY_merged');
worldX_ak = double(G_ak.worldX_merged);
worldY_ak = double(G_ak.worldY_merged);

B_ak      = load(ak.blSweepFile, 'blSweep');
bl_x      = B_ak.blSweep.x_mm;
bl_d99    = B_ak.blSweep.delta99_hybrid_mm;
valid_bl  = isfinite(bl_x) & isfinite(bl_d99);
bl_x_ak   = bl_x(valid_bl);
bl_d99_ak = bl_d99(valid_bl);
if any(diff(bl_x_ak) == 0)
    [bl_x_ak, ~, ic] = unique(bl_x_ak, 'sorted');
    bl_d99_ak = accumarray(ic, bl_d99_ak, [], @mean);
    fprintf('AK: duplicates removed -> %d unique x points\n', numel(bl_x_ak));
end
delta_ak = interp1(bl_x_ak, bl_d99_ak, ak.xref_noms, 'linear', NaN);

fprintf('\nAK delta99 at each xref:\n');
for i = 1:nXref
    fprintf('  x=%4.0f mm -> delta99 = %.2f mm\n', ak.xref_noms(i), delta_ak(i));
end

%% ===================== LOAD SW GRIDS AND delta99 ========================

delta_sw  = nan(1, nXref);
worldX_sw = cell(1, nXref);
worldY_sw = cell(1, nXref);

fprintf('\nSW delta99 at each xref:\n');
for iX = 1:nXref
    B      = load(sw.blSweepFiles{iX}, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;
    valid  = isfinite(bl_x) & isfinite(bl_d99);
    bl_xc  = bl_x(valid);
    bl_dc  = bl_d99(valid);
    if any(diff(bl_xc) == 0)
        [bl_xc, ~, ic] = unique(bl_xc, 'sorted');
        bl_dc = accumarray(ic, bl_dc, [], @mean);
    end
    delta_sw(iX) = interp1(bl_xc, bl_dc, sw.xref_noms(iX), 'linear', NaN);
    fprintf('  x=%5.0f mm -> delta99 = %.2f mm\n', sw.xref_noms(iX), delta_sw(iX));

    g1 = fullfile(sw.covFolders{iX}, 'grid.mat');
    g2 = fullfile(fileparts(sw.covFolders{iX}), 'grid.mat');
    if isfile(g1)
        G = load(g1, 'worldX_merged', 'worldY_merged');
    elseif isfile(g2)
        G = load(g2, 'worldX_merged', 'worldY_merged');
    else
        error('grid.mat not found for SW xref %d.\n  %s\n  %s', ...
            sw.xref_noms(iX), g1, g2);
    end
    worldX_sw{iX} = double(G.worldX_merged);
    worldY_sw{iX} = double(G.worldY_merged);
end

%% ===================== CACHE CONTOURS & ELLIPSES ========================

nTargets = numel(ydelta_targets);
nLevels  = numel(rho_levels);

bnd_ak = cell(nXref, nTargets, nLevels);
bnd_sw = cell(nXref, nTargets, nLevels);
ell_ak = cell(nXref, nTargets, nLevels);
ell_sw = cell(nXref, nTargets, nLevels);

fprintf('\nFetching contours and fitting ellipses...\n');
for iX = 1:nXref
    for iT = 1:nTargets
        yd_target = ydelta_targets(iT);
        for iL = 1:nLevels

            delta_ak_x = delta_ak(iX);
            if isfinite(delta_ak_x)
                yr_ak = yd_target * delta_ak_x;
                [bdx, bdy, ell] = get_contour_2(ak.covFolders{iX}, corrType, ...
                    ak.xref_noms(iX), yr_ak, ...
                    worldX_ak, worldY_ak, dx_back, dx_fwd, rho_levels(iL));
                if ~isempty(bdx)
                    bnd_ak{iX,iT,iL} = [bdx(:)/delta_ak_x, bdy(:)/delta_ak_x];
                end
                ell_ak{iX,iT,iL} = ell;
            end

            delta_sw_x = delta_sw(iX);
            if isfinite(delta_sw_x)
                yr_sw = yd_target * delta_sw_x;
                [bdx, bdy, ell] = get_contour_2(sw.covFolders{iX}, corrType, ...
                    sw.xref_noms(iX), yr_sw, ...
                    worldX_sw{iX}, worldY_sw{iX}, dx_back, dx_fwd, rho_levels(iL));
                if ~isempty(bdx)
                    bnd_sw{iX,iT,iL} = [bdx(:)/delta_sw_x, bdy(:)/delta_sw_x];
                end
                ell_sw{iX,iT,iL} = ell;
            end

        end
    end
    fprintf('  xref %d/%d done\n', iX, nXref);
end
fprintf('Done.\n\n');

%% ===================== AXIS LIMITS ======================================

y_row_lim = nan(nTargets, 2);
x_all_min =  Inf;
x_all_max = -Inf;

for iT = 1:nTargets
    y_min_row =  Inf;
    y_max_row = -Inf;
    for iX = 1:nXref
        for iL = 1:nLevels
            for src = 1:2
                if src == 1; b = bnd_ak{iX,iT,iL}; else; b = bnd_sw{iX,iT,iL}; end
                if isempty(b), continue; end
                y_min_row = min(y_min_row, min(b(:,2)));
                y_max_row = max(y_max_row, max(b(:,2)));
                x_all_min = min(x_all_min, min(b(:,1)));
                x_all_max = max(x_all_max, max(b(:,1)));
            end
        end
    end
    if isfinite(y_min_row)
        margin = max((y_max_row - y_min_row) * 0.20, 0.05);
        y_row_lim(iT,:) = [y_min_row - margin, y_max_row + margin];
    end
end

if isfinite(x_all_min)
    x_margin  = (x_all_max - x_all_min) * 0.05;
    x_lim_all = [x_all_min - x_margin, x_all_max + x_margin];
else
    x_lim_all = [-5 5];
end

t_ell = linspace(0, 2*pi, 361);

%% ===================== FIGURE 1: RAW CONTOURS ===========================

fig1 = figure('Color', 'w', ...
    'Position', [50 30 260*nXref 220*nTargets], ...
    'Name', sprintf('R_{%s} isocontours: AK vs SW', corrType));

axC = gobjects(nTargets, nXref);
for iT = 1:nTargets
    for iX = 1:nXref
        axC(iT,iX) = subplot(nTargets, nXref, (iT-1)*nXref + iX);
        hold on; box on; grid on;
        set(axC(iT,iX), 'FontSize', fontSize);
    end
end

for iX = 1:nXref
    col_ak = ak_colors{ak.regime_idx(iX)};
    col_sw = sw_colors{sw.regime_idx(iX)};
    for iT = 1:nTargets
        ax = axC(iT,iX);
        for iL = 1:nLevels
            if ~isempty(bnd_ak{iX,iT,iL})
                b = bnd_ak{iX,iT,iL};
                plot(ax, b(:,1), b(:,2), rho_linestyle{iL}, 'Color', col_ak, ...
                    'LineWidth', rho_linewidth(iL), 'HandleVisibility', 'off');
            end
            if ~isempty(bnd_sw{iX,iT,iL})
                b = bnd_sw{iX,iT,iL};
                plot(ax, b(:,1), b(:,2), rho_linestyle{iL}, 'Color', col_sw, ...
                    'LineWidth', rho_linewidth(iL), 'HandleVisibility', 'off');
            end
        end
        xline(ax, 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        yline(ax, 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        if iT == 1
            title(ax, regimeLabels{iX}, 'Interpreter', 'latex', ...
                'FontWeight', 'bold', 'FontSize', fontSize);
        end
        if iT == nTargets
            xlabel(ax, '$\Delta x / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
        else
            set(ax, 'XTickLabel', []);
        end
        if iX == 1
            ylabel(ax, '$\Delta y / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
        else
            set(ax, 'YTickLabel', []);
        end
    end
end

apply_limits(axC, nTargets, nXref, y_row_lim, x_lim_all);
drawnow;
add_panel_labels_and_row_annotations(fig1, axC, nTargets, nXref, ydelta_targets, fontSize);
add_figure_legend(fig1, axC(nTargets, ceil(nXref/2)), ...
    ak_colors, sw_colors, rho_levels, rho_linewidth, rho_linestyle, fontSize);

%% ===================== FIGURE 2: FITTED ELLIPSES ========================

fig2 = figure('Color', 'w', ...
    'Position', [100 80 260*nXref 220*nTargets], ...
    'Name', sprintf('R_{%s} ellipses: AK vs SW', corrType));

axE = gobjects(nTargets, nXref);
for iT = 1:nTargets
    for iX = 1:nXref
        axE(iT,iX) = subplot(nTargets, nXref, (iT-1)*nXref + iX);
        hold on; box on; grid on;
        set(axE(iT,iX), 'FontSize', fontSize);
    end
end

for iX = 1:nXref
    col_ak     = ak_colors{ak.regime_idx(iX)};
    col_sw     = sw_colors{sw.regime_idx(iX)};
    delta_ak_x = delta_ak(iX);
    delta_sw_x = delta_sw(iX);
    for iT = 1:nTargets
        ax = axE(iT,iX);
        for iL = 1:nLevels
            e = ell_ak{iX,iT,iL};
            if ~isempty(e) && e.valid && isfinite(delta_ak_x)
                a=e.a/delta_ak_x; b=e.b/delta_ak_x;
                cx=e.cx/delta_ak_x; cy=e.cy/delta_ak_x; th=e.theta;
                xe=cx+a.*cos(t_ell).*cosd(th)-b.*sin(t_ell).*sind(th);
                ye=cy+a.*cos(t_ell).*sind(th)+b.*sin(t_ell).*cosd(th);
                plot(ax,xe,ye,rho_linestyle{iL},'Color',col_ak, ...
                    'LineWidth',rho_linewidth(iL),'HandleVisibility','off');
            end
            e = ell_sw{iX,iT,iL};
            if ~isempty(e) && e.valid && isfinite(delta_sw_x)
                a=e.a/delta_sw_x; b=e.b/delta_sw_x;
                cx=e.cx/delta_sw_x; cy=e.cy/delta_sw_x; th=e.theta;
                xe=cx+a.*cos(t_ell).*cosd(th)-b.*sin(t_ell).*sind(th);
                ye=cy+a.*cos(t_ell).*sind(th)+b.*sin(t_ell).*cosd(th);
                plot(ax,xe,ye,rho_linestyle{iL},'Color',col_sw, ...
                    'LineWidth',rho_linewidth(iL),'HandleVisibility','off');
            end
        end
        xline(ax, 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        yline(ax, 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        if iT == 1
            title(ax, regimeLabels{iX}, 'Interpreter', 'latex', ...
                'FontWeight', 'bold', 'FontSize', fontSize);
        end
        if iT == nTargets
            xlabel(ax, '$\Delta x / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
        else
            set(ax, 'XTickLabel', []);
        end
        if iX == 1
            ylabel(ax, '$\Delta y / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
        else
            set(ax, 'YTickLabel', []);
        end
    end
end

apply_limits(axE, nTargets, nXref, y_row_lim, x_lim_all);
drawnow;
add_panel_labels_and_row_annotations(fig2, axE, nTargets, nXref, ydelta_targets, fontSize);
add_figure_legend(fig2, axE(nTargets, ceil(nXref/2)), ...
    ak_colors, sw_colors, rho_levels, rho_linewidth, rho_linestyle, fontSize);

%% ===================== FIGURE 3: CONTOURS + ELLIPSES OVERLAID ===========

fig3 = figure('Color', 'w', ...
    'Position', [150 130 260*nXref 220*nTargets], ...
    'Name', sprintf('R_{%s} contours + ellipses: AK vs SW', corrType));

axO = gobjects(nTargets, nXref);
for iT = 1:nTargets
    for iX = 1:nXref
        axO(iT,iX) = subplot(nTargets, nXref, (iT-1)*nXref + iX);
        hold on; box on; grid on;
        set(axO(iT,iX), 'FontSize', fontSize);
    end
end

for iX = 1:nXref
    col_ak     = ak_colors{ak.regime_idx(iX)};
    col_sw     = sw_colors{sw.regime_idx(iX)};
    col_ak_lt  = col_ak * 0.25 + 0.75;
    col_sw_lt  = col_sw * 0.25 + 0.75;
    delta_ak_x = delta_ak(iX);
    delta_sw_x = delta_sw(iX);
    for iT = 1:nTargets
        ax = axO(iT,iX);
        for iL = 1:nLevels
            if ~isempty(bnd_ak{iX,iT,iL})
                b = bnd_ak{iX,iT,iL};
                plot(ax, b(:,1), b(:,2), rho_linestyle{iL}, ...
                    'Color', col_ak_lt, 'LineWidth', rho_linewidth(iL), ...
                    'HandleVisibility', 'off');
            end
            if ~isempty(bnd_sw{iX,iT,iL})
                b = bnd_sw{iX,iT,iL};
                plot(ax, b(:,1), b(:,2), rho_linestyle{iL}, ...
                    'Color', col_sw_lt, 'LineWidth', rho_linewidth(iL), ...
                    'HandleVisibility', 'off');
            end
            e = ell_ak{iX,iT,iL};
            if ~isempty(e) && e.valid && isfinite(delta_ak_x)
                a=e.a/delta_ak_x; b=e.b/delta_ak_x;
                cx=e.cx/delta_ak_x; cy=e.cy/delta_ak_x; th=e.theta;
                xe=cx+a.*cos(t_ell).*cosd(th)-b.*sin(t_ell).*sind(th);
                ye=cy+a.*cos(t_ell).*sind(th)+b.*sin(t_ell).*cosd(th);
                plot(ax,xe,ye,rho_linestyle{iL},'Color',col_ak, ...
                    'LineWidth',rho_linewidth(iL),'HandleVisibility','off');
            end
            e = ell_sw{iX,iT,iL};
            if ~isempty(e) && e.valid && isfinite(delta_sw_x)
                a=e.a/delta_sw_x; b=e.b/delta_sw_x;
                cx=e.cx/delta_sw_x; cy=e.cy/delta_sw_x; th=e.theta;
                xe=cx+a.*cos(t_ell).*cosd(th)-b.*sin(t_ell).*sind(th);
                ye=cy+a.*cos(t_ell).*sind(th)+b.*sin(t_ell).*cosd(th);
                plot(ax,xe,ye,rho_linestyle{iL},'Color',col_sw, ...
                    'LineWidth',rho_linewidth(iL),'HandleVisibility','off');
            end
        end
        xline(ax, 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        yline(ax, 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
        if iT == 1
            title(ax, regimeLabels{iX}, 'Interpreter', 'latex', ...
                'FontWeight', 'bold', 'FontSize', fontSize);
        end
        if iT == nTargets
            xlabel(ax, '$\Delta x / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
        else
            set(ax, 'XTickLabel', []);
        end
        if iX == 1
            ylabel(ax, '$\Delta y / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
        else
            set(ax, 'YTickLabel', []);
        end
    end
end

apply_limits(axO, nTargets, nXref, y_row_lim, x_lim_all);
drawnow;
add_panel_labels_and_row_annotations(fig3, axO, nTargets, nXref, ydelta_targets, fontSize);
add_figure_legend(fig3, axO(nTargets, ceil(nXref/2)), ...
    ak_colors, sw_colors, rho_levels, rho_linewidth, rho_linestyle, fontSize);

%% ===================== LOCAL FUNCTIONS ==================================

function apply_limits(axArr, nTargets, nXref, y_row_lim, x_lim_all)
for iT = 1:nTargets
    row_axes = axArr(iT, :);
    set(row_axes, 'XLim', x_lim_all);
    if all(isfinite(y_row_lim(iT,:)))
        set(row_axes, 'YLim', y_row_lim(iT,:));
    end
    linkaxes(row_axes, 'xy');
end
end

function add_figure_legend(fig, ax_anchor, ak_colors, sw_colors, ...
    rho_levels, rho_linewidth, rho_linestyle, fontSize)

h_ak  = plot(ax_anchor, NaN, NaN, '-',              ...
    'Color', ak_colors{1}, 'LineWidth', 1.8, ...
    'DisplayName', 'Case 2');
h_sw  = plot(ax_anchor, NaN, NaN, '-',              ...
    'Color', sw_colors{1}, 'LineWidth', 1.8, ...
    'DisplayName', 'VP2025');
h_lw1 = plot(ax_anchor, NaN, NaN, rho_linestyle{1}, ...
    'Color', 'k', 'LineWidth', rho_linewidth(1), ...
    'DisplayName', ['$\rho = ' sprintf('%.2f', rho_levels(1)) '$']);
h_lw2 = plot(ax_anchor, NaN, NaN, rho_linestyle{2}, ...
    'Color', 'k', 'LineWidth', rho_linewidth(2), ...
    'DisplayName', ['$\rho = ' sprintf('%.2f', rho_levels(2)) '$']);

lgd = legend(ax_anchor, [h_ak, h_sw, h_lw1, h_lw2], ...
    'Interpreter',  'latex', ...
    'Orientation',  'horizontal', ...
    'FontSize',     fontSize);

drawnow;
lgd.Units    = 'normalized';
lgd_h        = lgd.Position(4);
legend_space = lgd_h + 0.06;

% Compress and shift all axes into space above legend
all_ax = findall(fig, 'Type', 'axes');
for k = 1:numel(all_ax)
    p          = all_ax(k).Position;
    new_bottom = p(2) * (1 - legend_space) + legend_space;
    new_height = p(4) * (1 - legend_space);
    all_ax(k).Position = [p(1), new_bottom, p(3), new_height];
end

drawnow;
lgd.Position = [0.5 - lgd.Position(3)/2, 0.01, ...
                lgd.Position(3),          lgd.Position(4)];
end


function add_panel_labels_and_row_annotations(fig, axArr, nTargets, nXref, ...
    ydelta_targets, fontSize)

panelLetters = 'abcdefghijklmnopqrstuvwxyz';

for iT = 1:nTargets
    for iX = 1:nXref
        ax  = axArr(iT, iX);
        idx = (iT-1)*nXref + iX;
        pos = ax.Position;

        annotation(fig, 'textbox', ...
            [pos(1), pos(2)+pos(4)-0.04, 0.04, 0.04], ...
            'String',             ['(' panelLetters(idx) ')'], ...
            'EdgeColor',          'none', ...
            'FontSize',           fontSize, ...
            'FontWeight',         'bold', ...
            'Interpreter',        'none', ...
            'VerticalAlignment',  'top', ...
            'HorizontalAlignment','left');
    end

    ax_last = axArr(iT, nXref);
    pos     = ax_last.Position;

    x_tag = pos(1) + pos(3) + 0.005;
    y_tag = pos(2);
    w_tag = 0.04;
    h_tag = pos(4);

    ax_tag = axes(fig, 'Position', [x_tag, y_tag, w_tag, h_tag], ...
        'Visible', 'off');

    if ydelta_targets(iT) == 0.05
        rowStr = '$y_{\mathrm{ref}}/\delta = 0.05$';
    else
        rowStr = ['$y_{\mathrm{ref}}/\delta = ' sprintf('%g', ydelta_targets(iT)) '$'];
    end

    text(ax_tag, 0.5, 0.5, rowStr, ...
        'Units',               'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'middle', ...
        'Interpreter',         'latex', ...
        'FontSize',            fontSize, ...
        'FontWeight',          'bold', ...
        'Rotation',            90);
end
end

%% ===================== EXPORT ===========================================

fig_width_in  = 8;
fig_height_in = 4.5;

save_dir  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\figures\';
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

for iFig = 1:3
    switch iFig
        case 1; fig_out = fig1; tag = 'rawcontours';
        case 2; fig_out = fig2; tag = 'ellipses';
        case 3; fig_out = fig3; tag = 'contours_ellipses';
    end

    drawnow;
    set(fig_out, 'Units',         'inches', ...
                 'Position',      [1, 1, fig_width_in, fig_height_in]);
    set(fig_out, 'PaperUnits',    'inches', ...
                 'PaperSize',     [fig_width_in, fig_height_in], ...
                 'PaperPosition', [0, 0, fig_width_in, fig_height_in]);
    % 
    % filename = fullfile(save_dir, ...
    %     sprintf('Ruu_isocontour_compareforTSFP_%s_%s.pdf', tag, timestamp));
    % 
    % exportgraphics(fig_out, filename, ...
    %     'ContentType',    'vector', ...
    %     'BackgroundColor','white');

    fprintf('Saved: %s\n', filename);
end
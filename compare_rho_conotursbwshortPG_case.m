% =========================================================================
% compare_isocontours_multiRho.m
%
% Compares R_uu isocontour geometry across cases from the SAME dataset at
% a single streamwise position. Two rho levels are overlaid per panel:
%   solid  line -> rho = 0.6
%   dashed line -> rho = 0.3
%
% Raw isocontours are plotted pale; fitted ellipses (via get_contour_2)
% are overlaid at full opacity — matching the style of
% CompareTomCase2_isocontourgeometries.m (Figure 3).
%
% Layout: 2 rows x 1 col
%   Row 1 -> y_ref / delta = 0.05
%   Row 2 -> y_ref / delta = 0.50
%
% Cases are distinguished by colour (red gradient). No legend.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS ======================================

% --- Regime / streamwise position ---
regime_request = 'zpg';   % 'fpg' | 'cross' | 'apg' | 'zpg'

corrType = 'uu';

% --- Two rho levels ---
rho_levels    = [0.3, 0.6];
rho_linestyle = {'--', '-'};     % dashed for 0.3, solid for 0.6
rho_linewidth = [2.0,  2.0];

% --- Wall-normal reference heights ---
ydelta_targets = [0.05, 0.5];

% --- Contour extraction window (mm either side of xref) ---
dx_back  = 300;
dx_fwd   = 600;

fontSize = 12;

% --- Case definitions ---
cases(1).label    = 'Case 1';
cases(1).cov_dir  = 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038';
cases(1).blFile   = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case1_PIVresults\blSweep_20260419_193106.mat';
cases(1).gridFile = 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\grid.mat';

cases(2).label    = 'Case 2';
cases(2).cov_dir  = 'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244';
cases(2).blFile   = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';
cases(2).gridFile = 'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244\grid.mat';

cases(3).label    = 'Case 3';
cases(3).cov_dir  = 'D:\FULLPROCESSEDY235AOAN11AOAFN11PIVDATA\two_point_covariance_20260421_083251';
cases(3).blFile   = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case6_PIVresults\blSweep_20260421_082649.mat';
cases(3).gridFile = 'D:\FULLPROCESSEDY235AOAN11AOAFN11PIVDATA\two_point_covariance_20260421_083251\grid.mat';

%% ===================== REGIME -> XREF MAP ================================

regimeMap = struct('fpg', 350, 'cross', 550, 'apg', 711, 'zpg', 1000);
if ~isfield(regimeMap, regime_request)
    error('regime_request must be one of: fpg, cross, apg, zpg');
end
xref_nom = regimeMap.(regime_request);

fprintf('Regime: %s  |  xref = %.0f mm\n\n', regime_request, xref_nom);

%% ===================== COLOUR SCHEME ====================================

nCases = numel(cases);

red_dark   = [0.6 0.0 0.0];
red_bright = [1.0 0.4 0.4];
caseColours = [linspace(red_dark(1), red_bright(1), nCases)', ...
               linspace(red_dark(2), red_bright(2), nCases)', ...
               linspace(red_dark(3), red_bright(3), nCases)'];

% Pale version for raw contours: 25% colour, 75% white
% (matches col_ak_lt = col_ak * 0.25 + 0.75 in CompareTomCase2)
paleFactor  = 0.25;
caseColPale = caseColours * paleFactor + (1 - paleFactor);

% Parametric angle vector for ellipse rendering
t_ell = linspace(0, 2*pi, 361);

%% ===================== LOAD GRID ========================================

G      = load(cases(1).gridFile, 'worldX_merged', 'worldY_merged');
worldX = double(G.worldX_merged);
worldY = double(G.worldY_merged);

%% ===================== LOAD delta99 FOR EACH CASE =======================

for iC = 1:nCases
    cases(iC).delta    = NaN;
    cases(iC).xref_dir = '';
end

for iC = 1:nCases
    fprintf('Loading blSweep: %s...\n', cases(iC).label);

    B      = load(cases(iC).blFile, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;
    valid  = isfinite(bl_x) & isfinite(bl_d99);
    bl_xv  = bl_x(valid);
    bl_dv  = bl_d99(valid);

    % Guard against duplicate x values (interp1 requirement)
    if any(diff(bl_xv) == 0)
        [bl_xv, ~, ic] = unique(bl_xv, 'sorted');
        bl_dv = accumarray(ic, bl_dv, [], @mean);
        fprintf('  Duplicates removed -> %d unique x points\n', numel(bl_xv));
    end

    delta = interp1(bl_xv, bl_dv, xref_nom, 'linear', NaN);

    if ~isfinite(delta)
        warning('  delta99 not available at xref=%.0f mm -- skipping\n', xref_nom);
        continue;
    end
    cases(iC).delta = delta;
    fprintf('  delta99 = %.2f mm\n', delta);

    % Locate xref subfolder
    subPattern = sprintf('R_uu_xref%g*', xref_nom);
    d = dir(fullfile(cases(iC).cov_dir, subPattern));
    d = d([d.isdir]);

    if isempty(d)
        warning('  No subfolder matching %s -- skipping', subPattern);
        continue;
    end
    cases(iC).xref_dir = fullfile(cases(iC).cov_dir, d(1).name);
    fprintf('  xref dir: %s\n\n', cases(iC).xref_dir);
end

%% ===================== CACHE CONTOURS & ELLIPSES ========================
% File-finding logic copied verbatim from run_piv_structure_comparison Fig1
% (the version that worked). Ellipse fit appended inline after bwboundaries.
%
% bnd_cache{iC,iT,iL} = [dx/delta, dy/delta]  normalised boundary coords
% ell_cache{iC,iT,iL} = struct: valid, a, b, cx, cy, theta(deg)  in mm

nTargets  = numel(ydelta_targets);
nLevels   = numel(rho_levels);
bnd_cache = cell(nCases, nTargets, nLevels);
ell_cache = cell(nCases, nTargets, nLevels);

x_min_domain = min(worldX(:));
x_max_domain = max(worldX(:));

t_ell = linspace(0, 2*pi, 361);

fprintf('Extracting contours and fitting ellipses...\n');

for iC = 1:nCases

    if isempty(cases(iC).xref_dir) || ~isfinite(cases(iC).delta)
        continue;
    end

    delta = cases(iC).delta;

    % ---- exact file-finding from run_piv_structure_comparison Fig1 ----
    files = dir(fullfile(cases(iC).xref_dir, 'R_uu_xref*_yref*.mat'));
    if isempty(files), continue; end

    yr_available = nan(numel(files), 1);
    xr_available = nan(numel(files), 1);
    for f = 1:numel(files)
        tok = regexp(files(f).name, ...
            'R_uu_xref([-+]?\d*\.?\d+)_yref([-+]?\d*\.?\d+)', ...
            'tokens', 'once');
        if ~isempty(tok)
            xr_available(f) = str2double(tok{1});
            yr_available(f) = str2double(tok{2});
        end
    end

    xr_ok = abs(xr_available - xref_nom) < 5;
    if ~any(xr_ok)
        [~, iBest] = min(abs(xr_available - xref_nom));
        xr_ok(iBest) = true;
    end
    % -------------------------------------------------------------------

    fprintf('\n  %s  (delta99 = %.2f mm)\n', cases(iC).label, delta);

    for iT = 1:nTargets
        yr_target = ydelta_targets(iT) * delta;

        yr_sub         = yr_available;
        yr_sub(~xr_ok) = Inf;
        [~, iClosest]  = min(abs(yr_sub - yr_target));

        S  = load(fullfile(cases(iC).xref_dir, files(iClosest).name), 'R_s', 'xr', 'yr');
        R  = double(S.R_s);
        xr = double(S.xr);
        yr = double(S.yr);

        fprintf('    iT=%d | y/d target=%.2f | actual=%.3f\n', ...
            iT, ydelta_targets(iT), yr/delta);

        % Separation grids (verbatim from run_piv_structure_comparison)
        dX = worldX - xr;
        dY = worldY - yr;

        x_lo = max(xr - dx_back, x_min_domain);
        x_hi = min(xr + dx_fwd,  x_max_domain);

        boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
        R_box           = R;
        R_box(~boxMask) = NaN;
        R_box(R_box >  1.0) = 1.0;
        R_box(R_box < -1.0) = NaN;

        dX_vec = dX(1, :);
        dY_vec = dY(:, 1);

        dist2           = dX.^2 + dY.^2;
        dist2(~boxMask) = Inf;
        [~, i0]         = min(dist2(:));

        for iL = 1:nLevels
            rho  = rho_levels(iL);

            ell_cache{iC, iT, iL} = struct('valid', false, ...
                'a', NaN, 'b', NaN, 'cx', NaN, 'cy', NaN, 'theta', NaN);

            mask = R_box >= rho;
            if sum(mask(:)) < 5, continue; end

            CC = bwconncomp(mask);
            if CC.NumObjects == 0, continue; end

            regionIdx = 0;
            for r = 1:CC.NumObjects
                if any(CC.PixelIdxList{r} == i0)
                    regionIdx = r; break;
                end
            end
            if regionIdx == 0
                [~, regionIdx] = max(cellfun(@numel, CC.PixelIdxList));
            end

            mask_origin = false(size(mask));
            mask_origin(CC.PixelIdxList{regionIdx}) = true;

            B_contour = bwboundaries(mask_origin);
            if isempty(B_contour), continue; end

            bnd    = B_contour{1};
            bnd_dx = dX_vec(bnd(:,2));   % physical mm
            bnd_dy = dY_vec(bnd(:,1));   % physical mm

            % Store normalised boundary
            bnd_cache{iC, iT, iL} = [bnd_dx(:)/delta, bnd_dy(:)/delta];

            % ---- Ellipse fit: verbatim tip-based method from get_contour_2 ----
            dx_spacing = abs(dX_vec(2) - dX_vec(1));
            dy_spacing = abs(dY_vec(2) - dY_vec(1));

            up_mask  = bnd_dx <= (min(bnd_dx) + dx_spacing);
            dn_mask  = bnd_dx >= (max(bnd_dx) - dx_spacing);
            top_mask = bnd_dy >= (max(bnd_dy) - dy_spacing);
            bot_mask = bnd_dy <= (min(bnd_dy) + dy_spacing);

            x_up  = mean(bnd_dx(up_mask));   y_up  = mean(bnd_dy(up_mask));
            x_dn  = mean(bnd_dx(dn_mask));   y_dn  = mean(bnd_dy(dn_mask));
            x_top = mean(bnd_dx(top_mask));  y_top = mean(bnd_dy(top_mask));
            x_bot = mean(bnd_dx(bot_mask));  y_bot = mean(bnd_dy(bot_mask));

            ell.theta = atan2d(y_dn - y_up, x_dn - x_up);
            ell.a     = sqrt((x_dn - x_up)^2 + (y_dn - y_up)^2) / 2;
            ell.cx    = (x_up  + x_dn)  / 2;
            ell.cy    = (y_top + y_bot)  / 2;

            th_rad   = deg2rad(ell.theta);
            perp_x   = -sin(th_rad);
            perp_y   =  cos(th_rad);
            proj_top = (x_top - ell.cx)*perp_x + (y_top - ell.cy)*perp_y;
            proj_bot = (x_bot - ell.cx)*perp_x + (y_bot - ell.cy)*perp_y;
            b_raw    = abs(proj_top - proj_bot) / 2;

            % Wall-clamp: ellipse lower extent must not go below min(bnd_dy)
            b_lower_tip_y = ell.cy - b_raw * cos(th_rad);
            if b_lower_tip_y < min(bnd_dy) && abs(cos(th_rad)) > 0.1
                b_clamped = (ell.cy - min(bnd_dy)) / abs(cos(th_rad));
                ell.b = min(b_raw, b_clamped);
            else
                ell.b = b_raw;
            end

            ell.valid = true;
            ell_cache{iC, iT, iL} = ell;
        end
    end
end
fprintf('\nDone.\n\n');

%% ===================== FIGURE: CONTOURS + ELLIPSES OVERLAID =============
% Two-pass draw order per panel:
%   Pass 1 -- pale raw isocontours (background)
%   Pass 2 -- full-opacity fitted ellipses (foreground)

fig = figure('Color', 'w', ...
    'Position', [80 80 480 260*nTargets], ...
    'Name', sprintf('R_{uu} contours + ellipses  |  %s  (x_{ref} = %.0f mm)', ...
    upper(regime_request), xref_nom));

axC = gobjects(nTargets, 1);

%% ===================== AXIS LIMITS ======================================

y_row_lim = nan(nTargets, 2);
x_all_min =  Inf;
x_all_max = -Inf;

for iT = 1:nTargets
    y_min_row =  Inf;
    y_max_row = -Inf;
    for iC = 1:nCases
        for iL = 1:nLevels
            b = bnd_cache{iC, iT, iL};
            if isempty(b), continue; end
            y_min_row = min(y_min_row, min(b(:,2)));
            y_max_row = max(y_max_row, max(b(:,2)));
            x_all_min = min(x_all_min, min(b(:,1)));
            x_all_max = max(x_all_max, max(b(:,1)));
        end
    end
    if isfinite(y_min_row)
        margin = max((y_max_row - y_min_row) * 0.20, 0.05);
        y_row_lim(iT, :) = [y_min_row - margin, y_max_row + margin];
    end
end

if isfinite(x_all_min)
    x_margin  = (x_all_max - x_all_min) * 0.05;
    x_lim_all = [x_all_min - x_margin, x_all_max + x_margin];
else
    x_lim_all = [-5 5];
end

for iT = 1:nTargets
    axC(iT) = subplot(nTargets, 1, iT);
    hold on; box on; grid on;
    set(axC(iT), 'FontSize', fontSize);

    % Pass 1: pale raw isocontours
    for iC = 1:nCases
        col_pale = caseColPale(iC, :);
        for iL = 1:nLevels
            b = bnd_cache{iC, iT, iL};
            if isempty(b), continue; end
            plot(axC(iT), b(:,1), b(:,2), rho_linestyle{iL}, ...
                'Color', col_pale, 'LineWidth', rho_linewidth(iL), ...
                'HandleVisibility', 'off');
        end
    end

    % Pass 2: full-opacity fitted ellipses
    for iC = 1:nCases
        col   = caseColours(iC, :);
        delta = cases(iC).delta;
        for iL = 1:nLevels
            e = ell_cache{iC, iT, iL};
            if ~e.valid || ~isfinite(delta), continue; end
            a  = e.a  / delta;
            eb = e.b  / delta;
            cx = e.cx / delta;
            cy = e.cy / delta;
            th = e.theta;
            xe = cx + a .*cos(t_ell).*cosd(th) - eb.*sin(t_ell).*sind(th);
            ye = cy + a .*cos(t_ell).*sind(th) + eb.*sin(t_ell).*cosd(th);
            plot(axC(iT), xe, ye, rho_linestyle{iL}, ...
                'Color', col, 'LineWidth', rho_linewidth(iL), ...
                'HandleVisibility', 'off');
        end
    end

    xline(axC(iT), 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
    yline(axC(iT), 0, 'k:', 'LineWidth', 0.7, 'HandleVisibility', 'off');

    ylabel(axC(iT), '$\Delta y / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
    if iT == nTargets
        xlabel(axC(iT), '$\Delta x / \delta$', 'Interpreter', 'latex', 'FontSize', fontSize);
    else
        set(axC(iT), 'XTickLabel', []);
    end

    xlim(axC(iT), x_lim_all);
    if all(isfinite(y_row_lim(iT, :)))
        ylim(axC(iT), y_row_lim(iT, :));
    end
end

linkaxes(axC, 'x');

%% ===================== ROW ANNOTATIONS & PANEL LABELS ===================

drawnow;
panelLetters = 'abcdefghijklmnopqrstuvwxyz';

for iT = 1:nTargets
    ax  = axC(iT);
    pos = ax.Position;

    annotation(fig, 'textbox', ...
        [pos(1), pos(2) + pos(4) - 0.04, 0.04, 0.04], ...
        'String',              ['(' panelLetters(iT) ')'], ...
        'EdgeColor',           'none', ...
        'FontSize',            fontSize, ...
        'FontWeight',          'bold', ...
        'Interpreter',         'none', ...
        'VerticalAlignment',   'top', ...
        'HorizontalAlignment', 'left');

    x_tag = pos(1) + pos(3) + 0.005;
    ax_tag = axes(fig, 'Position', [x_tag, pos(2), 0.05, pos(4)], 'Visible', 'off');

    if ydelta_targets(iT) == 0.05
        rowStr = '$y_{\mathrm{ref}}/\delta = 0.05$';
    else
        rowStr = ['$y_{\mathrm{ref}}/\delta = ' sprintf('%g', ydelta_targets(iT)) '$'];
    end

    text(ax_tag, 0.5, 0.5, rowStr, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Interpreter', 'latex', ...
        'FontSize', fontSize, 'FontWeight', 'bold', 'Rotation', 90);
end

%% ===================== EXPORT ===========================================

save_dir  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\figures\';
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

fig_width_cm  = 9.0;
fig_height_cm = 10;   % ~4 cm per row

drawnow;
fig.Units        = 'centimeters';
fig.Position(3:4) = [fig_width_cm, fig_height_cm];

filename = fullfile(save_dir, ...
    sprintf('Ruu_isocontour_multiRho_%s_%s.pdf', upper(regime_request), timestamp));

exportgraphics(fig, filename, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Saved: %s\n', filename);
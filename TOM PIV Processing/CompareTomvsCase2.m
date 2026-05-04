% =========================================================================
% CompareTomvsCase2.m
%
% Compares bounding box descriptors (Lx/delta, Ly/delta, Lx/Ly, theta)
% between Case AK and Case SW across matched flow regimes.
%
% Colour scheme (shared with CompareTomvsCase2.m):
%   AK:  grayscale gradient dark->light  (solid lines, square markers)
%   SW:  blue    gradient dark->light  (dashed lines, circle markers)
%   Both gradients run Inlet (dark) -> ZPG recovery (light)
%
% Layout: 4 figures, each with 1 x nLevels subplots (one per rho level)
%   - x axis: y_ref / delta99  (log scale)
%   - lines:  one per xref, colour = regime shade
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

rho_levels_plot = [0.3 , 0.6];
rho_labels      = {'\rho = 0.3', '\rho = 0.6'};

x_lim = [];   % x axis limits [y/delta] — [] for auto

% =========================================================================
% SHARED COLOUR SCHEME
%   6 regimes: Inlet | FPG max | Crossover | APG max | TE-relative | ZPG rec.
%   AK:  grayscale  dark -> light   (solid lines,  square markers)
%   SW:  blue       dark -> light   (dashed lines, circle markers)
%   KEEP THIS BLOCK IDENTICAL in both scripts.
% =========================================================================
t_vals = linspace(0, 1, 6);   % 0 = darkest, 1 = lightest

% AK: near-black -> light grey
ak_dark   = 	[1 0 0]; %[0.8500, 0.3250, 0.0980];
ak_dark = [0.65, 0.00, 0.00];
ak_light  = [0.95, 0.70, 0.70];
ak_colors = cell(1, 6);
for k = 1:6
    ak_colors{k} = (1 - t_vals(k)) * ak_dark + t_vals(k) * ak_light;
end


% sw_dark   = [0.05, 0.20, 0.45];
% sw_light  = [0.45, 0.72, 0.95];

sw_dark  = [0 0 1]; %[0, 0.4470, 0.7410];   % MATLAB default blue
sw_dark = [0.00, 0.00, 0.65];

sw_light = [0.68, 0.85, 0.95];    % light sky blue


sw_colors = cell(1, 6);







for k = 1:6
    sw_colors{k} = (1 - t_vals(k)) * sw_dark + t_vals(k) * sw_light;
end

% {1} = Inlet (darkest)  ...  {6} = ZPG recovery (lightest)

% --- Case AK ---
ak.blSweepFile  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';
ak.bboxDir      = 'G:\two_point_covariance_20260426_120726\bbox_analysis\';
ak.xref_noms    = [54, 550, 1000];   % mm  (350, 711 commented out)
ak.regimeLabels = {'Inlet', 'Crossover', '\Delta(x_{TE}-x)/\delta=3'};
% Regime index (into colour arrays) for each AK xref:
%   54=Inlet(1), 550=Crossover(3), 1000=TE-rel(5)
ak.regime_idx   = [1, 3, 5];

% --- Case SW (per-xref blSweep) ---
sw.blSweepFiles = { ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215729_Pos1.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215751_Pos2.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215811_Pos3.mat', ...
};
% sw.blSweepFiles = { ...
%     'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h400mm\blSweep_20260422_233218_Pos1.mat', ...
%     'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h400mm\blSweep_20260422_233236_Pos2.mat', ...
%     'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h400mm\blSweep_20260422_233256_Pos3.mat', ...
% };

sw.bboxDirs = { ...
    'G:\SW 500mm Mean Flow Fields\Pos_1\two_point_covariance_20260424_225929\bbox_analysis\', ...
    'G:\SW 500mm Mean Flow Fields\Pos_2\two_point_covariance_20260424_230729\bbox_analysis\', ...
    'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260424_231444\bbox_analysis\'  ...
};

% sw.bboxDirs = { ...
%     'G:\SW 400mm Processed Data Single Snapshots\-8\Pos_1\two_point_covariance_20260423_122842\bbox_analysis\', ...
%     'G:\SW 400mm Processed Data Single Snapshots\-8\Pos_2\two_point_covariance_20260423_114635\bbox_analysis\', ...
%     'G:\SW 400mm Processed Data Single Snapshots\-8\Pos_3\two_point_covariance_20260423_113156\bbox_analysis\', ...
% };


sw.xref_noms  = [5866, 6921, 8000];   % (6480, 7320, 8573 commented out)
% Regime index (into colour arrays) for each SW xref — must match ak.regime_idx:
%   5866=Inlet(1), 6921=Crossover(3), 8000=TE-rel(5)
sw.regime_idx   = [1, 3, 5];
sw.regimeLabels = ak.regimeLabels;

%% ===================== LOAD AK delta99 =====================

nXref = numel(ak.xref_noms);

B_ak     = load(ak.blSweepFile, 'blSweep');
bl_x     = B_ak.blSweep.x_mm;
bl_d99   = B_ak.blSweep.delta99_hybrid_mm;
valid    = isfinite(bl_x) & isfinite(bl_d99);
bl_x_ak  = bl_x(valid);  bl_d99_ak = bl_d99(valid);
if any(diff(bl_x_ak) == 0)
    [bl_x_ak, ~, ic] = unique(bl_x_ak, 'sorted');
    bl_d99_ak = accumarray(ic, bl_d99_ak, [], @mean);
end
delta_ak = interp1(bl_x_ak, bl_d99_ak, ak.xref_noms, 'linear', NaN);

fprintf('AK delta99:\n');
for i = 1:nXref
    fprintf('  x=%5.0f mm -> %.2f mm\n', ak.xref_noms(i), delta_ak(i));
end

%% ===================== LOAD SW delta99 =====================

delta_sw = nan(1, nXref);
for iX = 1:nXref
    B      = load(sw.blSweepFiles{iX}, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;
    valid  = isfinite(bl_x) & isfinite(bl_d99);
    bl_xc  = bl_x(valid);  bl_dc = bl_d99(valid);
    if any(diff(bl_xc) == 0)
        [bl_xc, ~, ic] = unique(bl_xc, 'sorted');
        bl_dc = accumarray(ic, bl_dc, [], @mean);
    end
    delta_sw(iX) = interp1(bl_xc, bl_dc, sw.xref_noms(iX), 'linear', NaN);
end

fprintf('\nSW delta99:\n');
for i = 1:nXref
    fprintf('  x=%5.0f mm -> %.2f mm\n', sw.xref_noms(i), delta_sw(i));
end

%% ===================== LOAD BBOX DATA =====================

ak_data = cell(1, nXref);
sw_data = cell(1, nXref);

for iX = 1:nXref
    ak_data{iX} = load_bbox(ak.bboxDir,       ak.xref_noms(iX));
    sw_data{iX} = load_bbox(sw.bboxDirs{iX},  sw.xref_noms(iX));
end

%% ===================== FIGURES =====================

nLevels  = numel(rho_levels_plot);
figNames = {'a / \delta', 'b / \delta', 'a / b', '\theta  (deg)'};
for iFig = 1:4

    figure('Color', 'w', ...
        'Position', [50 50 450*nLevels 500], ...
        'Name', figNames{iFig});

    axArr = gobjects(1, nLevels);

    for iL = 1:nLevels

        axArr(iL) = subplot(1, nLevels, iL);
        hold on; box on; grid on;

        leg_h   = gobjects(0);
        leg_lbl = {};

        for iX = 1:nXref

            col_ak = ak_colors{ak.regime_idx(iX)};
            col_sw = sw_colors{sw.regime_idx(iX)};

            % ---- AK: solid line + filled square markers (grayscale) ----
            D_ak = ak_data{iX};
            if ~isempty(D_ak) && isfinite(delta_ak(iX))
                [xax, yax] = extract_series(D_ak, iFig, iL, ...
                    rho_levels_plot, delta_ak(iX));
                if ~isempty(xax)
                    [xax, ui] = unique(xax, 'sorted');
                    yax = yax(ui);
                    h = plot(axArr(iL), xax, yax, '-', ...
                        'Color', col_ak, 'LineWidth', 1.4, ...
                        'DisplayName', ak.regimeLabels{iX});
                    valid_pts = isfinite(yax);
                    plot(axArr(iL), xax(valid_pts), yax(valid_pts), 's', ...
                        'Color', col_ak, 'MarkerFaceColor', col_ak, ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
                        'HandleVisibility', 'off');
                    leg_h(end+1)   = h;                    %#ok<AGROW>
                    leg_lbl{end+1} = ak.regimeLabels{iX}; %#ok<AGROW>
                end
            end

            % ---- SW: dashed line + filled circle markers (blue) ----
            D_sw = sw_data{iX};
            if ~isempty(D_sw) && isfinite(delta_sw(iX))
                [xax, yax] = extract_series(D_sw, iFig, iL, ...
                    rho_levels_plot, delta_sw(iX));
                if ~isempty(xax)
                    [xax, ui] = unique(xax, 'sorted');
                    yax = yax(ui);
                    plot(axArr(iL), xax, yax, '--', ...
                        'Color', col_sw, 'LineWidth', 1.4, ...
                        'HandleVisibility', 'off');
                    valid_pts = isfinite(yax);
                    plot(axArr(iL), xax(valid_pts), yax(valid_pts), 'o', ...
                        'Color', col_sw, 'MarkerFaceColor', col_sw, ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
                        'HandleVisibility', 'off');
                end
            end

        end % xref loop

        set(axArr(iL), 'XScale', 'log');
        xlabel(axArr(iL), '$y / \delta$', 'Interpreter', 'latex', 'FontSize', 12);

        if iL == 1
            ylabel(axArr(iL), ['$' get_ylabel_latex(iFig) '$'], ...
                'Interpreter', 'latex', 'FontSize', 12);
        end

        title(axArr(iL), rho_labels{iL}, 'Interpreter', 'tex', 'FontSize', 12);
        if ~isempty(x_lim), xlim(axArr(iL), x_lim); end

        % ---- Legend: first subplot only ----
        if iL == 1
            % Case style proxies
            h_ak_style = plot(axArr(iL), NaN, NaN, '-s', ...
                'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.4, ...
                'DisplayName', 'Case 2');
            h_sw_style = plot(axArr(iL), NaN, NaN, '--o', ...
                'Color', [0 0 1], 'MarkerFaceColor', [0 0 1], ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.4, ...
                'DisplayName', 'VP2025');
            if ~isempty(leg_h)
                legend(axArr(iL), [h_ak_style; h_sw_style; leg_h(:)], ...
                    [{'AK', 'VP2025'}, leg_lbl], ...
                    'Location', 'best', 'FontSize', 10, 'Interpreter', 'tex');
            end
        end

    end % level loop

    linkaxes(axArr, 'x');
    sgtitle(['$' get_ylabel_latex(iFig) '$'], ...
        'Interpreter', 'latex', 'FontSize', 13);

end % figure loop

%% ===================== LOCAL FUNCTIONS =====================

function D = load_bbox(bboxDir, xref_nom)
D = [];
fname = fullfile(bboxDir, sprintf('bbox_results_xref_%g.mat', xref_nom));
fprintf('  [load_bbox] Looking for: %s\n', fname);
if ~isfile(fname)
    fname2 = fullfile(bboxDir, sprintf('bbox_results_xref_%d.mat', round(xref_nom)));
    fprintf('  [load_bbox] Trying alt:  %s\n', fname2);
    if ~isfile(fname2)
        warning('bbox file not found: %s', fname);
        return;
    end
    fname = fname2;
end
fprintf('  [load_bbox] Loaded: %s\n', fname);
D = load(fname);
fprintf('  [load_bbox] Fields: %s\n', strjoin(fieldnames(D)', ', '));
fprintf('  [load_bbox] yr_vec: [%.2f, %.2f] mm, %d points\n', ...
    min(D.yr_vec), max(D.yr_vec), numel(D.yr_vec));
fprintf('  [load_bbox] rho_levels: %s\n', mat2str(D.rho_levels, 3));

[~, idx] = sort(D.yr_vec);
D.yr_vec          = D.yr_vec(idx);
D.Lx_mat          = D.Lx_mat(idx,:);
D.Ly_mat          = D.Ly_mat(idx,:);
D.aspect_mat      = D.aspect_mat(idx,:);
D.a_ellipse_mat   = D.a_ellipse_mat(idx,:);
D.b_ellipse_mat   = D.b_ellipse_mat(idx,:);

if isfield(D, 'theta_mat') && ~isempty(D.theta_mat)
    D.theta_mat = D.theta_mat(idx,:);
end
end
% -------------------------------------------------------------------------

function [xax, yax] = extract_series(D, iFig, iL, rho_levels_plot, delta)
xax = [];
yax = [];
if isempty(D), return; end
rho_target = rho_levels_plot(iL);
[~, col_idx] = min(abs(D.rho_levels - rho_target));
fprintf('  [extract] iFig=%d iL=%d rho_target=%.3f -> col_idx=%d (actual rho=%.3f)\n', ...
    iFig, iL, rho_target, col_idx, D.rho_levels(col_idx));
yr_vec = D.yr_vec;
xax    = yr_vec / delta;
switch iFig
    case 1;  yax = D.a_ellipse_mat(:, col_idx) / delta;
case 2;  yax = D.b_ellipse_mat(:, col_idx) / delta;
case 3;  yax = D.a_ellipse_mat(:, col_idx) ./ D.b_ellipse_mat(:, col_idx);
    case 4
        if isfield(D, 'theta_mat') && ~isempty(D.theta_mat)
            yax = D.theta_mat(:, col_idx);
        else
            xax = []; return;
        end
end
fprintf('    [extract] %d / %d valid points\n', sum(isfinite(yax)), numel(yax));
end

% -------------------------------------------------------------------------

function lbl = get_ylabel_latex(iFig)
switch iFig
    case 1; lbl = 'a / \delta';
    case 2; lbl = 'b / \delta';
    case 3; lbl = 'a / b';
    case 4; lbl = '\theta \ \mathrm{(deg)}';
end
end
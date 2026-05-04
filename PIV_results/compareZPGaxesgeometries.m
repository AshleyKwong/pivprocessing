% =========================================================================
% CompareTomvsCase2.m
%
% Compares bounding box descriptors (a/delta, b/delta, a/b, theta)
% between Case AK and Case VP across matched flow regimes.
%
% Colour scheme matches plot_ellipse_compare_2x2.m exactly:
%   AK:  red   family, light -> dark  (solid lines,  square markers)
%   VP:  blue  family, light -> dark  (dashed lines, circle markers)
%   Both gradients run Inlet (light) -> TE-relative (dark)
%
% Layout: 4 figures, each with 1 x nLevels subplots (one per rho level)
%   - x axis: y_ref / delta99  (log scale)
%   - lines:  one per xref, colour = station shade
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;


fs = 14;   % font size — applied to all text in the figure
% --- Export dimensions (cm) ----------------------------------------------
fig_width_cm  = 17.0;
fig_height_cm =  8.0;
save_dir      = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\figures';
rho_levels_plot = [0.3, 0.6];
rho_labels      = {'$\rho_{uu} = 0.3$', '$\rho_{uu} = 0.6$'};

x_lim = [];   % x axis limits [y/delta] — [] for auto

% =========================================================================
% COLOUR SCHEME — matches plot_ellipse_compare_2x2.m exactly
%   3 stations: Inlet (light) -> Crossover (mid) -> TE-relative (dark)
% =========================================================================
nXref  = 3;
t_vals = [0.0, 0.5, 1.0];   % pacing matches ellipse script

% AK: rose -> crimson
ak_light = [0.90, 0.55, 0.55];
ak_dark  = [0.65, 0.00, 0.00];

% VP: sky -> navy
vp_light = [0.30, 0.55, 0.85];
vp_dark  = [0.00, 0.00, 0.65];

ak_colors = cell(1, nXref);
vp_colors = cell(1, nXref);
for k = 1:nXref
    ak_colors{k} = (1 - t_vals(k)) * ak_light + t_vals(k) * ak_dark;
    vp_colors{k} = (1 - t_vals(k)) * vp_light + t_vals(k) * vp_dark;
end

% =========================================================================
% CASE AK
% =========================================================================
ak.label        = 'Case 2';
ak.blSweepFile  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';
ak.bboxDir      = 'G:\two_point_covariance_20260426_120726\bbox_analysis\';
ak.xref_noms    = [54, 550, 1000];   % mm
ak.stationLabels = {'Inlet', 'Crossover', '$\Delta(x_{TE}-x)/\delta = 3$'};

% =========================================================================
% CASE VP
% =========================================================================
vp.label        = 'VP2025';
vp.blSweepFiles = { ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215729_Pos1.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215751_Pos2.mat', ...
    'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215811_Pos3.mat', ...
};
vp.bboxDirs = { ...
    'G:\SW 500mm Mean Flow Fields\Pos_1\two_point_covariance_20260424_225929\bbox_analysis\', ...
    'G:\SW 500mm Mean Flow Fields\Pos_2\two_point_covariance_20260424_230729\bbox_analysis\', ...
    'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260424_231444\bbox_analysis\', ...
};
vp.xref_noms    = [5866, 6921, 8000];
vp.stationLabels = ak.stationLabels;   % same regime labels

%% ===================== LOAD AK delta99 =====================

B_ak     = load(ak.blSweepFile, 'blSweep');
bl_x     = B_ak.blSweep.x_mm;
bl_d99   = B_ak.blSweep.delta99_hybrid_mm;
valid    = isfinite(bl_x) & isfinite(bl_d99);
bl_x_ak  = bl_x(valid);
bl_d99_ak = bl_d99(valid);
if any(diff(bl_x_ak) == 0)
    [bl_x_ak, ~, ic] = unique(bl_x_ak, 'sorted');
    bl_d99_ak = accumarray(ic, bl_d99_ak, [], @mean);
end
delta_ak = interp1(bl_x_ak, bl_d99_ak, ak.xref_noms, 'linear', NaN);

fprintf('AK delta99:\n');
for i = 1:nXref
    fprintf('  x = %5.0f mm  ->  %.2f mm\n', ak.xref_noms(i), delta_ak(i));
end

%% ===================== LOAD VP delta99 =====================

delta_vp = nan(1, nXref);
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
end

fprintf('\nVP delta99:\n');
for i = 1:nXref
    fprintf('  x = %5.0f mm  ->  %.2f mm\n', vp.xref_noms(i), delta_vp(i));
end

%% ===================== LOAD BBOX DATA =====================

ak_data = cell(1, nXref);
vp_data = cell(1, nXref);

for iX = 1:nXref
    ak_data{iX} = load_bbox(ak.bboxDir,      ak.xref_noms(iX));
    vp_data{iX} = load_bbox(vp.bboxDirs{iX}, vp.xref_noms(iX));
end

%% ===================== FIGURES =====================
timestamp     = datestr(now, 'yyyymmdd_HHMMSS');
file_tags     = {'a_over_delta_NOLEGEND', 'b_over_delta_NOLEGEND', 'a_over_b_NOLEGEND', 'theta_deg_NOLEGEND'};
cm2in         = 1 / 2.54;
fig_width_in  = fig_width_cm  * cm2in;
fig_height_in = fig_height_cm * cm2in;

nLevels  = numel(rho_levels_plot);
figNames = {'$a \,/\, \delta$', '$b \,/\, \delta$', '$a \,/\, b$', '$\theta$ (deg)'};

for iFig = 1:4

    fig = figure('Color', 'w', ...
        'Position', [50 50 420*nLevels 480], ...
        'Name', figNames{iFig});

    axArr = gobjects(1, nLevels);

    for iL = 1:nLevels

        axArr(iL) = subplot(1, nLevels, iL);
        hold on; box on; grid on;
        set(axArr(iL), 'GridAlpha', 0.15, 'GridLineStyle', ':');

        for iX = 1:nXref

            col_ak = ak_colors{iX};
            col_vp = vp_colors{iX};

            % ---- AK: solid line + filled square markers -----------------
            D_ak = ak_data{iX};
            if ~isempty(D_ak) && isfinite(delta_ak(iX))
                [xax, yax] = extract_series(D_ak, iFig, iL, ...
                    rho_levels_plot, delta_ak(iX));
                if ~isempty(xax)
                    [xax, ui] = unique(xax, 'sorted');
                    yax = yax(ui);
                    plot(axArr(iL), xax, yax, '-', ...
                        'Color', col_ak, 'LineWidth', 1.4, ...
                        'HandleVisibility', 'off');
                    valid_pts = isfinite(yax);
                    plot(axArr(iL), xax(valid_pts), yax(valid_pts), 's', ...
                        'Color', col_ak, 'MarkerFaceColor', col_ak, ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
                        'HandleVisibility', 'off');
                end
            end

            % ---- VP: dashed line + filled circle markers ----------------
            D_vp = vp_data{iX};
            if ~isempty(D_vp) && isfinite(delta_vp(iX))
                [xax, yax] = extract_series(D_vp, iFig, iL, ...
                    rho_levels_plot, delta_vp(iX));
                if ~isempty(xax)
                    [xax, ui] = unique(xax, 'sorted');
                    yax = yax(ui);
                    plot(axArr(iL), xax, yax, '--', ...
                        'Color', col_vp, 'LineWidth', 1.4, ...
                        'HandleVisibility', 'off');
                    valid_pts = isfinite(yax);
                    plot(axArr(iL), xax(valid_pts), yax(valid_pts), 'o', ...
                        'Color', col_vp, 'MarkerFaceColor', col_vp, ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
                        'HandleVisibility', 'off');
                end
            end

        end % xref loop

        set(axArr(iL), 'XScale', 'log', ...
            'TickLabelInterpreter', 'latex', 'FontSize', fs);
        xlabel(axArr(iL), '$y_{ref} \,/\, \delta$', ...
            'Interpreter', 'latex', 'FontSize', fs);

        if iL == 1
            ylabel(axArr(iL), ['$' get_ylabel_latex(iFig) '$'], ...
                'Interpreter', 'latex', 'FontSize', fs);
        end

        title(axArr(iL), rho_labels{iL}, ...
            'Interpreter', 'latex', 'FontSize', fs, 'FontWeight', 'bold');
        % Panel label (a), (b), ... top-left corner in axes-normalised units
        % panel_labels = 'abcdefghij';
        panel_labels = 'ababababab';

        text(axArr(iL), 0.03, 0.97, ...
            sprintf('(%s)', panel_labels((iFig-1)*nLevels + iL)), ...
            'Units',               'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment',   'top', ...
            'Interpreter',         'latex', ...
            'FontSize',            fs, ...
            'FontWeight',          'bold');
        if ~isempty(x_lim)
            xlim(axArr(iL), x_lim);
        end

    end % level loop

    % ---- Legend: top-left panel only ------------------------------------
    a_leg = axArr(1);

    % Case style proxies
    h_ak_style = plot(a_leg, NaN, NaN, '-s', ...
        'Color', ak_dark, 'MarkerFaceColor', ak_dark, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.4, ...
        'DisplayName', ak.label);
    h_vp_style = plot(a_leg, NaN, NaN, '--o', ...
        'Color', vp_dark, 'MarkerFaceColor', vp_dark, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.4, ...
        'DisplayName', vp.label);

    % Station colour swatches
    % h_st = gobjects(nXref, 1);
    % for iX = 1:nXref
    %     col_mid = 0.5 * (ak_colors{iX} + vp_colors{iX});
    %     h_st(iX) = plot(a_leg, NaN, NaN, 's', ...
    %         'MarkerFaceColor', col_mid, ...
    %         'MarkerEdgeColor', col_mid * 0.6, ...
    %         'MarkerSize', 8, ...
    %         'DisplayName', ak.stationLabels{iX});
    % end

    % legend(a_leg, [h_ak_style; h_vp_style], ...
    %     'Interpreter', 'latex', ...
    %     'Location',    'best', ...
    %     'FontSize',    fs - 2, ...
    %     'Box',         'on');

    linkaxes(axArr, 'x');

    sgtitle(figNames{iFig}, 'Interpreter', 'latex', 'FontSize', fs);
    drawnow;
    set(fig, 'Units', 'inches', 'Position', [1, 1, fig_width_in, fig_height_in]);
    set(fig, 'PaperUnits', 'inches', ...
        'PaperSize',     [fig_width_in, fig_height_in], ...
        'PaperPosition', [0, 0, fig_width_in, fig_height_in]);
    filename = fullfile(save_dir, ...
        sprintf('bbox_compare_%s_%s.pdf', file_tags{iFig}, timestamp));
    exportgraphics(fig, filename, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Saved: %s\n', filename);
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
D.yr_vec        = D.yr_vec(idx);
D.Lx_mat        = D.Lx_mat(idx,:);
D.Ly_mat        = D.Ly_mat(idx,:);
D.aspect_mat    = D.aspect_mat(idx,:);
D.a_ellipse_mat = D.a_ellipse_mat(idx,:);
D.b_ellipse_mat = D.b_ellipse_mat(idx,:);
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
    case 1; lbl = 'a \,/\, \delta';
    case 2; lbl = 'b \,/\, \delta';
    case 3; lbl = 'a \,/\, b';
    case 4; lbl = '\theta \ \mathrm{(deg)}';
end
end
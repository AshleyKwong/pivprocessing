% =========================================================================
% plot_bbox_results.m
%
% Loads bounding box descriptor results produced by analyse_Ruu_maps.m
% and produces summary figures.
%
% Two modes:
%   'sweep' — fixed yref, varying xref → x axis = xr_vec, lines = yref
%   'point' — fixed xref, varying yref → x axis = yr_vec, lines = xref
%
% Figures produced:
%   1-7  — Lx, Ly, aspect, Lxu, Lxd, Lyt, Lyb vs independent variable
%   8-10 — Threshold sensitivity: Lx, Ly, aspect
%   11   — Ly vs Lx scatter
%   12   — theta vs independent variable (point mode only)
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

corrMode   = 'point';

resultsDir = 'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\bbox_analysis\';

fileLabels = {
    'xref_1150'
};
lineVals   = [1150];   % xref values in mm
lineName   = 'x_{ref}';
xaxisLabel = 'y_{ref}  (mm)';
xaxisField = 'yr';

% --- For point mode: normalise by delta99 ---
normalise_by_delta = true;

blSweepFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat'];

outlierWindow = 1;
x_lim = [];

%% ===================== LOAD delta99 =====================

if strcmp(corrMode, 'point') && normalise_by_delta

    B      = load(blSweepFile, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;

    % Remove non-finite points
    valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
    bl_x_clean   = bl_x(valid_bl);
    bl_d99_clean = bl_d99(valid_bl);

    % Deduplicate — average delta99 at repeated x locations
    % (duplicates arise from mag_factor window overlap in bl_sweep)

    if any(diff(bl_x_clean) == 0)
        [bl_x_clean, ~, ic] = unique(bl_x_clean, 'sorted');
        bl_d99_clean        = accumarray(ic, bl_d99_clean, [], @mean);
        fprintf('Duplicates found and removed: %d -> %d unique points\n', ...
            numel(bl_x), numel(bl_x_clean));
    end

    % [bl_x_clean, ~, ic] = unique(bl_x_clean, 'sorted');
    % bl_d99_clean        = accumarray(ic, bl_d99_clean, [], @mean);
    % 
    % fprintf('bl_x unique points after deduplication: %d\n', numel(bl_x_clean));

    % delta99 at each xref location
    delta_vals = interp1(bl_x_clean, bl_d99_clean, lineVals, 'linear', NaN);

    fprintf('delta99 at each xref:\n');
    for i = 1:numel(lineVals)
        fprintf('  x = %.1f mm  ->  delta99 = %.2f mm\n', ...
            lineVals(i), delta_vals(i));
    end

    xaxisLabel_norm = 'y_{ref} / \delta_{99}';
    Lscale_label    = '/ \delta_{99}';

else
    delta_vals          = ones(1, numel(lineVals));
    normalise_by_delta  = false;
    xaxisLabel_norm     = xaxisLabel;
    Lscale_label        = '  (mm)';
end

%% ===================== LOAD RESULTS =====================

nLines  = numel(fileLabels);
allData = struct();

for iY = 1:nLines
    fname = fullfile(resultsDir, ...
        sprintf('bbox_results_%s.mat', fileLabels{iY}));

    if ~isfile(fname)
        warning('File not found: %s — skipping', fname);
        allData(iY).loaded = false;
        continue;
    end

    D = load(fname);
    allData(iY).loaded  = true;
    allData(iY).xr      = D.xr_vec;
    allData(iY).yr      = D.yr_vec;
    allData(iY).Lx      = D.Lx_mat;
    allData(iY).Ly      = D.Ly_mat;
    allData(iY).Lxu     = D.Lxu_mat;
    allData(iY).Lxd     = D.Lxd_mat;
    allData(iY).Lyt     = D.Lyt_mat;
    allData(iY).Lyb     = D.Lyb_mat;
    allData(iY).aspect  = D.aspect_mat;
    allData(iY).rho     = D.rho_levels;
    allData(iY).yLabel  = fileLabels{iY};

    allData(iY).touch_upstream   = D.touch_upstream_mat;
    allData(iY).touch_downstream = D.touch_downstream_mat;
    allData(iY).touch_upper      = D.touch_upper_mat;
    allData(iY).touch_lower      = D.touch_lower_mat;

    % Sort by correct independent variable
    if strcmp(corrMode, 'point')
        [~, sortIdx] = sort(D.yr_vec);
    else
        [~, sortIdx] = sort(D.xr_vec);
    end

    allData(iY).xr     = D.xr_vec(sortIdx);
    allData(iY).yr     = D.yr_vec(sortIdx);
    allData(iY).Lx     = D.Lx_mat(sortIdx,:);
    allData(iY).Ly     = D.Ly_mat(sortIdx,:);
    allData(iY).Lxu    = D.Lxu_mat(sortIdx,:);
    allData(iY).Lxd    = D.Lxd_mat(sortIdx,:);
    allData(iY).Lyt    = D.Lyt_mat(sortIdx,:);
    allData(iY).Lyb    = D.Lyb_mat(sortIdx,:);
    allData(iY).aspect = D.aspect_mat(sortIdx,:);
    allData(iY).touch_upstream   = D.touch_upstream_mat(sortIdx,:);
    allData(iY).touch_downstream = D.touch_downstream_mat(sortIdx,:);
    allData(iY).touch_upper      = D.touch_upper_mat(sortIdx,:);
    allData(iY).touch_lower      = D.touch_lower_mat(sortIdx,:);

    if D.compute_theta && isfield(D, 'theta_mat')
        theta = D.theta_mat(sortIdx,:);
        theta(theta < 0)  = theta(theta < 0)  + 180;
        theta(theta > 90) = 180 - theta(theta > 90);
        allData(iY).theta = theta;
    else
        allData(iY).theta = [];
    end

    fprintf('Loaded: %s  (%d points)\n', fileLabels{iY}, ...
        numel(D.(sprintf('%s_vec', xaxisField(1:2)))));
end

rho_levels = allData(find([allData.loaded], 1)).rho;
nLevels    = numel(rho_levels);

%% ===================== COLOUR SCHEMES =====================

cmapLines = zeros(nLines, 3);
for i = 1:nLines
    t = (i-1) / max(nLines-1, 1);
    cmapLines(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
end

cmapRho = [0.85 0.33 0.10;
           0.47 0.67 0.19;
           0.13 0.47 0.71];

rhoLabels = arrayfun(@(r) sprintf('\\rho = %.2f', r), rho_levels, ...
    'UniformOutput', false);

%% ===================== HELPER: FILTER OUTLIERS =====================

    function y_out = filterVec(y_in, win)
        y_out = y_in;
        valid = isfinite(y_in);
        if sum(valid) < 3, return; end
        isOut = false(size(y_in));
        isOut(valid) = isoutlier(y_in(valid), 'movmedian', win);
        y_out(isOut) = NaN;
    end

%% ===================== FIGURES 1-7: OVERVIEW AT 1/e =====================

iL_main = 2;

descriptors = {'Lx','Ly','aspect','Lxu','Lxd','Lyt','Lyb'};

if normalise_by_delta && strcmp(corrMode, 'point')
    ylabels_str = {'L_x / \delta', 'L_y / \delta', 'L_x / L_y', ...
                   'L_x^u / \delta', 'L_x^d / \delta', ...
                   'L_y^t / \delta', 'L_y^b / \delta'};
else
    ylabels_str = {'L_x  (mm)', 'L_y  (mm)', 'L_x / L_y', ...
                   'L_x^u  (mm)', 'L_x^d  (mm)', ...
                   'L_y^t  (mm)', 'L_y^b  (mm)'};
end

titles_str = {'Total streamwise extent L_x', ...
              'Total wall-normal extent L_y', ...
              'Aspect ratio L_x / L_y', ...
              'Upstream extent L_x^u', ...
              'Downstream extent L_x^d', ...
              'Wall-normal top extent L_y^t', ...
              'Wall-normal bottom extent L_y^b'};

isLengthScale = [true, true, false, true, true, true, true];

for iD = 1:numel(descriptors)

    figure('Color','w', 'Position',[100 80 700 550], ...
        'Name', titles_str{iD});
    ax = axes; hold on; box on;

    for iY = 1:nLines
        if ~allData(iY).loaded, continue; end

        yr_raw = allData(iY).(xaxisField);
        if normalise_by_delta && strcmp(corrMode, 'point') ...
                && isfinite(delta_vals(iY))
            xaxis = yr_raw / delta_vals(iY);
        else
            xaxis = yr_raw;
        end

        dat = filterVec(allData(iY).(descriptors{iD})(:, iL_main), outlierWindow);
        if normalise_by_delta && strcmp(corrMode, 'point') ...
                && isLengthScale(iD) && isfinite(delta_vals(iY))
            dat = dat / delta_vals(iY);
        end

        plot(ax, xaxis, dat, 'o-', ...
            'Color',           cmapLines(iY,:), ...
            'MarkerFaceColor', cmapLines(iY,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize',      5, ...
            'LineWidth',       1.5, ...
            'DisplayName', sprintf('%s = %.0f mm', lineName, lineVals(iY)));
    end

    if strcmp(corrMode, 'point')
        set(ax, 'XScale', 'log');
    end

    xlabel(ax, xaxisLabel_norm);
    ylabel(ax, ylabels_str{iD});
    title(ax,  sprintf('%s  [\\rho = 1/e]', titles_str{iD}));
    legend(ax, 'Location', 'best');
    grid(ax,   'on');
    if ~isempty(x_lim), xlim(ax, x_lim); end

end

%% ===================== FIGURES 8-10: THRESHOLD SENSITIVITY =====================

sens_descriptors = {'Lx', 'Ly', 'aspect'};
if normalise_by_delta && strcmp(corrMode, 'point')
    sens_ylabels = {'L_x / \delta', 'L_y / \delta', 'L_x / L_y'};
else
    sens_ylabels = {'L_x  (mm)', 'L_y  (mm)', 'L_x / L_y'};
end
sens_titles  = {'Total streamwise extent — threshold sensitivity', ...
                'Total wall-normal extent — threshold sensitivity', ...
                'Aspect ratio — threshold sensitivity'};
sens_isLength = [true, true, false];

for iD = 1:numel(sens_descriptors)

    figure('Color','w', ...
        'Position', [150 50 700 250*nLines], ...
        'Name', sens_titles{iD});

    axAll = gobjects(nLines, 1);

    for iY = 1:nLines

        axAll(iY) = subplot(nLines, 1, iY);
        hold on; box on;

        if ~allData(iY).loaded, continue; end

        yr_raw = allData(iY).(xaxisField);
        if normalise_by_delta && strcmp(corrMode, 'point') ...
                && isfinite(delta_vals(iY))
            xaxis = yr_raw / delta_vals(iY);
        else
            xaxis = yr_raw;
        end

        for iL = 1:nLevels
            dat = filterVec(allData(iY).(sens_descriptors{iD})(:, iL), outlierWindow);
            if normalise_by_delta && strcmp(corrMode, 'point') ...
                    && sens_isLength(iD) && isfinite(delta_vals(iY))
                dat = dat / delta_vals(iY);
            end

            plot(axAll(iY), xaxis, dat, 'o-', ...
                'Color',           cmapRho(iL,:), ...
                'MarkerFaceColor', cmapRho(iL,:), ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize',      4, ...
                'LineWidth',       1.2, ...
                'DisplayName', rhoLabels{iL});
        end

        if strcmp(corrMode, 'point')
            set(axAll(iY), 'XScale', 'log');
        end

        ylabel(axAll(iY), sens_ylabels{iD});
        title(axAll(iY), ...
            sprintf('%s = %.0f mm', lineName, lineVals(iY)), ...
            'FontWeight', 'normal');
        grid(axAll(iY), 'on');

        if iY == 1
            legend(axAll(iY), 'Location', 'best');
        end
        if iY == nLines
            xlabel(axAll(iY), xaxisLabel_norm);
        end
        if ~isempty(x_lim)
            xlim(axAll(iY), x_lim);
        end

    end

    linkaxes(axAll, 'xy');

end

%% ===================== FIGURE 11: Ly vs Lx SCATTER =====================

Lx_all = [];
Ly_all = [];
for iY = 1:nLines
    if ~allData(iY).loaded, continue; end
    Lx_tmp = filterVec(allData(iY).Lx(:, iL_main), outlierWindow);
    Ly_tmp = filterVec(allData(iY).Ly(:, iL_main), outlierWindow);
    Lx_all = [Lx_all; Lx_tmp(isfinite(Lx_tmp))];
    Ly_all = [Ly_all; Ly_tmp(isfinite(Ly_tmp))];
end
Lx_lim = [min(Lx_all) max(Lx_all)];
Ly_lim = [min(Ly_all) max(Ly_all)];

figure('Color','w', 'Position',[100 80 1400 900], ...
    'Name', 'L_y vs L_x — bounding box scatter');

for iY = 1:nLines

    ax = subplot(2, ceil(nLines/2), iY);
    hold on; box on;

    if ~allData(iY).loaded, continue; end

    xaxis = allData(iY).(xaxisField);
    Lx    = filterVec(allData(iY).Lx(:, iL_main), outlierWindow);
    Ly    = filterVec(allData(iY).Ly(:, iL_main), outlierWindow);

    valid = isfinite(Lx) & isfinite(Ly);
    scatter(ax, Lx(valid), Ly(valid), 18, xaxis(valid), 'filled', ...
        'MarkerFaceAlpha', 0.6);

    plot(ax, Lx_lim, Lx_lim, 'k--', 'LineWidth', 1.0, ...
        'DisplayName', 'L_x/L_y = 1');
    for ratio = [2, 5, 10]
        plot(ax, Lx_lim, Lx_lim/ratio, '-', ...
            'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, ...
            'DisplayName', sprintf('L_x/L_y = %d', ratio));
    end

    cb = colorbar(ax);
    cb.Label.String = xaxisLabel;
    colormap(ax, turbo);

    if iY == 1
        legend(ax, 'Location', 'northwest', 'FontSize', 8);
    end

    xlim(ax, Lx_lim);
    ylim(ax, Ly_lim);
    xlabel(ax, 'L_x  (mm)');
    ylabel(ax, 'L_y  (mm)');
    title(ax, sprintf('%s = %.1f mm', lineName, lineVals(iY)), ...
        'FontWeight', 'normal');
    grid(ax, 'on');

end

sgtitle(sprintf('L_y vs L_x  [\\rho = 1/e]  —  coloured by %s', xaxisLabel), ...
    'FontSize', 13);

%% ===================== FIGURE 12: THETA (point mode only) =====================

if strcmp(corrMode, 'point')

    hasTheta = any(cellfun(@(d) ~isempty(d.theta), ...
        num2cell(allData(find([allData.loaded])))));

    if hasTheta
        figure('Color','w', 'Position',[100 80 1300 500], ...
            'Name', 'Inclination angle \theta vs y_{ref}');
        ax = axes; hold on; box on;

        for iY = 1:nLines
            if ~allData(iY).loaded || isempty(allData(iY).theta)
                continue;
            end

            yr_raw = allData(iY).(xaxisField);
            if normalise_by_delta && strcmp(corrMode, 'point') ...
                    && isfinite(delta_vals(iY))
                xaxis = yr_raw / delta_vals(iY);
            else
                xaxis = yr_raw;
            end

            dat = filterVec(allData(iY).theta(:, iL_main), outlierWindow);

            plot(ax, xaxis, dat, 'o-', ...
                'Color',           cmapLines(iY,:), ...
                'MarkerFaceColor', cmapLines(iY,:), ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize',      5, ...
                'LineWidth',       1.5, ...
                'DisplayName', sprintf('%s = %.1f mm', lineName, lineVals(iY)));
        end

        yline(ax, 13, 'k--', 'LineWidth', 1.0, ...
            'Label', '13° (ZPG log region)', ...
            'LabelHorizontalAlignment', 'left', 'DisplayName', 'Marusic Heuer 2007');
        yline(ax, 16, 'k--', 'LineWidth', 1.0, ...
            'Label', '16° (ZPG log region)', ...
            'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

        xlabel(ax, xaxisLabel_norm);
        ylabel(ax, '\theta  (deg)');
        set(ax, 'XScale', 'log');
        title(ax, 'Structural inclination angle \theta  [\rho = 1/e]');
        legend(ax, 'Location', 'best');
        grid(ax, 'on');
        if ~isempty(x_lim), xlim(ax, x_lim); end
    end

end
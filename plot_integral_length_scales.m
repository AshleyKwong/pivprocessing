% =========================================================================
% plot_integral_length_scales.m
%
% Description:
%   Loads and visualises integral length scales (L_uu) and Taylor
%   microscales (lambda_T) computed from two-point spatial correlation
%   planes at multiple reference positions.
%
% Two modes:
%   'sweep' — fixed yref, varying xref → x axis = xref, lines = yref
%   'point' — fixed xref, varying yref → x axis = yref/delta, lines = xref
%             log scale, normalised by delta99
%
%   For each folder, produces a 5-panel diagnostic figure:
%     (1) Normalised correlation decay curves with threshold levels
%     (2) Integral length scale L_uu (raw + outlier-filtered)
%     (3) Taylor microscale lambda_T (raw + outlier-filtered)
%     (4) Threshold-based length scales vs L_uu
%     (5) Ratio L_uu / L_threshold at targetRho
%
%   Final summary figure overlays filtered L_uu and lambda_T across
%   all folders, coloured from lightest (first) to darkest (last).
%
% Inputs (set by user):
%   corrMode        - 'sweep' or 'point'
%   baseDir         - Root directory containing result subfolders
%   folders         - String array of subfolder names to process
%   lineVals        - Numerical values for each folder (yref or xref in mm)
%   targetRho       - Threshold level used for ratio panel (panel 5)
%   normalise_by_delta - true/false (point mode only)
%   blSweepFile     - Path to blSweep mat file (point mode only)
%   outlierWindow   - Window for moving-median outlier detection (sweep only)
%   useSmoothCurves - Toggle smoothed vs raw decay curves
%
% Outputs:
%   - One 5-panel diagnostic figure per folder
%   - One combined overlay figure (L_uu and lambda_T)
%
% Dependencies:
%   integral_length_summary_dx.mat  (from compute_integral_length_scales_hpc.m)
%   decay_curves_dx.mat             (from compute_integral_length_scales_hpc.m)
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

corrMode = 'point';   % 'sweep' or 'point'
% Toggle Taylor microscale panels on/off
showMicroscale = false;

baseDir = 'G:\SW 500mm Mean Flow Fields\Pos_3\two_point_covariance_20260424_231444\';

if strcmp(corrMode, 'sweep')

    folders = [
        "sweep_x_yref2.5"
        "sweep_x_yref6.0"
        "sweep_x_yref14.0"
        "sweep_x_yref89.1"
    ];
    folderTitles = [
        "y_{ref} = 2.5 mm"
        "y_{ref} = 6.0 mm"
        "y_{ref} = 14.0 mm"
        "y_{ref} = 89.1 mm"
    ];
    lineVals   = [2.5, 6.0, 14.0, 89.1];
    lineName   = 'y_{ref}';
    xaxisLabel = 'x_{ref}  (mm)';
    xaxisField = 'xr_exact';

else   % point

    folders = [
    "R_uu_xref7320"
    "R_uu_xref8000"
    % "R_uu_xref711"
    % "R_uu_xref1000"
    ];
    folderTitles = [
        "x_{ref} = 7320 mm"
        "x_{ref} = 8000 mm"
        % "x_{ref} = 711 mm"
        % "x_{ref} = 1000 mm"
    ];
    lineVals   = [7320 8000];
    lineName   = 'x_{ref}';
    xaxisLabel = 'y/\delta';
    xaxisField = 'yr_exact';

end

% Threshold for ratio panel (panel 5)
targetRho = 1/exp(1);

% Normalise by delta99 (point mode only)
normalise_by_delta = true;

blSweepFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\blSweep_20260424_215811_Pos3.mat'];

% Outlier filtering — sweep mode only
outlierWindow   = 100;
useSmoothCurves = true;

%% ===================== LOAD delta99 =====================

if strcmp(corrMode, 'point') && normalise_by_delta

    B      = load(blSweepFile, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;

    valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
    bl_x_clean   = bl_x(valid_bl);
    bl_d99_clean = bl_d99(valid_bl);

    delta_vals = interp1(bl_x_clean, bl_d99_clean, lineVals, ...
        'linear', NaN);

    fprintf('delta99 at each xref:\n');
    for i = 1:numel(lineVals)
        fprintf('  x = %.1f mm  ->  delta99 = %.2f mm\n', ...
            lineVals(i), delta_vals(i));
    end

else
    delta_vals = ones(1, numel(lineVals));
    normalise_by_delta = false;
end

%% ===================== Y AXIS LABELS =====================

if strcmp(corrMode, 'point') && normalise_by_delta
    Lscale_ylabel = 'L_{uu} / \delta';
    lam_ylabel    = '\lambda_T / \delta';
    thr_ylabel    = 'Length scale / \delta';
    ratio_ylabel  = 'L_{uu} / L_{\rho}  (normalised)';
else
    Lscale_ylabel = 'L_{uu}  (mm)';
    lam_ylabel    = '\lambda_T  (mm)';
    thr_ylabel    = 'Length scale  (mm)';
    ratio_ylabel  = 'L_{uu} / L_{\rho}';
end

%% ===================== COLOUR MAP =====================

nFolders  = numel(folders);
cmapLines = zeros(nFolders, 3);
for i = 1:nFolders
    t = (i-1) / max(nFolders-1, 1);
    cmapLines(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
end

cInt     = [0.00 0.00 0.00];
cInt_raw = [0.70 0.70 0.70];
cLam     = [0.10 0.45 0.85];
cLam_raw = 0.65*cLam + 0.35*[1 1 1];

%% ===================== LOAD & STORE =====================

all_xaxis    = cell(nFolders, 1);
all_L_raw    = cell(nFolders, 1);
all_L_filt   = cell(nFolders, 1);
all_lam_raw  = cell(nFolders, 1);
all_lam_filt = cell(nFolders, 1);
all_W        = cell(nFolders, 1);
all_thr      = cell(nFolders, 1);
all_ratio    = cell(nFolders, 1);
all_rhoChosen= zeros(nFolders, 1);
all_curves_dx  = cell(nFolders, 1);
all_curves_rho = cell(nFolders, 1);
all_xzero    = cell(nFolders, 1);
all_isOut_L  = cell(nFolders, 1);
all_isOut_lam= cell(nFolders, 1);

for iF = 1:nFolders

    thisDir     = fullfile(baseDir, folders(iF));
    summaryFile = fullfile(thisDir, 'integral_length_summary_dx.mat');
    curvesFile  = fullfile(thisDir, 'decay_curves_dx.mat');

    if ~isfile(summaryFile) || ~isfile(curvesFile)
        warning('Skipping %s (missing files)', thisDir);
        continue;
    end

    S = load(summaryFile);
    D = load(curvesFile);

    tbl    = S.results.table;
    params = S.results.params;

    % Independent variable
    indepVar = tbl.(xaxisField);

    % Normalise x axis
    if strcmp(corrMode, 'point') && normalise_by_delta ...
            && isfinite(delta_vals(iF))
        xaxis = indepVar / delta_vals(iF);
    else
        xaxis = indepVar;
    end

    % Sort ascending
    [xaxis, sortIdx] = sort(xaxis);

    L_int    = tbl.L_int(sortIdx);
    lambda_T = tbl.lambda_T(sortIdx);
    x_zero   = tbl.zero_crossing(sortIdx);

    % Normalise length scales
    if strcmp(corrMode, 'point') && normalise_by_delta ...
            && isfinite(delta_vals(iF))
        L_int    = L_int    / delta_vals(iF);
        lambda_T = lambda_T / delta_vals(iF);
        x_zero   = x_zero   / delta_vals(iF);
    end

    % Threshold scales
    thr = params.thresholds(:).';
    W   = nan(height(tbl), numel(thr));
    for j = 1:numel(thr)
        varName = matlab.lang.makeValidName(sprintf('threshold_%g', thr(j)));
        if ismember(varName, tbl.Properties.VariableNames)
            col = tbl.(varName)(sortIdx);
            if strcmp(corrMode, 'point') && normalise_by_delta ...
                    && isfinite(delta_vals(iF))
                col = col / delta_vals(iF);
            end
            W(:,j) = col;
        end
    end

    % Target rho column for ratio panel
    [~, idxRho]  = min(abs(thr - targetRho));
    rhoChosen    = thr(idxRho);
    W_rho        = W(:, idxRho);

    % Outlier filtering — sweep only
    if strcmp(corrMode, 'point')
        L_filt   = L_int;
        lam_filt = lambda_T;
        isOut_L   = false(size(L_int));
        isOut_lam = false(size(lambda_T));
    else
        validL   = isfinite(xaxis) & isfinite(L_int);
        validLam = isfinite(xaxis) & isfinite(lambda_T);
        isOut_L   = false(size(L_int));
        isOut_lam = false(size(lambda_T));
        isOut_L(validL)    = isoutlier(L_int(validL),     'movmedian', outlierWindow);
        isOut_lam(validLam)= isoutlier(lambda_T(validLam),'movmedian', outlierWindow);
        L_filt   = L_int;    L_filt(isOut_L)    = NaN;
        lam_filt = lambda_T; lam_filt(isOut_lam)= NaN;
    end

    % Ratio L_uu / W_rho
    ratio_L_thr = L_filt ./ W_rho;
    ratio_L_thr(~isfinite(W_rho)) = NaN;

    % Decay curves
    if useSmoothCurves
        dx_cells  = D.curve_coord_smooth;
        rho_cells = D.curve_rho_smooth;
    else
        dx_cells  = D.curve_coord_raw;
        rho_cells = D.curve_rho_raw;
    end
    dx_cells  = dx_cells(sortIdx);
    rho_cells = rho_cells(sortIdx);

    % Store
    all_xaxis{iF}     = xaxis;
    all_L_raw{iF}     = L_int;
    all_L_filt{iF}    = L_filt;
    all_lam_raw{iF}   = lambda_T;
    all_lam_filt{iF}  = lam_filt;
    all_W{iF}         = W;
    all_thr{iF}       = thr;
    all_ratio{iF}     = ratio_L_thr;
    all_rhoChosen(iF) = rhoChosen;
    all_curves_dx{iF} = dx_cells;
    all_curves_rho{iF}= rho_cells;
    all_xzero{iF}     = x_zero;
    all_isOut_L{iF}   = isOut_L;
    all_isOut_lam{iF} = isOut_lam;

    fprintf('Loaded: %s  (%d points)\n', folders(iF), numel(xaxis));
end

%% ===================== 5-PANEL DIAGNOSTIC PER FOLDER =====================

for iF = 1:nFolders

    if isempty(all_xaxis{iF}), continue; end

    xaxis    = all_xaxis{iF};
    L_int    = all_L_raw{iF};
    L_filt   = all_L_filt{iF};
    lam_raw  = all_lam_raw{iF};
    lam_filt = all_lam_filt{iF};
    W        = all_W{iF};
    thr      = all_thr{iF};
    ratio    = all_ratio{iF};
    rhoChosen= all_rhoChosen(iF);
    x_zero   = all_xzero{iF};
    isOut_L  = all_isOut_L{iF};
    isOut_lam= all_isOut_lam{iF};
    dx_cells  = all_curves_dx{iF};
    rho_cells = all_curves_rho{iF};
    Nref      = numel(dx_cells);
    nPanels = 4 + showMicroscale;   % 5 with microscale, 4 without

    figure('Color','w','Position',[100 80 1400 220*nPanels], ...
        'Name', sprintf('Length scales — %s', folders(iF)));

    % --- Panel 1: Decay curves ---
    ax1 = subplot(nPanels, 1, 1);
    hold(ax1,'on'); box(ax1,'on');
    cmap = turbo(max(Nref,1));
    for k = 1:Nref
        dxk  = dx_cells{k};
        rhok = rho_cells{k};
        if isempty(dxk) || isempty(rhok), continue; end
        plot(ax1, dxk, rhok, 'Color', cmap(k,:), 'LineWidth', 0.8);
        if isfinite(x_zero(k)) && x_zero(k) > 0
            xz_mm = x_zero(k);
            if strcmp(corrMode,'point') && normalise_by_delta ...
                    && isfinite(delta_vals(iF))
                xz_mm = xz_mm * delta_vals(iF);
            end
            plot(ax1, xz_mm, 0, 'o', ...
                'MarkerFaceColor', cmap(k,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 3);
        end
    end
    yline(ax1, 0, 'k-', 'LineWidth', 1);
    for j = 1:numel(thr)
        yline(ax1, thr(j), '--', 'LineWidth', 1, ...
            'Label', sprintf('\\rho = %.3f', thr(j)), ...
            'LabelHorizontalAlignment', 'left');
    end
    xlabel(ax1, '\Deltax  (mm)');
    ylabel(ax1, '\rho_{uu}(\Deltax, \Deltay=0)');
    title(ax1, sprintf('Decay curves — %s', folderTitles(iF)));
    grid(ax1, 'on');
    set(ax1, 'Layer', 'top');

    % --- Panel 2: L_int ---
    ax2 = subplot(nPanels, 1, 2);
    hold(ax2,'on'); box(ax2,'on');

    if strcmp(corrMode, 'sweep')
        % Sweep: show raw, filtered, and rejected separately
        plot(ax2, xaxis, L_int, 'o-', ...
            'Color', cInt_raw, 'MarkerFaceColor', cInt_raw, ...
            'MarkerEdgeColor', cInt_raw, 'LineWidth', 1.0, ...
            'MarkerSize', 5, 'DisplayName', 'Raw');
    end

    
    plot(ax2, xaxis, L_filt, 'o-', ...
        'Color', cInt, 'MarkerFaceColor', cInt, ...
        'MarkerEdgeColor', cInt, 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'L_{uu} \rho =0');
    if strcmp(corrMode, 'sweep') && any(isOut_L)
        plot(ax2, xaxis(isOut_L), L_int(isOut_L), 'rx', ...
            'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Rejected');
    end
    if strcmp(corrMode, 'point'), set(ax2, 'XScale', 'log'); end
    xlabel(ax2, xaxisLabel);
    ylabel(ax2, Lscale_ylabel);
    title(ax2, 'Integral length scale L_{uu}');
    legend(ax2, 'Location', 'best');
    grid(ax2, 'on');

    % --- Panel 3: lambda_T ---
    if showMicroscale
    ax3 = subplot(nPanels, 1, 3);
    hold(ax3,'on'); box(ax3,'on');
    plot(ax3, xaxis, lam_raw, 's-', ...
        'Color', cLam_raw, 'MarkerFaceColor', cLam_raw, ...
        'MarkerEdgeColor', cLam_raw, 'LineWidth', 1.0, ...
        'MarkerSize', 5, 'DisplayName', 'Raw');
    plot(ax3, xaxis, lam_filt, 's-', ...
        'Color', cLam, 'MarkerFaceColor', cLam, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'L_{uu} \rho =0');
    if any(isOut_lam)
        plot(ax3, xaxis(isOut_lam), lam_raw(isOut_lam), 'rx', ...
            'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Rejected');
    end
    if strcmp(corrMode, 'point'), set(ax3, 'XScale', 'log'); end
    xlabel(ax3, xaxisLabel);
    ylabel(ax3, lam_ylabel);
    title(ax3, 'Taylor microscale \lambda_T');
    legend(ax3, 'Location', 'best');
    grid(ax3, 'on');
    end

    % --- Panel 4: Threshold comparison ---
    ax4 = subplot(nPanels, 1, 3 + showMicroscale);
    hold(ax4,'on'); box(ax4,'on');
    colors = lines(size(W,2));
    for j = 1:size(W,2)
        plot(ax4, xaxis, W(:,j), 'o-', ...
            'LineWidth', 1.0, 'Color', colors(j,:), ...
            'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 5, ...
            'DisplayName', sprintf('\\rho = %.3f', thr(j)));
    end
    if strcmp(corrMode, 'sweep')
        plot(ax4, xaxis, L_int, 'o-', ...
            'Color', cInt_raw, 'MarkerFaceColor', cInt_raw, ...
            'MarkerEdgeColor', cInt_raw, 'LineWidth', 1.0, ...
            'MarkerSize', 4, 'DisplayName', 'L_{uu} raw');
    end
    plot(ax4, xaxis, L_filt, 'ko-', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 2.0, 'MarkerSize', 5, ...
        'DisplayName', 'L_{uu} \rho =0');
    if strcmp(corrMode, 'point'), set(ax4, 'XScale', 'log'); end
    xlabel(ax4, xaxisLabel);
    ylabel(ax4, thr_ylabel);
    title(ax4, 'Threshold widths vs L_{uu}');
    legend(ax4, 'Location', 'best');
    grid(ax4, 'on');

    % --- Panel 5: Ratio L_uu / W_rho ---
    ax5 = subplot(nPanels, 1, 4 + showMicroscale);
    hold(ax5,'on'); box(ax5,'on');
    plot(ax5, xaxis, ratio, 'o-', ...
        'Color',           [0.2 0.2 0.2], ...
        'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5, 'MarkerSize', 5);
    yline(ax5, 1.0, 'r--', 'LineWidth', 1.2, ...
        'Label', 'ratio = 1', ...
        'LabelHorizontalAlignment', 'left');
    if strcmp(corrMode, 'point'), set(ax5, 'XScale', 'log'); end
    xlabel(ax5, xaxisLabel);
    ylabel(ax5, sprintf('L_{uu} / L_{\\rho=%.3f}', rhoChosen));
    title(ax5, sprintf('Ratio L_{uu} to threshold width at \\rho = %.3f', rhoChosen));
    grid(ax5, 'on');

    if showMicroscale
        linkaxes([ax2 ax3 ax4 ax5], 'x');
    else
        linkaxes([ax2 ax4 ax5], 'x');
    end
end

%% ===================== OVERLAY FIGURE =====================

figure('Color','w','Position',[200 80 900 700], ...
    'Name', 'Length scales — overlay all folders');
if showMicroscale
    axL   = subplot(2,1,1); hold(axL,'on');   box(axL,'on');
    axLam = subplot(2,1,2); hold(axLam,'on'); box(axLam,'on');
else
    axL   = axes; hold(axL,'on'); box(axL,'on');
end

for iF = 1:nFolders
    if isempty(all_xaxis{iF}), continue; end
    col   = cmapLines(iF,:);
    label = sprintf('%s = %.0f mm', lineName, lineVals(iF));

    plot(axL, all_xaxis{iF}, all_L_filt{iF}, 'o-', ...
        'Color', col, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 5, ...
        'LineWidth', 1.8, 'DisplayName', label);

    if showMicroscale
        plot(axLam, all_xaxis{iF}, all_lam_filt{iF}, 's-', ...
            'Color', col, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 5, ...
            'LineWidth', 1.8, 'DisplayName', label);
    end
end

if strcmp(corrMode, 'point')
    set(axL, 'XScale', 'log');
    if showMicroscale, set(axLam, 'XScale', 'log'); end
end

xlabel(axL, xaxisLabel); ylabel(axL, Lscale_ylabel);
title(axL, sprintf('L_{uu} — all %s', lineName));
legend(axL, 'Location', 'best'); grid(axL, 'on');

if showMicroscale
    xlabel(axLam, xaxisLabel); ylabel(axLam, lam_ylabel);
    title(axLam, sprintf('\\lambda_T — all %s', lineName));
    legend(axLam, 'Location', 'best'); grid(axLam, 'on');
    linkaxes([axL axLam], 'x');
end
% axL   = subplot(2,1,1); hold(axL,'on');   box(axL,'on');
% axLam = subplot(2,1,2); hold(axLam,'on'); box(axLam,'on');
% 
% for iF = 1:nFolders
%     if isempty(all_xaxis{iF}), continue; end
% 
%     col   = cmapLines(iF,:);
%     label = sprintf('%s = %.0f mm', lineName, lineVals(iF));
% 
%     plot(axL, all_xaxis{iF}, all_L_filt{iF}, 'o-', ...
%         'Color', col, 'MarkerFaceColor', col, ...
%         'MarkerEdgeColor', 'k', 'MarkerSize', 5, ...
%         'LineWidth', 1.8, 'DisplayName', label);
% 
%     plot(axLam, all_xaxis{iF}, all_lam_filt{iF}, 's-', ...
%         'Color', col, 'MarkerFaceColor', col, ...
%         'MarkerEdgeColor', 'k', 'MarkerSize', 5, ...
%         'LineWidth', 1.8, 'DisplayName', label);
% end
% 
% if strcmp(corrMode, 'point')
%     set(axL,   'XScale', 'log');
%     set(axLam, 'XScale', 'log');
% end
% 
% xlabel(axL,   xaxisLabel);  ylabel(axL,   Lscale_ylabel);
% title(axL,    sprintf('L_{uu} — all %s', lineName));
% legend(axL,   'Location', 'best');  grid(axL, 'on');
% 
% xlabel(axLam, xaxisLabel);  ylabel(axLam, lam_ylabel);
% title(axLam,  sprintf('\\lambda_T — all %s', lineName));
% legend(axLam, 'Location', 'best');  grid(axLam, 'on');
% 
% linkaxes([axL axLam], 'x');
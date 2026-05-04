% =========================================================================
% run_piv_structure_comparison.m
%
% Compares turbulent structure across all cases at a chosen streamwise
% regime location. Produces the following figures:
%
%   Figure 1 — R_uu isocontour overlay
%              n panels (one per y/delta target), cases as coloured lines
%
%   Figure 2 — Integral length scale profiles
%              L_uu/delta vs y/delta, cases as coloured lines
%
%   Figure 3 — Threshold sensitivity (optional)
%              L_uu (solid) vs W_rho (dashed) per case
%
%   Figure 4 — Per-case diagnostics (optional, one figure per case)
%              Decay curves, L_uu profile, threshold widths, ratio
%
%   Figure 5 — Streamwise extent Lx/delta vs y/delta
%   Figure 6 — Wall-normal extent Ly/delta vs y/delta
%   Figure 7 — Aspect ratio Lx/Ly vs y/delta
%   Figure 8 — Inclination angle theta vs y/delta (if theta available)
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

% --- Regime selection ---
% Choose ONE of: 'fpg' | 'cross' | 'apg' | 'zpg'
%   fpg   -> xref = 350  mm
%   cross -> xref = 550  mm
%   apg   -> xref = 711  mm
%   zpg   -> xref = 1000 mm
regime_request = 'zpg';

% --- Correlation threshold (used for BOTH contour and length scale plots) ---
%   1/exp(1)  ~ 0.368  canonical 1/e decay
%   0.5               mid-level, robust to noise
%   0.2               outer contour, large-scale footprint
rho_level = 0.3;

% --- Wall-normal targets (y/delta) shown in contour figure ---
% ydelta_targets = [0.02, 0.04, 0.05, 0.1, 0.5, 0.8];
ydelta_targets = [ 0.05,   0.5, ];
% --- Contour window (mm either side of xref) ---
dx_back = 300;
dx_fwd  = 600;

% --- Optional figures ---
showThresholdSensitivity = false;   % Figure 3: L_uu vs W_rho overlay
showDiagnostics          = false;  % Figure 4: per-case 5-panel inspection
useSmoothCurves          = false;   % diagnostics only: smooth vs raw decay curves

% --- Case definitions ---
% Add/uncomment cases as needed. Each case requires:
%   label    : short name for legend
%   cov_dir  : root two_point_covariance directory for this case
%   blFile   : path to blSweep .mat for this case
%   gridFile : path to grid.mat (cases sharing a grid can use the same file)

cases(1).label   = 'Case 1';
cases(1).cov_dir = 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038';
cases(1).blFile  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case1_PIVresults\blSweep_20260419_193106.mat';
cases(1).gridFile= 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\grid.mat';

cases(2).label   = 'Case 2';
cases(2).cov_dir = 'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244';
cases(2).blFile  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';
cases(2).gridFile= 'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244\grid.mat';

cases(3).label   = 'Case 3';
cases(3).cov_dir = 'D:\FULLPROCESSEDY235AOAN11AOAFN11PIVDATA\two_point_covariance_20260421_083251';
cases(3).blFile  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case6_PIVresults\blSweep_20260421_082649.mat';
cases(3).gridFile= 'D:\FULLPROCESSEDY235AOAN11AOAFN11PIVDATA\two_point_covariance_20260421_083251\grid.mat';

% cases(4).label   = 'Case 4';
% cases(4).cov_dir = 'D:\...';
% cases(4).blFile  = 'C:\...';
% cases(4).gridFile= 'D:\...';

% cases(5).label   = 'Case 5';
% cases(5).cov_dir = 'D:\...';
% cases(5).blFile  = 'C:\...';
% cases(5).gridFile= 'D:\...';

%% ===================== REGIME -> XREF MAP =====================

regimeMap = struct('fpg', 350, 'cross', 550, 'apg', 711, 'zpg', 1000);

if ~isfield(regimeMap, regime_request)
    error('regime_request must be one of: fpg, cross, apg, zpg');
end
xref_nom = regimeMap.(regime_request);

fprintf('Regime: %s  |  xref = %.0f mm  |  rho = %.3f\n\n', ...
    regime_request, xref_nom, rho_level);

%% ===================== COLOUR AND MARKER SCHEME =====================

nCases = numel(cases);

red_dark   = [0.6 0.0 0.0];
red_bright = [1.0 0.4 0.4];
red_grad   = [linspace(red_dark(1), red_bright(1), nCases)', ...
              linspace(red_dark(2), red_bright(2), nCases)', ...
              linspace(red_dark(3), red_bright(3), nCases)'];
% Fixed colours per case (greyscale + blue for case 5)
% caseColours = [
% 
%     0.35  0.35  0.35;   % Case 2 — dark grey
%     0.55  0.55  0.55;   % Case 3 — mid grey
%     0.75  0.75  0.75;   % Case 4 — light grey
%     0.08  0.40  0.74;   % Case 5 — blue
%     0.00  0.00  0.00;   % Case 1 — black
% ];
caseColours = red_grad; 
% Pad with blue if somehow more than 5 cases are defined
if nCases > size(caseColours, 1)
    caseColours(end+1:nCases, :) = repmat([0.08 0.40 0.74], nCases-size(caseColours,1), 1);
end

% One distinct marker per case — redundant encoding alongside colour
% Case 1: diamond  Case 2: square  Case 3: star  Case 4: triangle  Case 5: circle
markerList = {'s','d', 'p', '^', 'o'};
if nCases > numel(markerList)
    markerList(end+1:nCases) = {'o'};
end

%% ===================== LOAD GRID =====================

G      = load(cases(1).gridFile, 'worldX_merged', 'worldY_merged');
worldX = double(G.worldX_merged);
worldY = double(G.worldY_merged);

x_min_domain = min(worldX(:));
x_max_domain = max(worldX(:));

%% ===================== LOAD ALL CASE DATA =====================

% Pre-allocate storage
for iC = 1:nCases
    cases(iC).delta    = NaN;
    cases(iC).xref_dir = '';
    % length scale fields
    cases(iC).ydelta   = [];
    cases(iC).L_filt   = [];
    cases(iC).L_raw    = [];
    cases(iC).W_thr    = [];
    cases(iC).W_all    = [];
    cases(iC).thr      = [];
    cases(iC).ratio    = [];
    cases(iC).rhoChosen= NaN;
    cases(iC).curves_dx = {};
    cases(iC).curves_rho= {};
    cases(iC).xzero    = [];
    % bbox / spatial extent fields (loaded from bbox_analysis subfolder)
    cases(iC).bbox_ydelta = [];
    cases(iC).bbox_Lx     = [];
    cases(iC).bbox_Ly     = [];
    cases(iC).bbox_aspect = [];
    cases(iC).bbox_theta  = [];
    % bbox / spatial extent fields
    cases(iC).bbox_ydelta = [];
    cases(iC).bbox_Lx     = [];
    cases(iC).bbox_Ly     = [];
    cases(iC).bbox_aspect = [];
    cases(iC).bbox_theta  = [];
end

for iC = 1:nCases

    fprintf('Loading %s...\n', cases(iC).label);

    % --- delta99 at xref_nom ---
    B      = load(cases(iC).blFile, 'blSweep');
    bl_x   = B.blSweep.x_mm;
    bl_d99 = B.blSweep.delta99_hybrid_mm;
    valid  = isfinite(bl_x) & isfinite(bl_d99);
    delta  = interp1(bl_x(valid), bl_d99(valid), xref_nom, 'linear', NaN);

    if ~isfinite(delta)
        warning('  delta99 not available at xref=%.0f mm — skipping\n', xref_nom);
        continue;
    end
    cases(iC).delta = delta;
    fprintf('  delta99 = %.2f mm\n', delta);

    % --- Locate xref subfolder ---
    subPattern = sprintf('R_uu_xref%g*', xref_nom);
    d = dir(fullfile(cases(iC).cov_dir, subPattern));
    d = d([d.isdir]);

    if isempty(d)
        warning('  No subfolder matching %s — skipping', subPattern);
        continue;
    end
    cases(iC).xref_dir = fullfile(cases(iC).cov_dir, d(1).name);
    fprintf('  xref dir: %s\n', cases(iC).xref_dir);

    % --- Load integral length scale summary ---
    summaryFile = fullfile(cases(iC).xref_dir, 'integral_length_summary_dx.mat');
    curvesFile  = fullfile(cases(iC).xref_dir, 'decay_curves_dx.mat');

    if ~isfile(summaryFile) || ~isfile(curvesFile)
        warning('  Missing summary or curves file — skipping length scale load');
    else
        S = load(summaryFile);
        D = load(curvesFile);

        tbl    = S.results.table;
        params = S.results.params;

        [ydelta, sortIdx] = sort(tbl.yr_exact / delta);

        L_int    = tbl.L_int(sortIdx)    / delta;
        lambda_T = tbl.lambda_T(sortIdx) / delta;
        x_zero   = tbl.zero_crossing(sortIdx) / delta;

        thr = params.thresholds(:).';
        W   = nan(height(tbl), numel(thr));
        for j = 1:numel(thr)
            varName = matlab.lang.makeValidName(sprintf('threshold_%g', thr(j)));
            if ismember(varName, tbl.Properties.VariableNames)
                W(:,j) = tbl.(varName)(sortIdx) / delta;
            end
        end

        [~, idxRho] = min(abs(thr - rho_level));
        rhoChosen   = thr(idxRho);
        W_rho       = W(:, idxRho);

        ratio = L_int ./ W_rho;
        ratio(~isfinite(W_rho)) = NaN;

        if useSmoothCurves && isfield(D, 'curve_coord_smooth')
            dx_cells  = D.curve_coord_smooth(sortIdx);
            rho_cells = D.curve_rho_smooth(sortIdx);
        else
            dx_cells  = D.curve_coord_raw(sortIdx);
            rho_cells = D.curve_rho_raw(sortIdx);
        end

        cases(iC).ydelta    = ydelta;
        cases(iC).L_filt    = L_int;
        cases(iC).L_raw     = L_int;
        cases(iC).W_thr     = W_rho;
        cases(iC).W_all     = W;
        cases(iC).thr       = thr;
        cases(iC).ratio     = ratio;
        cases(iC).rhoChosen = rhoChosen;
        cases(iC).curves_dx = dx_cells;
        cases(iC).curves_rho= rho_cells;
        cases(iC).xzero     = x_zero;

        fprintf('  Length scales: %d wall-normal points loaded\n', numel(ydelta));
    end

    % --- Load bbox / spatial extent results ---
    % Expects: <cov_dir>/bbox_analysis/bbox_results_xref_<xref_nom>.mat
    bboxFile = fullfile(cases(iC).cov_dir, 'bbox_analysis', ...
        sprintf('bbox_results_xref_%.0f.mat', xref_nom));

    if ~isfile(bboxFile)
        warning('  No bbox file found at %s — skipping spatial extents', bboxFile);
    else
        BD = load(bboxFile);

        % [~, bsortIdx] = sort(BD.yr_vec(:));   % force column — safe matrix row indexing
        % yr_sorted  = BD.yr_vec(bsortIdx);
        % yd_bbox    = yr_sorted / delta;
        % 
        % % Pick the column matching rho_level (closest available threshold)
        % [~, iRhoB] = min(abs(BD.rho_levels - rho_level))
        % 
        % Lx_raw  = BD.Lx_mat(bsortIdx, iRhoB)     / delta;
        % Ly_raw  = BD.Ly_mat(bsortIdx, iRhoB)     / delta;
        % asp_raw = BD.aspect_mat(bsortIdx, iRhoB);   % dimensionless
        yr_col  = BD.yr_vec(:);
        yd_bbox = yr_col / delta;

        [~, iRhoB] = min(abs(BD.rho_levels - rho_level));

        Lx_raw  = BD.Lx_mat(:, iRhoB)     / delta;
        Ly_raw  = BD.Ly_mat(:, iRhoB)     / delta;
        asp_raw = BD.aspect_mat(:, iRhoB);

        % Theta — map to [0, 90] following original convention
        % if BD.compute_theta && isfield(BD, 'theta_mat')
        %     th = BD.theta_mat(bsortIdx, iRhoB);
        %     th(th < 0)  = th(th < 0)  + 180;
        %     th(th > 90) = 180 - th(th > 90);
        % else
        %     th = nan(size(Lx_raw));
        % end
        if BD.compute_theta && isfield(BD, 'theta_mat')
            th = BD.theta_mat(:, iRhoB);
            th(th < 0)  = th(th < 0)  + 180;
            th(th > 90) = 180 - th(th > 90);
        else
            th = nan(size(Lx_raw));
        end

        cases(iC).bbox_ydelta = yd_bbox;
        cases(iC).bbox_Lx     = Lx_raw;
        cases(iC).bbox_Ly     = Ly_raw;
        cases(iC).bbox_aspect = asp_raw;
        cases(iC).bbox_theta  = th;

        fprintf('  Bbox extents: %d wall-normal points loaded\n', numel(yd_bbox));
    end
fprintf('\n');
end 
%% ===================== FIGURE 1: R_uu CONTOUR OVERLAY =====================
% n panels (one per y/delta), cases as coloured lines

nTargets = numel(ydelta_targets);

figure('Color','w', ...
    'Position', [50 50 750 230*nTargets], ...
    'Name', sprintf('Fig 1 — R_uu contours  |  %s  rho=%.3f', ...
    upper(regime_request), rho_level));

axC      = gobjects(nTargets, 1);
leg_h    = gobjects(nCases, 1);
leg_lbl  = cell(nCases, 1);

for iT = 1:nTargets

    yd_target = ydelta_targets(iT);
    axC(iT)   = subplot(nTargets, 1, iT);
    hold on; box on; grid on;

    for iC = 1:nCases

        if isempty(cases(iC).xref_dir) || ~isfinite(cases(iC).delta)
            continue;
        end

        delta     = cases(iC).delta;
        col       = caseColours(iC,:);
        yr_target = yd_target * delta;

        % Find closest file by yref
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
        yr_sub            = yr_available;
        yr_sub(~xr_ok)    = Inf;
        [~, iClosest]     = min(abs(yr_sub - yr_target));

        S  = load(fullfile(cases(iC).xref_dir, files(iClosest).name), 'R_s','xr','yr');
        R  = double(S.R_s);
        xr = double(S.xr);
        yr = double(S.yr);

        fprintf('  Fig1 | %s | y/d target=%.2f | actual=%.3f\n', ...
            cases(iC).label, yd_target, yr/delta);

        % Separation grids
        dX = worldX - xr;
        dY = worldY - yr;

        x_lo = max(xr - dx_back, x_min_domain);
        x_hi = min(xr + dx_fwd,  x_max_domain);

        boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
        R_box           = R;
        R_box(~boxMask) = NaN;
        R_box(R_box >  1.0) = 1.0;
        R_box(R_box < -1.0) = NaN;

        dX_vec = dX(1,:);
        dY_vec = dY(:,1);

        % Contour extraction
        mask = R_box >= rho_level;
        if sum(mask(:)) < 5, continue; end

        dist2           = dX.^2 + dY.^2;
        dist2(~boxMask) = Inf;
        [~, i0]         = min(dist2(:));

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
        bnd_dx = dX_vec(bnd(:,2)) / delta;
        bnd_dy = dY_vec(bnd(:,1)) / delta;

        h = plot(axC(iT), bnd_dx, bnd_dy, '-', ...
            'Color', col, 'LineWidth', 1.8);
        plot(axC(iT), 0, 0, '+', ...
            'Color', col, 'MarkerSize', 8, 'LineWidth', 1.5);

        if iT == 1
            leg_h(iC)   = h;
            leg_lbl{iC} = cases(iC).label;
        end

    end % cases

    xline(axC(iT), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
    yline(axC(iT), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
    ylabel(axC(iT), '\Deltay / \delta');
    title(axC(iT), sprintf('y/\\delta \\approx %.2f', yd_target), ...
        'FontWeight','normal','FontSize',10);
    if iT == nTargets
        xlabel(axC(iT), '\Deltax / \delta');
    end

end % targets

linkaxes(axC, 'xy');

validH = leg_h(isgraphics(leg_h));
validL = leg_lbl(isgraphics(leg_h));
if ~isempty(validH)
    legend(axC(1), validH, validL, 'Location','best','FontSize',10);
end

sgtitle(sprintf('R_{uu} isocontour  [\\rho = %.3f]  |  %s  (x_{ref} = %.0f mm)', ...
    rho_level, upper(regime_request), xref_nom), 'FontSize', 13);

%% ===================== FIGURE 2: INTEGRAL LENGTH SCALE PROFILES =====================
% y/delta (log) on x-axis, L_uu/delta on y-axis, cases as coloured lines

figure('Color','w','Position',[850 50 700 500], ...
    'Name', sprintf('Fig 2 — L_uu profiles  |  %s', upper(regime_request)));

axL = axes; hold(axL,'on'); box(axL,'on'); grid(axL,'on');

for iC = 1:nCases
    if isempty(cases(iC).ydelta), continue; end
    col = caseColours(iC,:);
    plot(axL, cases(iC).ydelta, cases(iC).L_filt, ...
        [markerList{iC} '-'], ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
        'MarkerSize', 6, 'LineWidth', 1.8, ...
        'DisplayName', cases(iC).label);
end

set(axL,'XScale','log');
xlabel(axL,'y / \delta');
ylabel(axL,'L_{uu} / \delta');
title(axL, sprintf('Integral length scale profiles  |  %s  (x_{ref} = %.0f mm)', ...
    upper(regime_request), xref_nom), 'FontSize', 11);
legend(axL,'Location','best','FontSize',10);

%% ===================== FIGURE 3: THRESHOLD SENSITIVITY (optional) =====================

if showThresholdSensitivity

    figure('Color','w','Position',[850 600 700 500], ...
        'Name', sprintf('Fig 3 — Threshold sensitivity rho=%.3f  |  %s', ...
        rho_level, upper(regime_request)));

    axT = axes; hold(axT,'on'); box(axT,'on'); grid(axT,'on');

    for iC = 1:nCases
        if isempty(cases(iC).ydelta), continue; end
        col = caseColours(iC,:);
        plot(axT, cases(iC).ydelta, cases(iC).L_filt, ...
            [markerList{iC} '-'], ...
            'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
            'MarkerSize', 6, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s  L_{uu}', cases(iC).label));
        plot(axT, cases(iC).ydelta, cases(iC).W_thr, ...
            [markerList{iC} '--'], ...
            'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
            'MarkerSize', 4, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s  W_{\\rho=%.3f}', cases(iC).label, rho_level));
    end

    set(axT,'XScale','log');
    xlabel(axT,'y / \delta');
    ylabel(axT,'Length scale / \delta');
    title(axT, sprintf('L_{uu} (solid) vs W_{\\rho=%.3f} (dashed)  |  %s', ...
        rho_level, upper(regime_request)), 'FontSize', 11);
    legend(axT,'Location','best','FontSize',9);

end

%% ===================== FIGURE 4: PER-CASE DIAGNOSTICS (optional) =====================

if ~showDiagnostics, return; end

for iC = 1:nCases

    if isempty(cases(iC).ydelta), continue; end

    ydelta    = cases(iC).ydelta;
    L_filt    = cases(iC).L_filt;
    W         = cases(iC).W_all;
    thr       = cases(iC).thr;
    ratio     = cases(iC).ratio;
    rhoChosen = cases(iC).rhoChosen;
    dx_cells  = cases(iC).curves_dx;
    rho_cells = cases(iC).curves_rho;
    Nref      = numel(dx_cells);

    figure('Color','w','Position',[100 80 1200 900], ...
        'Name', sprintf('Fig 4 — Diagnostics  |  %s  |  %s', ...
        cases(iC).label, upper(regime_request)));

    ax1 = subplot(4,1,1); hold(ax1,'on'); box(ax1,'on');
    cmap = turbo(max(Nref,1));
    for k = 1:Nref
        dxk  = dx_cells{k}; rhok = rho_cells{k};
        if isempty(dxk), continue; end
        plot(ax1, dxk, rhok, 'Color', cmap(k,:), 'LineWidth', 0.8);
    end
    yline(ax1, 0, 'k-', 'LineWidth', 1);
    for j = 1:numel(thr)
        yline(ax1, thr(j), '--', 'LineWidth', 1, ...
            'Label', sprintf('\\rho=%.3f', thr(j)), ...
            'LabelHorizontalAlignment','left');
    end
    xlabel(ax1,'\Deltax  (mm)'); ylabel(ax1,'\rho_{uu}');
    title(ax1, sprintf('Decay curves — %s  |  %s', cases(iC).label, upper(regime_request)));
    grid(ax1,'on');

    col = caseColours(iC,:);
    mkr = markerList{iC};

    ax2 = subplot(4,1,2); hold(ax2,'on'); box(ax2,'on');
    plot(ax2, ydelta, L_filt, [mkr '-'], ...
        'Color', col, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor','k','MarkerSize',5,'LineWidth',1.8);
    set(ax2,'XScale','log');
    xlabel(ax2,'y / \delta'); ylabel(ax2,'L_{uu} / \delta');
    title(ax2,'Integral length scale profile'); grid(ax2,'on');

    ax3 = subplot(4,1,3); hold(ax3,'on'); box(ax3,'on');
    colors = lines(size(W,2));
    for j = 1:size(W,2)
        plot(ax3, ydelta, W(:,j), 'o-', ...
            'Color', colors(j,:), 'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1.0, ...
            'DisplayName', sprintf('\\rho=%.3f', thr(j)));
    end
    plot(ax3, ydelta, L_filt, [mkr '-'], ...
        'Color', col, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor','k','MarkerSize',5,'LineWidth',2.0,'DisplayName','L_{uu}');
    set(ax3,'XScale','log');
    xlabel(ax3,'y / \delta'); ylabel(ax3,'Length scale / \delta');
    title(ax3,'Threshold widths vs L_{uu}');
    legend(ax3,'Location','best'); grid(ax3,'on');

    ax4 = subplot(4,1,4); hold(ax4,'on'); box(ax4,'on');
    plot(ax4, ydelta, ratio, [mkr '-'], ...
        'Color', col, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor','k','MarkerSize',5,'LineWidth',1.5);
    yline(ax4, 1.0, 'r--','LineWidth',1.2, ...
        'Label','ratio=1','LabelHorizontalAlignment','left');
    set(ax4,'XScale','log');
    xlabel(ax4,'y / \delta');
    ylabel(ax4, sprintf('L_{uu} / W_{\\rho=%.3f}', rhoChosen));
    title(ax4,'Ratio L_{uu} to threshold width'); grid(ax4,'on');

    linkaxes([ax2 ax3 ax4],'x');

end

%% ===================== FIGURE 5: Lx/delta vs y/delta =====================

figure('Color','w','Position',[50 100 700 500], ...
    'Name', sprintf('Fig 5 — L_x/delta  |  %s', upper(regime_request)));

ax5 = axes; hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');

for iC = 1:nCases
    if isempty(cases(iC).bbox_ydelta), continue; end
    col = caseColours(iC,:);
    plot(ax5, cases(iC).bbox_ydelta, cases(iC).bbox_Lx, ...
        [markerList{iC} '-'], ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
        'MarkerSize', 6, 'LineWidth', 1.8, ...
        'DisplayName', cases(iC).label);
end

set(ax5,'XScale','log');
xlabel(ax5,'y / \delta');
ylabel(ax5,'L_x / \delta');
title(ax5, sprintf('Streamwise extent L_x  [\\rho = %.3f]  |  %s  (x_{ref} = %.0f mm)', ...
    rho_level, upper(regime_request), xref_nom), 'FontSize', 11);
legend(ax5,'Location','best','FontSize',10);

%% ===================== FIGURE 6: Ly/delta vs y/delta =====================

figure('Color','w','Position',[800 100 700 500], ...
    'Name', sprintf('Fig 6 — L_y/delta  |  %s', upper(regime_request)));

ax6 = axes; hold(ax6,'on'); box(ax6,'on'); grid(ax6,'on');

for iC = 1:nCases
    if isempty(cases(iC).bbox_ydelta), continue; end
    col = caseColours(iC,:);
    plot(ax6, cases(iC).bbox_ydelta, cases(iC).bbox_Ly, ...
        [markerList{iC} '-'], ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
        'MarkerSize', 6, 'LineWidth', 1.8, ...
        'DisplayName', cases(iC).label);
end

set(ax6,'XScale','log');
xlabel(ax6,'y / \delta');
ylabel(ax6,'L_y / \delta');
title(ax6, sprintf('Wall-normal extent L_y  [\\rho = %.3f]  |  %s  (x_{ref} = %.0f mm)', ...
    rho_level, upper(regime_request), xref_nom), 'FontSize', 11);
legend(ax6,'Location','best','FontSize',10);

%% ===================== FIGURE 7: Lx/Ly (aspect ratio) vs y/delta =====================

figure('Color','w','Position',[50 650 700 500], ...
    'Name', sprintf('Fig 7 — L_x/L_y aspect  |  %s', upper(regime_request)));

ax7 = axes; hold(ax7,'on'); box(ax7,'on'); grid(ax7,'on');

for iC = 1:nCases
    if isempty(cases(iC).bbox_ydelta), continue; end
    col = caseColours(iC,:);
    plot(ax7, cases(iC).bbox_ydelta, cases(iC).bbox_aspect, ...
        [markerList{iC} '-'], ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
        'MarkerSize', 6, 'LineWidth', 1.8, ...
        'DisplayName', cases(iC).label);
end

yline(ax7, 1.0, 'k:', 'LineWidth', 1.0, ...
    'Label', 'isotropic', 'LabelHorizontalAlignment','right', ...
    'HandleVisibility','off');

set(ax7,'XScale','log');
xlabel(ax7,'y / \delta');
ylabel(ax7,'L_x / L_y');
title(ax7, sprintf('Aspect ratio L_x/L_y  [\\rho = %.3f]  |  %s  (x_{ref} = %.0f mm)', ...
    rho_level, upper(regime_request), xref_nom), 'FontSize', 11);
legend(ax7,'Location','best','FontSize',10);

%% ===================== FIGURE 8: theta vs y/delta =====================

hasTheta = any(arrayfun(@(c) ~isempty(c.bbox_theta) && any(isfinite(c.bbox_theta)), cases));

if hasTheta

    figure('Color','w','Position',[800 650 700 500], ...
        'Name', sprintf('Fig 8 — theta  |  %s', upper(regime_request)));

    ax8 = axes; hold(ax8,'on'); box(ax8,'on'); grid(ax8,'on');

    for iC = 1:nCases
        if isempty(cases(iC).bbox_theta), continue; end
        col = caseColours(iC,:);
        plot(ax8, cases(iC).bbox_ydelta, cases(iC).bbox_theta, ...
            [markerList{iC} '-'], ...
            'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
            'MarkerSize', 6, 'LineWidth', 1.8, ...
            'DisplayName', cases(iC).label);
    end

    % Literature reference lines
    yline(ax8, 13, 'k--', 'LineWidth', 1.0, ...
        'Label', '13°', 'LabelHorizontalAlignment','left', ...
        'HandleVisibility','off');
    yline(ax8, 16, 'k--', 'LineWidth', 1.0, ...
        'Label', '16°', 'LabelHorizontalAlignment','left', ...
        'HandleVisibility','off');

    set(ax8,'XScale','log');
    xlabel(ax8,'y / \delta');
    ylabel(ax8,'\theta  (deg)');
    title(ax8, sprintf('Inclination angle \\theta  [\\rho = %.3f]  |  %s  (x_{ref} = %.0f mm)', ...
        rho_level, upper(regime_request), xref_nom), 'FontSize', 11);
    legend(ax8,'Location','best','FontSize',10);

end
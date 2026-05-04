% =========================================================================
% plot_Ruu_yref_comparison.m
%
% For a set of user-specified x_ref positions, loads the R_uu correlation
% map from each y_ref folder whose x_ref is closest to the requested value,
% and overlays the isocontour at a chosen rho level for all y_ref values
% on a single axes — one figure per requested x_ref position.
%
% This lets you directly compare how the correlation structure shape
% changes with wall-normal position at a fixed streamwise location.
%
% Author:  Ashley Kwong 
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

% Full paths to each sweep folder — one per y_ref
sweepFolders = {
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\sweep_x_yref2.0'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260414_120505\sweep_x_yref2.6'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260414_120904\sweep_x_yref4.6'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260414_140922\sweep_x_yref8.0'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260414_141202\sweep_x_yref14.0'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260414_152812\sweep_x_yref24.5'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260414_170430\sweep_x_yref42.7'
    'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\sweep_x_yref89.1'
};

% Numerical y_ref values matching folder order (mm)
y_vals = [2.0, 2.6, 4.6, 8.0, 14.0, 24.5, 42.7 89.1];

% x_ref positions to inspect (mm) — script finds closest available file
x_targets = [350, 711, 1000];

% Correlation threshold for isocontour
rho_level = 1/exp(1);

% Search box (mm) — must match what was used in analyse_Ruu_maps.m
dx_back = 300;   % upstream
dx_fwd  = 600;   % downstream

%% ===================== COLOUR SCHEME =====================

nY    = numel(sweepFolders);
nX    = numel(x_targets);

% y_ref gradient: light blue (near wall) -> dark navy (outer)
cmapY = zeros(nY, 3);
for i = 1:nY
    t = (i-1) / max(nY-1, 1);
    cmapY(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
end

%% ===================== PARSE y_ref FROM FOLDER NAMES =====================

% Extract y_ref from folder name for file search pattern
yref_strs = cell(nY, 1);
for iY = 1:nY
    tok = regexp(sweepFolders{iY}, 'sweep_x_yref([\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        yref_strs{iY} = tok{1};
    else
        yref_strs{iY} = sprintf('%.1f', y_vals(iY));
    end
end

%% ===================== MAIN LOOP OVER x_targets =====================

for iX = 1:nX

    x_target = x_targets(iX);
    fprintf('\n========================================\n');
    fprintf('x_target = %.1f mm\n', x_target);

    figure('Color','w', ...
        'Position', [50 + (iX-1)*30, 80 + (iX-1)*30, 900, 600], ...
        'Name', sprintf('R_{uu} contours at x_{ref} \\approx %.0f mm', x_target));
    ax = axes; hold on; box on;

    x_found = nan(nY, 1);   % actual x_ref used for each y_ref

    for iY = 1:nY

        thisDir = sweepFolders{iY};

        % --- Load grid ---
        gridFile = fullfile(thisDir, 'grid.mat');
        if ~isfile(gridFile)
            warning('No grid.mat in %s — skipping', thisDir);
            continue;
        end
        G      = load(gridFile, 'worldX_merged', 'worldY_merged');
        worldX = double(G.worldX_merged);
        worldY = double(G.worldY_merged);

        x_min_domain = min(worldX(:));
        x_max_domain = max(worldX(:));
        y_min_domain = min(worldY(:));
        y_max_domain = max(worldY(:));

        % --- Find file with closest x_ref ---
        pattern = sprintf('R_uu_xref*_yref%s.mat', yref_strs{iY});
        files   = dir(fullfile(thisDir, pattern));
        if isempty(files)
            warning('No files matching %s in %s — skipping', pattern, thisDir);
            continue;
        end

        % Parse x_ref from each filename
        xr_available = nan(numel(files), 1);
        for f = 1:numel(files)
            tok = regexp(files(f).name, ...
                'R_uu_xref([-+]?\d*\.?\d+)_yref', 'tokens', 'once');
            if ~isempty(tok)
                xr_available(f) = str2double(tok{1});
            end
        end

        % Find closest
        [~, iClosest] = min(abs(xr_available - x_target));
        chosenFile    = fullfile(thisDir, files(iClosest).name);
        xr_actual     = xr_available(iClosest);
        x_found(iY)   = xr_actual;

        fprintf('  y_ref = %5.1f mm | using x_ref = %.2f mm  (%s)\n', ...
            y_vals(iY), xr_actual, files(iClosest).name);

        % --- Load correlation map ---
        S  = load(chosenFile, 'R_s', 'xr', 'yr');
        R  = double(S.R_s);
        xr = double(S.xr);
        yr = double(S.yr);

        % --- Delta grids ---
        dX = worldX - xr;
        dY = worldY - yr;

        % --- Clamped search box ---
        x_lo = max(xr - dx_back, x_min_domain);
        x_hi = min(xr + dx_fwd,  x_max_domain);

        boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
        R_box           = R;
        R_box(~boxMask) = NaN;
        R_box(R_box >  1.0) = 1.0;
        R_box(R_box < -1.0) = NaN;

        dX_box           = dX;  dX_box(~boxMask) = NaN;
        dY_box           = dY;  dY_box(~boxMask) = NaN;

        dx_vec = dX(1,:);
        dy_vec = dY(:,1);

        % --- Extract isocontour of origin-connected region ---
        mask = R_box >= rho_level;
        if sum(mask(:)) < 5
            fprintf('    WARNING: fewer than 5 points above threshold — skipping\n');
            continue;
        end

        % Find origin pixel
        dist2           = dX.^2 + dY.^2;
        dist2(~boxMask) = Inf;
        [~, i0]         = min(dist2(:));

        % Connected components
        CC = bwconncomp(mask);
        if CC.NumObjects == 0, continue; end

        regionIdx = 0;
        for r = 1:CC.NumObjects
            if any(CC.PixelIdxList{r} == i0)
                regionIdx = r;
                break;
            end
        end
        if regionIdx == 0
            regionSizes = cellfun(@numel, CC.PixelIdxList);
            [~, regionIdx] = max(regionSizes);
        end

        mask_origin = false(size(mask));
        mask_origin(CC.PixelIdxList{regionIdx}) = true;

        % Extract boundary of origin region
        B = bwboundaries(mask_origin);

        for b = 1:numel(B)
            bnd    = B{b};
            bnd_dx = dx_vec(bnd(:,2));
            bnd_dy = dy_vec(bnd(:,1));
            plot(ax, bnd_dx, bnd_dy, '-', ...
                'Color',       cmapY(iY,:), ...
                'LineWidth',   2.0, ...
                'DisplayName', sprintf('y_{ref} = %.1f mm', y_vals(iY)));
        end

        % Bounding box
        dX_pts = dX_box(mask_origin);
        dY_pts = dY_box(mask_origin);
        Lxu = abs(min(dX_pts));
        Lxd = max(dX_pts);
        Lyt = max(dY_pts);
        Lyb = abs(min(dY_pts));

        bx = [-Lxu,  Lxd,  Lxd, -Lxu, -Lxu];
        by = [-Lyb, -Lyb,  Lyt,  Lyt, -Lyb];
        plot(ax, bx, by, '--', ...
            'Color',            cmapY(iY,:), ...
            'LineWidth',        0.8, ...
            'HandleVisibility', 'off');

    end % y_ref loop

    % --- Reference point marker and zero lines ---
    plot(ax, 0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    xline(ax, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(ax, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Deduplicate legend entries (one per y_ref)
    h = get(ax, 'Children');
    seen   = false(nY, 1);
    hKeep  = gobjects(0);
    for k = 1:numel(h)
        dn = get(h(k), 'DisplayName');
        for iY = 1:nY
            lbl = sprintf('y_{ref} = %.1f mm', y_vals(iY));
            if strcmp(dn, lbl) && ~seen(iY)
                seen(iY)  = true;
                hKeep(end+1) = h(k); %#ok<AGROW>
            end
        end
    end
    legend(ax, hKeep, 'Location', 'best', 'FontSize', 9);

    % Report actual x_ref values used
    x_found_str = strjoin(arrayfun(@(v) sprintf('%.1f', v), ...
        x_found(isfinite(x_found)), 'UniformOutput', false), ', ');

    xlabel(ax, '\Deltax  (mm)');
    ylabel(ax, '\Deltay  (mm)');
    title(ax, sprintf(['R_{uu} isocontour  [\\rho = 1/e]  |  ' ...
        'x_{target} = %.0f mm  (actual x_{ref}: %s mm)'], ...
        x_target, x_found_str), 'FontSize', 10);
    grid(ax, 'on');
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');

end % x_target loop

%% ===================== COLOUR SCHEME =====================

nY = numel(sweepFolders);
nX = numel(x_targets);

% x_ref colour: light (upstream) -> dark (downstream)
cmapX = zeros(nX, 3);
for i = 1:nX
    t = (i-1) / max(nX-1, 1);
    cmapX(i,:) = (1-t)*[0.95 0.70 0.30] + t*[0.50 0.05 0.05];
end

%% ===================== PARSE y_ref STRINGS FROM FOLDER NAMES =====================

yref_strs = cell(nY, 1);
for iY = 1:nY
    tok = regexp(sweepFolders{iY}, 'sweep_x_yref([\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        yref_strs{iY} = tok{1};
    else
        yref_strs{iY} = sprintf('%.1f', y_vals(iY));
    end
end

%% ===================== FIGURE: nY x 1 SUBPLOTS =====================

figure('Color','w', ...
    'Position', [50 30 800 220*nY], ...
    'Name',     'R_{uu} isocontour evolution in x_{ref} — per y_{ref}');

axAll = gobjects(nY, 1);

for iY = 1:nY

    axAll(iY) = subplot(nY, 1, iY);
    hold on; box on;

    thisDir = sweepFolders{iY};

    % --- Load grid ---
    gridFile = fullfile(thisDir, 'grid.mat');
    if ~isfile(gridFile)
        warning('No grid.mat in %s — skipping', thisDir);
        title(axAll(iY), sprintf('y_{ref} = %.1f mm — NO GRID', y_vals(iY)));
        continue;
    end
    G      = load(gridFile, 'worldX_merged', 'worldY_merged');
    worldX = double(G.worldX_merged);
    worldY = double(G.worldY_merged);

    x_min_domain = min(worldX(:));
    x_max_domain = max(worldX(:));

    % --- Find all available files and their x_ref values ---
    pattern = sprintf('R_uu_xref*_yref%s.mat', yref_strs{iY});
    files   = dir(fullfile(thisDir, pattern));
    if isempty(files)
        warning('No files matching %s — skipping', pattern);
        title(axAll(iY), sprintf('y_{ref} = %.1f mm — NO FILES', y_vals(iY)));
        continue;
    end

    xr_available = nan(numel(files), 1);
    for f = 1:numel(files)
        tok = regexp(files(f).name, ...
            'R_uu_xref([-+]?\d*\.?\d+)_yref', 'tokens', 'once');
        if ~isempty(tok)
            xr_available(f) = str2double(tok{1});
        end
    end

    fprintf('\ny_ref = %.1f mm\n', y_vals(iY));

    x_found = nan(nX, 1);

    for iX = 1:nX

        x_target = x_targets(iX);

        % Find closest file
        [~, iClosest] = min(abs(xr_available - x_target));
        chosenFile    = fullfile(thisDir, files(iClosest).name);
        xr_actual     = xr_available(iClosest);
        x_found(iX)   = xr_actual;

        fprintf('  x_target = %6.1f mm | using x_ref = %.2f mm\n', ...
            x_target, xr_actual);

        % --- Load ---
        S  = load(chosenFile, 'R_s', 'xr', 'yr');
        R  = double(S.R_s);
        xr = double(S.xr);
        yr = double(S.yr);

        % --- Delta grids ---
        dX = worldX - xr;
        dY = worldY - yr;

        % --- Clamped search box ---
        x_lo = max(xr - dx_back, x_min_domain);
        x_hi = min(xr + dx_fwd,  x_max_domain);

        boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
        R_box           = R;
        R_box(~boxMask) = NaN;
        R_box(R_box >  1.0) = 1.0;
        R_box(R_box < -1.0) = NaN;

        dX_box           = dX;  dX_box(~boxMask) = NaN;
        dY_box           = dY;  dY_box(~boxMask) = NaN;

        dx_vec = dX(1,:);
        dy_vec = dY(:,1);

        % --- Threshold mask ---
        mask = R_box >= rho_level;
        if sum(mask(:)) < 5
            fprintf('    WARNING: fewer than 5 points above threshold\n');
            continue;
        end

        % --- Find origin pixel ---
        dist2           = dX.^2 + dY.^2;
        dist2(~boxMask) = Inf;
        [~, i0]         = min(dist2(:));

        % --- Connected components ---
        CC = bwconncomp(mask);
        if CC.NumObjects == 0, continue; end

        regionIdx = 0;
        for r = 1:CC.NumObjects
            if any(CC.PixelIdxList{r} == i0)
                regionIdx = r;
                break;
            end
        end
        if regionIdx == 0
            regionSizes = cellfun(@numel, CC.PixelIdxList);
            [~, regionIdx] = max(regionSizes);
        end

        mask_origin = false(size(mask));
        mask_origin(CC.PixelIdxList{regionIdx}) = true;

        % --- Draw isocontour ---
        B = bwboundaries(mask_origin);
        for b = 1:numel(B)
            bnd    = B{b};
            bnd_dx = dx_vec(bnd(:,2));
            bnd_dy = dy_vec(bnd(:,1));
            plot(axAll(iY), bnd_dx, bnd_dy, '-', ...
                'Color',       cmapX(iX,:), ...
                'LineWidth',   1.8, ...
                'DisplayName', sprintf('x_{ref} = %.0f mm', xr_actual));
        end

        % --- Bounding box (dashed) ---
        dX_pts = dX_box(mask_origin);
        dY_pts = dY_box(mask_origin);
        if numel(dX_pts) >= 5
            Lxu = abs(min(dX_pts));
            Lxd = max(dX_pts);
            Lyt = max(dY_pts);
            Lyb = abs(min(dY_pts));
            bx  = [-Lxu,  Lxd,  Lxd, -Lxu, -Lxu];
            by  = [-Lyb, -Lyb,  Lyt,  Lyt, -Lyb];
            plot(axAll(iY), bx, by, '--', ...
                'Color',            cmapX(iX,:), ...
                'LineWidth',        0.7, ...
                'HandleVisibility', 'off');
        end

    end % x_target loop

    % Reference point marker and zero lines
    plot(axAll(iY), 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    xline(axAll(iY), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(axAll(iY), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Deduplicate legend — one entry per x_target
    h    = get(axAll(iY), 'Children');
    seen = false(nX, 1);
    hKeep = gobjects(0);
    for k = 1:numel(h)
        dn = get(h(k), 'DisplayName');
        for iX = 1:nX
            lbl = sprintf('x_{ref} = %.0f mm', x_found(iX));
            if strcmp(dn, lbl) && ~seen(iX)
                seen(iX)     = true;
                hKeep(end+1) = h(k); %#ok<AGROW>
            end
        end
    end
    legend(axAll(iY), hKeep, 'Location', 'eastoutside', 'FontSize', 8);

    ylabel(axAll(iY), '\Deltay  (mm)');
    title(axAll(iY), ...
        sprintf('y_{ref} = %.1f mm', y_vals(iY)), ...
        'FontWeight', 'normal', 'FontSize', 10);
    grid(axAll(iY), 'on');
    set(axAll(iY), 'YDir', 'normal');
    axis(axAll(iY), 'equal');

    if iY == nY
        xlabel(axAll(iY), '\Deltax  (mm)');
    end

end % y_ref loop

% Link x and y axes across all subplots for direct comparison
linkaxes(axAll, 'xy');

sgtitle(sprintf('R_{uu} isocontour evolution  [\\rho = 1/e]  —  x_{ref} = %s mm', ...
    strjoin(arrayfun(@(v) sprintf('%.0f', v), x_targets, ...
    'UniformOutput', false), ', ')), ...
    'FontSize', 11);

%% ===================== FIGURE: Lx vs y/delta =====================
% For each x_ref target, plot Lx vs y/delta
% Lines coloured by x_ref position
% Requires delta (boundary layer thickness) at each x_ref — user input

% --- User input: delta values at each x_target (mm) ---
% These should come from your boundary layer sweep results
% Order must match x_targets
% ===================== LOAD delta FROM blSweep =====================

blSweepFile = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';

B       = load(blSweepFile, 'blSweep');
bl_x    = B.blSweep.x_mm;
bl_d99  = B.blSweep.delta99_hybrid_mm;

% Remove NaNs before interpolation
valid_bl       = isfinite(bl_x) & isfinite(bl_d99);
bl_x_clean     = bl_x(valid_bl);
bl_d99_clean   = bl_d99(valid_bl);

% Interpolate delta99 at each x_target
% Use 'linear' with NaN for any x_target outside the bl_x range
delta_vals = interp1(bl_x_clean, bl_d99_clean, x_targets, 'linear', NaN);
figure('Color','w', 'Position',[100 80 700 500], ...
    'Name', 'L_x vs y/delta');
axLx = axes; hold on; box on;


% Report
fprintf('\ndelta99 at each x_target:\n');
for iX = 1:nX
    if isfinite(delta_vals(iX))
        fprintf('  x_target = %6.1f mm  ->  delta99 = %.2f mm\n', ...
            x_targets(iX), delta_vals(iX));
    else
        fprintf('  x_target = %6.1f mm  ->  delta99 = NaN (outside bl_x range — skipping)\n', ...
            x_targets(iX));
    end
end
% -------------------------------------------------------


for iX = 1:nX

    % Skip if delta unavailable at this x_target
    if ~isfinite(delta_vals(iX))
        fprintf('Skipping x_target = %.1f mm — no delta99 available\n', ...
            x_targets(iX));
        continue;
    end

    x_target = x_targets(iX);
    delta    = delta_vals(iX);


    Lx_prof  = nan(nY, 1);
    Lx_norm  = nan(nY,1) ; 
    y_norm   = nan(nY, 1);

    for iY = 1:nY

        thisDir = sweepFolders{iY};

        % Load grid
        gridFile = fullfile(thisDir, 'grid.mat');
        if ~isfile(gridFile), continue; end
        G      = load(gridFile, 'worldX_merged', 'worldY_merged');
        worldX = double(G.worldX_merged);
        worldY = double(G.worldY_merged);

        x_min_domain = min(worldX(:));
        x_max_domain = max(worldX(:));

        % Find closest file
        pattern = sprintf('R_uu_xref*_yref%s.mat', yref_strs{iY});
        files   = dir(fullfile(thisDir, pattern));
        if isempty(files), continue; end

        xr_available = nan(numel(files), 1);
        for f = 1:numel(files)
            tok = regexp(files(f).name, ...
                'R_uu_xref([-+]?\d*\.?\d+)_yref', 'tokens', 'once');
            if ~isempty(tok)
                xr_available(f) = str2double(tok{1});
            end
        end

        [~, iClosest] = min(abs(xr_available - x_target));
        chosenFile    = fullfile(thisDir, files(iClosest).name);

        % Load
        S  = load(chosenFile, 'R_s', 'xr', 'yr');
        R  = double(S.R_s);
        xr = double(S.xr);
        yr = double(S.yr);

        % Delta grids and search box
        dX = worldX - xr;
        dY = worldY - yr;

        x_lo = max(xr - dx_back, x_min_domain);
        x_hi = min(xr + dx_fwd,  x_max_domain);

        boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
        R_box           = R;
        R_box(~boxMask) = NaN;
        R_box(R_box >  1.0) = 1.0;
        R_box(R_box < -1.0) = NaN;

        dX_box           = dX;  dX_box(~boxMask) = NaN;
        dY_box           = dY;  dY_box(~boxMask) = NaN;

        % Threshold and connected region
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
                regionIdx = r;
                break;
            end
        end
        if regionIdx == 0
            regionSizes = cellfun(@numel, CC.PixelIdxList);
            [~, regionIdx] = max(regionSizes);
        end

        mask_origin = false(size(mask));
        mask_origin(CC.PixelIdxList{regionIdx}) = true;

        dX_pts = dX_box(mask_origin);
        dY_pts = dY_box(mask_origin);
        if numel(dX_pts) < 5, continue; end

        % Bounding box
        Lxu = abs(min(dX_pts));
        Lxd = max(dX_pts);

        Lx_prof(iY) = Lxu + Lxd;
        Lx_norm(iY) = (Lxu + Lxd) /delta; 
        y_norm(iY)  = yr / delta;

    end % y_ref loop

    % Plot profile for this x_target
    valid = isfinite(Lx_norm) & isfinite(y_norm);
    plot(axLx, y_norm(valid), Lx_norm(valid), 'o-', ...
        'Color',           cmapX(iX,:), ...
        'LineWidth',       1.8, ...
        'MarkerFaceColor', cmapX(iX,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize',      6, ...
        'DisplayName',     sprintf('x_{ref} \\approx %.0f mm', x_target));

end % x_target loop

xlabel(axLx, 'y / \delta');
ylabel(axLx, 'L_x / \delta');
title(axLx,  sprintf('Streamwise extent vs wall-normal position  [\\rho = 1/e]'));
legend(axLx, 'Location', 'best');
grid(axLx,   'on');
set(axLx,    'XScale', 'log');   % log scale on y/delta is standard in TBL lit
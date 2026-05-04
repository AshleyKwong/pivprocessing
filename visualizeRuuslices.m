% =========================================================================
% plot_Ruu_3D_stack.m
%
% Two modes:
%   'sweep' — fixed xref, varying yref from sweep folders
%             (yref in mm, not normalised — delta varies with x)
%   'point' — fixed yref targets in y/delta, varying xref
%             Each xref folder contains files at the correct yref for
%             that location's delta. Ribbons are at consistent y/delta.
%
% Axes:
%   X  — Deltax / delta
%   Y  — y_ref / delta  (depth axis)
%   Z  — Deltay / delta
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

% 'sweep' — original mode, fixed xref, yref sweep folders
% 'point' — new mode, fixed yref/delta targets, xref folders
corrMode = 'point';

blSweepFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case1_PIVresults\blSweep_20260419_193106.mat'];

rho_level = 1/exp(1);
dx_back   = 300;
dx_fwd    = 600;

% Shared grid file
gridFile = 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\grid.mat';

if strcmp(corrMode, 'sweep')
    sweepFolders = {
        'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244\sweep_x_yref2.0'
        'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244\sweep_x_yref6.0'
        'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244\sweep_x_yref14.0'
        'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\two_point_covariance_20260417_160244\sweep_x_yref89.1'
    };
    y_vals   = [2.0, 6.0, 14.0, 89.1];   % mm
    x_target = 550;                        % mm — single xref to inspect

else   % point

    xrefFolders = {
        % 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\R_uu_xref350'
        % 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\R_uu_xref550'
        % 'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\R_uu_xref711'
        'D:\FULLYPROCESSEDY250AOAN04AOAFN04PIVDATA\two_point_covariance_20260419_201038\R_uu_xref1000'
    };
    xref_vals = [ 1000];   % mm 350, 550, 711,

    % Target y/delta values — ribbons shown at these wall-normal positions
    ydelta_targets = [0.02, 0.04, 0.05, 0.1, 0.5, 0.8];
end

%% ===================== LOAD delta99 =====================

B      = load(blSweepFile, 'blSweep');
bl_x   = B.blSweep.x_mm;
bl_d99 = B.blSweep.delta99_hybrid_mm;

% Remove non-finite points
valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
bl_x_clean   = bl_x(valid_bl);
bl_d99_clean = bl_d99(valid_bl);
if any(diff(bl_x_clean) == 0)
    [bl_x_clean, ~, ic] = unique(bl_x_clean, 'sorted');
    bl_d99_clean        = accumarray(ic, bl_d99_clean, [], @mean);
    fprintf('Duplicates found and removed: %d -> %d unique points\n', ...
        numel(bl_x), numel(bl_x_clean));
end

if strcmp(corrMode, 'sweep')
    delta = interp1(bl_x_clean, bl_d99_clean, x_target, 'linear', NaN);
    assert(isfinite(delta), 'delta99 not available at x_target = %.1f mm', x_target);
    fprintf('x_target = %.1f mm  |  delta99 = %.2f mm\n', x_target, delta);
else
    delta_at_xref = interp1(bl_x_clean, bl_d99_clean, xref_vals, 'linear', NaN);
    fprintf('delta99 at each xref:\n');
    for i = 1:numel(xref_vals)
        fprintf('  x = %.1f mm  ->  delta99 = %.2f mm\n', ...
            xref_vals(i), delta_at_xref(i));
    end
end

%% ===================== LOAD GRID =====================

G      = load(gridFile, 'worldX_merged', 'worldY_merged');
worldX = double(G.worldX_merged);
worldY = double(G.worldY_merged);

x_min_domain = min(worldX(:));
x_max_domain = max(worldX(:));
dx_vec_full  = worldX(1,:);
dy_vec_full  = worldY(:,1);

%% ===================== COLOUR SCHEME =====================

if strcmp(corrMode, 'sweep')
    nLines = numel(sweepFolders);
else
    nLines = numel(ydelta_targets);
end

cmapLines = zeros(nLines, 3);
for i = 1:nLines
    t = (i-1) / max(nLines-1, 1);
    cmapLines(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
end

%% ===================== FIGURE SETUP =====================

if strcmp(corrMode, 'sweep')
    figTitle = sprintf('3D R_{uu} stack  [\\rho=1/e]  |  x_{ref} = %.0f mm', x_target);
else
    figTitle = '3D R_{uu} stack  [\rho=1/e]  |  point mode — consistent y/\delta';
end

figure('Color','w','Position',[100 80 900 700], 'Name', figTitle);
ax3d = axes;
hold(ax3d, 'on');
view(ax3d, -35, 25);
grid(ax3d, 'on');
box(ax3d,  'on');

xlabel(ax3d, '\Deltax / \delta');
ylabel(ax3d, 'y_{ref} / \delta');
zlabel(ax3d, '\Deltay / \delta');
title(ax3d, figTitle, 'FontSize', 10);

legend_handles = gobjects(0);
legend_labels  = {};

%% ===================== HELPER: EXTRACT AND PLOT CONTOUR =====================

    function plot_contour(ax3d, R, xr, yr, delta, worldX, worldY, ...
            dx_back, dx_fwd, x_min_domain, x_max_domain, ...
            rho_level, col, label)

        dX = worldX - xr;
        dY = worldY - yr;

        x_lo = max(xr - dx_back, x_min_domain);
        x_hi = min(xr + dx_fwd,  x_max_domain);

        boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
        R_box           = R;
        R_box(~boxMask) = NaN;
        R_box(R_box >  1.0) = 1.0;
        R_box(R_box < -1.0) = NaN;

        dx_vec = worldX(1,:);
        dy_vec = worldY(:,1);

        mask = R_box >= rho_level;
        if sum(mask(:)) < 5, return; end

        dist2           = dX.^2 + dY.^2;
        dist2(~boxMask) = Inf;
        [~, i0]         = min(dist2(:));

        CC = bwconncomp(mask);
        if CC.NumObjects == 0, return; end

        regionIdx = 0;
        for r = 1:CC.NumObjects
            if any(CC.PixelIdxList{r} == i0)
                regionIdx = r; break;
            end
        end
        if regionIdx == 0
            regionSizes = cellfun(@numel, CC.PixelIdxList);
            [~, regionIdx] = max(regionSizes);
        end

        mask_origin = false(size(mask));
        mask_origin(CC.PixelIdxList{regionIdx}) = true;

        B_contour = bwboundaries(mask_origin);
        if isempty(B_contour), return; end

        ydelta = yr / delta;

        for b = 1:numel(B_contour)
            bnd    = B_contour{b};
            bnd_dx = dx_vec(bnd(:,2)) / delta;
            bnd_dy = dy_vec(bnd(:,1)) / delta;
            bnd_z  = repmat(ydelta, size(bnd_dx));

            h = plot3(ax3d, bnd_dx, bnd_z, bnd_dy, '-', ...
                'Color', col, 'LineWidth', 2.0);

            if b == 1 && ~isempty(label)
                assignin('caller', 'h_out', h);
                assignin('caller', 'label_out', label);
            end
        end

        % Filled translucent patch
        bnd    = B_contour{1};
        bnd_dx = dx_vec(bnd(:,2)) / delta;
        bnd_dy = dy_vec(bnd(:,1)) / delta;
        bnd_z  = repmat(ydelta, size(bnd_dx));
        fill3(ax3d, bnd_dx, bnd_z, bnd_dy, col, ...
            'FaceAlpha', 0.08, 'EdgeColor', 'none');

        % Reference point marker
        plot3(ax3d, 0, ydelta, 0, '+', ...
            'Color', col, 'MarkerSize', 8, 'LineWidth', 1.5);
    end

%% ===================== SWEEP MODE =====================

if strcmp(corrMode, 'sweep')

    for iY = 1:numel(sweepFolders)

        thisDir = sweepFolders{iY};
        col     = cmapLines(iY,:);

        % Parse yref from folder name
        tok = regexp(thisDir, 'sweep_x_yref([\d.]+)', 'tokens', 'once');
        if isempty(tok), continue; end
        yref_str = tok{1};

        pattern = sprintf('R_uu_xref*_yref%s.mat', yref_str);
        files   = dir(fullfile(thisDir, pattern));
        if isempty(files), continue; end

        % Find closest file to x_target
        xr_available = nan(numel(files), 1);
        for f = 1:numel(files)
            t2 = regexp(files(f).name, ...
                'R_uu_xref([-+]?\d*\.?\d+)_yref', 'tokens', 'once');
            if ~isempty(t2), xr_available(f) = str2double(t2{1}); end
        end
        [~, iClosest] = min(abs(xr_available - x_target));
        S  = load(fullfile(thisDir, files(iClosest).name), 'R_s', 'xr', 'yr');
        R  = double(S.R_s);
        xr = double(S.xr);
        yr = double(S.yr);

        fprintf('  sweep | y_ref = %.1f mm | y/delta = %.3f\n', yr, yr/delta);

        label = sprintf('y/\\delta = %.2f', yr/delta);
        h_out = []; label_out = '';
        plot_contour(ax3d, R, xr, yr, delta, worldX, worldY, ...
            dx_back, dx_fwd, x_min_domain, x_max_domain, ...
            rho_level, col, label);
        if ~isempty(h_out)
            legend_handles(end+1) = h_out;
            legend_labels{end+1}  = label_out;
        end
    end

%% ===================== POINT MODE =====================

else

    nXref    = numel(xrefFolders);
    nTargets = numel(ydelta_targets);

    % For each y/delta target, plot one ribbon per xref
    % Each ribbon is at the consistent y/delta level, just at a different x
    % The depth axis (Y) is y/delta — same for all xref at this target

    for iT = 1:nTargets

        yd_target = ydelta_targets(iT);
        col       = cmapLines(iT,:);

        for iX = 1:nXref

            thisDir = xrefFolders{iX};
            xr_nom  = xref_vals(iX);
            delta   = delta_at_xref(iX);

            if ~isfinite(delta), continue; end

            % Target yref in mm for this xref and y/delta target
            yr_target = yd_target * delta;

            % Find all files in this folder
            files = dir(fullfile(thisDir, 'R_uu_xref*_yref*.mat'));
            if isempty(files), continue; end

            % Parse yref from each filename
            yr_available = nan(numel(files), 1);
            xr_available = nan(numel(files), 1);
            for f = 1:numel(files)
                tok1 = regexp(files(f).name, ...
                    'R_uu_xref([-+]?\d*\.?\d+)_yref([-+]?\d*\.?\d+)', ...
                    'tokens', 'once');
                if ~isempty(tok1)
                    xr_available(f) = str2double(tok1{1});
                    yr_available(f) = str2double(tok1{2});
                end
            end

            % Find file with yref closest to yr_target
            % and xref closest to xr_nom
            xr_ok = abs(xr_available - xr_nom) < 5;   % within 5mm
            if ~any(xr_ok)
                [~, iBest] = min(abs(xr_available - xr_nom));
                xr_ok(iBest) = true;
            end

            yr_sub = yr_available;
            yr_sub(~xr_ok) = Inf;
            [~, iClosest] = min(abs(yr_sub - yr_target));

            chosenFile = fullfile(thisDir, files(iClosest).name);
            S  = load(chosenFile, 'R_s', 'xr', 'yr');
            R  = double(S.R_s);
            xr = double(S.xr);
            yr = double(S.yr);

            actual_ydelta = yr / delta;

            fprintf('  point | xref = %.0f mm | y/delta target = %.2f | actual = %.3f | yr = %.1f mm\n', ...
                xr_nom, yd_target, actual_ydelta, yr);

            % Only add to legend on first xref for this y/delta target
            if iX == 1
                label = sprintf('y/\\delta \\approx %.2f', yd_target);
            else
                label = '';
            end

            h_out = []; label_out = '';
            plot_contour(ax3d, R, xr, yr, delta, worldX, worldY, ...
                dx_back, dx_fwd, x_min_domain, x_max_domain, ...
                rho_level, col, label);
            if ~isempty(h_out) && ~isempty(label_out)
                legend_handles(end+1) = h_out; %#ok<AGROW>
                legend_labels{end+1}  = label_out; %#ok<AGROW>
            end

        end % xref loop

    end % ydelta target loop

end

% ===================== FINISHING TOUCHES =====================

% Anchor line through reference points
yl = ylim(ax3d);
plot3(ax3d, [0 0], yl, [0 0], 'k:', 'LineWidth', 1.0);

% Zero planes
xl = xlim(ax3d);
zl = zlim(ax3d);
yl = ylim(ax3d);

[Yp, Zp] = meshgrid(yl, zl);
surf(ax3d, zeros(size(Yp)), Yp, Zp, ...
    'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.05, 'EdgeColor', 'none');

[Xp, Yp] = meshgrid(xl, yl);
surf(ax3d, Xp, Yp, zeros(size(Xp)), ...
    'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.05, 'EdgeColor', 'none');

if ~isempty(legend_handles)
    legend(ax3d, legend_handles, legend_labels, ...
        'Location', 'best', 'FontSize', 8);
end

rotate3d(ax3d, 'on');
fprintf('\nTip: click and drag to rotate the 3D view\n');
%% ===================== FIGURE 2: 2D OVERLAY — ONE ROW PER y/delta =====================

if strcmp(corrMode, 'point')

    nTargets = numel(ydelta_targets);
    nXref    = numel(xrefFolders);

    % xref colour scheme
    cmapX = zeros(nXref, 3);
    for i = 1:nXref
        t = (i-1) / max(nXref-1, 1);
        cmapX(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
    end

    figure('Color','w', ...
        'Position', [200 50 700 250*nTargets], ...
        'Name', 'R_uu isocontour overlay — fixed y/delta, varying xref');

    axArr = gobjects(nTargets, 1);

    for iT = 1:nTargets

        yd_target = ydelta_targets(iT);
        axArr(iT) = subplot(nTargets, 1, iT);
        hold on; box on;

        leg_h   = gobjects(0);
        leg_lbl = {};

        for iX = 1:nXref

            thisDir = xrefFolders{iX};
            xr_nom  = xref_vals(iX);
            delta   = delta_at_xref(iX);
            col     = cmapX(iX,:);

            if ~isfinite(delta), continue; end

            yr_target = yd_target * delta;

            % Find files
            files = dir(fullfile(thisDir, 'R_uu_xref*_yref*.mat'));
            if isempty(files), continue; end

            yr_available = nan(numel(files), 1);
            xr_available = nan(numel(files), 1);
            for f = 1:numel(files)
                tok1 = regexp(files(f).name, ...
                    'R_uu_xref([-+]?\d*\.?\d+)_yref([-+]?\d*\.?\d+)', ...
                    'tokens', 'once');
                if ~isempty(tok1)
                    xr_available(f) = str2double(tok1{1});
                    yr_available(f) = str2double(tok1{2});
                end
            end

            % Find closest xref then closest yref
            xr_ok = abs(xr_available - xr_nom) < 5;
            if ~any(xr_ok)
                [~, iBest] = min(abs(xr_available - xr_nom));
                xr_ok(iBest) = true;
            end
            yr_sub = yr_available;
            yr_sub(~xr_ok) = Inf;
            [~, iClosest] = min(abs(yr_sub - yr_target));

            S  = load(fullfile(thisDir, files(iClosest).name), 'R_s', 'xr', 'yr');
            R  = double(S.R_s);
            xr = double(S.xr);
            yr = double(S.yr);

            actual_ydelta = yr / delta;

            % Delta grids — separation from reference point
            dX = worldX - xr;   % [Ny x Nx] separation in x
            dY = worldY - yr;   % [Ny x Nx] separation in y

            x_lo = max(xr - dx_back, x_min_domain);
            x_hi = min(xr + dx_fwd,  x_max_domain);

            boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
            R_box           = R;
            R_box(~boxMask) = NaN;
            R_box(R_box >  1.0) = 1.0;
            R_box(R_box < -1.0) = NaN;

            % 1D separation vectors for boundary indexing
            dX_vec = dX(1,:);   % [1 x Nx] — row of separations
            dY_vec = dY(:,1);   % [Ny x 1] — column of separations

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
                regionSizes = cellfun(@numel, CC.PixelIdxList);
                [~, regionIdx] = max(regionSizes);
            end

            mask_origin = false(size(mask));
            mask_origin(CC.PixelIdxList{regionIdx}) = true;

            B_contour = bwboundaries(mask_origin);
            if isempty(B_contour), continue; end

            % Plot — Deltax/delta on x, Deltay/delta on y
            % Both centred at 0 by construction
            bnd    = B_contour{1};
            bnd_dx = dX_vec(bnd(:,2)) / delta;   % Deltax/delta
            bnd_dy = dY_vec(bnd(:,1)) / delta;   % Deltay/delta

            h = plot(axArr(iT), bnd_dx, bnd_dy, '-', ...
                'Color',     col, ...
                'LineWidth', 1.8);

            % Mark reference point
            plot(axArr(iT), 0, 0, '+', ...
                'Color',      col, ...
                'MarkerSize', 8, ...
                'LineWidth',  1.5);

            leg_h(end+1)   = h; %#ok<AGROW>
            leg_lbl{end+1} = sprintf('x_{ref}=%.0f mm', ...
                xr_nom); %#ok<AGROW>

        end % xref loop

        % Reference lines through origin
        xline(axArr(iT), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(axArr(iT), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        ylabel(axArr(iT), '\Deltay / \delta');
        title(axArr(iT), sprintf('y/\\delta \\approx %.2f', yd_target), ...
            'FontWeight', 'normal');
        grid(axArr(iT), 'on');
        % axis(axArr(iT), 'equal');

        if iT == nTargets
            xlabel(axArr(iT), '\Deltax / \delta');
        end
        if iT == 1
            legend(axArr(iT), leg_h, leg_lbl, ...
                'Location', 'best', 'FontSize', 12);
        end

    end % target loop

    linkaxes(axArr, 'xy');

    sgtitle([sprintf('R_{uu} isocontour evolution [\\rho = %.2f]', rho_level)], 'FontSize', 14);

end
%%
% =========================================================================
% plot_Ruu_2D_overlay.m
%
% Plots R_uu isocontours (2D) at consistent y/delta wall-normal positions,
% for one or more xref locations. Each subplot shows one y/delta level
% with all xref locations overlaid.
%
% Axes:  Deltax/delta (x)  vs  Deltay/delta (y)
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc; close all;

%% ===================== USER INPUTS =====================

blSweepFile = ['C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat'];

rho_level = 0.5;
dx_back   = 300;   % mm upstream of xref to show
dx_fwd    = 600;   % mm downstream of xref to show

% grid.mat from the covariance output folder for this position
gridFile = 'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\grid.mat';

% Covariance result folders — one per xref location
xrefFolders = {
    'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\R_uu_xref54\'
    'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\R_uu_xref350\'
    'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\R_uu_xref550\'
    'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\R_uu_xref711\'
    'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\R_uu_xref1000\'
    'G:\SW 400mm Processed Data Single Snapshots\-8\two_point_covariance_20260423_175816_compTom\R_uu_xref1150\'
};
xref_vals      = [54 350 550 711 1000 1150];        % mm — must match folder list above
ydelta_targets = [0.02, 0.04, 0.05, 0.1, 0.5, 0.8];

%% ===================== LOAD delta99 =====================

B      = load(blSweepFile, 'blSweep');
bl_x   = B.blSweep.x_mm;
bl_d99 = B.blSweep.delta99_hybrid_mm;

valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
bl_x_clean   = bl_x(valid_bl);
bl_d99_clean = bl_d99(valid_bl);

if any(diff(bl_x_clean) == 0)
    [bl_x_clean, ~, ic] = unique(bl_x_clean, 'sorted');
    bl_d99_clean        = accumarray(ic, bl_d99_clean, [], @mean);
    fprintf('Duplicates found and removed: %d -> %d unique points\n', ...
        numel(bl_x), numel(bl_x_clean));
end

delta_at_xref = interp1(bl_x_clean, bl_d99_clean, xref_vals, 'linear', NaN);
fprintf('delta99 at each xref:\n');
for i = 1:numel(xref_vals)
    fprintf('  x = %.1f mm  ->  delta99 = %.2f mm\n', ...
        xref_vals(i), delta_at_xref(i));
end

%% ===================== LOAD GRID =====================

% grid.mat is saved as single — cast to double for all subsequent arithmetic
G      = load(gridFile, 'worldX_merged', 'worldY_merged');
worldX = double(G.worldX_merged);   % [Ny x Nx] mm — per-position native grid
worldY = double(G.worldY_merged);   % [Ny x Nx] mm

x_min_domain = min(worldX(:));
x_max_domain = max(worldX(:));

%% ===================== COLOUR SCHEME =====================

nLines = numel(ydelta_targets);
cmapLines = zeros(nLines, 3);
for i = 1:nLines
    t = (i-1) / max(nLines-1, 1);
    cmapLines(i,:) = (1-t)*[0.75 0.88 1.0] + t*[0.0 0.05 0.25];
end

%% ===================== FIGURE: 2D OVERLAY =====================

    nTargets = numel(ydelta_targets);
    nXref    = numel(xrefFolders);

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

            dX = worldX - xr;
            dY = worldY - yr;

            x_lo = max(xr - dx_back, x_min_domain);
            x_hi = min(xr + dx_fwd,  x_max_domain);

            boxMask = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));

            % Cast to double before masking
            R_box           = double(R);
            R_box(~boxMask) = NaN;
            R_box(R_box >  1.0) = 1.0;
            R_box(R_box < -1.0) = NaN;

            dX_vec = dX(1,:);
            dY_vec = dY(:,1);

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

            bnd    = B_contour{1};
            bnd_dx = dX_vec(bnd(:,2)) / delta;
            bnd_dy = dY_vec(bnd(:,1)) / delta;

            h = plot(axArr(iT), bnd_dx, bnd_dy, '-', ...
                'Color',     col, ...
                'LineWidth', 1.8);

            plot(axArr(iT), 0, 0, '+', ...
                'Color',      col, ...
                'MarkerSize', 8, ...
                'LineWidth',  1.5);

            leg_h(end+1)   = h; %#ok<AGROW>
            leg_lbl{end+1} = sprintf('x_{ref}=%.0f mm  (y/\\delta=%.3f)', ...
                xr_nom, actual_ydelta); %#ok<AGROW>

        end

        xline(axArr(iT), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        yline(axArr(iT), 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        ylabel(axArr(iT), '\Deltay / \delta');
        title(axArr(iT), sprintf('y/\\delta \\approx %.2f', yd_target), ...
            'FontWeight', 'normal');
        grid(axArr(iT), 'on');

        if iT == nTargets
            xlabel(axArr(iT), '\Deltax / \delta');
        end
        if iT == 1
            legend(axArr(iT), leg_h, leg_lbl, ...
                'Location', 'best', 'FontSize', 12);
        end

    end

    linkaxes(axArr, 'xy');
    sgtitle(sprintf('R_{uu} isocontour evolution [\\rho = %.2f]', rho_level), ...
        'FontSize', 14);
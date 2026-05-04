% make_Ruu_mp4_sync_pressure.m
% MP4 + GIF of correlation planes with synchronized dCp/dx subplot
%
% Pressure input is supplied in x/delta, then converted to physical x:
%   x_phys = x_over_delta * delta_0 + x_offset
%
% Normalization by BL_0 is done only for plotting.

clc; close all; clear   

%% ===== USER INPUTS =====
baseDir   = 'D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\sweep_x_yref89.1\';
gridFile  = fullfile(baseDir, 'grid.mat');
csvFile   = fullfile(baseDir, 'integral_length_summary_dx.csv');

corrPrefix = 'R_uu_xref';     % e.g. R_uu_xref / R_uv_xref / R_vv_xref
fixedYref  = 89.1;             % this folder's y_ref

gifName    = 'Ruu_sweep_sync_yref89p1_nooutlier.gif';
mp4Name    = 'Ruu_sweep_sync_yref89p1_nooutlier.mp4';
frameDelay = 0.12;            % seconds per frame for GIF
safeDelay  = 0.08;            % used for MP4 frame rate
x_LE_mm = 7650; 

% --- Pressure-trace inputs ---
load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\Pressure\pressuredata_wingcases_revisedNOZPGZERO.mat'); 
load('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat'); 
[~, x0] = min(abs(blSweep.x_mm - 0)); 
[~, xf] = min(abs(blSweep.x_mm - 40)); 
BL_0_mm = mean(blSweep.delta99_hybrid_mm(x0:xf)) ; % taking first 10 cm

caseNo = 2; 
% --- pressure data ---
xdelta    = case_pdata( caseNo).xloc_full_xdelta;
deltaBL_m = case_pdata( caseNo).refBLmm .* 0.01;  % in 
mean_Cp   = mean(case_pdata( caseNo).Cp_full, 2, 'omitnan');

x_pdata_m  = xdelta .* deltaBL_m + 7.65;   % recover global x [mm]: undo original non-dim + LE offset

[x_p_sorted, p_sortidx] = sort(x_pdata_m, 'ascend');
Cp_sorted    = mean_Cp(p_sortidx);

xdelta_sorted = x_pdata_m(p_sortidx) ./ (BL_0_mm/1000);   % x/delta_actual: global x / 116mm, no LE subtraction

valid_idx    = ~isnan(Cp_sorted) & ~isnan(xdelta_sorted) & ~isinf(xdelta_sorted);
xdelta_valid = xdelta_sorted(valid_idx);
Cp_valid     = Cp_sorted(valid_idx);
x_p_valid_sc = x_p_sorted(valid_idx);

[xdelta_unique, ~, ic] = unique(xdelta_valid, 'sorted');
Cp_unique              = accumarray(ic, Cp_valid,     [], @mean);
x_p_unique_sc          = accumarray(ic, x_p_valid_sc, [], @mean);

dCp_dxdelta = gradient(Cp_unique, xdelta_unique);
sg_order  = 3;
sg_win_dcp = 11;

sg_win_dcp  = min(sg_win_dcp, 2*floor(length(dCp_dxdelta)/2)-1);
sg_win_dcp  = max(sg_win_dcp, sg_order+2);
if mod(sg_win_dcp, 2)==0, sg_win_dcp = sg_win_dcp+1; end
dcpdx_deltanorm = sgolayfilt(dCp_dxdelta, sg_order, sg_win_dcp);

% x_delta = (case_pdata(caseNo).xloc_testsection_xdelta.*case_pdata(caseNo).refBLmm)./ (BL_0_mm/10); % in cm
% dcpdx_deltanorm = case_pdata(caseNo).mean_dcpdxdelta.* (((BL_0_mm/10)/case_pdata(caseNo).refBLmm));
% filter_dcpdxdelta = movmean(dcpdx_deltanorm(~isoutlier(dcpdx_deltanorm)),2);
% filter_dcpdxdelta = movmean(dcpdx_deltanorm,2);
% filter_dcpdxdelta(7) = [];
% x_delta = x_delta(~isoutlier(dcpdx_deltanorm));
% x_delta(7) = [];

x_over_delta = xdelta_unique - x_LE_mm/BL_0_mm;    % vector of x/delta values for dCp/dx zeroed at xLE. 
dCpdx        = dcpdx_deltanorm;     % vector of dCp/dx values, same length as x_over_delta
delta_0      = BL_0_mm;               % converts x/delta to physical x  
x_offset     = 0;           % physical x-offset to match correlation coordinates

% --- Plot normalization (post only) ---
BL_0         = delta_0;               % use 1 for physical axes, or e.g. delta_0 for x/BL_0 and y/BL_0

% --- Plot settings ---
corrClim     = [0 1];
corrLabel    = 'R_{uu}';
lineColor    = [1 0.2 0.2];
figPos       = [100 100 1500 900];

% --- Fixed dCp/dx limits ---
dCpdxYLim    = [-0.15 0.1];

% --- Camera boundary / overlap overlay ---
cameraBoundaryFile = 'C:\Users\ak1u24\Downloads\PG_fixedmergandcal\windowCenterCameras_mm.mat';
showCameraEdges    = true;     % plots every camera min/max boundary
showOverlapRegions = true;     % shades overlap regions between adjacent cameras

cameraEdgeColor    = [1.0 0.9 0.2];
cameraEdgeWidth    = 1.1;

overlapColor       = [0.2 1.0 0.2];
overlapFaceAlpha   = 0.10;
% ===============================

%% ===== INPUT CHECKS =====
if isempty(x_over_delta) || isempty(dCpdx)
    error('Provide x_over_delta and dCpdx as user inputs.');
end

if isempty(delta_0)
    error('Provide delta_0 as user input.');
end

if isempty(BL_0) || ~isscalar(BL_0) || BL_0 == 0
    error('BL_0 must be a nonzero scalar.');
end

if numel(x_over_delta) ~= numel(dCpdx)
    error('x_over_delta and dCpdx must have the same length.');
end

x_over_delta = x_over_delta(:);
dCpdx        = dCpdx(:);

%% ===== LOAD GRID =====
G = load(gridFile, 'worldX_merged', 'worldY_merged');
X = G.worldX_merged;
Y = G.worldY_merged;

xMin = min(X(:));
xMax = max(X(:));
yMin = min(Y(:));
yMax = max(Y(:));
yBottom = yMin;

%% ===== READ CSV =====
Traw = readtable(csvFile);
T = Traw;

if any(strcmpi(T.Properties.VariableNames, 'status_code'))
    ok = (T.status_code == 2) | (T.status_code == 1);
    T = T(ok, :);
end

if isempty(T)
    error('No successful rows (status_code == 1 or 2) in %s', csvFile);
end

fprintf('Total rows in CSV: %d\n', height(Traw));
fprintf('Rows kept for sweep: %d\n', height(T));

dXR = diff(T.xr_exact);
if ~isempty(dXR)
    badIdx = find(dXR > 1.5 * median(dXR));
    if ~isempty(badIdx)
        fprintf('Gaps in xr_exact at rows:\n');
        disp([badIdx, T.xr_exact(badIdx), T.xr_exact(badIdx+1)]);
    end
end

[~, idx] = sort(T.xr_exact);
T = T(idx, :);
nFrames = height(T);

xrFrame = T.xr_exact(:);
yrFrame = T.yr_exact(:);

%% ===== CONVERT PRESSURE X/DELTA TO PHYSICAL X =====

xCp_phys = x_over_delta .* delta_0 + x_offset;

[xCp_phys, idxSort] = sort(xCp_phys);
dCpdx = dCpdx(idxSort);

[xCp_phys, iu] = unique(xCp_phys, 'stable');
dCpdx = dCpdx(iu);

if numel(xCp_phys) < 2
    error('Need at least two unique pressure x-locations after conversion.');
end

%% ===== INTERPOLATE PRESSURE DATA TO CORRELATION FRAMES =====
xrFrame_global = xrFrame + 7100; % so its in the global space.
dCpdxFrame = interp1(xCp_phys, dCpdx, xrFrame_global, 'linear', 'extrap');

% %% ===== LOAD CAMERA SPANS =====
% camBounds_phys = [];
% camNames = strings(0);
% overlapBounds_phys = [];
% 
% if showCameraEdges || showOverlapRegions
%     C = load(cameraBoundaryFile, 'windowCenterCameras_mm');
%     wc = C.windowCenterCameras_mm;
% 
%     nCams = numel(wc.x1_mm);
%     camBounds_phys = nan(nCams, 2);   % [xmin xmax]
% 
%     for k = 1:nCams
%         x1 = wc.x1_mm{k};
%         camBounds_phys(k,1) = min(x1(:));
%         camBounds_phys(k,2) = max(x1(:));
%     end
% 
%     % sort by left edge in case cell order is not strictly left-to-right
%     [~, camOrder] = sort(camBounds_phys(:,1), 'ascend');
%     camBounds_phys = camBounds_phys(camOrder,:);
% 
%     if isfield(wc, 'cameras')
%         camNames = string(wc.cameras);
%         camNames = camNames(camOrder);
%     else
%         camNames = "cam" + string(1:nCams);
%         camNames = camNames(camOrder);
%     end
% 
%     fprintf('\nCamera x spans (physical mm):\n');
%     for k = 1:nCams
%         fprintf('%s : [%.3f, %.3f]\n', camNames(k), camBounds_phys(k,1), camBounds_phys(k,2));
%     end
% 
%     % Adjacent overlap regions
%     overlapList = [];
%     for k = 1:(nCams-1)
%         leftBound  = max(camBounds_phys(k,1),   camBounds_phys(k+1,1));
%         rightBound = min(camBounds_phys(k,2),   camBounds_phys(k+1,2));
% 
%         if rightBound > leftBound
%             overlapList = [overlapList; leftBound rightBound]; %#ok<AGROW>
%             fprintf('Overlap %s-%s : [%.3f, %.3f]\n', ...
%                 camNames(k), camNames(k+1), leftBound, rightBound);
%         else
%             fprintf('No overlap %s-%s\n', camNames(k), camNames(k+1));
%         end
%     end
% 
%     overlapBounds_phys = overlapList;
% end
%% ===== LOAD CAMERA SPANS =====
camBounds_phys = [];
camBounds_plot = [];
camNames_native = strings(0);
camNames_sorted = strings(0);
camOrder = [];
overlapBounds_phys = [];
overlapPairs = [];

if showCameraEdges || showOverlapRegions
    C = load(cameraBoundaryFile, 'windowCenterCameras_mm');
    wc = C.windowCenterCameras_mm;

    nCams = numel(wc.x1_mm);

    % --- Native/original order from the struct ---
    camBounds_native = nan(nCams, 2);   % [xmin xmax] in original cell order
    camMid_native    = nan(nCams, 1);

    for k = 1:nCams
        x1 = wc.x1_mm{k};
        camBounds_native(k,1) = min(x1(:));
        camBounds_native(k,2) = max(x1(:));
        camMid_native(k) = mean(camBounds_native(k,:));
    end

    if isfield(wc, 'cameras')
        camNames_native = string(wc.cameras(:));
    else
        camNames_native = "cam" + string((1:nCams)');
    end

    % --- Spatial order used ONLY for overlap logic ---
    [~, camOrder] = sort(camBounds_native(:,1), 'ascend');
    camBounds_sorted = camBounds_native(camOrder,:);
    camNames_sorted  = camNames_native(camOrder);

    fprintf('\nCamera x spans in ORIGINAL struct order:\n');
    for k = 1:nCams
        fprintf('cell %d (%s) : [%.3f, %.3f] mm | xmid = %.3f mm\n', ...
            k, camNames_native(k), ...
            camBounds_native(k,1), camBounds_native(k,2), camMid_native(k));
    end

    fprintf('\nCamera x spans in SPATIAL left-to-right order:\n');
    for k = 1:nCams
        origIdx = camOrder(k);
        fprintf('x-order %d = cell %d (%s) : [%.3f, %.3f] mm\n', ...
            k, origIdx, camNames_sorted(k), ...
            camBounds_sorted(k,1), camBounds_sorted(k,2));
    end

    % Adjacent overlaps are computed in spatial order only
    overlapList = [];
    overlapPairs = strings(0,1);

    fprintf('\nAdjacent overlaps based on SPATIAL order:\n');
    for k = 1:(nCams-1)
        leftBound  = max(camBounds_sorted(k,1),   camBounds_sorted(k+1,1));
        rightBound = min(camBounds_sorted(k,2),   camBounds_sorted(k+1,2));

        idxA = camOrder(k);
        idxB = camOrder(k+1);

        if rightBound > leftBound
            overlapList = [overlapList; leftBound rightBound]; %#ok<AGROW>
            overlapPairs(end+1,1) = sprintf('cell %d (%s) <-> cell %d (%s)', ...
                idxA, camNames_native(idxA), idxB, camNames_native(idxB));
            fprintf('%s : [%.3f, %.3f] mm\n', overlapPairs(end), leftBound, rightBound);
        else
            fprintf('No overlap: cell %d (%s) <-> cell %d (%s)\n', ...
                idxA, camNames_native(idxA), idxB, camNames_native(idxB));
        end
    end

    % For plotting the edge lines, use ORIGINAL order so identity is preserved
    camBounds_phys = camBounds_native;
    overlapBounds_phys = overlapList;
end
%% ===== NORMALIZATION FOR PLOTTING ONLY =====
Xplot     = (X + (7100 - 7650) ) ./ BL_0;
Yplot     = Y ./ BL_0;
xrPlot    = (xrFrame + (7100 - 7650) ) ./ BL_0;
yrPlot    = yrFrame ./ BL_0;
xCpPlot   = (xCp_phys) ./ BL_0;

xMinPlot  = min(Xplot(:));
xMaxPlot  = max(Xplot(:));
yMinPlot  = min(Yplot(:));
yMaxPlot  = max(Yplot(:));
yBotPlot  = yMinPlot;

if ~isempty(camBounds_phys)
    camBounds_plot = (camBounds_phys + (7100 - 7650))./ BL_0;
else
    camBounds_plot = [];
end

if ~isempty(overlapBounds_phys)
    overlapBounds_plot = (overlapBounds_phys+ (7100 - 7650)) ./ BL_0;
else
    overlapBounds_plot = [];
end

%% ===== FIGURE / VIDEO SETUP =====
fig = figure('Color', 'w', ...
    'Units', 'pixels', ...
    'Position', figPos, ...
    'Resize', 'off', ...
    'Renderer', 'opengl');
set(fig, 'ToolBar', 'none', 'MenuBar', 'none');

gifFile = fullfile(baseDir, gifName);
mp4File = fullfile(baseDir, mp4Name);

v = VideoWriter(mp4File, 'MPEG-4');
v.FrameRate = round(1 / safeDelay);
open(v);

missingList = [];

%% ===== FRAME LOOP =====
for i = 1:nFrames
    clf(fig);

    xr   = xrFrame(i);   % physical x used for file lookup
    yr   = yrFrame(i);   % physical y
    xrN  = xrPlot(i);    % normalized for plotting
    yrN  = yrPlot(i);

    matName = sprintf('%s%.1f_yref%.1f.mat', corrPrefix, xr, fixedYref);
    matPath = fullfile(baseDir, matName);

    if ~isfile(matPath)
        warning('Missing MAT file for frame %d (x_ref = %.3f). Skipping.', i, xr);
        missingList(end+1,:) = [i xr]; %#ok<AGROW>
        continue;
    end

    S = load(matPath, 'R_s');
    Rplane = S.R_s;

    % ---------------------------
    % Top subplot: correlation
    % ---------------------------
    ax1 = subplot(2,1,1, 'Parent', fig);
    hold(ax1, 'on');
    disableDefaultInteractivity(ax1);

    imagesc(ax1, Xplot(1,:), Yplot(:,1), Rplane);
    set(ax1, 'YDir', 'normal');
    clim(ax1, corrClim);
    colormap(ax1, redblue(10));

    cb = colorbar(ax1);
    cb.Color = 'k';
    cb.Label.String = corrLabel;
    cb.Label.Color = 'k';

    axis(ax1, 'image');
    xlim(ax1, [xMinPlot xMaxPlot]);
    ylim(ax1, [yMinPlot yMaxPlot]);

    % --- Overlap regions (shaded) ---
    if showOverlapRegions && ~isempty(overlapBounds_plot)
        for k = 1:size(overlapBounds_plot,1)
            xL = overlapBounds_plot(k,1);
            xR = overlapBounds_plot(k,2);

            patch(ax1, [xL xR xR xL], [yMinPlot yMinPlot yMaxPlot yMaxPlot], ...
                overlapColor, ...
                'FaceAlpha', overlapFaceAlpha, ...
                'EdgeColor', 'none');
        end
    end

    % Re-draw image on top if you prefer less tinting:
    % imagesc(ax1, Xplot(1,:), Yplot(:,1), Rplane); set(ax1, 'YDir', 'normal');

    % --- Camera edges ---
    if showCameraEdges && ~isempty(camBounds_plot)
        for k = 1:size(camBounds_plot,1)
            xline(ax1, camBounds_plot(k,1), '--', ...
                'Color', cameraEdgeColor, ...
                'LineWidth', cameraEdgeWidth);

            xline(ax1, camBounds_plot(k,2), '--', ...
                'Color', cameraEdgeColor, ...
                'LineWidth', cameraEdgeWidth);
        end
    end

    % Vertical tracker line at current x_ref
    plot(ax1, [xrN xrN], [yrN yBotPlot], '--', ...
        'Color', lineColor, ...
        'LineWidth', 1.5);

    % Horizontal floor line
    yline(ax1, 0 / BL_0, '--', 'Color', 'k', 'LineWidth', 1);

    xlabel(ax1, '$x / \delta$', 'Color', 'k', 'Interpreter','latex');
    ylabel(ax1, '$y / \delta$', 'Color', 'k', 'Interpreter','latex');
    title(ax1, sprintf('%s | x_{ref}/delta = %.3f, y_{ref}/delta = %.3f', corrLabel, xrN, yrN), ...
        'Color', 'k', 'FontSize', 11);

    ax1.Color = 'k';
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 13;
    ax1.LineWidth = 0.8;
    ax1.TickDir = 'out';

    % ---------------------------
    % Bottom subplot: dCp/dx
    % ---------------------------
    ax2 = subplot(2,1,2, 'Parent', fig);
    hold(ax2, 'on');
    disableDefaultInteractivity(ax2);

    plot(ax2, xCpPlot, dCpdx, 'w-o', 'LineWidth', 1.8);
    xline(ax2, xrN, '-', 'Color', lineColor, 'LineWidth', 1.2);

    % plot(ax2, xrN, dCpdxFrame(i), 'o', ...
    %     'MarkerSize', 7, ...
    %     'MarkerFaceColor', lineColor, ...
    %     'MarkerEdgeColor', 'w', ...
    %     'LineWidth', 1.0);

    xlabel(ax2, '$x / \delta$', 'Color', 'k', 'Interpreter','latex');
    ylabel(ax2, '$ (dC_p/dx) \delta$', 'Color', 'k', 'Interpreter','latex');
    title(ax2, 'Pressure-gradient evolution', 'Color', 'k', 'FontSize', 11);

    xlim(ax2,  [xMinPlot xMaxPlot]);
    ylim(ax2, dCpdxYLim);

    grid(ax2, 'on');
    ax2.GridColor = [1 1 1] * 0.25;
    ax2.Color = 'k';
    ax2.XColor = 'k';
    ax2.YColor = 'k';
    ax2.FontName = 'Times New Roman';
    ax2.FontSize = 13;
    ax2.LineWidth = 0.8;
    ax2.TickDir = 'out';

    drawnow;

    frame = getframe(fig);
    writeVideo(v, frame);

    [im, map] = rgb2ind(frame.cdata, 256);
    if i == 1
        imwrite(im, map, gifFile, 'gif', ...
            'LoopCount', inf, 'DelayTime', frameDelay);
    else
        imwrite(im, map, gifFile, 'gif', ...
            'WriteMode', 'append', 'DelayTime', frameDelay);
    end

    if mod(i,10) == 0
        fprintf('Frame %d / %d written.\n', i, nFrames);
    end
end

%% ===== CLEAN UP =====
if ~isempty(missingList)
    fprintf('Total missing MAT files: %d\n', size(missingList,1));
    disp(array2table(missingList, 'VariableNames', {'rowIndex','xr_exact'}));
end

close(v);
close(fig);

fprintf('GIF saved to:\n%s\n', gifFile);
fprintf('MP4 saved to:\n%s\n', mp4File);
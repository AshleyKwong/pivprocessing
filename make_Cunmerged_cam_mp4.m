% make_Cunmerged_cam_mp4.m
% GIF + MP4 of unmerged camera-native covariance planes:
% C_unmerged_cam<cam>_<corrType>_xrefXXX.X_yrefYYY.Y.mat

clear; clc; close all;

%% ===== USER INPUTS ====================================================
baseDir = 'C:\Users\ak1u24\Downloads\sweep_x_yref4.0';
camNum = 5;
corrType = 'uu';              % 'uu', 'vv', 'uv', or 'vu'
fixedYref = 4.0;              % y_ref used in the sweep folder

gifName = 'Cunmerged_cam4_sweep_yref4p0.gif';
mp4Name = 'Cunmerged_cam4_sweep_yref4p0.mp4';
frameDelay = 0.12;            % GIF delay [s]
safeDelay = 0.08;             % MP4 frame delay equivalent

% Color scaling options
useManualCLim = false;
manualCLim = [-0.5 0.5];      % used only if useManualCLim = true

% Plot options
showWallLine = true;
showRefLine = true;
showRefMarker = true;
figPos = [100 100 1600 700];
%% ======================================================================

pattern = sprintf('C_unmerged_cam%d_%s_xref*_yref%.1f.mat', camNum, corrType, fixedYref);
F = dir(fullfile(baseDir, pattern));
if isempty(F)
    error('No files found matching pattern:\n%s', fullfile(baseDir, pattern));
end

% Parse xref values from filenames and sort
xrefs = nan(numel(F),1);
for k = 1:numel(F)
    tok = regexp(F(k).name, 'xref([-+]?\d*\.?\d+)_yref', 'tokens', 'once');
    if ~isempty(tok)
        xrefs(k) = str2double(tok{1});
    end
end
[~, idx] = sort(xrefs);
F = F(idx);
xrefs = xrefs(idx);
nFrames = numel(F);

fprintf('Found %d unmerged covariance files.\n', nFrames);
fprintf('x_ref range: [%.3f, %.3f] mm\n', min(xrefs), max(xrefs));

% Determine global color limits if requested
if useManualCLim
    climUse = manualCLim;
else
    maxAbs = 0;
    for k = 1:nFrames
        S = load(fullfile(baseDir, F(k).name), 'C_cam_debug');
        C = double(S.C_cam_debug);
        maxAbs = max(maxAbs, max(abs(C(:)), [], 'omitnan'));
    end
    climUse = [-maxAbs, maxAbs];
end
fprintf('Using color limits: [%.4g, %.4g]\n', climUse(1), climUse(2));

% Load first file for grid extents
S0 = load(fullfile(baseDir, F(1).name), 'C_cam_debug', 'x_cam_debug', 'y_cam_debug', 'xr', 'yr');
X0 = double(S0.x_cam_debug);
Y0 = double(S0.y_cam_debug);
xMin = min(X0(:));
xMax = max(X0(:));
yMin = min(Y0(:));
yMax = max(Y0(:));

fig = figure('Color','k', ...
    'Units','pixels', ...
    'Position', figPos, ...
    'Resize','off', ...
    'Renderer','opengl');
set(fig, 'ToolBar','none', 'MenuBar','none');

gifFile = fullfile(baseDir, gifName);
mp4File = fullfile(baseDir, mp4Name);

v = VideoWriter(mp4File, 'MPEG-4');
v.FrameRate = round(1/safeDelay);
open(v);

for i = 1:nFrames
    clf(fig);
    ax = axes('Parent', fig, 'Color','k', 'XColor','w', 'YColor','w');
    hold(ax, 'on');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;
    disableDefaultInteractivity(ax);

    S = load(fullfile(baseDir, F(i).name), 'C_cam_debug', 'x_cam_debug', 'y_cam_debug', 'xr', 'yr');
    Cplane = double(S.C_cam_debug);
    X = double(S.x_cam_debug);
    Y = double(S.y_cam_debug);
    xr = double(S.xr);
    yr = double(S.yr);

    imagesc(ax, X(1,:), Y(:,1), Cplane);
    set(ax, 'YDir', 'normal');
    clim(ax, climUse);
    colormap(ax, redblue(32));
    cb = colorbar(ax);
    cb.Label.String = sprintf('C_{%s} (unmerged, Cam %d)', corrType, camNum);
    cb.Label.Color = 'w';
    cb.Color = 'w';

    axis(ax, 'image');
    xlim(ax, [xMin xMax]);
    ylim(ax, [yMin yMax]);

    if showRefLine
        plot(ax, [xr xr], [yr yMin], '--', 'Color', [1 0.2 0.2], 'LineWidth', 1.5);
    end
    if showRefMarker
        plot(ax, xr, yr, 'wo', 'MarkerFaceColor', [1 0.2 0.2], 'MarkerSize', 5);
    end
    if showWallLine
        yline(ax, 0, '--', 'Color', 'k', 'LineWidth', 1);
    end

    xlabel(ax, 'x (mm)', 'Color','w');
    ylabel(ax, 'y (mm)', 'Color','w');
    title(ax, sprintf('Unmerged Cam %d covariance | C_{%s} | x_{ref}=%.2f, y_{ref}=%.2f', ...
        camNum, corrType, xr, yr), 'Color','w', 'FontSize', 12);

    ax.TickDir = 'out';
    ax.LineWidth = 0.8;
    ax.Position = [0.08 0.12 0.78 0.82];

    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);

    [im, map] = rgb2ind(frame.cdata, 256);
    if i == 1
        imwrite(im, map, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay);
    else
        imwrite(im, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
    end

    if mod(i,10) == 0 || i == nFrames
        fprintf('Frame %d / %d written.\n', i, nFrames);
    end
end

close(v);
close(fig);
fprintf('GIF saved to:\n%s\n', gifFile);
fprintf('MP4 saved to:\n%s\n', mp4File);

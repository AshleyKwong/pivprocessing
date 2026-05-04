% make_Ruu_gif_yref25p7.m
% GIF + MP4 of R_uu correlation planes with red core tracker

clear; clc; close all;

% ===== USER INPUTS =====
baseDir    = 'C:\Users\ak1u24\Downloads\sweep_x_yref1.5';
gridFile   = fullfile(baseDir, 'grid.mat');
csvFile    = fullfile(baseDir, 'integral_length_summary_dx.csv');
gifName    = 'Ruu_sweep_yref1p5.gif';
frameDelay = 0.12;   % seconds per frame for GIF
corrType = "R_{uu}"; 
corrPrefix = 'R_uu_xref';   % or R_uv_xref / R_vv_xref
fixedYref  = 1.5;          % this folder's y_ref
mp4Name    = 'Ruu_sweep_yref2p6.mp4';
safeDelay  = 0.08;          % used for MP4 frame rate

delta_HWA1 = 130; %140; %mm this is at the point of measurement for HWA
delta_0 = 110; % mm this is at the inlet

% =======================

% --- Load grid ---
G = load(gridFile, 'worldX_merged', 'worldY_merged');
X = G.worldX_merged;
Y = G.worldY_merged;

xMin = min(X(:));
xMax = max(X(:));
yMin = min(Y(:));
yMax = max(Y(:));
yBottom = yMin;

% % --- Read CSV and keep successful rows (status_code == 2) ---
% T = readtable(csvFile);
% if any(strcmpi(T.Properties.VariableNames,'status_code'))
%     ok = (T.status_code == 2) | (T.status_code == 1);
%     T  = T(ok, :);
% end
% 
% if isempty(T)
%     error('No successful rows (status_code==2) in %s', csvFile);
% end
% fprintf('Total rows in CSV: %d\n', height(readtable(csvFile)));
% fprintf('Rows with status_code==2: %d\n', height(T));
% 
% % Check for a hole in the xr_exact sequence
% dXR = diff(T.xr_exact);
% badIdx = find(dXR > 1.5 * median(dXR));   % heuristic gap detector
% if ~isempty(badIdx)
%     fprintf('Gaps in xr_exact at rows:\n');
%     disp([badIdx, T.xr_exact(badIdx), T.xr_exact(badIdx+1)]);
% end
% 
% 
% % Sort rows by xr_exact to make the sweep orderly
% [~, idx] = sort(T.xr_exact);
% T        = T(idx, :);
% nFrames  = height(T);

% --- Read CSV (optional) ---
if isfile(csvFile)
    T = readtable(csvFile);
    if any(strcmpi(T.Properties.VariableNames, 'status_code'))
        ok = (T.status_code == 2) | (T.status_code == 1);
        T = T(ok, :);
    end
    if isempty(T)
        warning('No successful rows (status_code==2) in %s. Scanning all MAT files instead.', csvFile);
        T = [];
    else
        fprintf('Total rows in CSV: %d\n', height(readtable(csvFile)));
        fprintf('Rows with status_code==2: %d\n', height(T));
        % Check for gaps
        dXR = diff(T.xr_exact);
        badIdx = find(dXR > 1.5 * median(dXR));
        if ~isempty(badIdx)
            fprintf('Gaps in xr_exact at rows:\n');
            disp([badIdx, T.xr_exact(badIdx), T.xr_exact(badIdx+1)]);
        end
        [~, idx] = sort(T.xr_exact);
        T = T(idx, :);
    end
else
    warning('CSV not found: %s\nFalling back to directory scan.', csvFile);
    T = [];
end

% --- Fallback: scan directory for matching MAT files ---
if isempty(T)
    pattern = fullfile(baseDir, sprintf('%s*_yref%.1f.mat', corrPrefix, fixedYref));
    files = dir(pattern);
    if isempty(files)
        error('No MAT files found matching pattern: %s', pattern);
    end
    xr_vals = zeros(numel(files), 1);
    for k = 1:numel(files)
        % Parse xr from filename, e.g. R_vv_xref123.4_yref25.7.mat
        % tok = regexp(files(k).name, [corrPrefix, '([\d.]+)_yref'], 'tokens');
        % With this (escape the prefix for use in regex):
        safePfx = regexptranslate('escape', corrPrefix);
        tok = regexp(files(k).name, [safePfx, '([-\d.]+)_yref'], 'tokens');
        fprintf('File: %s  →  xr parsed: %.4f\n', files(k).name, xr_vals(k));
        if ~isempty(tok)
            xr_vals(k) = str2double(tok{1}{1});
        end
    end
    [xr_vals, sidx] = sort(xr_vals);
    files = files(sidx);
    % Build a minimal table with xr_exact and yr_exact = fixedYref
    T = table(xr_vals, repmat(fixedYref, numel(xr_vals), 1), ...
        'VariableNames', {'xr_exact', 'yr_exact'});
    fprintf('CSV absent — found %d MAT files via directory scan.\n', height(T));
end

nFrames = height(T);

% --- Prepare figure (fixed pixel size) ---
fig = figure('Color','k', ...
             'Units','pixels', ...
             'Position',[100 100 1600 700], ...
             'Resize','off', ...
             'Renderer','opengl');
set(fig, 'ToolBar','none', 'MenuBar','none');

gifFile = fullfile(baseDir, gifName);
mp4File = fullfile(baseDir, mp4Name);

% --- Prepare MP4 writer ---
v = VideoWriter(mp4File, 'MPEG-4');
v.FrameRate = round(1/safeDelay);    % e.g. 12 fps
open(v);

missingList = [];   % place before the loop
gifInitialised = false;   % Add this BEFORE the loop

for i = 1:nFrames
    clf(fig);
    ax = axes('Parent', fig, 'Color','k', ...
              'XColor','w', 'YColor','w');
    hold(ax, 'on');
    ax.FontName  = 'Times New Roman';
    ax.FontSize  = 14;
    disableDefaultInteractivity(ax);
    xr = T.xr_exact(i) ;
    yr = T.yr_exact(i) ;      % physical ref point

    % Build MAT filename: R_uu_xref<xr>_yref25.7.mat
    matName = sprintf('%s%.1f_yref%.1f.mat', corrPrefix, xr, fixedYref);
    matPath = fullfile(baseDir, matName);
    if ~isfile(matPath)
        warning('Missing MAT file for frame %d (x_ref=%.3f). Skipping.', i, xr);
        missingList(end+1,:) = [i xr]; %#ok<AGROW>
        continue;
    end

    S      = load(matPath, 'R_s');    % 1789x3873
    Rplane = S.R_s;

    % Plot correlation plane
    imagesc(ax, X(1,:)./delta_HWA1, Y(:,1)./delta_HWA1, Rplane);
    set(ax, 'YDir', 'normal');
    clim(ax, [0 1]);
    colormap(ax, redblue(10));
    cb = colorbar(ax);
    cb.Label.String = 'R_{uu}';
    cb.Label.Color  = 'w';
    cb.Color        = 'w';
    axis(ax, 'image');

    % Axes limits & aspect based on worldY_merged only
    xlim(ax, [xMin xMax]./ delta_HWA1);
    ylim(ax, [yMin yMax]./ delta_HWA1);

    % Dashed red line from y_ref down to bottom of domain
    plot(ax, [xr xr], [yr yBottom], '--', ...
        'Color', [1 0.2 0.2], ...
        'LineWidth', 1.5);

    % Horizontal black line at y = 0 (floor)
    yline(ax, 0, '--', 'Color','k', 'LineWidth', 1);

    % xlabel(ax, 'x (mm)', 'Color','w');
    % ylabel(ax, 'y (mm)', 'Color','w');
    xlabel(ax, 'x / \delta (-)', 'Color','w');
    ylabel(ax, 'y / \delta (-)', 'Color','w');
    title(ax, sprintf('%s | x_{ref}=%.2f, y_{ref}=%.2f', corrType, xr, yr), ...
        'Color','w', 'FontSize', 10);

    % Tight framing, minimal whitespace
    ax.TickDir   = 'out';
    ax.LineWidth = 0.8;
    ax.Position  = [0.08 0.12 0.78 0.82];

        drawnow;

    % Capture frame once from the fixed-size figure
    frame = getframe(fig);          % struct with .cdata (h×w×3)
    writeVideo(v, frame);           % MP4

    % Convert same frame to indexed image for GIF
    [im, map] = rgb2ind(frame.cdata, 256);

    % if i == 1
    %     imwrite(im, map, gifFile, 'gif', ...
    %         'LoopCount', inf, 'DelayTime', frameDelay);
    % else
    %     imwrite(im, map, gifFile, 'gif', ...
    %         'WriteMode', 'append', 'DelayTime', frameDelay);
    % end
    % Inside the loop, replace the imwrite block:
    if ~gifInitialised
        imwrite(im, map, gifFile, 'gif', ...
            'LoopCount', inf, 'DelayTime', frameDelay);
        gifInitialised = true;
    else
        imwrite(im, map, gifFile, 'gif', ...
            'WriteMode', 'append', 'DelayTime', frameDelay);
    end

    if mod(i,10) == 0
        fprintf('Frame %d / %d written.\n', i, nFrames);
    end
end

if ~isempty(missingList)
    fprintf('Total missing MAT files: %d\n', size(missingList,1));
    disp(array2table(missingList, ...
        'VariableNames', {'rowIndex','xr_exact'}));
end
close(v);
close(fig);
fprintf('GIF saved to:\n%s\n', gifFile);
fprintf('MP4 saved to:\n%s\n', mp4File);
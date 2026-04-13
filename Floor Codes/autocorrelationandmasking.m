%% crop_calibration_im7.m
% Reads .im7 calibration images from each camera subfolder,
% centre-crops to the target sensor size, and saves as 16-bit TIFs.
%
% Requirements: readimx (LaVision DaVis toolbox) must be on the MATLAB path.
addpath(genpath('C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\calib_tools')) % this has the im7 reader
%% ── Paths ──────────────────────────────────────────────────────────────────
srcPath = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\';
outPath = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\y250_aoan04_aoafn04\';
%% ── Target crop sizes [rows, cols] — one row per camera ────────────────────
% targetSize = [3528, 5312;   % Camera 1
%     3472, 5312;   % Camera 2
%     3512, 5312;   % Camera 3
%     3536, 5312;   % Camera 4
%     3520, 5312];  % Camera 5
% targetSize = [4600, 5312;   % Camera 1
%     4600, 5312;   % Camera 2
%     4600, 5312;   % Camera 3
%     4600, 5312;   % Camera 4
%     4600, 5312];  % Camera 5
targetSize = [3536, 5312;   % Camera 1
    3536, 5312;   % Camera 2
    3536, 5312;   % Camera 3
    3536, 5312;   % Camera 4
    3536, 5312];  % Camera 5

%% ── Quick size check: raw PNGs in outPath ───────────────────────────────────
fprintf('\n===== RAW IMAGE SIZE CHECK =====\n');
fprintf('%-10s  %-30s  %s\n', 'Camera', 'File', 'Size [rows x cols]');
fprintf('%s\n', repmat('-',1,60));

camDirs = dir(outPath);
camDirs = camDirs([camDirs.isdir]);
camDirs = camDirs(~ismember({camDirs.name},{'.','..'}));
camDirs = camDirs(arrayfun(@(x) ~isempty(regexp(x.name,'^[Cc]am(era)?[1-5]$','once')), camDirs));

for i = 1:numel(camDirs)
    camDir  = fullfile(outPath, camDirs(i).name);
    pngs    = dir(fullfile(camDir, '*.png'));
    if isempty(pngs)
        fprintf('%-10s  No PNG files found\n', camDirs(i).name);
        continue
    end
    % Read just the first PNG — imfinfo is faster than imread (no pixel load)
    info = imfinfo(fullfile(camDir, pngs(1).name));
    fprintf('%-10s  %-30s  %d x %d\n', ...
        camDirs(i).name, pngs(1).name, info.Height, info.Width);
end
fprintf('================================\n\n');
% ── Auto-set targetSize from actual image dimensions ────────────────────────
targetSize = zeros(numel(camDirs), 2);

for i = 1:numel(camDirs)
    camDir = fullfile(outPath, camDirs(i).name);
    pngs   = dir(fullfile(camDir, '*.png'));
    if isempty(pngs)
        error('No PNG found in %s — cannot set targetSize for camera %d.', camDir, i);
    end
    info = imfinfo(fullfile(camDir, pngs(1).name));
    targetSize(i,:) = [info.Height, info.Width];
end

fprintf('targetSize set automatically:\n');
for i = 1:numel(camDirs)
    fprintf('  Camera %d:  [%d x %d]\n', i, targetSize(i,1), targetSize(i,2));
end
fprintf('\n');

%% ── Discover camera subdirectories ─────────────────────────────────────────
d = dir(srcPath);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));

% Keep only folders named 'camera1' through 'camera5'
isCamFolder = arrayfun(@(x) ~isempty(regexp(x.name, '^camera[1-5]$', 'once')), d);
d = d(isCamFolder);

% Keep only folders that actually contain .im7 files
hasIm7 = arrayfun(@(x) ~isempty(dir(fullfile(srcPath, x.name, '*.im7'))), d);
d = d(hasIm7);

fprintf('Found %d camera folder(s) with .im7 files.\n', numel(d));
% d = dir(srcPath);
% d = d([d.isdir]);
% d = d(~ismember({d.name}, {'.', '..'}));
%
% % Keep only folders that actually contain .im7 files
% hasIm7 = arrayfun(@(x) ~isempty(dir(fullfile(srcPath, x.name, '*.im7'))), d);
% d = d(hasIm7);
%
% fprintf('Found %d camera folder(s) with .im7 files.\n', numel(d));
%
% assert(numel(d) == size(targetSize, 1), ...
%     'Found %d folders with .im7 files but %d target sizes defined — check targetSize.', ...
%     numel(d), size(targetSize, 1));

%% ── Create top-level output folder ─────────────────────────────────────────
if ~exist(outPath, 'dir'), mkdir(outPath); end

%% ── Main loop ───────────────────────────────────────────────────────────────
for i = 1:numel(d)

    camSrc = fullfile(srcPath, d(i).name);
    camOut = fullfile(outPath, d(i).name);
    if ~exist(camOut, 'dir'), mkdir(camOut); end

    tRows = targetSize(i, 1);
    tCols = targetSize(i, 2);

    imFiles = dir(fullfile(camSrc, '*.im7'));
    if isempty(imFiles)
        fprintf('[SKIP] No .im7 files in %s\n', camSrc);
        continue
    end

    fprintf('[Cam %d]  %s  →  target [%d × %d]  (%d file(s))\n', ...
        i, d(i).name, tRows, tCols, numel(imFiles));

    for j = 1:numel(imFiles)

        % ── Read & convert ────────────────────────────────────────────────
        A   = readimx(fullfile(camSrc, imFiles(j).name));
        img = im7ToIm(A, 1);          % floating-point intensity
        img = uint16(img');            % transpose + cast (matches your snippet)

        % ── Centre-crop to target size ────────────────────────────────────
        [imgR, imgC] = size(img);

        % If the image is already smaller than the target on any axis, warn.
        if imgR < tRows || imgC < tCols
            warning('Cam %d, file %s: image [%d×%d] is SMALLER than target [%d×%d]. Skipping crop.', ...
                i, imFiles(j).name, imgR, imgC, tRows, tCols);
        end

        r1 = 1;
        r2 = tRows;
        c1 = max(1, floor((imgC - tCols)/2) + 1);   % columns still centred
        c2 = c1 + tCols - 1;

        imgCropped = img(r1:r2, c1:c2);

        % ── Save as 16-bit TIFF ───────────────────────────────────────────
        [~, baseName] = fileparts(imFiles(j).name);
        outFile = fullfile(camOut, [baseName, '.tif']);
        imwrite(imgCropped, outFile, 'tif');
        % ── Visualise crop for every file ────────────────────────────────
        figure(100 + i); clf;

        subplot(1,2,1);
        imagesc_clim(gca, img, 1, 99);
        colormap gray; axis image on;
        xlabel('Column (px)'); ylabel('Row (px)');
        hold on;

        % Shade the cropped-away bottom region
        x_shade = [1, imgC, imgC, 1, 1];
        y_shade  = [tRows, tRows, imgR, imgR, tRows];
        fill(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', ...
            'LineWidth', 1.5, 'DisplayName', 'Cropped away');

        % Shade left/right column margins if any
        if c1 > 1
            fill([1, c1, c1, 1, 1], [1, 1, imgR, imgR, 1], ...
                'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
        if c2 < imgC
            fill([c2, imgC, imgC, c2, c2], [1, 1, imgR, imgR, 1], ...
                'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end

        % Draw the crop boundary
        rectangle('Position', [c1, r1, tCols, tRows], ...
            'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '--');

        legend('show', 'Location', 'northwest', 'FontSize', 8);
        title(sprintf('Cam %d Frame %d — Original [%d × %d]', i, j, imgR, imgC), 'FontSize', 10);
        ax1 = gca;
        ax1.XLim = [1, imgC];
        ax1.YLim = [1, imgR];
        ax1.XTick = 0:500:imgC;
        ax1.YTick = 0:500:imgR;
        hold off;

        subplot(1,2,2);
        imagesc_clim(gca, imgCropped, 1, 99);
        colormap gray; axis image on;
        xlabel('Column (px)'); ylabel('Row (px)');
        title(sprintf('Cam %d Frame %d — Cropped [%d × %d]', i, j, tRows, tCols), 'FontSize', 10);
        ax2 = gca;
        ax2.XLim = [1, tCols];
        ax2.YLim = [1, tRows];
        ax2.XTick = 0:500:tCols;
        ax2.YTick = 0:500:tRows;

        sgtitle(sprintf('Cam %d Frame %d — keeping rows %d:%d, cols %d:%d', ...
            i, j, r1, r2, c1, c2), 'FontWeight', 'bold', 'FontSize', 11);

        % Save
        crop_fig_name = sprintf('Cam%d_frame%04d_crop_preview.png', i, j);
        exportgraphics(figure(100+i), fullfile(outPath, crop_fig_name), 'Resolution', 300);
        fprintf('  Crop preview saved: %s\n', crop_fig_name);

        fprintf('  [%d/%d] Saved: %s\n', j, numel(imFiles), outFile);
    end
end

fprintf('\nAll done. Cropped TIFs saved to:\n  %s\n', outPath);


%% ── Visualise cropped images ────────────────────────────────────────────────
figure('Name', 'Cropped Calibration Images', 'Units', 'normalized', ...
    'Position', [0.05 0.1 0.9 0.8]);

for i = 1:numel(d)

    camOut = fullfile(outPath, d(i).name);
    tifFiles = dir(fullfile(camOut, '*.tif'));

    % subplot(5, 1, i);
    figure(i);
    if isempty(tifFiles)
        text(0.5, 0.5, sprintf('No TIF\nCam %d', i), ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        axis off;
        continue
    end

    % Just display the first TIF from each camera folder
    img = imread(fullfile(camOut, tifFiles(1).name));
    imagesc(img);
    colormap gray;
    axis image off;
    title(sprintf('Cam %d\n[%d × %d]', i, size(img,1), size(img,2)), ...
        'FontSize', 9);
    colorbar;
end

sgtitle('Cropped Calibration Images', 'FontSize', 12, 'FontWeight', 'bold');


%% CAMERA 1 — RAW PIVTools vs Cropped Calibration
camera1_rawpivtools = imread(fullfile(outPath, 'camera1\Cam1_frame00001_A_raw.png'));
cropped_img_cal     = imread(fullfile(outPath, 'camera1\B00001.tif'));

% Build title strings with sizes
raw_sz  = size(camera1_rawpivtools);
crop_sz = size(cropped_img_cal);
raw_title  = sprintf('RAW PIVTools — Cam 1  [%d × %d]',  raw_sz(1),  raw_sz(2));
crop_title = sprintf('Cropped TIF  — Cam 1  [%d × %d]', crop_sz(1), crop_sz(2));

figure('Units','normalized','Position',[0.05 0.15 0.9 0.7]);

subplot(1,2,1);
imagesc(camera1_rawpivtools);
colormap gray; axis image on; colorbar;
gamma = 0.3;   % < 1 brightens, > 1 darkens — try 0.3–0.5 for particles
camera1_gamma = imadjust(camera1_rawpivtools, [], [], gamma);
imagesc(camera1_gamma);

title(raw_title, 'FontSize', 10);

subplot(1,2,2);
imagesc(cropped_img_cal);
colormap gray; axis image off; colorbar;
clim([0 100]);           % ← same scale so the two panels are comparable
title(crop_title, 'FontSize', 10);

sgtitle('RAW PIVTools vs Cropped Calibration — Camera 1', ...
    'FontSize', 13, 'FontWeight', 'bold');


%% autocorrelate to find the floor
% look in the outPath
d = dir(outPath);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));

% Pre-allocate results struct
floorResults = struct('cam', {}, 'roi_row', {}, 'roi_col', {}, ...
    'floor_y_A', {}, 'floor_y_B', {});

for i = 1:numel(d)
    camName  = d(i).name;                   % e.g. 'cam1'
    camDir   = fullfile(outPath, camName);

    % ── Load A and B frames ───────────────────────────────────────────
    fnA = fullfile(camDir, sprintf('Cam%d_frame00001_A_raw.png', i));
    fnB = fullfile(camDir, sprintf('Cam%d_frame00001_B_raw.png', i));

    if ~isfile(fnA) || ~isfile(fnB)
        warning('Missing images for %s — skipping.', camName);
        continue
    end

    imgA = imread(fnA);
    imgB = imread(fnB);

    % Convert to grayscale if RGB
    if size(imgA, 3) == 3; imgA = rgb2gray(imgA); end
    if size(imgB, 3) == 3; imgB = rgb2gray(imgB); end

    % ── Figure: 1×2 subplot ──────────────────────────────────────────
    hFig = figure('Name', sprintf('Cam %d — Select ROI', i), ...
        'NumberTitle', 'off');

    ax1 = subplot(1,2,1);
    % imagesc(imgA);
    imagesc_clim(ax1, imgA, 1, 99);
    colormap gray; axis image off;
    title(sprintf('Cam%d  Frame A', i));

    ax2 = subplot(1,2,2);
    % imagesc(imgB);
    imagesc_clim(ax2, imgB, 1, 99);
    colormap gray; axis image off;
    title(sprintf('Cam%d  Frame B', i));

    sgtitle(sprintf('Cam %d — Draw ROI box on Frame A to select floor region', i));

    % ── User draws ROI on Frame A ─────────────────────────────────────
    % Select axes A for drawing
    axes(ax1); %#ok<LAXES>
    disp(['Cam' num2str(i) ': Draw a rectangle around the floor reflection region on Frame A, then double-click to confirm.']);
    hROI = drawrectangle(ax1, 'Color', 'r', 'LineWidth', 1.5, ...
        'Label', 'Floor ROI');
    % Block until user double-clicks to confirm
    wait(hROI);

    % ROI position: [xmin ymin width height] in axes (image pixel) coords
    pos = hROI.Position;                    % [xmin, ymin, w, h]
    xmin = round(pos(1));
    ymin = round(pos(2));
    w    = round(pos(3));
    h    = round(pos(4));

    xmax = xmin + w - 1;
    ymax = ymin + h - 1;

    % Clamp to image bounds
    [imgH, imgW] = size(imgA);
    xmin = max(1, xmin);  xmax = min(imgW, xmax);
    ymin = max(1, ymin);  ymax = min(imgH, ymax);

    fprintf('Cam%d ROI (original pixel coords): rows [%d:%d], cols [%d:%d]\n', ...
        i, ymin, ymax, xmin, xmax);

    % Draw the same box on Frame B axes for visual reference
    drawrectangle(ax2, 'Position', pos, 'Color', 'cyan', ...
        'LineWidth', 1.5, 'Label', 'Same ROI', 'InteractionsAllowed', 'none');




% ── Interactive floor line selection ──────────────────────────────────
    fprintf('\nCam%d: Use zoom/pan freely, then place and drag the two floor points.\n', i);
    fprintf('  Double-click each point to confirm its position.\n');

    floor_confirmed = false;
    while ~floor_confirmed

        figure(200 + i); clf;
        imagesc_clim(gca, imgA, 1, 99);
        colormap gray; axis image off;
        title(sprintf('Cam %d — Zoom/pan, then place LEFT floor point (double-click to confirm)', i));

        % Enable zoom/pan toolbar so user can navigate before placing points
        zoom on;
        disp('  Zoom/pan to the LEFT end of the floor, then press any key to place point...');
        waitforbuttonpress;
        zoom off;

        % Place left point — fully draggable before double-click confirm
        hL = drawpoint(gca, 'Color', 'r', 'Label', 'Floor L', 'LabelVisible', 'hover');
        wait(hL);
        pt_L = hL.Position;

        title(sprintf('Cam %d — Zoom/pan, then place RIGHT floor point (double-click to confirm)', i));
        zoom on;
        disp('  Zoom/pan to the RIGHT end of the floor, then press any key to place point...');
        waitforbuttonpress;
        zoom off;

        % Place right point
        hR = drawpoint(gca, 'Color', 'b', 'Label', 'Floor R', 'LabelVisible', 'hover');
        wait(hR);
        pt_R = hR.Position;

        % Sort left to right by column
        pts = sortrows([pt_L; pt_R], 1);
        px_clicks = pts(:,1)';
        py_clicks = pts(:,2)';

        floor_left_input  = py_clicks(1);
        floor_right_input = py_clicks(2);
        col_left          = px_clicks(1);
        col_right         = px_clicks(2);

        % Draw the extrapolated floor line across the full image width
        x_extrap = [1, imgW];
        y_extrap = interp1(px_clicks, py_clicks, x_extrap, 'linear', 'extrap');

        hold on;
        plot(x_extrap, y_extrap, 'r--', 'LineWidth', 2, 'DisplayName', 'Floor estimate');
        hold off;
        title(sprintf('Cam %d — Floor line preview', i));

        reply = input('  Happy with floor line? y to proceed, n to redo: ', 's');
        if isempty(reply) || strcmpi(reply, 'y')
            floor_confirmed = true;
        else
            % Clear points and redraw for another attempt
            delete(hL); delete(hR);
        end
    end

    floor_search_tol = 2;   % ± px search window around expected floor in autocorr

    % Linear interpolation for expected floor row at any column
    floor_expected = @(col) interp1(px_clicks, py_clicks, col, 'linear', 'extrap');

    floorResults(i).floor_left_input  = floor_left_input;
    floorResults(i).floor_right_input = floor_right_input;
    floorResults(i).floor_col_left    = col_left;
    floorResults(i).floor_col_right   = col_right;

    fprintf('Cam%d: floor line confirmed — L=(%.0f, %.0f)  R=(%.0f, %.0f)\n', ...
        i, col_left, floor_left_input, col_right, floor_right_input);

    % ── Autocorrelation on the ROI crop — multiple strips ─────────────────
    cropA = double(imgA(ymin:ymax, xmin:xmax));
    cropB = double(imgB(ymin:ymax, xmin:xmax));

    nStrips   = floor(5312/6);
    stripW    = floor((xmax - xmin + 1) / nStrips);
    h         = ymax - ymin + 1;

    strip_col_centres       = zeros(1, nStrips);
    floor_y_A_strips        = nan(1, nStrips);
    floor_y_B_strips        = nan(1, nStrips);
    floor_y_combined_strips = nan(1, nStrips);
    strip_weights           = zeros(1, nStrips);

    for s = 1:nStrips
        c1s = (s-1)*stripW + 1;
        c2s = s*stripW;

        sA = cropA(:, c1s:c2s);
        sB = cropB(:, c1s:c2s);

        strip_col_centres(s) = xmin + (c1s + c2s)/2 - 1;

        % Skip flat strips
        if std(sA(:)) < 1 || std(sB(:)) < 1
            fprintf('  Cam%d strip %d: skipped (flat/no texture)\n', i, s);
            continue
        end

        % Expected floor row for this strip, converted to expected lag index
        floor_y_expected_s   = floor_expected(strip_col_centres(s));
        floor_y_roi_expected = floor_y_expected_s - ymin + 1;       % ROI coords
        lag_expected         = 2 * (floor_y_roi_expected - h/2);    % lag in pixels

        % ── Autocorrelations and cross-correlation ─────────────────────────
        acA  = normxcorr2(sA, sA);
        acB  = normxcorr2(sB, sB);
        acAB = normxcorr2(sA, sB);

        ac1d_A  = mean(acA,  2);
        ac1d_B  = mean(acB,  2);
        ac1d_AB = mean(acAB, 2);

        nR      = length(ac1d_A);
        zeroLag = ceil(nR / 2);

        % Mask central peak — adaptive to ROI height
        halfMask = max(5, round(0.05 * h));
        maskIdx  = max(1, zeroLag-halfMask) : min(nR, zeroLag+halfMask);

        % Search window from user constraint
        lag_expected_idx = zeroLag + round(lag_expected);
        searchMin = max(1,  lag_expected_idx - floor_search_tol);
        searchMax = min(nR, lag_expected_idx + floor_search_tol);

        % ── Constrained max after masking, then subpixel fit ──────────────
        [pkA,  idxA]  = find_reflection_max(ac1d_A,  maskIdx, searchMin, searchMax);
        [pkB,  idxB]  = find_reflection_max(ac1d_B,  maskIdx, searchMin, searchMax);
        [pkAB, idxAB] = find_reflection_max(ac1d_AB, maskIdx, searchMin, searchMax);

        lagA  = subpixel_peak(ac1d_A,  idxA,  5) - zeroLag;
        lagB  = subpixel_peak(ac1d_B,  idxB,  5) - zeroLag;
        lagAB = subpixel_peak(ac1d_AB, idxAB, 5) - zeroLag;

        % ── Weighted combination ───────────────────────────────────────────
        weights = [pkA, pkB, pkAB];
        lags    = [lagA, lagB, lagAB];

        valid_w = isfinite(lags) & isfinite(weights) & weights > 0;
        if sum(valid_w) == 0
            fprintf('  Cam%d strip %d: no valid peak found (expected lag idx=%d, window=[%d %d])\n', ...
                i, s, lag_expected_idx, searchMin, searchMax);
            continue
        end

        lag_combined = sum(weights(valid_w) .* lags(valid_w)) / sum(weights(valid_w));

        floor_y_A_strips(s)        = ymin + (h/2 + lagA/2)         - 1;
        floor_y_B_strips(s)        = ymin + (h/2 + lagB/2)         - 1;
        floor_y_combined_strips(s) = ymin + (h/2 + lag_combined/2) - 1;
        strip_weights(s)           = mean(weights(valid_w));

        fprintf('  Cam%d strip %d: expected floor=%.0f px  found=%.1f px  (lagA=%.2f lagB=%.2f lagAB=%.2f)\n', ...
            i, s, floor_y_expected_s, floor_y_combined_strips(s), lagA, lagB, lagAB);
    end

    % ── Consistency check A vs B, downweight bad strips ───────────────────
    diff_AB    = abs(floor_y_A_strips - floor_y_B_strips);
    bad_strips = diff_AB > 8;
    if any(bad_strips)
        warning('Cam%d: strips %s have A/B disagreement > 8 px — downweighted.', ...
            i, num2str(find(bad_strips)));
        strip_weights(bad_strips) = strip_weights(bad_strips) * 0.2;
    end

    good = isfinite(floor_y_combined_strips) & strip_weights > 0;

    % ── Weighted least-squares line fit ───────────────────────────────────
    if sum(good) >= 2
        W  = diag(strip_weights(good));
        X  = [strip_col_centres(good)', ones(sum(good), 1)];
        y  = floor_y_combined_strips(good)';
        pC = ((X'*W*X) \ (X'*W*y))';
    elseif sum(good) == 1
        warning('Cam%d: only 1 good strip, using horizontal floor.', i);
        pC = [0, floor_y_combined_strips(good)];
    else
        warning('Cam%d: no good strips found!', i);
        pC = [0, NaN];
    end

    pA = polyfit_safe(strip_col_centres, floor_y_A_strips, isfinite(floor_y_A_strips));
    pB = polyfit_safe(strip_col_centres, floor_y_B_strips, isfinite(floor_y_B_strips));

    floor_y_img_A = polyval(pA, (xmin+xmax)/2);
    floor_y_img_B = polyval(pB, (xmin+xmax)/2);
    floor_y_img_C = polyval(pC, (xmin+xmax)/2);

    fprintf('Cam%d  floor at ROI centre: A=%.1f  B=%.1f  combined=%.1f px\n', ...
        i, floor_y_img_A, floor_y_img_B, floor_y_img_C);

    % ── Store results ──────────────────────────────────────────────────────
    floorResults(i).cam                     = camName;
    floorResults(i).roi_row                 = [ymin, ymax];
    floorResults(i).roi_col                 = [xmin, xmax];
    floorResults(i).floor_y_A              = floor_y_img_A;
    floorResults(i).floor_y_B              = floor_y_img_B;
    floorResults(i).floor_y_combined       = floor_y_img_C;
    floorResults(i).floor_line_A           = pA;
    floorResults(i).floor_line_B           = pB;
    floorResults(i).floor_line_combined    = pC;
    floorResults(i).strip_col_centres      = strip_col_centres;
    floorResults(i).floor_y_A_strips       = floor_y_A_strips;
    floorResults(i).floor_y_B_strips       = floor_y_B_strips;
    floorResults(i).floor_y_combined_strips = floor_y_combined_strips;
    floorResults(i).strip_weights          = strip_weights;

    % ── Overlay on image subplots ──────────────────────────────────────────
    x_full = [1, imgW];

    hold(ax1, 'on');
    plot(ax1, strip_col_centres, floor_y_A_strips, 'r-+', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Strip estimates A');
    plot(ax1, x_full, polyval(pA, x_full), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Fit A');
    plot(ax1, x_full, polyval(pC, x_full), 'y-',  'LineWidth', 2,   'DisplayName', 'Floor combined');
    legend(ax1, 'show', 'Location', 'northwest');
    hold(ax1, 'off');

    hold(ax2, 'on');
    plot(ax2, strip_col_centres, floor_y_B_strips, 'c-+', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Strip estimates B');
    plot(ax2, x_full, polyval(pB, x_full), 'c--', 'LineWidth', 1.5, 'DisplayName', 'Fit B');
    plot(ax2, x_full, polyval(pC, x_full), 'y-',  'LineWidth', 2,   'DisplayName', 'Floor combined');
    legend(ax2, 'show', 'Location', 'northwest');
    hold(ax2, 'off');

    % ── Diagnostic figure ─────────────────────────────────────────────────
    figure('Name', sprintf('Cam %d — Floor strip estimates', i), 'NumberTitle', 'off');

    subplot(1,2,1);
    hold on;
    plot(strip_col_centres, floor_y_A_strips,        'r-+', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Frame A');
    plot(strip_col_centres, floor_y_B_strips,        'b-+', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Frame B');
    plot(strip_col_centres, floor_y_combined_strips, 'k-+', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Combined');
    if any(bad_strips)
        plot(strip_col_centres(bad_strips), floor_y_combined_strips(bad_strips), ...
            'kx', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', 'Downweighted');
    end
    x_plot = linspace(xmin, xmax, 200);
    plot(x_plot, polyval(pA, x_plot), 'r-',  'LineWidth', 1,   'DisplayName', 'Fit A');
    plot(x_plot, polyval(pB, x_plot), 'b--', 'LineWidth', 1,   'DisplayName', 'Fit B');
    plot(x_plot, polyval(pC, x_plot), 'k-',  'LineWidth', 2,   'DisplayName', 'Fit combined');
    % Show user-defined floor line for reference
    plot(x_extrap, y_extrap, 'm--', 'LineWidth', 1.5, 'DisplayName', 'User floor hint');
    set(gca, 'YDir', 'reverse');
    xlabel('Column (image px)');  ylabel('Floor row (image px)');
    title(sprintf('Cam %d — floor row per strip', i));
    legend('show', 'Location', 'best');
    grid on;

    subplot(1,2,2);
    bar(strip_col_centres, diff_AB, 'FaceColor', [0.8 0.3 0.3]);
    yline(8, 'k--', 'LineWidth', 1.5, 'Label', '8 px threshold');
    xlabel('Strip column centre (px)');
    ylabel('|floor\_A - floor\_B|  (px)');
    title('A vs B consistency per strip');
    grid on;

    sgtitle(sprintf('Cam %d — Multi-strip floor detection', i), 'FontWeight', 'bold');

    %     % ── Autocorrelation on the ROI crop (Frame A) ─────────────────────
    %     cropA = double(imgA(ymin:ymax, xmin:xmax));
    %     cropB = double(imgB(ymin:ymax, xmin:xmax));
    %     % 2-D normalised autocorrelation (xcorr2 via FFT)
    %     acorrA = normxcorr2(cropA, cropA);   % already in [-1, 1], no further scaling needed
    %     acorrB = normxcorr2(cropB, cropB);
    %
    %     % ── Find floor reflection lag from autocorrelation ────────────────
    %     % The floor reflection creates a secondary peak symmetric about zero lag.
    %     % Collapse to 1-D by averaging across columns (x direction).
    %     acorrA_1d = mean(acorrA, 2);   % [2*h-1 x 1]
    %     acorrB_1d = mean(acorrB, 2);
    %
    %     nRows    = size(acorrA, 1);    % = 2*h - 1
    %     zeroLag  = ceil(nRows / 2);    % index of zero-lag peak
    %
    %     % Mask out the central peak region (±5% of height) to find secondary peak
    %     halfMask = max(3, round(0.05 * h));
    %     maskIdx  = (zeroLag - halfMask):(zeroLag + halfMask);
    %
    %     acorrA_masked = acorrA_1d;
    %     acorrA_masked(maskIdx) = -Inf;
    %
    %     [~, peakIdx_A] = max(acorrA_masked);
    %     lag_A = peakIdx_A - zeroLag;   % signed lag in pixels
    %
    %     % The floor is at lag/2 below the ROI top (reflection symmetry)
    %     floor_y_roi_A = h / 2 + lag_A / 2;           % within ROI
    %     floor_y_img_A = ymin + floor_y_roi_A - 1;     % in original image
    %
    %     % Repeat for B
    %     acorrB_masked = mean(acorrB, 2);
    %     acorrB_masked(maskIdx) = -Inf;
    %     [~, peakIdx_B] = max(acorrB_masked);
    %     lag_B = peakIdx_B - zeroLag;
    %     floor_y_roi_B = h / 2 + lag_B / 2;
    %     floor_y_img_B = ymin + floor_y_roi_B - 1;
    %
    %     fprintf('Cam%d  floor estimate → Frame A: y = %.1f px,  Frame B: y = %.1f px (original image)\n', ...
    %         i, floor_y_img_A, floor_y_img_B);
    %
    %     % ── Store results ─────────────────────────────────────────────────
    %     floorResults(i).cam        = camName;
    %     floorResults(i).roi_row    = [ymin, ymax];   % ← preserved original-image y range
    %     floorResults(i).roi_col    = [xmin, xmax];
    %     floorResults(i).floor_y_A  = floor_y_img_A;
    %     floorResults(i).floor_y_B  = floor_y_img_B;
    %
    %     % ── Overlay floor line on both subplots ───────────────────────────
    %     yline(ax1, floor_y_img_A, 'r--', 'LineWidth', 1.5, ...
    %         'Label', sprintf('Floor %.0f px', floor_y_img_A));
    %     yline(ax2, floor_y_img_B, 'c--', 'LineWidth', 1.5, ...
    %         'Label', sprintf('Floor %.0f px', floor_y_img_B));
    %
    %
    %
    %     %── Visualise autocorrelation results ────────────────────────────────────────
    %     hAC = figure('Name', sprintf('Cam %d — Autocorrelation', i), 'NumberTitle', 'off');
    %
    %     nRows   = length(acorrA_1d);
    %     lagAxis = (-(nRows-1)/2 : (nRows-1)/2);   % symmetric lag axis in pixels
    %
    %     % ── 1D column-averaged autocorrelation ───────────────────────────────────────
    %     subplot(2, 3, [1 2]);
    %     plot(lagAxis, acorrA_1d, 'r',  'LineWidth', 1.5); hold on;
    %     plot(lagAxis, acorrB_1d, 'b--','LineWidth', 1.5);
    %     xline(0, 'k:', 'LineWidth', 1);
    %     xline( lag_A, 'r:', 'Label', sprintf('lag_A = %d px', lag_A), 'LineWidth', 1);
    %     xline( lag_B, 'b:', 'Label', sprintf('lag_B = %d px', lag_B), 'LineWidth', 1);
    %     xlabel('Row lag  \tau_y  (pixels)');
    %     ylabel('Normalised autocorrelation');
    %     title(sprintf('Cam%d — Column-averaged 1D autocorrelation', i));
    %     legend('Frame A', 'Frame B', 'Location', 'northeast');
    %     grid on;
    %
    %     % ── 2D autocorrelation maps ───────────────────────────────────────────────────
    %     subplot(2, 3, 4);
    %     imagesc(acorrA);
    %     colormap(gca, 'parula'); colorbar;
    %     %axis image;
    %     hold on;
    %     % Mark zero-lag centre
    %     [cy, cx] = deal(ceil(size(acorrA,1)/2), ceil(size(acorrA,2)/2));
    %     plot(cx, cy,        'w+', 'MarkerSize', 12, 'LineWidth', 2);         % zero lag
    %     plot(cx, cy + lag_A,'r+', 'MarkerSize', 12, 'LineWidth', 2);         % secondary peak
    %     title('Frame A — 2D autocorr');
    %     xlabel('\tau_x'); ylabel('\tau_y');
    %
    %     subplot(2, 3, 5);
    %     imagesc(acorrB);
    %     colormap(gca, 'parula'); colorbar;
    %     % axis image;
    %     hold on;
    %     plot(cx, cy,        'w+', 'MarkerSize', 12, 'LineWidth', 2);
    %     plot(cx, cy + lag_B,'b+', 'MarkerSize', 12, 'LineWidth', 2);
    %     title('Frame B — 2D autocorr');
    %     xlabel('\tau_x'); ylabel('\tau_y');
    %
    %     % ── Original crop with floor line overlaid ────────────────────────────────────
    %     subplot(2, 3, 3);
    %     imagesc(cropA); colormap(gca, 'gray'); %axis image off;
    %     hold on;
    %     yline(floor_y_roi_A, 'r--', 'LineWidth', 1.5, ...
    %         'Label', sprintf('Floor ~%.0f px (ROI)', floor_y_roi_A));
    %     title(sprintf('Frame A crop  [rows %d:%d]', ymin, ymax));
    %
    %     subplot(2, 3, 6);
    %     imagesc(cropB); colormap(gca, 'gray'); %axis image off;
    %     hold on;
    %     yline(floor_y_roi_B, 'b--', 'LineWidth', 1.5, ...
    %         'Label', sprintf('Floor ~%.0f px (ROI)', floor_y_roi_B));
    %     title(sprintf('Frame B crop  [rows %d:%d]', ymin, ymax));
    %
    %     sgtitle(sprintf('Cam %d — Autocorrelation diagnostics  (floor y_{img} A=%.0f, B=%.0f px)', ...
    %         i, floor_y_img_A, floor_y_img_B));
end

disp('Done. Floor results stored in floorResults struct.');
% ── Save floorResults struct ─────────────────────────────────────────────────
floor_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
floor_save_path = fullfile(outPath, ['floorResults_', floor_timestamp, '.mat']);
save(floor_save_path, 'floorResults');
fprintf('\nfloorResults saved to:\n  %s\n', floor_save_path);

%% ── Overlay autocorrelation floor estimates on calibration images ────────────
for i = 1:numel(d)

    % Load calibration TIF for this camera
    camOut   = fullfile(outPath, d(i).name);
    tifFiles = dir(fullfile(camOut, '*.tif'));

    figure(300+i); clf;

    if isempty(tifFiles)
        text(0.5, 0.5, sprintf('No calibration TIF found — Cam %d', i), ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        axis off;
        continue
    end

    calImg = imread(fullfile(camOut, tifFiles(1).name));
    [calH, calW] = size(calImg);

    ax = gca;
    imagesc_clim(ax, calImg, 1, 99);
    colormap(ax, 'gray');
    axis(ax, 'image'); axis(ax, 'on');
    xlabel(ax, 'Column (px)', 'FontSize', 9);
    ylabel(ax, 'Row (px)', 'FontSize', 9);
    ax.XLim = [1, calW];
    ax.YLim = [1, calH];
    ax.XTick = 0:500:calW;
    ax.YTick = 0:500:calH;
    hold(ax, 'on');

    % Floor line from autocorrelation (combined fit)
    x_full  = [1, calW];
    pC      = floorResults(i).floor_line_combined;
    pA      = floorResults(i).floor_line_A;
    pB      = floorResults(i).floor_line_B;

    plot(ax, x_full, polyval(pC, x_full), 'y-',  'LineWidth', 2,   'DisplayName', 'Floor combined (autocorr)');
    plot(ax, x_full, polyval(pA, x_full), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Floor frame A');
    plot(ax, x_full, polyval(pB, x_full), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Floor frame B');

    % User floor hint
    px_hint      = [floorResults(i).floor_col_left,   floorResults(i).floor_col_right];
    py_hint      = [floorResults(i).floor_left_input,  floorResults(i).floor_right_input];
    x_extrap_hint = [1, calW];
    y_extrap_hint = interp1(px_hint, py_hint, x_extrap_hint, 'linear', 'extrap');
    % plot(ax, x_extrap_hint, y_extrap_hint, 'm--', 'LineWidth', 1.5, 'DisplayName', 'User floor hint');

    % Strip estimates
    good = isfinite(floorResults(i).floor_y_combined_strips) & ...
           floorResults(i).strip_weights > 0;
    plot(ax, floorResults(i).strip_col_centres(good), ...
             floorResults(i).floor_y_combined_strips(good), ...
        'm+', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'Strip estimates');

    legend(ax, 'show', 'Location', 'northwest', 'FontSize', 8);
    hold(ax, 'off');

    title(sprintf('Cam %d — calibration image  |  floor combined=%.1f px  A=%.1f px  B=%.1f px', ...
        i, polyval(pC, calW/2), polyval(pA, calW/2), polyval(pB, calW/2)), 'FontSize', 10);
    % Save figure as PNG
    fig_save_name = sprintf('Cam%d_floor_calibration_overlay_%s.png', i, floor_timestamp);
    fig_save_path = fullfile(outPath, fig_save_name);
    exportgraphics(figure(i), fig_save_path, 'Resolution', 300);
    fprintf('Saved: %s\n', fig_save_path);
end

%% drawing masks.
max_height = targetSize(i, 1);  % Camera 4 is already max, no padding needed
img_width  = 5312;

for i = 1:5
    orig_height = targetSize(i, 1);
    pad_top     = (max_height - orig_height) / 2;
    floorResults(i).pad_top = pad_top;

    % Floor y at left and right edges, in original image coordinates
    floor_y_L = floorAnchors(i).y_mean(1);   % left  edge (col 1)
    floor_y_R = floorAnchors(i).y_mean(2);   % right edge (col 5312)

    % Shift into padded image space
    adjusted_y_L = floor_y_L + pad_top;
    adjusted_y_R = floor_y_R + pad_top;

    floorResults(i).adjusted_y_L = adjusted_y_L;
    floorResults(i).adjusted_y_R = adjusted_y_R;

    % Per-column floor row: linear interpolation from left to right edge
    cols          = 1:img_width;
    floor_row_col = adjusted_y_L + (adjusted_y_R - adjusted_y_L) * (cols - 1) / (img_width - 1);
    floor_row_col = round(floor_row_col);                       % integer rows
    floor_row_col = min(max(floor_row_col, 1), max_height);     % clamp to image bounds

    floorResults(i).floor_row_col = floor_row_col;

    fprintf("\tCam %d | floor L = %.1f, R = %.1f (padded) | tilt = %+.1f px over %d cols\n", ...
        i, adjusted_y_L, adjusted_y_R, adjusted_y_R - adjusted_y_L, img_width);

    % Build mask: 1 above the tilted floor line, 0 at and below it.
    % Vectorised: compare each row index to the per-column threshold.
    [~, R] = meshgrid(cols, 1:max_height);                 % R(r,c) = r
    adjusted_mask = double(R < floor_row_col);             % 1 strictly above, 0 at/below
    floorResults(i).adjusted_mask = adjusted_mask;
end

% Plotting loop
for i = 1:numel(d)
    camName = d(i).name;
    camDir  = fullfile(outPath, camName);

    % Load frame A
    fnA = fullfile(camDir, sprintf('Cam%d_frame00001_A_raw.png', i));
    fnB = fullfile(camDir, sprintf('Cam%d_frame00001_B_raw.png', i));

    if ~isfile(fnA) || ~isfile(fnB)
        warning('Missing images for %s — skipping.', camName);
        continue
    end

    imgA = imread(fnA);

    % Masked image using adjusted (tilted) mask
    masked_adj = double(imgA) .* floorResults(i).adjusted_mask;

    figure(i);
    ax1 = subplot(1,2,1);
    imagesc_clim(ax1, imgA, 1, 99);
    colormap(ax1, 'gray');
    hold(ax1, 'on');

    % Draw the tilted floor line between the two anchor points
    plot(ax1, [1 size(imgA,2)], ...
              [floorResults(i).adjusted_y_L, floorResults(i).adjusted_y_R], ...
              'b--', 'LineWidth', 2, 'DisplayName', ...
              sprintf('Tilted floor (L=%.1f, R=%.1f)', ...
                      floorResults(i).adjusted_y_L, floorResults(i).adjusted_y_R));
    legend(ax1, 'show', 'Location', 'northwest');
    title(ax1, sprintf('Camera %d — Original + tilted floor line', i));
    axis(ax1, 'image'); axis(ax1, 'off');
    hold(ax1, 'off');

    ax2 = subplot(1,2,2);
    imagesc_clim(ax2, masked_adj, 1, 99);
    colormap(ax2, 'gray');
    title(ax2, sprintf('Camera %d — Masked (pad top = %d px)', i, floorResults(i).pad_top));
    axis(ax2, 'image'); axis(ax2, 'off');

    sgtitle(sprintf('Camera %d  |  orig height: %d  |  pad: %d px each side  |  tilt: %+.1f px', ...
        i, targetSize(i,1), floorResults(i).pad_top, ...
        floorResults(i).adjusted_y_R - floorResults(i).adjusted_y_L));
end
%%
for i = 1:numel(d)
    camName = d(i).name;
    camDir  = fullfile(outPath, camName);

    % Load frame A
    fnA = fullfile(camDir, sprintf('Cam%d_frame00001_A_raw.png', i));
    fnB = fullfile(camDir, sprintf('Cam%d_frame00001_B_raw.png', i));
    if ~isfile(fnA) || ~isfile(fnB)
        warning('Missing images for %s — skipping.', camName);
        continue
    end
    imgA = imread(fnA);

    % Masked image using adjusted mask
    masked_adj = double(imgA) .* floorResults(i).adjusted_mask;

    figure(i);

    ax1 = subplot(1,2,1);
    imagesc_clim(ax1, imgA, 1, 99);
    colormap(ax1, 'gray');

    title(sprintf('Cam%d  Frame A', i));
    hold(ax1, 'on');
    yline(ax1, floorResults(i).adjusted_row, ...
        'b--', 'LineWidth', 2, 'DisplayName', ...
        sprintf('Adjusted floor (row %.1f)', floorResults(i).adjusted_row));

    legend(ax1, 'show', 'Location', 'northwest');
    title(ax1, sprintf('Camera %d — Original + floor line', i));
    axis(ax1, 'image'); axis(ax1, 'off');
    hold(ax1, 'off');

    ax2 = subplot(1,2,2);
    imagesc_clim(ax2, masked_adj, 1, 99);
    colormap(ax2, 'gray');
    title(ax2, sprintf('Camera %d — Masked (pad top = %d px)', i, floorResults(i).pad_top));
    axis(ax2, 'image'); axis(ax2, 'off');

    sgtitle(sprintf('Camera %d  |  orig height: %d  |  pad: %d px each side', ...
        i, targetSize(i,1), floorResults(i).pad_top));
end


%% interactive polygon mask preview (example: Camera 2)
camToAdjust = 5;

camName = d(camToAdjust).name;
camDir  = fullfile(outPath, camName);
fnA     = fullfile(camDir, sprintf('Cam%d_frame00001_A_raw.png', camToAdjust));

if ~isfile(fnA)
    warning('Missing image for %s — cannot draw polygon mask.', camName);
else
    imgA = imread(fnA);
    if size(imgA,3) == 3
        imgA = rgb2gray(imgA);
    end

    [imgH, imgW] = size(imgA);
    donePoly = false;

    while ~donePoly
        figure(300); clf;
        imagesc_clim(gca, imgA, 1, 99);
        colormap gray;
        axis image off;
        title(sprintf(['Camera %d — draw 4-point polygon to MASK.\n' ...
            'Double-click to finish polygon.'], camToAdjust));

        hPoly = drawpolygon(gca, 'Color', 'r', 'LineWidth', 1.5);
        wait(hPoly);

        pts = round(hPoly.Position);

        % enforce exactly 4 points
        if size(pts,1) ~= 4
            warning('Please select exactly 4 points. You selected %d.', size(pts,1));
            continue
        end

        % keep original clicked points for checking
        pts_raw = pts;

        % snap out-of-bounds clicks to nearest valid image edge
        pts(:,1) = min(max(pts(:,1), 1), imgW);   % x in [1, imgW]
        pts(:,2) = min(max(pts(:,2), 1), imgH);   % y in [1, imgH]

        if any(pts(:) ~= pts_raw(:))
            fprintf('Some polygon points were outside the image and were snapped to the nearest edge.\n');
        end

        % build logical mask: true inside polygon = masked region
        mask = poly2mask(pts(:,1), pts(:,2), imgH, imgW);

        % preview masked image
        masked_img = double(imgA) .* double(~mask);

        figure(301); clf;

        ax1 = subplot(1,2,1);
        imagesc_clim(ax1, imgA, 1, 99);
        colormap(ax1, 'gray');
        axis(ax1, 'image'); axis(ax1, 'off');
        hold(ax1, 'on');
        plot(ax1, [pts(:,1); pts(1,1)], [pts(:,2); pts(1,2)], ...
            'r-', 'LineWidth', 2);
        plot(ax1, pts(:,1), pts(:,2), 'yo', ...
            'MarkerFaceColor', 'y', 'MarkerSize', 6);
        title(ax1, sprintf('Camera %d — Raw + polygon', camToAdjust));
        hold(ax1, 'off');

        ax2 = subplot(1,2,2);
        imagesc_clim(ax2, masked_img, 1, 99);
        colormap(ax2, 'gray');
        axis(ax2, 'image'); axis(ax2, 'off');
        title(ax2, sprintf('Camera %d — Polygon-masked preview', camToAdjust));

        sgtitle(sprintf('Camera %d polygon mask preview', camToAdjust));

        reply = input('Accept this polygon mask? y/n [y]: ', 's');
        if isempty(reply) || strcmpi(reply, 'y')
            donePoly = true;
        end
    end

    floorResults(camToAdjust).manual_polygon_points = double(pts);
    floorResults(camToAdjust).manual_mask = logical(mask);
    floorResults(camToAdjust).mask_mode = 'polygon';

    fprintf('Stored manual polygon mask for Camera %d.\n', camToAdjust);
end
%% saving the masks in the right mat file
% maskHeight = 3536;
% maskWidth  = 5312;
%% saving the masks in the right mat file
maskHeight = targetSize(1,1);
maskWidth  = targetSize(1,2);

for i = 1:numel(d)

    if isfield(floorResults(i), 'mask_mode') && strcmp(floorResults(i).mask_mode, 'polygon')
        pts = floorResults(i).manual_polygon_points;

        % safety clamp in case dimensions changed
        pts(:,1) = min(max(round(pts(:,1)), 1), maskWidth);
        pts(:,2) = min(max(round(pts(:,2)), 1), maskHeight);

        mask = poly2mask(pts(:,1), pts(:,2), maskHeight, maskWidth);

        polyStruct.index = int64(0);
        polyStruct.name  = 'Polygon 1';
        polyStruct.points = double(pts);
        polygons = {polyStruct};

        fprintf('Using manual polygon mask for Cam %d.\n', i);

    else
        mask = false(maskHeight, maskWidth);

        floor_row = round(floorResults(i).floor_y_A) + 5;
        floor_row = max(1, min(maskHeight, floor_row));

        mask(floor_row:end, :) = true;

        polyStruct.index = int64(0);
        polyStruct.name  = 'Polygon 1';
        polyStruct.points = [0,         double(floor_row); ...
            maskWidth, double(floor_row); ...
            maskWidth, double(maskHeight); ...
            0,         double(maskHeight)];
        polygons = {polyStruct};

        fprintf('Using horizontal floor mask for Cam %d.\n', i);
    end

    saveName = fullfile(outPath, sprintf('mask_loop=0_Cam%d.mat', i));
    save(saveName, 'mask', 'polygons');
    fprintf('Saved: %s | nnz(mask) = %d\n', saveName, nnz(mask));

    % preview saved mask on raw image
    camName = d(i).name;
    camDir  = fullfile(outPath, camName);
    fnA     = fullfile(camDir, sprintf('Cam%d_frame00001_A_raw.png', i));

    if ~isfile(fnA)
        warning('Missing image for %s — skipping plot.', camName);
        continue
    end

    imgA = imread(fnA);
    if size(imgA,3) == 3
        imgA = rgb2gray(imgA);
    end

    masked_img = double(imgA) .* double(~mask);

    figure(400 + i); clf;

    ax1 = subplot(1,2,1);
    imagesc_clim(ax1, imgA, 1, 99);
    colormap(ax1, 'gray');
    axis(ax1, 'image'); axis(ax1, 'off');
    hold(ax1, 'on');

    if isfield(floorResults(i), 'mask_mode') && strcmp(floorResults(i).mask_mode, 'polygon')
        pts = floorResults(i).manual_polygon_points;
        plot(ax1, [pts(:,1); pts(1,1)], [pts(:,2); pts(1,2)], ...
            'r-', 'LineWidth', 2);
        plot(ax1, pts(:,1), pts(:,2), 'yo', ...
            'MarkerFaceColor', 'y', 'MarkerSize', 6);
        title(ax1, sprintf('Camera %d — Raw + polygon mask', i));
    else
        yline(ax1, floor_row, 'r--', 'LineWidth', 2, ...
            'DisplayName', sprintf('Floor row %d', floor_row));
        legend(ax1, 'show', 'Location', 'northwest');
        title(ax1, sprintf('Camera %d — Raw + floor line', i));
    end
    hold(ax1, 'off');

    ax2 = subplot(1,2,2);
    imagesc_clim(ax2, masked_img, 1, 99);
    colormap(ax2, 'gray');
    axis(ax2, 'image'); axis(ax2, 'off');
    title(ax2, sprintf('Camera %d — Saved mask preview', i));

    sgtitle(sprintf('Camera %d saved mask preview', i));
end

%% copy one loop's mask template to other loops
sourceLoop = 0;          % loop to copy FROM
targetLoops = 1:13;      % loops to copy TO
camsToCopy = 1:5;        % cameras to copy

sourceFolder = fullfile(outPath, sprintf('loop=%d_data', sourceLoop));

if ~exist(sourceFolder, 'dir')
    error('Source folder does not exist: %s', sourceFolder);
end

for loop_id = targetLoops

    destFolder = fullfile(outPath, sprintf('loop=%d_data', loop_id));
    if ~exist(destFolder, 'dir')
        mkdir(destFolder);
    end

    for cam = camsToCopy

        sourceFile = fullfile(sourceFolder, ...
            sprintf('mask_loop=%d_Cam%d.mat', sourceLoop, cam));

        if ~isfile(sourceFile)
            warning('Missing source file: %s (skipping)', sourceFile);
            continue
        end

        destFile = fullfile(destFolder, ...
            sprintf('mask_loop=%d_Cam%d.mat', loop_id, cam));

        copyfile(sourceFile, destFile);

        fprintf('Copied: %s -> %s\n', sourceFile, destFile);
    end
end


function [pk_val, pk_idx] = find_reflection_max(ac1d, maskIdx, searchMin, searchMax)
    ac_search = -Inf(size(ac1d));
    ac_search(searchMin:searchMax) = ac1d(searchMin:searchMax);
    ac_search(maskIdx) = -Inf;
    [pk_val, pk_idx] = max(ac_search);
    if ~isfinite(pk_val)
        pk_val = NaN;  pk_idx = NaN;
    end
end

function pk_sub = subpixel_peak(ac1d, pk_int, npts)
    if isnan(pk_int)
        pk_sub = NaN;  return
    end
    half = floor(npts/2);
    idx  = (pk_int - half):(pk_int + half);
    idx  = idx(idx >= 1 & idx <= length(ac1d));
    if length(idx) < 3
        pk_sub = pk_int;  return
    end
    x     = idx(:);
    y     = ac1d(idx);
    valid = isfinite(y);
    if sum(valid) < 3
        pk_sub = pk_int;  return
    end
    p      = polyfit(x(valid), y(valid), 2);
    pk_sub = -p(2) / (2*p(1));
    if pk_sub < idx(1) || pk_sub > idx(end)
        pk_sub = pk_int;
    end
end

function p = polyfit_safe(x, y, good)
    if sum(good) >= 2
        p = polyfit(x(good), y(good), 1);
    elseif sum(good) == 1
        p = [0, y(good)];
    else
        p = [0, NaN];
    end
end
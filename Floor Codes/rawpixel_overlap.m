% search for autocorrelation between images to ensure that the overlaps are
% where they say they are.

%% ── User settings ────────────────────────────────────────────────────────────
clear; clc; close all;
nFrame   = 100;
baseDir  = 'C:\Users\ak1u24\Downloads\PG_fixedmergandcal\loop=06\uncalibrated_piv\150\';
camNames = {'Cam1','Cam2','Cam3','Cam4','Cam5'};
nCams    = numel(camNames);

% Calibration file paths — point to each case's .mat file
calib_path_A = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\calib_caseA_20260412_114104.mat';
calib_path_B = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\calib_caseB_20260412_114104.mat';
calib_path_C = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\calib_caseC_20260412_114104.mat';

%% ── Load calibrations ────────────────────────────────────────────────────────
tmp = load(calib_path_A);  calib_A = tmp.calib_A;
tmp = load(calib_path_B);  calib_B = tmp.calib_B;
tmp = load(calib_path_C);  calib_C = tmp.calib_C;

calibs      = {calib_B};
calib_names = {'Case 1'};
nCalibs     = numel(calibs);

fprintf('Loaded %d calibration cases.\n', nCalibs);

%% ── Load PIV fields for requested frame ──────────────────────────────────────
fileName = sprintf('%05d.mat', nFrame);
piv = struct();

for i = 1:nCams
    filePath = fullfile(baseDir, camNames{i}, 'instantaneous', fileName);

    if ~isfile(filePath)
        error('Frame %d not found for %s:\n  %s', nFrame, camNames{i}, filePath);
    end

    raw = load(filePath);
    p   = raw.piv_result(end);   % always take final pass

    piv(i).ux          = p.ux;
    piv(i).uy          = p.uy;
    piv(i).x_px        = p.win_ctrs_x;     % [ny x nx] grid of x pixel coords
    piv(i).y_px        = p.win_ctrs_y;     % [ny x nx] grid of y pixel coords
    piv(i).b_mask      = p.b_mask;
    piv(i).nan_mask    = p.nan_mask;
    piv(i).n_windows   = p.n_windows;
    piv(i).window_size = p.window_size;
    piv(i).cam         = camNames{i};
    piv(i).valid       = ~p.b_mask & ~p.nan_mask;

    fprintf('Loaded %s — frame %05d  |  grid: %d x %d  |  window: %d x %d px\n', ...
        camNames{i}, nFrame, p.n_windows(1), p.n_windows(2), ...
        p.window_size(1), p.window_size(2));
end

fprintf('\nAll cameras and calibrations loaded. Ready for coordinate transform.\n');
%% ── Transform pixel coords to mm for each camera and each calibration ────────
% polyvaln expects [n x 2] input: [x_px(:), y_px(:)]
% fit_cx1 maps pixel -> mm in x1 direction
% fit_cx2 maps pixel -> mm in x2 direction

for c = 1:nCalibs

    calib = calibs{c};

    fprintf('\nTransforming coordinates — %s\n', calib_names{c});

    % ── Transform PIV window centres from pixel to mm ────────────────────────────
    for i = 1:nCams

        % Build full 2D grids from 1D position vectors
        % x_px (1x882)  — horizontal, column direction
        % y_px (1x1765) — vertical,   row direction
        [x_grid, y_grid] = meshgrid(piv(i).x_px, piv(i).y_px);   % [1765 x 882]

        piv(i).x_grid = x_grid;
        piv(i).y_grid = y_grid;

        % Flatten to [n x 2] in confirmed order [col_px, row_px]
        px_in   = [x_grid(:), y_grid(:)];
        grid_sz = size(x_grid);   % [1765 x 882]

        for c = 1:nCalibs

            calib = calibs{c};

            x1_flat = polyvaln(calib.fit_cx1{i}.p, px_in);
            x2_flat = polyvaln(calib.fit_cx2{i}.p, px_in);

            piv(i).calib(c).x1_mm = reshape(x1_flat, grid_sz);
            piv(i).calib(c).x2_mm = reshape(x2_flat, grid_sz);
            piv(i).calib(c).name  = calib_names{c};

            fprintf('%s — %s: x1=[%.1f, %.1f] mm  x2=[%.1f, %.1f] mm\n', ...
                camNames{i}, calib_names{c}, ...
                min(x1_flat), max(x1_flat), ...
                min(x2_flat), max(x2_flat));
        end

        fprintf('\nTransform complete.\n');
    end

end
%% ── Find overlap regions between adjacent cameras ────────────────────────────
for c = 1:nCalibs

    fprintf('\n=== %s ===\n', calib_names{c});

    for i = 1:nCams-1

        % x1 extent of each camera
        x1_max_L = max(piv(i).calib(c).x1_mm(:));
        x1_min_R = min(piv(i+1).calib(c).x1_mm(:));
        x1_min_L = min(piv(i).calib(c).x1_mm(:));
        x1_max_R = max(piv(i+1).calib(c).x1_mm(:));

        % Overlap bounds in mm
        x1_overlap_min = x1_min_R;
        x1_overlap_max = x1_max_L;

        if x1_overlap_min >= x1_overlap_max
            warning('%s — %s vs %s: no overlap found!', ...
                calib_names{c}, camNames{i}, camNames{i+1});
            continue
        end

        overlap_width = x1_overlap_max - x1_overlap_min;
        fprintf('  %s <-> %s: overlap x1=[%.2f, %.2f] mm  (%.2f mm wide)\n', ...
            camNames{i}, camNames{i+1}, x1_overlap_min, x1_overlap_max, overlap_width);

        % Extract PIV vectors in overlap region for left camera
        mask_L = piv(i).calib(c).x1_mm   >= x1_overlap_min & ...
            piv(i).calib(c).x1_mm   <= x1_overlap_max & ...
            piv(i).valid;

        % Extract PIV vectors in overlap region for right camera
        mask_R = piv(i+1).calib(c).x1_mm >= x1_overlap_min & ...
            piv(i+1).calib(c).x1_mm <= x1_overlap_max & ...
            piv(i+1).valid;

        % Store overlap info
        piv(i).calib(c).overlap(i).x1_bounds    = [x1_overlap_min, x1_overlap_max];
        piv(i).calib(c).overlap(i).mask_L        = mask_L;
        piv(i).calib(c).overlap(i).mask_R        = mask_R;
        piv(i).calib(c).overlap(i).x1_mm_L       = piv(i).calib(c).x1_mm(mask_L);
        piv(i).calib(c).overlap(i).x2_mm_L       = piv(i).calib(c).x2_mm(mask_L);
        piv(i).calib(c).overlap(i).ux_L          = piv(i).ux(mask_L);
        piv(i).calib(c).overlap(i).uy_L          = piv(i).uy(mask_L);
        piv(i).calib(c).overlap(i).x1_mm_R       = piv(i+1).calib(c).x1_mm(mask_R);
        piv(i).calib(c).overlap(i).x2_mm_R       = piv(i+1).calib(c).x2_mm(mask_R);
        piv(i).calib(c).overlap(i).ux_R          = piv(i+1).ux(mask_R);
        piv(i).calib(c).overlap(i).uy_R          = piv(i+1).uy(mask_R);

        fprintf('    Points in overlap — %s: %d  |  %s: %d\n', ...
            camNames{i},   nnz(mask_L), ...
            camNames{i+1}, nnz(mask_R));
    end
end

fprintf('\nOverlap extraction complete.\n');

%% ── Find overlap column indices in pixel-space PIV grid ─────────────────────
% Use mm only to identify which columns overlap, then work in pixel space

for c = 1:nCalibs

    fprintf('\n=== %s ===\n', calib_names{c});

    for i = 1:nCams-1

        ov = piv(i).calib(c).overlap(i);
        x1_lo = ov.x1_bounds(1);
        x1_hi = ov.x1_bounds(2);

        % ── Find overlap column indices for left camera ────────────────────
        % x1_mm is [ny x nx] — take mean along rows to get 1D x1 per column
        x1_cols_L = mean(piv(i).calib(c).x1_mm, 1);       % [1 x nx_L]
        col_idx_L = find(x1_cols_L >= x1_lo & x1_cols_L <= x1_hi);

        % ── Find overlap column indices for right camera ───────────────────
        x1_cols_R = mean(piv(i+1).calib(c).x1_mm, 1);     % [1 x nx_R]
        col_idx_R = find(x1_cols_R >= x1_lo & x1_cols_R <= x1_hi);

        fprintf('  %s <-> %s\n', camNames{i}, camNames{i+1});
        fprintf('    %s overlap cols: %d:%d  (%d columns)\n', ...
            camNames{i},   col_idx_L(1), col_idx_L(end), numel(col_idx_L));
        fprintf('    %s overlap cols: %d:%d  (%d columns)\n', ...
            camNames{i+1}, col_idx_R(1), col_idx_R(end), numel(col_idx_R));

        % ── Extract overlap sub-grids in pixel space ───────────────────────
        ux_L_ov = piv(i).ux(:,   col_idx_L);    % [ny x n_cols_L]
        uy_L_ov = piv(i).uy(:,   col_idx_L);
        ux_R_ov = piv(i+1).ux(:, col_idx_R);    % [ny x n_cols_R]
        uy_R_ov = piv(i+1).uy(:, col_idx_R);

        % Pixel coords of overlap windows
        x_px_L_ov = piv(i).x_grid(:,   col_idx_L);  % [ny x n_cols_L]
        y_px_L_ov = piv(i).y_grid(:,   col_idx_L);
        x_px_R_ov = piv(i+1).x_grid(:, col_idx_R);  % [ny x n_cols_R]
        y_px_R_ov = piv(i+1).y_grid(:, col_idx_R);

        % Corresponding mm coords for reference
        x1_L_ov = piv(i).calib(c).x1_mm(:,   col_idx_L);
        x1_R_ov = piv(i+1).calib(c).x1_mm(:, col_idx_R);
        x2_L_ov = piv(i).calib(c).x2_mm(:,   col_idx_L);
        x2_R_ov = piv(i+1).calib(c).x2_mm(:, col_idx_R);

        fprintf('    Sub-grid size — L: [%d x %d]  R: [%d x %d]\n', ...
            size(ux_L_ov,1), size(ux_L_ov,2), ...
            size(ux_R_ov,1), size(ux_R_ov,2));

        % ── Store ─────────────────────────────────────────────────────────
        piv(i).calib(c).overlap(i).col_idx_L  = col_idx_L;
        piv(i).calib(c).overlap(i).col_idx_R  = col_idx_R;
        piv(i).calib(c).overlap(i).ux_L_ov    = ux_L_ov;
        piv(i).calib(c).overlap(i).uy_L_ov    = uy_L_ov;
        piv(i).calib(c).overlap(i).ux_R_ov    = ux_R_ov;
        piv(i).calib(c).overlap(i).uy_R_ov    = uy_R_ov;
        piv(i).calib(c).overlap(i).x_px_L_ov  = x_px_L_ov;
        piv(i).calib(c).overlap(i).y_px_L_ov  = y_px_L_ov;
        piv(i).calib(c).overlap(i).x_px_R_ov  = x_px_R_ov;
        piv(i).calib(c).overlap(i).y_px_R_ov  = y_px_R_ov;
        piv(i).calib(c).overlap(i).x1_L_ov    = x1_L_ov;
        piv(i).calib(c).overlap(i).x1_R_ov    = x1_R_ov;
        piv(i).calib(c).overlap(i).x2_L_ov    = x2_L_ov;
        piv(i).calib(c).overlap(i).x2_R_ov    = x2_R_ov;
    end
end

fprintf('\nPixel-space overlap extraction complete.\n');

%% ── Compare overlap vector fields in pixel space ─────────────────────────────
for c = 1:nCalibs

    fprintf('\n=== %s ===\n', calib_names{c});

    for i = 1:nCams-1

        ov = piv(i).calib(c).overlap(i);

        % ── Trim to same number of columns ────────────────────────────────
        nCols = min(size(ov.ux_L_ov, 2), size(ov.ux_R_ov, 2));
        nRows = size(ov.ux_L_ov, 1);

        % Left camera: take rightmost columns (highest x1, closest to right cam)
        ux_L = ov.ux_L_ov(:, end-nCols+1:end);
        uy_L = ov.uy_L_ov(:, end-nCols+1:end);

        % Right camera: take leftmost columns (lowest x1, closest to left cam)
        ux_R = ov.ux_R_ov(:, 1:nCols);
        uy_R = ov.uy_R_ov(:, 1:nCols);

        % ── Valid mask — neither camera has NaN or masked points ──────────
        valid = isfinite(ux_L) & isfinite(ux_R) & ...
            isfinite(uy_L) & isfinite(uy_R);

        fprintf('  %s <-> %s  [%d x %d grid]  valid: %d / %d (%.1f%%)\n', ...
            camNames{i}, camNames{i+1}, nRows, nCols, ...
            nnz(valid), numel(valid), 100*nnz(valid)/numel(valid));

        % ── Point-wise difference statistics ──────────────────────────────
        dux = ux_L(valid) - ux_R(valid);
        duy = uy_L(valid) - uy_R(valid);

        fprintf('    ux diff — mean: %+.4f px  std: %.4f px  rms: %.4f px\n', ...
            mean(dux), std(dux), rms(dux));
        fprintf('    uy diff — mean: %+.4f px  std: %.4f px  rms: %.4f px\n', ...
            mean(duy), std(duy), rms(duy));

        % ── 2D cross-correlation of ux fields ─────────────────────────────
        % Zero out invalid points and remove mean before correlating
        ux_L_cc = ux_L - mean(ux_L(valid));  ux_L_cc(~valid) = 0;
        ux_R_cc = ux_R - mean(ux_R(valid));  ux_R_cc(~valid) = 0;

        uy_L_cc = uy_L - mean(uy_L(valid));  uy_L_cc(~valid) = 0;
        uy_R_cc = uy_R - mean(uy_R(valid));  uy_R_cc(~valid) = 0;

        xc_ux = normxcorr2(ux_L_cc, ux_R_cc);
        xc_uy = normxcorr2(uy_L_cc, uy_R_cc);

        % Peak location relative to zero lag
        [~, idx_ux] = max(xc_ux(:));
        [~, idx_uy] = max(xc_uy(:));

        [pr_ux, pc_ux] = ind2sub(size(xc_ux), idx_ux);
        [pr_uy, pc_uy] = ind2sub(size(xc_uy), idx_uy);

        % Zero lag is at [nRows, nCols] in normxcorr2 output
        lag_col_ux = pc_ux - nCols;   % columns = x direction
        lag_row_ux = pr_ux - nRows;   % rows    = y direction
        lag_col_uy = pc_uy - nCols;
        lag_row_uy = pr_uy - nRows;

        fprintf('    ux cross-corr peak lag — col: %+d px  row: %+d px\n', ...
            lag_col_ux, lag_row_ux);
        fprintf('    uy cross-corr peak lag — col: %+d px  row: %+d px\n', ...
            lag_col_uy, lag_row_uy);

        % Peak correlation value — how similar are the fields?
        fprintf('    ux peak corr: %.4f  |  uy peak corr: %.4f\n', ...
            max(xc_ux(:)), max(xc_uy(:)));

        % ── Store ─────────────────────────────────────────────────────────
        piv(i).calib(c).overlap(i).ux_L_trimmed  = ux_L;
        piv(i).calib(c).overlap(i).ux_R_trimmed  = ux_R;
        piv(i).calib(c).overlap(i).uy_L_trimmed  = uy_L;
        piv(i).calib(c).overlap(i).uy_R_trimmed  = uy_R;
        piv(i).calib(c).overlap(i).valid_trimmed  = valid;
        piv(i).calib(c).overlap(i).dux_mean       = mean(dux);
        piv(i).calib(c).overlap(i).dux_std        = std(dux);
        piv(i).calib(c).overlap(i).dux_rms        = rms(dux);
        piv(i).calib(c).overlap(i).duy_mean       = mean(duy);
        piv(i).calib(c).overlap(i).duy_std        = std(duy);
        piv(i).calib(c).overlap(i).duy_rms        = rms(duy);
        piv(i).calib(c).overlap(i).xcorr_ux_lag   = [lag_col_ux, lag_row_ux];
        piv(i).calib(c).overlap(i).xcorr_uy_lag   = [lag_col_uy, lag_row_uy];
        piv(i).calib(c).overlap(i).xcorr_ux_peak  = max(xc_ux(:));
        piv(i).calib(c).overlap(i).xcorr_uy_peak  = max(xc_uy(:));
    end
end

fprintf('\nOverlap comparison complete.\n');


%% ── Visualise overlap comparison in the pixel to pixel space──────────────────────────────── 
for c = 1:nCalibs

    fig = figure(200 + c);
    clf;
    set(fig, 'Name', sprintf('Overlap comparison — %s', calib_names{c}), ...
        'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.88]);

    % Use a custom layout: top row 5 wide, bottom row 4 wide
    % Achieved with subplot using a 2x5 grid but merging isn't needed —
    % just use tiledlayout for clean control
    tl = tiledlayout(2, max(nCams, nCams-1), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Overlap comparison — %s  |  Frame %05d', calib_names{c}, nFrame), ...
        'FontSize', 12, 'FontWeight', 'bold');

    % ── Top row: raw ux for each camera with overlap shaded ───────────────
    for i = 1:nCams

        ax = nexttile(tl, i);

        % Plot ux as image
        imagesc(ax, piv(i).x_px, piv(i).y_px, piv(i).ux);
        colormap(ax, 'redblue');
        colorbar(ax);
        axis(ax, 'image'); axis(ax, 'on');
        xlabel(ax, 'x (px)'); ylabel(ax, 'y (px)');
        title(ax, sprintf('%s', camNames{i}), 'FontSize', 9);

        % Clim symmetric around zero
        % Tighter clim — clip at 5th/95th percentile of valid data
        ux_valid = piv(i).ux(isfinite(piv(i).ux) & piv(i).valid);
        % clim(ax, [prctile(ux_valid, 5), prctile(ux_valid, 95)]);
        clim(ax, [-17 -10]);

        hold(ax, 'on');

        % ── Shade overlap with left neighbour ─────────────────────────────
        if i > 1
            ov_prev = piv(i-1).calib(c).overlap(i-1);
            col_idx = ov_prev.col_idx_R;
            x_lo    = piv(i).x_px(col_idx(1));
            x_hi    = piv(i).x_px(col_idx(end));
            y_lo    = piv(i).y_px(1);
            y_hi    = piv(i).y_px(end);
            fill(ax, [x_lo x_hi x_hi x_lo], [y_lo y_lo y_hi y_hi], ...
                'c', 'FaceAlpha', 0.25, 'EdgeColor', 'c', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Overlap w/ Cam%d', i-1));
        end

        % ── Shade overlap with right neighbour ────────────────────────────
        if i < nCams
            ov_next = piv(i).calib(c).overlap(i);
            col_idx = ov_next.col_idx_L;
            x_lo    = piv(i).x_px(col_idx(1));
            x_hi    = piv(i).x_px(col_idx(end));
            y_lo    = piv(i).y_px(1);
            y_hi    = piv(i).y_px(end);
            fill(ax, [x_lo x_hi x_hi x_lo], [y_lo y_lo y_hi y_hi], ...
                'y', 'FaceAlpha', 0.25, 'EdgeColor', 'y', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Overlap w/ Cam%d', i+1));
        end
        % ── Draw floor line from calibration in pixel space ───────────────
        fp_i = calibs{c}.var{i}.floor_px;
        plot(ax, [fp_i(1,1), fp_i(2,1)], [fp_i(1,2), fp_i(2,2)], ...
            'g-', 'LineWidth', 2, 'DisplayName', 'Floor (autocorr)');
        legend(ax, 'show', 'Location', 'northwest', 'FontSize', 6);
        hold(ax, 'off');
    end
    % ── Bottom row: cross-correlation field for each adjacent pair ────────
    for i = 1:nCams-1

        ax = nexttile(tl, nCams + i);

        ov    = piv(i).calib(c).overlap(i);
        calib = calibs{c};

        % ── Build floor mask in pixel space ───────────────────────────────
        % floor_px = [col_left, row_left; col_right, row_right]
        fp_L  = calib.var{i}.floor_px;
        fp_R  = calib.var{i+1}.floor_px;

        % Get pixel row coords of overlap sub-grid
        y_px_L = piv(i).y_px;     % [1 x ny] row positions
        y_px_R = piv(i+1).y_px;

        % Interpolate floor row at each PIV row position for left camera
        floor_row_L = interp1([fp_L(1,1), fp_L(2,1)], ...
            [fp_L(1,2), fp_L(2,2)], ...
            piv(i).x_px(ov.col_idx_L), ...
            'linear', 'extrap');
        % Use mean floor row across overlap columns
        floor_row_L_mean = mean(floor_row_L);

        floor_row_R = interp1([fp_R(1,1), fp_R(2,1)], ...
            [fp_R(1,2), fp_R(2,2)], ...
            piv(i+1).x_px(ov.col_idx_R), ...
            'linear', 'extrap');
        floor_row_R_mean = mean(floor_row_R);

        % Row mask: keep rows ABOVE the floor (row index < floor row)
        above_floor_L = y_px_L(:) < floor_row_L_mean;
        above_floor_R = y_px_R(:) < floor_row_R_mean;

        % ── Extract and trim overlap sub-grids ────────────────────────────
        nCols     = min(size(ov.ux_L_trimmed,2), size(ov.ux_R_trimmed,2));
        ux_L_full = ov.ux_L_trimmed(:, end-nCols+1:end);
        ux_R_full = ov.ux_R_trimmed(:, 1:nCols);
        val_full  = ov.valid_trimmed(:, 1:nCols);

        % Apply floor mask — zero out rows at or below floor
        ux_L_cc = ux_L_full;
        ux_R_cc = ux_R_full;
        ux_L_cc(~above_floor_L, :) = 0;
        ux_R_cc(~above_floor_R, :) = 0;

        % Also zero invalid points and remove mean from valid region
        valid_L = val_full & above_floor_L;
        valid_R = val_full & above_floor_R;
        ux_L_cc(~valid_L) = 0;
        ux_R_cc(~valid_R) = 0;

        if any(valid_L(:))
            ux_L_cc = ux_L_cc - mean(ux_L_cc(valid_L));
        end
        if any(valid_R(:))
            ux_R_cc = ux_R_cc - mean(ux_R_cc(valid_R));
        end

        xc_ux = normxcorr2(ux_L_cc, ux_R_cc);

        % Lag axes in PIV windows
        nR = size(ux_L_cc, 1);
        nC = size(ux_L_cc, 2);
        lag_col = (1:size(xc_ux,2)) - nC;
        lag_row = (1:size(xc_ux,1)) - nR;

        imagesc(ax, lag_col, lag_row, xc_ux);
        colormap(ax, 'parula'); colorbar(ax);
        axis(ax, 'on');
        xlabel(ax, 'Col lag (PIV windows)');
        ylabel(ax, 'Row lag (PIV windows)');

        zoom_win = 200;
        xlim(ax, [-zoom_win, zoom_win]);
        ylim(ax, [-zoom_win, zoom_win]);

        % Find peak
        [~, peak_idx]    = max(xc_ux(:));
        [pr, pc]         = ind2sub(size(xc_ux), peak_idx);
        lag_c            = lag_col(pc);
        lag_r            = lag_row(pr);
        peak_val         = max(xc_ux(:));

        hold(ax, 'on');
        plot(ax, 0,     0,     'w+', 'MarkerSize', 12, 'LineWidth', 2, ...
            'DisplayName', 'Zero lag');
        plot(ax, lag_c, lag_r, 'rx', 'MarkerSize', 12, 'LineWidth', 2, ...
            'DisplayName', sprintf('Peak [%+d, %+d]', lag_c, lag_r));
        legend(ax, 'show', 'Location', 'northeast', 'FontSize', 6);
        hold(ax, 'off');

        title(ax, sprintf('%s <-> %s  |  peak=%.3f  lag=[%+d,%+d] win', ...
            camNames{i}, camNames{i+1}, peak_val, lag_c, lag_r), ...
            'FontSize', 8);
    end

    nexttile(tl, nCams + nCams);
    axis off;

    % Leave last bottom tile empty since nCams-1 = 4 < nCams = 5
    nexttile(tl, nCams + nCams);
    axis(ax, 'on');

    % ── Save ──────────────────────────────────────────────────────────────
    comp_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save_name = sprintf('overlap_comparison_%s_%s.png', ...
        strrep(calib_names{c}, ' ', '_'), comp_timestamp);
    exportgraphics(fig, fullfile(baseDir, save_name), 'Resolution', 150);
    fprintf('Saved: %s\n', save_name);
end



%% ── Verify the pixel shift in C2,C3 do the autocorrelation in the physical axes whilst keeping the vector field in pix─────────── 

% For each overlap pair, match points by proximity in mm
for c = 1: nCalibs
    for i = 1:nCams-1
        ov = piv(i).calib(c).overlap(i);

        % Recompute valid masks locally — sized to match ov.ux_L_ov / ux_R_ov
        valid_L = isfinite(ov.ux_L_ov) & piv(i).valid(:,   ov.col_idx_L);
        valid_R = isfinite(ov.ux_R_ov) & piv(i+1).valid(:, ov.col_idx_R);

        % Flatten valid mm coords + velocities for each camera
        x1_L_v = ov.x1_L_ov(valid_L);  x2_L_v = ov.x2_L_ov(valid_L);  ux_L_v = ov.ux_L_ov(valid_L);
        x1_R_v = ov.x1_R_ov(valid_R);  x2_R_v = ov.x2_R_ov(valid_R);  ux_R_v = ov.ux_R_ov(valid_R);

        % For each L point, find nearest R point in mm space
        idx = knnsearch([x1_R_v, x2_R_v], [x1_L_v, x2_L_v]);

        % Only keep pairs within a distance threshold (e.g. 1 mm)
        dist = sqrt((x1_L_v - x1_R_v(idx)).^2 + (x2_L_v - x2_R_v(idx)).^2);
        closestPix = dist < 1.0;

        ux_L_matched = ux_L_v(closestPix);
        ux_R_matched = ux_R_v(idx(closestPix));

        r = corr(ux_L_matched, ux_R_matched);
        fprintf('%s <-> %s:  r = %.4f   mean diff = %+.4f px   rms = %.4f px\n', ...
            camNames{i}, camNames{i+1}, r, mean(ux_L_matched - ux_R_matched), ...
            rms(ux_L_matched - ux_R_matched));
    end
end

%% ── Cross-correlation in mm space (interpolated) ─────────────────────────────

for c = 1:nCalibs

    fig_mm = figure(300 + c);
    clf;
    set(fig_mm, ...
        'Name',     sprintf('Overlap xcorr — mm space — %s', calib_names{c}), ...
        'Units',    'normalized', ...
        'Position', [0.02 0.05 0.96 0.45]);

    tl_mm = tiledlayout(1, nCams-1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl_mm, ...
        sprintf('Cross-correlation in mm space — %s  |  Frame %05d', ...
                calib_names{c}, nFrame), ...
        'FontSize', 12, 'FontWeight', 'bold');

    fprintf('\n=== mm-space xcorr — %s ===\n', calib_names{c});

    for i = 1:nCams-1

        ov = piv(i).calib(c).overlap(i);

        % ── Pull overlap arrays ───────────────────────────────────────────
        x1_L = ov.x1_L_ov;   % [ny x nCols_L]
        x2_L = ov.x2_L_ov;
        x1_R = ov.x1_R_ov;   % [ny x nCols_R]
        x2_R = ov.x2_R_ov;
        ux_L = ov.ux_L_ov;
        ux_R = ov.ux_R_ov;

        % ── Recompute valid masks locally (sized to match ov arrays) ──────
        valid_L = isfinite(ux_L) & piv(i).valid(:,   ov.col_idx_L);
        valid_R = isfinite(ux_R) & piv(i+1).valid(:, ov.col_idx_R);

        if nnz(valid_L) < 10 || nnz(valid_R) < 10
            warning('Too few valid points for %s <-> %s — skipping.', ...
                camNames{i}, camNames{i+1});
            continue
        end

        % ── Common regular mm grid ────────────────────────────────────────
        x1_lo = max(min(x1_L(valid_L)), min(x1_R(valid_R)));
        x1_hi = min(max(x1_L(valid_L)), max(x1_R(valid_R)));
        x2_lo = max(min(x2_L(valid_L)), min(x2_R(valid_R)));
        x2_hi = min(max(x2_L(valid_L)), max(x2_R(valid_R)));

        if x1_lo >= x1_hi || x2_lo >= x2_hi
            warning('No valid mm intersection for %s <-> %s — skipping.', ...
                camNames{i}, camNames{i+1});
            continue
        end

        % Grid spacing from mean calibrated step size
        dx1 = abs(mean(diff(mean(x1_L, 1))));
        dx2 = abs(mean(diff(mean(x2_L, 2))));

        [x1g, x2g] = meshgrid(x1_lo:dx1:x1_hi, x2_lo:dx2:x2_hi);

        fprintf('  %s <-> %s  |  mm grid [%d x %d]  dx1=%.3f mm  dx2=%.3f mm\n', ...
            camNames{i}, camNames{i+1}, size(x1g,1), size(x1g,2), dx1, dx2);

        % ── Interpolate onto common grid ──────────────────────────────────
        F_L = scatteredInterpolant( ...
            x1_L(valid_L), x2_L(valid_L), ux_L(valid_L), 'linear', 'none');
        F_R = scatteredInterpolant( ...
            x1_R(valid_R), x2_R(valid_R), ux_R(valid_R), 'linear', 'none');

        ux_L_mm = F_L(x1g, x2g);
        ux_R_mm = F_R(x1g, x2g);

        valid_mm = isfinite(ux_L_mm) & isfinite(ux_R_mm);

        fprintf('    valid on mm grid: %d / %d  (%.1f%%)\n', ...
            nnz(valid_mm), numel(valid_mm), 100*nnz(valid_mm)/numel(valid_mm));

        if nnz(valid_mm) < 10
            warning('Too few valid grid points after interpolation for %s <-> %s', ...
                camNames{i}, camNames{i+1});
            continue
        end

        % ── Zero mean, mask invalid, cross-correlate ──────────────────────
        ux_L_cc = ux_L_mm;  ux_L_cc(~valid_mm) = 0;
        ux_R_cc = ux_R_mm;  ux_R_cc(~valid_mm) = 0;
        ux_L_cc = ux_L_cc - mean(ux_L_cc(valid_mm));
        ux_R_cc = ux_R_cc - mean(ux_R_cc(valid_mm));

        xc_ux = normxcorr2(ux_L_cc, ux_R_cc);

        % ── Lag axes in mm ────────────────────────────────────────────────
        nR_grid = size(ux_L_cc, 1);
        nC_grid = size(ux_L_cc, 2);
        lag_x1_mm = ((1:size(xc_ux, 2)) - nC_grid) * dx1;
        lag_x2_mm = ((1:size(xc_ux, 1)) - nR_grid) * dx2;

        % ── Peak ──────────────────────────────────────────────────────────
        [peak_val, peak_idx] = max(xc_ux(:));
        [pr, pc]             = ind2sub(size(xc_ux), peak_idx);
        lag_x1_peak          = lag_x1_mm(pc);
        lag_x2_peak          = lag_x2_mm(pr);

        dux_mm = ux_L_mm(valid_mm) - ux_R_mm(valid_mm);
        fprintf('    peak = %.4f   lag = [%+.3f mm, %+.3f mm]  (x1, x2)\n', ...
            peak_val, lag_x1_peak, lag_x2_peak);
        fprintf('    ux diff — mean: %+.4f px  std: %.4f px  rms: %.4f px\n', ...
            mean(dux_mm), std(dux_mm), rms(dux_mm));

        % ── Plot ──────────────────────────────────────────────────────────
        ax = nexttile(tl_mm, i);
        imagesc(ax, lag_x1_mm, lag_x2_mm, xc_ux);
        colormap(ax, 'parula');
        colorbar(ax);
        axis(ax, 'on');
        xlabel(ax, 'x_1 lag (mm)');
        ylabel(ax, 'x_2 lag (mm)');

        % Zoom to ±half the overlap width in each direction
        zoom_x1 = (x1_hi - x1_lo) / 2;
        zoom_x2 = (x2_hi - x2_lo) / 2;
        xlim(ax, [-zoom_x1, zoom_x1]);
        ylim(ax, [-zoom_x2, zoom_x2]);

        hold(ax, 'on');
        plot(ax, 0, 0, 'w+', 'MarkerSize', 14, 'LineWidth', 2, ...
            'DisplayName', 'Zero lag');
        plot(ax, lag_x1_peak, lag_x2_peak, 'rx', 'MarkerSize', 14, 'LineWidth', 2, ...
            'DisplayName', sprintf('Peak [%+.2f, %+.2f] mm', lag_x1_peak, lag_x2_peak));
        legend(ax, 'show', 'Location', 'northeast', 'FontSize', 7);
        hold(ax, 'off');

        title(ax, sprintf('%s \\leftrightarrow %s\npeak = %.3f   lag = [%+.2f, %+.2f] mm', ...
            camNames{i}, camNames{i+1}, peak_val, lag_x1_peak, lag_x2_peak), ...
            'FontSize', 8);

        % ── Store ─────────────────────────────────────────────────────────
        piv(i).calib(c).overlap(i).mm_xcorr.dx1         = dx1;
        piv(i).calib(c).overlap(i).mm_xcorr.dx2         = dx2;
        piv(i).calib(c).overlap(i).mm_xcorr.xc_ux       = xc_ux;
        piv(i).calib(c).overlap(i).mm_xcorr.lag_x1_mm   = lag_x1_mm;
        piv(i).calib(c).overlap(i).mm_xcorr.lag_x2_mm   = lag_x2_mm;
        piv(i).calib(c).overlap(i).mm_xcorr.peak_val    = peak_val;
        piv(i).calib(c).overlap(i).mm_xcorr.lag_x1_peak = lag_x1_peak;
        piv(i).calib(c).overlap(i).mm_xcorr.lag_x2_peak = lag_x2_peak;
        piv(i).calib(c).overlap(i).mm_xcorr.dux_mean    = mean(dux_mm);
        piv(i).calib(c).overlap(i).mm_xcorr.dux_std     = std(dux_mm);
        piv(i).calib(c).overlap(i).mm_xcorr.dux_rms     = rms(dux_mm);
    end

    % ── Save ──────────────────────────────────────────────────────────────
    mm_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save_name_mm = sprintf('overlap_xcorr_mm_%s_%s.png', ...
        strrep(calib_names{c}, ' ', '_'), mm_timestamp);
    exportgraphics(fig_mm, fullfile(baseDir, save_name_mm), 'Resolution', 150);
    fprintf('Saved: %s\n', save_name_mm);

end

fprintf('\nMM-space cross-correlation complete.\n');
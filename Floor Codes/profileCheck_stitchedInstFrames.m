% ── User settings ─────────────────────────────────────────────────────────────
dt       = 49.6e-6;   % s
baseDir  = 'C:\Users\ak1u24\Downloads\PG_fixedmergandcal\loop=06\uncalibrated_piv\150\';
camNames = {'Cam1','Cam2','Cam3','Cam4','Cam5'};
nCams    = numel(camNames);
spot_frames  = [10];
plot_prctile = [2, 98];
alpha_overlap = 0.5;
alpha_normal  = 1.0;
case_colors   = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.7 0.3]};

% ── Load calibrations ──────────────────────────────────────────────────────────
calib_path_A = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\calib_caseA_20260412_114104.mat';
calib_path_B = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\calib_caseB_20260412_114104.mat';
calib_path_C = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check\0524CAL\CALIBRATION0524\calib_caseC_20260412_114104.mat';

tmp = load(calib_path_A);  calib_A = tmp.calib_A;
tmp = load(calib_path_B);  calib_B = tmp.calib_B;
tmp = load(calib_path_C);  calib_C = tmp.calib_C;

calibs      = {calib_B};
calib_names = {''};
nCalibs     = numel(calibs);

fprintf('Loaded %d calibration cases.\n', nCalibs);

% ── Derive nFrames from folder contents ───────────────────────────────────────
refFiles = dir(fullfile(baseDir, camNames{1}, 'instantaneous', '*.mat'));
refFiles = refFiles(~strcmpi({refFiles.name}, 'coordinates.mat'));
[~, sortIdx] = sort({refFiles.name});
refFiles   = refFiles(sortIdx);
nFrames    = numel(refFiles);
frameNames = {refFiles.name};

fprintf('Found %d frames (coordinates.mat excluded)\n', nFrames);

% ── Load one frame to get PIV grid info per camera ────────────────────────────
fprintf('\nLoading PIV grid info...\n');

for i = 1:nCams
    filePath = fullfile(baseDir, camNames{i}, 'instantaneous', frameNames{1});
    raw      = load(filePath);
    p        = raw.piv_result(end);
    clear raw;
    piv(i).x_px        = p.win_ctrs_x;
    piv(i).y_px        = p.win_ctrs_y;
    piv(i).window_size = p.window_size;
    piv(i).n_windows   = p.n_windows;
    piv(i).cam         = camNames{i};
    clear p;
    fprintf('  %s: loaded\n', camNames{i});
end

% ── Build pixel grids and pre-compute window centre mm positions ───────────────
fprintf('\nPre-computing window centre mm positions...\n');

for i = 1:nCams

    [wincenter_px1, wincenter_px2] = meshgrid(piv(i).x_px, piv(i).y_px);
    piv(i).x_grid = wincenter_px1;
    piv(i).y_grid = wincenter_px2;

    px_in   = [wincenter_px1(:), wincenter_px2(:)];
    grid_sz = size(wincenter_px1);

    for c = 1:nCalibs
        calib = calibs{c};

        x1_ctr = reshape(polyvaln(calib.fit_cx1{i}.p, px_in), grid_sz);
        x2_ctr = reshape(polyvaln(calib.fit_cx2{i}.p, px_in), grid_sz);

        if c == 3
            x2_ctr = x2_ctr - (calib.global_affine.a * x1_ctr + calib.global_affine.b);
        end

        piv(i).calib(c).x1_ctr_mm = x1_ctr;
        piv(i).calib(c).x2_ctr_mm = x2_ctr;
        piv(i).calib(c).name      = calib_names{c};
    end
end

% ── Find overlap mm bounds ────────────────────────────────────────────────────
fprintf('\nFinding overlap mm bounds...\n');

for c = 1:nCalibs
    for i = 1:nCams-1

        x1_min_L = min(piv(i  ).calib(c).x1_ctr_mm(:));
        x1_max_L = max(piv(i  ).calib(c).x1_ctr_mm(:));
        x1_min_R = min(piv(i+1).calib(c).x1_ctr_mm(:));
        x1_max_R = max(piv(i+1).calib(c).x1_ctr_mm(:));

        x2_min_L = min(piv(i  ).calib(c).x2_ctr_mm(:));
        x2_max_L = max(piv(i  ).calib(c).x2_ctr_mm(:));
        x2_min_R = min(piv(i+1).calib(c).x2_ctr_mm(:));
        x2_max_R = max(piv(i+1).calib(c).x2_ctr_mm(:));

        overlap(i,c).x1_lo    = max(x1_min_L, x1_min_R);
        overlap(i,c).x1_hi    = min(x1_max_L, x1_max_R);
        overlap(i,c).x2_lo    = max(x2_min_L, x2_min_R);
        overlap(i,c).x2_hi    = min(x2_max_L, x2_max_R);
        overlap(i,c).x1_width = overlap(i,c).x1_hi - overlap(i,c).x1_lo;
        overlap(i,c).x2_width = overlap(i,c).x2_hi - overlap(i,c).x2_lo;

        fprintf('  %s<->%s  %s: x1=[%.1f, %.1f] mm (%.1f mm wide)\n', ...
            camNames{i}, camNames{i+1}, calib_names{c}, ...
            overlap(i,c).x1_lo, overlap(i,c).x1_hi, overlap(i,c).x1_width);
    end
end

% ── Pre-compute calibrated u fields for spot frames ───────────────────────────
fprintf('\nPre-computing calibrated u fields for spot frames...\n');

u_spot   = cell(numel(spot_frames), nCams, nCalibs);
vld_spot = cell(numel(spot_frames), nCams);

for sf = 1:numel(spot_frames)

    nFrame_spot   = spot_frames(sf);
    fileName_spot = frameNames{nFrame_spot};
    fprintf('  Spot frame %d (%s)\n', nFrame_spot, fileName_spot);

    for i = 1:nCams
        filePath = fullfile(baseDir, camNames{i}, 'instantaneous', fileName_spot);
        if ~isfile(filePath), continue; end
        raw  = load(filePath);
        p    = raw.piv_result(end);
        ux_i = double(p.ux);
        uy_i = double(p.uy);
        vld  = ~p.b_mask & ~p.nan_mask;
        clear raw p;

        vld_spot{sf,i} = vld;

        for c = 1:nCalibs
            calib = calibs{c};

            px1_i = piv(i).x_grid + uy_i;
            px2_i = piv(i).y_grid + ux_i;
            pd_i  = [px1_i(:), px2_i(:)];

            x1d = reshape(polyvaln(calib.fit_cx1{i}.p, pd_i), size(ux_i));
            x2d = reshape(polyvaln(calib.fit_cx2{i}.p, pd_i), size(ux_i));
            if c == 3
                x2d = x2d - (calib.global_affine.a * x1d + calib.global_affine.b);
            end

            u_i = (x2d - piv(i).calib(c).x2_ctr_mm) * 1e-3 / dt;
            u_i(~vld) = NaN;

            u_clamp = prctile(u_i(isfinite(u_i)), plot_prctile);
            u_i(u_i < u_clamp(1) | u_i > u_clamp(2)) = NaN;

            u_spot{sf,i,c} = u_i;
        end
    end
end

fprintf('Pre-computation complete.\n');

% ══════════════════════════════════════════════════════════════════════════════
% CHECK 1 — STITCHED INSTANTANEOUS VIEW
% ══════════════════════════════════════════════════════════════════════════════
fprintf('\nCheck 1 — Stitched instantaneous view...\n');

for sf = 1:numel(spot_frames)

    nFrame_spot = spot_frames(sf);
    fprintf('  Frame %d\n', nFrame_spot);

    fig = figure(900 + sf); clf;
    set(fig, 'Name', sprintf('Stitched view — Frame %d', nFrame_spot), ...
        'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.88]);

    tl = tiledlayout(nCalibs, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Stitched instantaneous U (m/s) — Frame %d  |  overlap = %.0f%% alpha', ...
        nFrame_spot, alpha_overlap*100), 'FontSize', 12, 'FontWeight', 'bold');

    for c = 1:nCalibs

        ax = nexttile(tl, c);
        hold(ax, 'on');

        all_u = [];
        for i = 1:nCams
            u_i = u_spot{sf,i,c};
            all_u = [all_u; u_i(isfinite(u_i))]; %#ok<AGROW>
        end
        clim_vals = prctile(all_u, plot_prctile);

        for i = 1:nCams

            u_i  = u_spot{sf,i,c};
            x1_i = piv(i).calib(c).x1_ctr_mm(1,:);
            x2_i = piv(i).calib(c).x2_ctr_mm(:,1);

            row_above = x2_i > 0;
            u_plot    = u_i(row_above, :);
            x2_plot   = x2_i(row_above);

            [x2_plot_s, si] = sort(x2_plot);
            u_plot_s        = u_plot(si, :);

            u_plot_s(u_plot_s < clim_vals(1)) = clim_vals(1);
            u_plot_s(u_plot_s > clim_vals(2)) = clim_vals(2);

            alpha_mask = alpha_normal * ones(size(u_plot_s));
            if i > 1
                ov_prev  = overlap(i-1, c);
                col_left = x1_i >= ov_prev.x1_lo & x1_i <= ov_prev.x1_hi;
                alpha_mask(:, col_left) = alpha_overlap;
            end
            if i < nCams
                ov_next   = overlap(i, c);
                col_right = x1_i >= ov_next.x1_lo & x1_i <= ov_next.x1_hi;
                alpha_mask(:, col_right) = alpha_overlap;
            end

            h = imagesc(ax, x1_i, x2_plot_s, u_plot_s);
            h.AlphaData = alpha_mask;
            clim(ax, clim_vals);
        end

        for i = 1:nCams-1
            ov = overlap(i,c);
            xline(ax, ov.x1_lo, 'w--', 'LineWidth', 1.5, ...
                'Label', sprintf('%s|%s', camNames{i}, camNames{i+1}), ...
                'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
            xline(ax, ov.x1_hi, 'w:', 'LineWidth', 1.0);
        end

        colormap(ax, 'turbo');
        colorbar(ax);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal');
        xlabel(ax, 'x_1 streamwise (mm)');
        ylabel(ax, 'x_2 wall-normal (mm)');
        title(ax, sprintf('%s — frame %d', calib_names{c}, nFrame_spot), 'FontSize', 10);
        hold(ax, 'off');
    end

    ts = datestr(now, 'yyyymmdd_HHMMSS');
    exportgraphics(fig, fullfile(baseDir, ...
        sprintf('stitched_frame%03d_%s.png', nFrame_spot, ts)), 'Resolution', 150);
    fprintf('  Saved: frame %d\n', nFrame_spot);
end
%%
% ══════════════════════════════════════════════════════════════════════════════
% CHECK 2 — WALL-NORMAL PROFILE COMPARISON IN OVERLAP REGIONS
% nCalibs x 3 subplots per figure (one figure per camera pair)
% ══════════════════════════════════════════════════════════════════════════════
fprintf('\nCheck 2 — Wall-normal profile comparison...\n');

ls_styles  = {'-', '--', ':'};
loc_labels = {'left', 'centre', 'right'};

for i = 1:nCams-1

    fig = figure(1100 + i); clf;
    set(fig, 'Name', sprintf('Profiles — %s <-> %s', camNames{i}, camNames{i+1}), ...
        'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.88]);

    tl = tiledlayout(nCalibs, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Wall-normal profiles in overlap — %s <-> %s\n(averaged over frames %s)', ...
        camNames{i}, camNames{i+1}, num2str(spot_frames)), ...
        'FontSize', 12, 'FontWeight', 'bold');

    for c = 1:nCalibs

        ov  = overlap(i,c);
        dx1 = abs(mean(diff(piv(i).calib(c).x1_ctr_mm(1,:))));

        x1_locs = linspace(ov.x1_lo, ov.x1_hi, 3);

        % x2 axis for left camera — sorted ascending
        x2_L_col = piv(i).calib(c).x2_ctr_mm(:,1);
        [x2_sorted, si_L] = sort(x2_L_col);
        above_floor = x2_sorted > 0;
        x2_axis     = x2_sorted(above_floor);

        % x2 axis for right camera — sorted ascending
        x2_R_col = piv(i+1).calib(c).x2_ctr_mm(:,1);
        [x2_R_sorted, si_R] = sort(x2_R_col);
        above_floor_R = x2_R_sorted > 0;
        x2_axis_R     = x2_R_sorted(above_floor_R);

        for px = 1:3

            % Tile index: row = calibration case, col = streamwise location
            ax = nexttile(tl, (c-1)*3 + px);
            hold(ax, 'on');

            x1_target  = x1_locs(px);
            prof_L_acc = zeros(sum(above_floor), 1);
            prof_R_acc = zeros(sum(above_floor_R), 1);
            n_acc      = 0;

            for sf = 1:numel(spot_frames)

                % Left camera
                x1_L_row   = piv(i).calib(c).x1_ctr_mm(1,:);
                col_mask_L = abs(x1_L_row - x1_target) <= dx1;
                if ~any(col_mask_L)
                    [~, nn] = min(abs(x1_L_row - x1_target));
                    col_mask_L(nn) = true;
                end
                fld_L  = u_spot{sf, i, c};
                prof_L = mean(fld_L(:, col_mask_L), 2, 'omitnan');
                prof_L = prof_L(si_L);
                prof_L_acc = prof_L_acc + prof_L(above_floor);

                % Right camera
                x1_R_row   = piv(i+1).calib(c).x1_ctr_mm(1,:);
                col_mask_R = abs(x1_R_row - x1_target) <= dx1;
                if ~any(col_mask_R)
                    [~, nn] = min(abs(x1_R_row - x1_target));
                    col_mask_R(nn) = true;
                end
                fld_R  = u_spot{sf, i+1, c};
                prof_R = mean(fld_R(:, col_mask_R), 2, 'omitnan');
                prof_R = prof_R(si_R);
                prof_R_acc = prof_R_acc + prof_R(above_floor_R);

                n_acc = n_acc + 1;
            end

            prof_L_avg = prof_L_acc / n_acc;
            prof_R_avg = prof_R_acc / n_acc;

            % Clamp outliers
            lims = prctile([prof_L_avg; prof_R_avg], plot_prctile);
            prof_L_avg(prof_L_avg < lims(1) | prof_L_avg > lims(2)) = NaN;
            prof_R_avg(prof_R_avg < lims(1) | prof_R_avg > lims(2)) = NaN;

            % Left camera — solid
            h_L = semilogx(ax, x2_axis(isfinite(prof_L_avg)), ...
                prof_L_avg(isfinite(prof_L_avg)), ...
                'Color', case_colors{c+1}, 'LineStyle', ':', ...
                'LineWidth', 1, 'DisplayName', sprintf('%s L', camNames{i}));
            h_L.Marker        = 'o';
            h_L.MarkerSize    = 5;
            h_L.MarkerIndices = 1:8:sum(valid_R);

            % Right camera — circles
            valid_R = isfinite(prof_R_avg);
            h_R = semilogx(ax, x2_axis_R(valid_R), prof_R_avg(valid_R), ...
                'Color', case_colors{c}, 'LineStyle', ':', 'LineWidth', 1);
            h_R.Marker        = 'o';
            h_R.MarkerSize    = 5;
            h_R.MarkerIndices = 1:8:sum(valid_R);
            h_R.DisplayName   = sprintf('%s R', camNames{i+1});

            set(ax, 'XScale', 'log');
            grid(ax, 'on');
            xlabel(ax, 'x_2 wall-normal (mm)');
            ylabel(ax, 'U (m/s)');
            title(ax, sprintf('%s  |  %s  |  x1=%.0f mm', ...
                calib_names{c}, loc_labels{px}, x1_target), 'FontSize', 8);
            legend(ax, 'show', 'Location', 'best', 'FontSize', 7);
            hold(ax, 'off');
        end
    end

    ts = datestr(now, 'yyyymmdd_HHMMSS');
    exportgraphics(fig, fullfile(baseDir, ...
        sprintf('profiles_%s_vs_%s_%s.png', camNames{i}, camNames{i+1}, ts)), ...
        'Resolution', 150);
    fprintf('  Saved: %s <-> %s\n', camNames{i}, camNames{i+1});
end

fprintf('\nDone.\n');


%% ── Compare floor results per day ────────────────────────────────────────────
clc;

rootPath = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check';
camNums  = 1:5;

% ── Discover date folders ─────────────────────────────────────────────────────
dateFolders = dir(rootPath);
dateFolders = dateFolders([dateFolders.isdir]);
dateFolders = dateFolders(~ismember({dateFolders.name}, {'.','..'}));

% ── Loop per day ──────────────────────────────────────────────────────────────
for df = 1:numel(dateFolders)

    dateFolder  = fullfile(rootPath, dateFolders(df).name);
    dayLabel    = dateFolders(df).name;

    fprintf('\n========== Processing day: %s ==========\n', dayLabel);

    % ── Find CALIBRATION subfolder for this day ───────────────────────────
    contents = dir(dateFolder);
    contents = contents([contents.isdir]);
    calIdx   = find(~cellfun(@isempty, strfind({contents.name}, 'CALIBRATION')), 1);

    if isempty(calIdx)
        fprintf('  No CALIBRATION folder found in %s — skipping day.\n', dateFolder);
        continue
    end
    calFolder = fullfile(dateFolder, contents(calIdx).name);
    fprintf('  Calibration folder: %s\n', calFolder);

    % ── Find case subfolders for this day ─────────────────────────────────
    caseDirs = dir(dateFolder);
    caseDirs = caseDirs([caseDirs.isdir]);
    caseDirs = caseDirs(~ismember({caseDirs.name}, {'.','..'}));
    caseDirs = caseDirs(cellfun(@isempty, strfind({caseDirs.name}, 'CALIBRATION')));

    % Collect cases for this day
    day_cases = {};   % {label, floorResults}
    for cd = 1:numel(caseDirs)

        caseFolder = fullfile(dateFolder, caseDirs(cd).name);
        caseLabel  = caseDirs(cd).name;

        matFiles = dir(fullfile(caseFolder, 'floorResults_*.mat'));
        if isempty(matFiles)
            fprintf('  No floorResults found in %s — skipping.\n', caseFolder);
            continue
        end

        % Most recent mat file
        [~, sortIdx] = sort([matFiles.datenum], 'descend');
        mostRecent   = fullfile(caseFolder, matFiles(sortIdx(1)).name);
        fprintf('  Loading: %s\n', mostRecent);

        loaded = load(mostRecent, 'floorResults');
        day_cases{end+1, 1} = caseLabel;
        day_cases{end,   2} = loaded.floorResults;
    end

    if isempty(day_cases)
        fprintf('  No valid cases found for day %s — skipping.\n', dayLabel);
        continue
    end

    nCases = size(day_cases, 1);
    cmap   = lines(nCases);
    fprintf('  Found %d case(s) for day %s.\n', nCases, dayLabel);
    % ── Load calibration TIFs once per day ───────────────────────────────
    calTifs = cell(1, numel(camNums));
    calSizes = zeros(numel(camNums), 2);   % [height, width] per camera

    for camIdx = 1:numel(camNums)
        calTifPath = fullfile(calFolder, sprintf('camera%d', camNums(camIdx)), 'B00001.tif');
        if isfile(calTifPath)
            calTifs{camIdx} = imread(calTifPath);
            [calSizes(camIdx,1), calSizes(camIdx,2)] = size(calTifs{camIdx});
            fprintf('  Cam%d: loaded %s\n', camNums(camIdx), calTifPath);
        else
            fprintf('  Cam%d: B00001.tif not found in %s\n', camNums(camIdx), ...
                fullfile(calFolder, sprintf('camera%d', camNums(camIdx))));
            calTifs{camIdx} = [];
            calSizes(camIdx,:) = [4600, 5312];   % fallback
        end
    end

    % ── Per camera figure for this day ────────────────────────────────────
    for camIdx = 1:numel(camNums)

        figure(100*df + camIdx); clf;
        ax = gca;
        % ── Use pre-loaded calibration TIF ────────────────────────────────
        if ~isempty(calTifs{camIdx})
            calImg = calTifs{camIdx};
            calH   = calSizes(camIdx, 1);
            calW   = calSizes(camIdx, 2);
            imagesc_clim(ax, calImg, 1, 99);
            colormap(ax, 'gray');
        else
            calH = calSizes(camIdx, 1);
            calW = calSizes(camIdx, 2);
            fprintf('  Cam%d: no calibration image — blank background\n', camNums(camIdx));
        end

        axis(ax, 'image'); axis(ax, 'on');
        xlabel(ax, 'Column (px)'); ylabel(ax, 'Row (px)');
        ax.XLim  = [1, calW];
        ax.YLim  = [1, calH];
        ax.XTick = 0:500:calW;
        ax.YTick = 0:500:calH;
        hold(ax, 'on');

        % ── Evaluate each case's line at common x points ──────────────────
        x_eval    = linspace(1, calW, 500);
        line_vals = nan(nCases, length(x_eval));   % one row per case

        for cc = 1:nCases
            fr = day_cases{cc, 2};
            if camIdx > numel(fr) || ~isfield(fr(camIdx), 'floor_line_combined')
                continue
            end
            line_vals(cc, :) = polyval(fr(camIdx).floor_line_combined, x_eval);
        end

        % Only use cases that have data
        valid_cases = ~all(isnan(line_vals), 2);

        % ── Plot each case ────────────────────────────────────────────────
        for cc = 1:nCases

            fr        = day_cases{cc, 2};
            caseLabel = day_cases{cc, 1};
            col       = cmap(cc, :);

            if camIdx > numel(fr) || ~isfield(fr(camIdx), 'floor_line_combined')
                fprintf('  Case %s: no data for Camera %d\n', caseLabel, camNums(camIdx));
                continue
            end

            pC = fr(camIdx).floor_line_combined;
            pA = fr(camIdx).floor_line_A;
            pB = fr(camIdx).floor_line_B;
            x_full = [1, calW];

            % Combined fit — labelled
            plot(ax, x_full, polyval(pC, x_full), '-', ...
                'Color', col, 'LineWidth', 2, 'DisplayName', caseLabel);

            % A and B fits — unlabelled, lighter
            plot(ax, x_full, polyval(pA, x_full), '--', ...
                'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(ax, x_full, polyval(pB, x_full), ':', ...
                'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');

            % Strip estimates
            good = isfinite(fr(camIdx).floor_y_combined_strips) & ...
                fr(camIdx).strip_weights > 0;
            plot(ax, fr(camIdx).strip_col_centres(good), ...
                fr(camIdx).floor_y_combined_strips(good), ...
                '+', 'Color', col, 'MarkerSize', 8, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end

        % ── Mean and std across lines — only if more than one case ────────
        if sum(valid_cases) > 1

            mean_line = mean(line_vals(valid_cases, :), 1);
            std_line  = std(line_vals(valid_cases, :), 0, 1);
            % ── Mean floor anchor points at image edges ───────────────────────────


            % ±1σ shaded band
            x_band  = [x_eval, fliplr(x_eval)];
            y_band  = [mean_line + std_line, fliplr(mean_line - std_line)];
            fill(ax, x_band, y_band, 'w', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
                'DisplayName', sprintf('\\pm1\\sigma = %.2f px (mean)', mean(std_line)));

            % Mean line
            plot(ax, x_eval, mean_line, 'w-', 'LineWidth', 2, ...
                'DisplayName', sprintf('Mean = %.1f px', mean(mean_line)));

            % Stats — mean std across the image width, and range between cases
            mean_std   = mean(std_line);
            max_spread = max(max(line_vals(valid_cases,:),[],1) - min(line_vals(valid_cases,:),[],1));
            mean_floor = mean(mean_line);

            fprintf('\n  Camera %d  |  Day %s\n',   camNums(camIdx), dayLabel);
            fprintf('    Mean floor : %.2f px\n',   mean_floor);
            fprintf('    Mean StDev : %.2f px\n',   mean_std);
            fprintf('    Max spread : %.2f px\n',   max_spread);

            stats_str = sprintf('Mean=%.1f px  |  \\sigma=%.2f px (mean across width)  |  Max spread (bw leftmost, rightmost pixel=%.1f px', ...
                mean_floor, mean_std, max_spread);
        else
            fprintf('\n  Camera %d  |  Day %s: only 1 case, no std computed.\n', ...
                camNums(camIdx), dayLabel);
            stats_str = sprintf('Single case — Mean=%.1f px', ...
                mean(line_vals(valid_cases,:)));
        end
        % ── Mean floor anchor points at image edges ───────────────────────────
        x_anchors    = [1, calW];
        y_at_anchors = nan(nCases, 2);

        for cc = 1:nCases
            fr = day_cases{cc, 2};
            if camIdx > numel(fr) || ~isfield(fr(camIdx), 'floor_line_combined')
                continue
            end
            y_at_anchors(cc, :) = polyval(fr(camIdx).floor_line_combined, x_anchors);
        end

        valid_rows = ~any(isnan(y_at_anchors), 2);
        y_mean     = mean(y_at_anchors(valid_rows, :), 1);
        y_std      = std( y_at_anchors(valid_rows, :), 0, 1);
        y_range    = range(y_at_anchors(valid_rows, :), 1);

        if sum(valid_rows) == 1
            fprintf('\n  Camera %d — Single case day: using directly.\n', camNums(camIdx));
            fprintf('    x=%-5d  →  y = %7.2f px  (no spread — single case)\n', x_anchors(1), y_mean(1));
            fprintf('    x=%-5d  →  y = %7.2f px  (no spread — single case)\n', x_anchors(2), y_mean(2));
        else
            fprintf('\n  Camera %d — Mean floor from calibration:\n', camNums(camIdx));
            fprintf('    x=%-5d  →  y_mean = %7.2f px  |  std = %.2f px  |  range = %.2f px\n', ...
                x_anchors(1), y_mean(1), y_std(1), y_range(1));
            fprintf('    x=%-5d  →  y_mean = %7.2f px  |  std = %.2f px  |  range = %.2f px\n', ...
                x_anchors(2), y_mean(2), y_std(2), y_range(2));
        end
        fprintf('    → Anchor: [%d %.2f %d %.2f]\n', x_anchors(1), y_mean(1), x_anchors(2), y_mean(2));

        % Store
        floorAnchors(camIdx).x       = x_anchors;
        floorAnchors(camIdx).y_mean  = y_mean;
        floorAnchors(camIdx).y_std   = y_std;
        floorAnchors(camIdx).y_range = y_range;
        floorAnchors(camIdx).nCases  = sum(valid_rows);
        % ── Plot mean floor anchor over calibration TIF ───────────────────────
        plot(ax, x_anchors, y_mean, 'r-o', ...ans*
            'LineWidth', 3, ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', 'r', ...
            'DisplayName', sprintf('Mean floor (std: L=%.1fpx R=%.1fpx)', y_std(1), y_std(2)));

        % ±1σ band around the anchor line
        x_band_anchor = [x_anchors(1), x_anchors(2), x_anchors(2), x_anchors(1)];
        y_band_anchor = [y_mean(1)+y_std(1), y_mean(2)+y_std(2), ...
            y_mean(2)-y_std(2), y_mean(1)-y_std(1)];
        fill(ax, x_band_anchor, y_band_anchor, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        legend(ax, 'show', 'Location', 'northwest', 'FontSize', 8, ...
            'TextColor', 'w', 'Color', [0.2 0.2 0.2]);
        hold(ax, 'off');

        title(ax, sprintf('[%s]  Camera %d — Floor estimates across cases\n%s', ...
            dayLabel, camNums(camIdx), stats_str), 'FontSize', 10);

        % ── Save into the day folder ───────────────────────────────────────
        comp_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        save_name      = sprintf('%s_Cam%d_floor_comparison_%s.png', ...
            dayLabel, camNums(camIdx), comp_timestamp);
        exportgraphics(figure(100*df + camIdx), fullfile(dateFolder, save_name), 'Resolution', 300);
        fprintf('  Saved: %s\n', save_name);
    end
    anchorSavePath = fullfile(dateFolder, sprintf('%s_floorAnchors.mat', dayLabel));
    save(anchorSavePath, 'floorAnchors');
    fprintf('\n  Saved floor anchors: %s\n', anchorSavePath);
end

fprintf('\nAll days processed.\n');


% compareFloorCalibVsAutocorr.m
% Compares the floor location from:
%   1. Calibration polynomial (calib_B) — where is x2=0 in pixel space
%   2. Autocorrelation (floorResults)   — where was the floor detected
%
% One figure per camera per day, saved into the day folder

clc; clear; close all;

% ── User settings ─────────────────────────────────────────────────────────────
rootPath = 'C:\Users\ak1u24\Downloads\PIV_Raw_Images_Check';
camNums  = 1:5;

% ── Discover date folders ─────────────────────────────────────────────────────
dateFolders = dir(rootPath);
dateFolders = dateFolders([dateFolders.isdir]);
dateFolders = dateFolders(~ismember({dateFolders.name}, {'.','..'}));

% ── Loop per day ──────────────────────────────────────────────────────────────
for df = 1:numel(dateFolders)

    dateFolder = fullfile(rootPath, dateFolders(df).name);
    dayLabel   = dateFolders(df).name;
    fprintf('\n========== Day: %s ==========\n', dayLabel);

    % Find CALIBRATION subfolder
    contents = dir(dateFolder);
    contents = contents([contents.isdir]);
    calIdx   = find(~cellfun(@isempty, strfind({contents.name}, 'CALIBRATION')), 1);
    if isempty(calIdx)
        fprintf('  No CALIBRATION folder — skipping.\n');
        continue
    end
    calFolder = fullfile(dateFolder, contents(calIdx).name);

    % Look for calib_caseA mat file in the CALIBRATION folder
    calibMats = dir(fullfile(calFolder, 'calib_caseB_*.mat'));
    if isempty(calibMats)
        fprintf('  No calib_caseB_*.mat found in %s — skipping day.\n', calFolder);
        continue
    end

    % Take most recent if multiple exist
    [~, si]       = sort([calibMats.datenum], 'descend');
    calibFilePath = fullfile(calFolder, calibMats(si(1)).name);
    fprintf('  Loaded calibration: %s\n', calibFilePath);
    tmp   = load(calibFilePath, 'calib_B');
    calib = tmp.calib_B;

    % Find case subfolders with floorResults
    caseDirs = dir(dateFolder);
    caseDirs = caseDirs([caseDirs.isdir]);
    caseDirs = caseDirs(~ismember({caseDirs.name}, {'.','..'}));
    caseDirs = caseDirs(cellfun(@isempty, strfind({caseDirs.name}, 'CALIBRATION')));

    day_cases = {};
    for cd = 1:numel(caseDirs)
        caseFolder = fullfile(dateFolder, caseDirs(cd).name);
        matFiles   = dir(fullfile(caseFolder, 'floorResults_*.mat'));
        if isempty(matFiles), continue; end
        [~, si]  = sort([matFiles.datenum], 'descend');
        loaded   = load(fullfile(caseFolder, matFiles(si(1)).name), 'floorResults');
        day_cases{end+1, 1} = caseDirs(cd).name;
        day_cases{end,   2} = loaded.floorResults;
        fprintf('  Loaded floorResults: %s\n', caseDirs(cd).name);
    end

    if isempty(day_cases)
        fprintf('  No floorResults found — skipping.\n');
        continue
    end

    nCases = size(day_cases, 1);
    cmap   = lines(nCases);

    % Load calibration TIFs once per day
    calTifs  = cell(1, numel(camNums));
    calSizes = zeros(numel(camNums), 2);
    for camIdx = 1:numel(camNums)
        tifPath = fullfile(calFolder, sprintf('camera%d', camNums(camIdx)), 'B00001.tif');
        if isfile(tifPath)
            calTifs{camIdx} = imread(tifPath);
            [calSizes(camIdx,1), calSizes(camIdx,2)] = size(calTifs{camIdx});
            fprintf('  Cam%d: loaded TIF\n', camNums(camIdx));
        else
            calTifs{camIdx} = [];
            calSizes(camIdx,:) = [4600, 5312];
            fprintf('  Cam%d: TIF not found — blank background\n', camNums(camIdx));
        end
    end

    % Per camera figure
    for camIdx = 1:numel(camNums)

        figure(100*df + camIdx); clf;
        ax = gca;

        if ~isempty(calTifs{camIdx})
            calImg = calTifs{camIdx};
            calH   = calSizes(camIdx, 1);
            calW   = calSizes(camIdx, 2);
            imagesc_clim(ax, calImg, 1, 99);
            colormap(ax, 'gray');
        else
            calH = calSizes(camIdx, 1);
            calW = calSizes(camIdx, 2);
        end

        axis(ax, 'image'); axis(ax, 'on');
        xlabel(ax, 'Column (px)'); ylabel(ax, 'Row (px)');
        ax.XLim  = [1, calW];
        ax.YLim  = [1, calH];
        ax.XTick = 0:500:calW;
        ax.YTick = 0:500:calH;
        hold(ax, 'on');

        % Calibration floor line: where is x2=0 in pixel space
        dot_c   = calib.dot_c{camIdx};
        x1_min  = min(dot_c(:,1));
        x1_max  = max(dot_c(:,1));
        x1_eval = linspace(x1_min, x1_max, 500)';

        cal_input    = [x1_eval, zeros(size(x1_eval))];
        floor_px_col = polyvaln(calib.fit_px1{camIdx}.p, cal_input);
        floor_px_row = polyvaln(calib.fit_px2{camIdx}.p, cal_input);



        % Autocorrelation floor lines per case
        for cc = 1:nCases

            fr        = day_cases{cc, 2};
            caseLabel = day_cases{cc, 1};
            col       = cmap(cc, :);

            if camIdx > numel(fr) || ~isfield(fr(camIdx), 'floor_line_combined')
                continue
            end

            pC     = fr(camIdx).floor_line_combined;
            x_full = [1, calW];

            plot(ax, x_full, polyval(pC, x_full), '-', ...
                'Color', col, 'LineWidth', 1, ...
                'DisplayName', sprintf('%s (autocorr)', caseLabel));

            % plot(ax, x_full, polyval(fr(camIdx).floor_line_A, x_full), '--', ...
            %     'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
            % plot(ax, x_full, polyval(fr(camIdx).floor_line_B, x_full), ':', ...
            %     'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');

            good = isfinite(fr(camIdx).floor_y_combined_strips) & ...
                fr(camIdx).strip_weights > 0;
            plot(ax, fr(camIdx).strip_col_centres(good), ...
                fr(camIdx).floor_y_combined_strips(good), ...
                'o', 'Color', col, 'MarkerFaceColor', col ,'MarkerSize', 5, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');

            % Print offset at image centre
            mid_col      = calW / 2;
            autocorr_row = polyval(pC, mid_col);
            calib_row    = interp1(floor_px_col, floor_px_row, mid_col, 'linear', 'extrap');
            offset_px    = autocorr_row - calib_row;
            fprintf('  Cam%d  %s: floor offset = %.1f px  (autocorr - calib)\n', ...
                camNums(camIdx), caseLabel, offset_px);
        end
        % ── Statistics: deviation of each case's floor from calibration ───
        if nCases > 0

            x_eval     = linspace(1, calW, 500);
            calib_eval = interp1(floor_px_col, floor_px_row, x_eval, 'linear', 'extrap');

            dev_all    = nan(nCases, length(x_eval));
            valid_case = false(nCases, 1);

            for cc = 1:nCases
                fr = day_cases{cc, 2};
                if camIdx > numel(fr) || ~isfield(fr(camIdx), 'floor_line_combined')
                    continue
                end
                autocorr_eval  = polyval(fr(camIdx).floor_line_combined, x_eval);
                dev_all(cc, :) = autocorr_eval - calib_eval;
                valid_case(cc) = true;
            end

            % Per-case console output
            fprintf('\n  Camera %d  |  Day %s — Floor deviation (autocorr - calib):\n', ...
                camNums(camIdx), dayLabel);
            for cc = 1:nCases
                if ~valid_case(cc), continue; end
                dev_cc = dev_all(cc, :);
                fprintf('    %-30s  mean=%+.2f px  std=%.2f px  max=%.2f px\n', ...
                    day_cases{cc,1}, mean(dev_cc), std(dev_cc), max(abs(dev_cc)));
            end

            all_devs = dev_all(valid_case, :);

            % std of each case's deviation along x → one value per case
            per_case_std = std(all_devs, 0, 2);       % (nValid × 1)

            % average of those per-case stds → single representative spread
            mean_case_std = mean(per_case_std);        % scalar

            % overall mean deviation (bias) and max
            mean_dev = mean(all_devs(:));
            max_dev  = max(abs(all_devs(:)));

            fprintf('    %-30s  mean=%+.2f px  avg_std=%.2f px  max=%.2f px\n', ...
                'ALL CASES COMBINED', mean_dev, mean_case_std, max_dev);

            % Constant-width ±1σ band centred on calibration line
            x_band = [x_eval, fliplr(x_eval)];
            y_band = [(calib_eval + mean_case_std), ...
                fliplr(calib_eval - mean_case_std)];

            fill(ax, x_band, y_band, 'c', 'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
                'DisplayName', sprintf('Calib \\pm1\\sigma = %.2f px', mean_case_std));

            stats_str = sprintf('mean dev=%+.2f px  |  avg case \\sigma=%.2f px  |  max=%.2f px', ...
                mean_dev, mean_case_std, max_dev);

            title(ax, sprintf('[%s]  Camera %d — Calibration floor vs autocorrelation floor\n%s', ...
                dayLabel, camNums(camIdx), stats_str), 'FontSize', 10);
        end
        plot(ax, floor_px_col, floor_px_row, 'c-', 'LineWidth', 1.5, ...
            'DisplayName', 'Calib x_2=0 (floor)');
        legend(ax, 'show', 'Location', 'northwest', 'FontSize', 8, ...
            'TextColor', 'w', 'Color', [0.2 0.2 0.2]);
        hold(ax, 'off');

        title(ax, sprintf('[%s]  Camera %d — Calibration floor vs autocorrelation floor', ...
            dayLabel, camNums(camIdx)), 'FontSize', 10);

        ts        = datestr(now, 'yyyymmdd_HHMMSS');
        save_name = sprintf('%s_Cam%d_floor_calib_vs_autocorr_%s.png', ...
            dayLabel, camNums(camIdx), ts);
        exportgraphics(figure(100*df + camIdx), fullfile(dateFolder, save_name), ...
            'Resolution', 300);
        fprintf('  Saved: %s\n', save_name);
    end
end

fprintf('\nDone.\n');

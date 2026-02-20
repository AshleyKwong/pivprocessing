 % =========================================================================
%  PIV_FloorCorrection.m
%  Floor (y=0) correction for multi-camera PIV data.
%
%  Workflow:
%   1. User inputs: savedPath, nImages, method (i/e), nCams, nPasses
%   2. For each camera: loads coordinates.mat, runs floor detection,
%      computes per-column y-shift from linear fit (average of Method A & C),
%      applies shift to all PIV passes, plots before/after.
%
%  Ashley Kwong — University of Southampton
%  February 2026
% =========================================================================

clc; clear; close all;

%% =========================================================================
%  SECTION 1: User Inputs
% =========================================================================

savedPath = input('Enter savedPath (e.g. E:\\ProcessedPIV_case2fullpipe): ', 's');
nImages   = input('Enter number of images (e.g. 150): ');
nPasses   = input('Enter number of PIV passes (e.g. 3 → uses pass 3): ');
nCams     = input('Enter number of cameras (e.g. 4): ');

methodChar = '';
while ~ismember(methodChar, {'i','e'})
    methodChar = lower(strtrim(input('Method: instantaneous (i) or ensemble (e)? ', 's')));
end
if strcmp(methodChar, 'i')
    methodStr = 'instantaneous';
else
    methodStr = 'ensemble';
end

fprintf('\nSettings:\n');
fprintf('  Path    : %s\n', savedPath);
fprintf('  Images  : %d\n', nImages);
fprintf('  Passes  : %d (using pass %d)\n', nPasses, nPasses);
fprintf('  Cameras : %d\n', nCams);
fprintf('  Method  : %s\n\n', methodStr);

%% =========================================================================
%  SECTION 2: Find loop folder
% =========================================================================

loopDir = dir(fullfile(savedPath, 'loop=*'));
if isempty(loopDir)
    error('No loop= folder found in: %s', savedPath);
end
if numel(loopDir) > 1
    warning('Multiple loop= folders found — using the first: %s', loopDir(1).name);
end
loopPath = fullfile(savedPath, loopDir(1).name);
fprintf('Using loop folder: %s\n', loopPath);

%% =========================================================================
%  SECTION 3: Loop over cameras
% =========================================================================

for camIdx = 1:nCams

    camName  = sprintf('Cam%d', camIdx);
    coordFile = fullfile(loopPath, 'calibrated_piv', num2str(nImages), ...
                          camName, methodStr, 'coordinates.mat');

    fprintf('\n========== %s ==========\n', camName);

    if ~isfile(coordFile)
        warning('coordinates.mat not found for %s:\n  %s\n  Skipping.', camName, coordFile);
        continue;
    end

    % ------------------------------------------------------------------
    %  Load coordinates
    % ------------------------------------------------------------------
    fprintf('Loading: %s\n', coordFile);
    raw = load(coordFile);

    % coordinates is a cell array — field name may vary; handle both cases
    if isfield(raw, 'coordinates')
        coordinates = raw.coordinates;
    else
        fnames = fieldnames(raw);
        coordinates = raw.(fnames{1});
        warning('%s: expected field ''coordinates'', found ''%s''.', camName, fnames{1});
    end

    if nPasses > numel(coordinates)
        warning('%s: requested pass %d but only %d passes in file. Using last pass.', ...
            camName, nPasses, numel(coordinates));
        finalPass = numel(coordinates);
    else
        finalPass = nPasses;
    end

    x_orig = coordinates{finalPass}.x;   % [nY x nX]  mm
    y_orig = coordinates{finalPass}.y;   % [nY x nX]  mm

    x_vec = x_orig(1, :);   % row vector of x positions
    y_vec = y_orig(:, 1);   % column vector of y positions

    % ------------------------------------------------------------------
    %  Prompt for mean U field (needed for floor detection)
    % ------------------------------------------------------------------
    fprintf('  For floor detection, a mean-U field is required.\n');
    uFile = input(sprintf('  Enter full path to mean-U .mat for %s\n  (or press Enter to skip floor detection for this camera): ', camName), 's');

    if isempty(strtrim(uFile))
        warning('%s: no U field provided — skipping floor detection and coordinate shift.', camName);
        continue;
    end

    if ~isfile(uFile)
        warning('%s: file not found: %s\n  Skipping.', camName, uFile);
        continue;
    end

    rawU = load(uFile);
    % Try to auto-detect the U field
    fnames = fieldnames(rawU);
    % Prefer a field that looks like U (case-insensitive 'u' or 'meanU' etc.)
    uFieldIdx = find(strcmpi(fnames, 'U') | strcmpi(fnames, 'meanU') | ...
                     strcmpi(fnames, 'Umean') | strcmpi(fnames, 'ptMeanU_sub'), 1);
    if isempty(uFieldIdx)
        fprintf('  Available fields in %s:\n', uFile);
        disp(fnames);
        uFieldName = input('  Enter the field name for mean-U: ', 's');
    else
        uFieldName = fnames{uFieldIdx};
    end
    U_field = single(rawU.(uFieldName));

    % ------------------------------------------------------------------
    %  Run floor detection (wrapped from PIV_floordeviation.m)
    % ------------------------------------------------------------------
    fprintf('  Running floor detection on %s...\n', camName);
    [p_A, p_C, floor_y_A, floor_y_C] = PIV_detectFloor(x_vec, y_vec, U_field);

    % ------------------------------------------------------------------
    %  Average the two linear fits to get floor(x) line
    % ------------------------------------------------------------------
    % floor_line(x) = mean slope * x + mean intercept
    p_floor = 0.5 * (p_A + p_C);   % average of Method A and Method C fits

    floor_at_x = polyval(p_floor, x_vec);   % [1 x nX] — detected floor y per column

    fprintf('  Floor fit (avg A+C): slope = %.6f mm/mm, intercept = %.4f mm\n', ...
        p_floor(1), p_floor(2));

    % ------------------------------------------------------------------
    %  Shift: subtract floor_at_x from y so that floor = 0
    %  Each column j gets shifted by -floor_at_x(j)
    % ------------------------------------------------------------------
    % y_new(i,j) = y_orig(i,j) - floor_at_x(j)
    shift_matrix = repmat(floor_at_x, size(y_orig, 1), 1);  % [nY x nX]

    y_new = y_orig - shift_matrix;

    % ------------------------------------------------------------------
    %  Apply the same shift to ALL passes
    % ------------------------------------------------------------------
    coordinates_corrected = coordinates;   % copy all passes

    for p = 1:numel(coordinates)
        yp = coordinates{p}.y;
        % The grid structure is the same across passes; apply same column shifts
        shift_p = repmat(floor_at_x, size(yp, 1), 1);
        coordinates_corrected{p}.y = yp - shift_p;
        % x is unchanged
    end

    % ------------------------------------------------------------------
    %  Plot: before vs after floor correction
    % ------------------------------------------------------------------
    figure('Name', sprintf('%s — Floor Correction', camName), ...
           'Position', [50 50 1600 900]);

    % --- Top row: old y field (imagesc of U overlaid with old floor line)
    ax1 = subplot(2,2,1);
    imagesc(x_vec, y_vec, U_field);
    set(gca,'YDir','normal'); colormap(gca, jet); colorbar;
    hold on;
    plot(x_vec, floor_y_A,             'r.', 'MarkerSize', 4, 'DisplayName', 'Floor A raw');
    plot(x_vec, floor_y_C,             'g.', 'MarkerSize', 4, 'DisplayName', 'Floor C raw');
    plot(x_vec, polyval(p_A, x_vec),   'r-', 'LineWidth', 2,  'DisplayName', 'Floor A fit');
    plot(x_vec, polyval(p_C, x_vec),   'g-', 'LineWidth', 2,  'DisplayName', 'Floor C fit');
    plot(x_vec, floor_at_x,            'w--','LineWidth', 2,   'DisplayName', 'Floor avg fit');
    yline(0, 'k--', 'LineWidth', 1.2);
    xlabel('X (mm)'); ylabel('Y (mm)');
    title(sprintf('%s — BEFORE correction', camName));
    legend('Location', 'northeast', 'FontSize', 7);

    % --- Top right: after correction (y_new)
    y_vec_new = y_new(:, 1);   % new y column vector

    ax2 = subplot(2,2,2);
    imagesc(x_vec, y_vec_new, U_field);
    set(gca,'YDir','normal'); colormap(gca, jet); colorbar;
    hold on;
    yline(0, 'w--', 'LineWidth', 2, 'DisplayName', 'y = 0');
    xlabel('X (mm)'); ylabel('Y (mm)');
    title(sprintf('%s — AFTER correction', camName));
    legend('Location', 'northeast', 'FontSize', 7);

    % --- Bottom left: old y_vec
    ax3 = subplot(2,2,3);
    plot(x_vec, y_orig(1,:),   'b-', 'LineWidth', 1.5, 'DisplayName', 'Old y_{min} (row 1)'); hold on;
    plot(x_vec, y_orig(end,:), 'b--','LineWidth', 1.5, 'DisplayName', 'Old y_{max} (last row)');
    plot(x_vec, floor_at_x,    'r-', 'LineWidth', 2,   'DisplayName', 'Detected floor');
    yline(0, 'k:', 'LineWidth', 1);
    xlabel('X (mm)'); ylabel('Y (mm)');
    title('Window centres — BEFORE'); legend; grid on;

    % --- Bottom right: new y_vec
    ax4 = subplot(2,2,4);
    plot(x_vec, y_new(1,:),   'b-', 'LineWidth', 1.5, 'DisplayName', 'New y_{min} (row 1)'); hold on;
    plot(x_vec, y_new(end,:), 'b--','LineWidth', 1.5, 'DisplayName', 'New y_{max} (last row)');
    yline(0, 'r--', 'LineWidth', 2, 'DisplayName', 'y = 0 floor');
    xlabel('X (mm)'); ylabel('Y (mm)');
    title('Window centres — AFTER'); legend; grid on;

    linkaxes([ax1 ax2], 'xy');
    linkaxes([ax3 ax4], 'xy');
    sgtitle(sprintf('%s — PIV Floor Correction (avg Method A+C)', camName), 'FontSize', 13);

    % ------------------------------------------------------------------
    %  Save corrected coordinates back to file (optional prompt)
    % ------------------------------------------------------------------
    doSave = lower(strtrim(input(sprintf('\n  Save corrected coordinates for %s? (y/n): ', camName), 's')));
    if strcmp(doSave, 'y')
        coordinates = coordinates_corrected; %#ok<NASGU>
        save(coordFile, 'coordinates', '-v7.3');
        fprintf('  Saved corrected coordinates → %s\n', coordFile);
    else
        fprintf('  Skipped saving for %s.\n', camName);
    end

end   % end camera loop

fprintf('\nAll cameras processed.\n');


%% =========================================================================
%  LOCAL FUNCTION: PIV_detectFloor
%  Extracted and wrapped from PIV_floordeviation.m
%
%  Inputs:
%    x_vec   [1 x nX]  X positions (mm)
%    y_vec   [nY x 1]  Y positions (mm)
%    U_pt    [nY x nX] Mean streamwise velocity (m/s)
%
%  Outputs:
%    p_A      [1x2] polyfit coefficients for Method A floor line
%    p_C      [1x2] polyfit coefficients for Method C floor line
%    floor_y_A [1 x nX] raw Method A detections (mm, NaN where invalid)
%    floor_y_C [1 x nX] raw Method C detections (mm, NaN where invalid)
% =========================================================================

function [p_A, p_C, floor_y_A, floor_y_C] = PIV_detectFloor(x_vec, y_vec, U_pt)

U_pt = single(U_pt);

% --- 1. Search band ------------------------------------------------------
y_lo = -1;
y_hi = 15;

band_mask = y_vec >= y_lo & y_vec <= y_hi;
band_rows = find(band_mask);
y_band    = y_vec(band_rows);
nX        = size(U_pt, 2);
U_band    = U_pt(band_rows, :);

% --- 2. Tuning parameters ------------------------------------------------
min_valid_pts   = 5;
flat_threshold  = 0.5;    % m/s
near_zero_thresh= 0.3;    % m/s
valid_frac_min  = 0.6;

% Method A
pos_thresh = 0.5;   % (m/s)/mm
sustained  = 2;

% Method C
alpha = 0.10;

% --- 3. Preallocate ------------------------------------------------------
floor_y_A = nan(1, nX);
floor_y_C = nan(1, nX);

% --- 4. Column-wise detection --------------------------------------------
for col = 1:nX

    u_col = double(U_band(:, col));
    y_col = double(y_band);
    valid = ~isnan(u_col);

    if sum(valid) < min_valid_pts,                      continue; end
    if sum(valid)/length(u_col) < valid_frac_min,       continue; end

    u_v = u_col(valid);
    y_v = y_col(valid);

    if mean(abs(u_v)) < near_zero_thresh,               continue; end
    if (max(u_v) - min(u_v)) < flat_threshold,          continue; end

    [~, idx_min_u] = min(u_v);
    if idx_min_u >= length(u_v) - sustained,            continue; end

    u_sorted = sort(u_v, 'descend');
    u_fs     = mean(u_sorted(1:max(1, round(0.1*length(u_sorted)))));

    % --- Method A: first sustained positive dU/dy above U minimum --------
    dudy       = gradient(u_v, y_v);
    dudy_above = dudy(idx_min_u:end);
    y_above    = y_v(idx_min_u:end);

    for r = 1:length(dudy_above) - sustained
        if all(dudy_above(r:r+sustained) > pos_thresh)
            if r == 1
                floor_y_A(col) = y_above(r);
            else
                floor_y_A(col) = interp1(dudy_above(r-1:r), y_above(r-1:r), ...
                                         pos_thresh, 'linear', 'extrap');
            end
            break;
        end
    end

    % --- Method C: U recovers to alpha fraction above minimum ------------
    u_thresh  = min(u_v) + alpha * (u_fs - min(u_v));
    u_above_c = u_v(idx_min_u:end);
    y_above_c = y_v(idx_min_u:end);

    for r = 2:length(u_above_c)
        if u_above_c(r) >= u_thresh
            floor_y_C(col) = interp1(u_above_c(r-1:r), y_above_c(r-1:r), ...
                                     u_thresh, 'linear', 'extrap');
            break;
        end
    end

end   % col loop

% --- 5. Linear fits ------------------------------------------------------
valid_A = ~isnan(floor_y_A);
if sum(valid_A) >= 2
    p_A = polyfit(x_vec(valid_A), floor_y_A(valid_A), 1);
else
    warning('PIV_detectFloor: too few valid Method A detections — returning zero fit.');
    p_A = [0, 0];
end

valid_C = ~isnan(floor_y_C);
if sum(valid_C) >= 2
    p_C = polyfit(x_vec(valid_C), floor_y_C(valid_C), 1);
else
    warning('PIV_detectFloor: too few valid Method C detections — returning zero fit.');
    p_C = [0, 0];
end

% --- 6. Console report ---------------------------------------------------
pitch_A_mm  = p_A(1) * (x_vec(end) - x_vec(1));
pitch_A_deg = atand(p_A(1));
pitch_C_mm  = p_C(1) * (x_vec(end) - x_vec(1));
pitch_C_deg = atand(p_C(1));

fprintf('\n--- Wall pitch: Method A (sustained dU/dy) ---\n');
fprintf('  Slope   : %.6f mm/mm\n', p_A(1));
fprintf('  Rise    : %.3f mm over %.1f mm span\n', pitch_A_mm, x_vec(end)-x_vec(1));
fprintf('  Angle   : %.5f deg\n', pitch_A_deg);
fprintf('  Valid   : %d / %d (%.1f%%)\n', sum(valid_A), nX, 100*sum(valid_A)/nX);

fprintf('\n--- Wall pitch: Method C (U recovery alpha=%.2f) ---\n', alpha);
fprintf('  Slope   : %.6f mm/mm\n', p_C(1));
fprintf('  Rise    : %.3f mm over %.1f mm span\n', pitch_C_mm, x_vec(end)-x_vec(1));
fprintf('  Angle   : %.5f deg\n', pitch_C_deg);
fprintf('  Valid   : %d / %d (%.1f%%)\n', sum(valid_C), nX, 100*sum(valid_C)/nX);

% --- 7. Diagnostic figure -----------------------------------------------
nX_plot   = size(U_pt, 2);
check_cols = round(linspace(1, nX_plot, 6));

figure('Name', 'Floor detection diagnostics', 'Position', [100 100 1400 600]);
for k = 1:numel(check_cols)
    col   = check_cols(k);
    u_col = double(U_band(:, col));
    valid = ~isnan(u_col);
    u_v   = u_col(valid);
    y_v   = double(y_band(valid));
    dudy  = gradient(u_v, y_v);

    subplot(2, numel(check_cols), k);
    plot(u_v, y_v, 'b-o', 'MarkerSize', 3); hold on;
    if ~isnan(floor_y_A(col)), yline(floor_y_A(col), 'r--', 'LineWidth', 1.5); end
    if ~isnan(floor_y_C(col)), yline(floor_y_C(col), 'g--', 'LineWidth', 1.5); end
    xlabel('U (m/s)');
    if k==1, ylabel('Y (mm)'); end
    title(sprintf('x=%.0f mm', x_vec(col)), 'FontSize', 9);
    ylim([y_lo-0.5, y_hi+0.5]); grid on;

    subplot(2, numel(check_cols), k + numel(check_cols));
    plot(dudy, y_v, 'r-o', 'MarkerSize', 3); hold on;
    xline(0, 'k--', 'LineWidth', 1.2);
    xline(pos_thresh, 'k:', 'LineWidth', 1.0);
    if ~isnan(floor_y_A(col)), yline(floor_y_A(col), 'r--', 'LineWidth', 1.5); end
    if ~isnan(floor_y_C(col)), yline(floor_y_C(col), 'g--', 'LineWidth', 1.5); end
    xlabel('dU/dy [(m/s)/mm]');
    if k==1, ylabel('Y (mm)'); end
    ylim([y_lo-0.5, y_hi+0.5]); grid on;
end
sgtitle('Floor detection diagnostics — red: Method A | green: Method C', 'FontSize', 11);

end   % PIV_detectFloor

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
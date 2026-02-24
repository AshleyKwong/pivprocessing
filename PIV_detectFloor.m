% =========================================================================
% PIV_detectFloor.m
%
% Inputs:
%   xvec  [1  x nX]  X coordinate vector (mm) — one value per column
%   yvec  [nY x 1 ]  Y coordinate vector (mm) — one value per row
%   Ux    [nY x nX]  Mean streamwise velocity field (m/s)
%
% Outputs:
%   p_A       [1x2]   polyfit coefficients for Method A floor line
%   p_C       [1x2]   polyfit coefficients for Method C floor line (primary)
%   floor_y_A [1 x nX] per-column Method A detections (mm, NaN where invalid)
%   floor_y_C [1 x nX] per-column Method C detections (mm, NaN where invalid)
%   rise_A_mm [1 x nX] per-column Method C delta y   (mm, NaN where invalid)
%   rise_C_mm [1 x nX] per-column Method A delta y   (mm, NaN where invalid)

% =========================================================================
function [p_A, p_C, floor_y_A, floor_y_C, rise_A_mm, rise_C_mm] = PIV_detectFloor(xvec, yvec, Ux, y_lo, y_hi)
U_pt = single(Ux);
nX   = size(U_pt, 2);

% --- 1. Search band -------------------------------------------------------
if nargin < 4 || isempty(y_lo), y_lo = -0.5; end   % default fallback
if nargin < 5 || isempty(y_hi), y_hi = 10;   end


% --- 2. Tuning parameters -------------------------------------------------
min_valid_pts   = 5;
flat_threshold  = 0.5;   % m/s
near_zero_thresh= 0.3;   % m/s
valid_frac_min  = 0.6;
pos_thresh      = 0.5;   % (m/s)/mm  — Method A
sustained       = 2;     % consecutive points
alpha           = 0.10;  % Method C recovery fraction

% --- 3. Preallocate -------------------------------------------------------
floor_y_A = nan(1, nX);
floor_y_C = nan(1, nX);
floor_x_A = nan(1, nX);
floor_x_C = nan(1, nX);

% --- 4. Column-wise detection ---------------------------------------------
y_full = double(yvec);                      % [nY x 1] — shared across all cols
band_mask = y_full >= y_lo & y_full <= y_hi;
band_rows = find(band_mask);
y_band    = y_full(band_rows);              % precompute once — same every col

for col = 1:nX

    u_col = double(Ux(band_rows, col));
    valid = ~isnan(u_col);

    % Data quality gates
    if sum(valid) < min_valid_pts,                       continue; end
    if sum(valid)/length(u_col) < valid_frac_min,        continue; end

    u_v = u_col(valid);
    y_v = y_band(valid);

    if mean(abs(u_v)) < near_zero_thresh,                continue; end
    if (max(u_v) - min(u_v)) < flat_threshold,           continue; end

    [min_u, idx_min_u] = min(u_v);
    if idx_min_u >= length(u_v) - 1,                     continue; end

    n_fs  = max(2, round(0.1 * length(u_v)));
    u_sorted = sort(u_v, 'descend');
    u_fs     = mean(u_sorted(1:max(1, round(0.1*length(u_sorted)))));

    % --- METHOD A: first sustained positive dU/dy above U minimum ---------
    if idx_min_u <= length(u_v) - sustained
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
                floor_x_A(col) = xvec(col);   % scalar — same for whole column
                break;
            end
        end
    end

    % --- METHOD C (primary): U recovers to alpha above minimum ------------
    u_thresh   = min_u + alpha * (u_fs - min_u);
    u_above_c  = u_v(idx_min_u:end);
    y_above_c  = y_v(idx_min_u:end);

    for r = 2:length(u_above_c)
        if u_above_c(r) >= u_thresh
            floor_y_C(col) = interp1(u_above_c(r-1:r), y_above_c(r-1:r), ...
                                     u_thresh, 'linear', 'extrap');
            floor_x_C(col) = xvec(col);        % scalar — same for whole column
            break;
        end
    end

end % col loop

% --- 5. Linear fits -------------------------------------------------------
valid_A = ~isnan(floor_y_A);
if sum(valid_A) >= 2
    p_A = polyfit(floor_x_A(valid_A), floor_y_A(valid_A), 1);
else
    warning('PIV_detectFloor: too few valid Method A detections — returning zero fit.');
    p_A = [0, 0];
end

valid_C = ~isnan(floor_y_C);
if sum(valid_C) >= 2
    p_C = polyfit(floor_x_C(valid_C), floor_y_C(valid_C), 1);
else
    warning('PIV_detectFloor: too few valid Method C detections — returning zero fit.');
    p_C = [0, 0];
end

% --- 6. Console report ----------------------------------------------------
x_span_A = floor_x_A(valid_A);
x_span_C = floor_x_C(valid_C);
rise_A_mm = p_A(1) * (max(x_span_A) - min(x_span_A));
rise_C_mm = p_C(1) * (max(x_span_C) - min(x_span_C));

fprintf('\n--- Wall pitch: Method A (sustained dU/dy) ---\n');
fprintf('  Slope : %.6f mm/mm\n',  p_A(1));
fprintf('  Rise  : %.3f mm over %.1f mm\n', rise_A_mm, max(x_span_A)-min(x_span_A));
fprintf('  Angle : %.5f deg\n',    atand(p_A(1)));
fprintf('  Valid : %d / %d (%.1f%%)\n', sum(valid_A), nX, 100*sum(valid_A)/nX);

fprintf('\n--- Wall pitch: Method C (U recovery alpha=%.2f) ---\n', alpha);
fprintf('  Slope : %.6f mm/mm\n',  p_C(1));
fprintf('  Rise  : %.3f mm over %.1f mm\n', rise_C_mm, max(x_span_C)-min(x_span_C));
fprintf('  Angle : %.5f deg\n',    atand(p_C(1)));
fprintf('  Valid : %d / %d (%.1f%%)\n', sum(valid_C), nX, 100*sum(valid_C)/nX);

% --- 7. Diagnostic figure (6 sample columns) ------------------------------
check_cols = round(linspace(1, nX, 6));
figure('Name','PIV_detectFloor diagnostics','Position',[100 100 1400 600]);

for k = 1:numel(check_cols)
    col = check_cols(k);

    u_col = double(Ux(band_rows, col));
    valid = ~isnan(u_col);
    u_v   = u_col(valid);
    y_v   = y_band(valid);
    dudy  = gradient(u_v, y_v);

    subplot(2, numel(check_cols), k);
    plot(u_v, y_v, 'b-o', 'MarkerSize', 3); hold on;
    if ~isnan(floor_y_A(col)), yline(floor_y_A(col), 'r--', 'LineWidth', 1.5); end
    if ~isnan(floor_y_C(col)), yline(floor_y_C(col), 'g--', 'LineWidth', 1.5); end
    xlabel('U (m/s)');
    if k == 1, ylabel('Y (mm)'); end
    title(sprintf('x=%.0f mm', xvec(col)), 'FontSize', 9);   % xvec(col) directly
    ylim([y_lo-0.5, y_hi+0.5]); grid on;

    subplot(2, numel(check_cols), k + numel(check_cols));
    plot(dudy, y_v, 'r-o', 'MarkerSize', 3); hold on;
    xline(0, 'k--', 'LineWidth', 1.2);
    xline(pos_thresh, 'k:', 'LineWidth', 1.0);
    if ~isnan(floor_y_A(col)), yline(floor_y_A(col), 'r--', 'LineWidth', 1.5); end
    if ~isnan(floor_y_C(col)), yline(floor_y_C(col), 'g--', 'LineWidth', 1.5); end
    xlabel('dU/dy [(m/s)/mm]');
    if k == 1, ylabel('Y (mm)'); end
    ylim([y_lo-0.5, y_hi+0.5]); grid on;
end

sgtitle('Floor detection: U (top) | dU/dy (bottom) — red: Method A | green: Method C', ...
    'FontSize', 11);

end % PIV_detectFloor

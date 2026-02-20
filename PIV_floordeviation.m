% =========================================================================
%% FLOOR DETECTION: top-of-wall via methods A and C
%% Method A: first sustained positive dU/dy above U minimum
%% Method C: U recovery to alpha fraction above minimum
% =========================================================================

pt_x = pivData.coordinates(final_pass_index).x;
pt_y = pivData.coordinates(final_pass_index).y;
U_pt = single(ptMeanU_sub);

x_vec_pt = pt_x(1,:);
y_vec    = pt_y(:,1);

% --- 1. Search band ------------------------------------------------------
y_lo = -1;
y_hi = 15;
band_mask = y_vec >= y_lo & y_vec <= y_hi;
band_rows = find(band_mask);
y_band    = y_vec(band_rows);
nX        = size(U_pt, 2);
U_band    = U_pt(band_rows, :);

% --- 2. Tuning parameters ------------------------------------------------
min_valid_pts    = 5;
flat_threshold   = 0.5;    % m/s  — reject flat velocity profiles
near_zero_thresh = 0.3;    % m/s  — reject near-zero columns
valid_frac_min   = 0.6;    % fraction of band rows that must be non-NaN

% Method A
pos_thresh = 0.5;          % (m/s)/mm — dU/dy must exceed this
sustained  = 2;            % consecutive points dU/dy must stay positive

% Method C
alpha      = 0.10;          % U recovery fraction (10% above minimum)
% in method C we find the min U for that region 

% --- 3. Preallocate ------------------------------------------------------
floor_y_A = nan(1, nX);
floor_y_C = nan(1, nX);

% --- 4. Column-wise detection --------------------------------------------
for col = 1:nX
    u_col = double(U_band(:, col));
    y_col = double(y_band);

    valid = ~isnan(u_col);

    % Data quality rejections
    if sum(valid) < min_valid_pts,                    continue; end
    if sum(valid)/length(u_col) < valid_frac_min,     continue; end

    u_v = u_col(valid);
    y_v = y_col(valid);

    if mean(abs(u_v)) < near_zero_thresh,             continue; end
    if (max(u_v) - min(u_v)) < flat_threshold,        continue; end

    % Locate U minimum (bottom of wall trough)
    [min_u, idx_min_u] = min(u_v);

    % Nothing to search if minimum is at the top of the band
    if idx_min_u >= length(u_v) - sustained,          continue; end

    % Freestream proxy: mean of top 10% of U values in band
    u_sorted = sort(u_v, 'descend');
    u_fs     = mean(u_sorted(1:max(1, round(0.1*length(u_sorted)))));

    % -----------------------------------------------------------------
    % METHOD A: first sustained positive dU/dy above U minimum
    % -----------------------------------------------------------------
    dudy         = gradient(u_v, y_v);
    dudy_above   = dudy(idx_min_u:end);
    y_above      = y_v(idx_min_u:end);

    for r = 1:length(dudy_above) - sustained
        if all(dudy_above(r:r+sustained) > pos_thresh)
            if r == 1
                floor_y_A(col) = y_above(r);
            else
                % Linear interpolation to exact threshold crossing
                floor_y_A(col) = interp1( ...
                    dudy_above(r-1:r), y_above(r-1:r), pos_thresh, 'linear', 'extrap');
            end
            break;
        end
    end

    % -----------------------------------------------------------------
    % METHOD C: U recovers to alpha fraction above minimum
    % -----------------------------------------------------------------
    u_thresh     = min_u + alpha * (u_fs - min_u);
    u_above      = u_v(idx_min_u:end);
    y_above_c    = y_v(idx_min_u:end);

    for r = 2:length(u_above)
        if u_above(r) >= u_thresh
            % Linear interpolation for sub-grid precision
            floor_y_C(col) = interp1( ...
                u_above(r-1:r), y_above_c(r-1:r), u_thresh, 'linear', 'extrap');
            break;
        end
    end
end

% --- 5. Linear fits ------------------------------------------------------
% Method A
valid_A  = ~isnan(floor_y_A);
p_A      = polyfit(x_vec_pt(valid_A), floor_y_A(valid_A), 1);
pitch_A_mm  = p_A(1) * (x_vec_pt(end) - x_vec_pt(1));
pitch_A_deg = atand(p_A(1));

% Method C
valid_C  = ~isnan(floor_y_C);
p_C      = polyfit(x_vec_pt(valid_C), floor_y_C(valid_C), 1);
pitch_C_mm  = p_C(1) * (x_vec_pt(end) - x_vec_pt(1));
pitch_C_deg = atand(p_C(1));

fprintf('\n--- Wall pitch: Method A (sustained dU/dy) ---\n');
fprintf('  Slope      : %.6f mm/mm\n', p_A(1));
fprintf('  Total rise : %.3f mm  over %.1f mm span\n', pitch_A_mm, x_vec_pt(end)-x_vec_pt(1));
fprintf('  Angle      : %.5f deg\n',   pitch_A_deg);
fprintf('  Valid cols : %d / %d  (%.1f%%)\n', sum(valid_A), nX, 100*sum(valid_A)/nX);

fprintf('\n--- Wall pitch: Method C (U recovery alpha=%.2f) ---\n', alpha);
fprintf('  Slope      : %.6f mm/mm\n', p_C(1));
fprintf('  Total rise : %.3f mm  over %.1f mm span\n', pitch_C_mm, x_vec_pt(end)-x_vec_pt(1));
fprintf('  Angle      : %.5f deg\n',   pitch_C_deg);
fprintf('  Valid cols : %d / %d  (%.1f%%)\n', sum(valid_C), nX, 100*sum(valid_C)/nX);

% --- 6. Diagnostic profiles (6 sample columns) --------------------------
check_cols = round(linspace(1, nX, 6));
figure('Position', [100 100 1400 600]);
for k = 1:length(check_cols)
    col   = check_cols(k);
    u_col = double(U_band(:, col));
    valid = ~isnan(u_col);
    u_v   = u_col(valid);
    y_v   = double(y_band(valid));
    dudy  = gradient(u_v, y_v);

    subplot(2, length(check_cols), k);
    plot(u_v, y_v, 'b-o', 'MarkerSize', 3); hold on;
    if ~isnan(floor_y_A(col))
        yline(floor_y_A(col), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Method A');
    end
    if ~isnan(floor_y_C(col))
        yline(floor_y_C(col), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Method C');
    end
    xlabel('U (m/s)');
    if k==1, ylabel('Y (mm)'); end
    title(sprintf('x=%.0f mm', x_vec_pt(col)), 'FontSize', 9);
    ylim([y_lo-0.5, y_hi+0.5]); grid on;

    subplot(2, length(check_cols), k + length(check_cols));
    plot(dudy, y_v, 'r-o', 'MarkerSize', 3); hold on;
    xline(0,         'k--', 'LineWidth', 1.2);
    xline(pos_thresh,'k:',  'LineWidth', 1.0);
    if ~isnan(floor_y_A(col))
        yline(floor_y_A(col), 'r--', 'LineWidth', 1.5);
    end
    if ~isnan(floor_y_C(col))
        yline(floor_y_C(col), 'g--', 'LineWidth', 1.5);
    end
    xlabel('dU/dy [(m/s)/mm]');
    if k==1, ylabel('Y (mm)'); end
    ylim([y_lo-0.5, y_hi+0.5]); grid on;
end
sgtitle('Diagnostic: U (top) and dU/dy (bottom) — red=Method A, green=Method C', ...
    'FontSize', 11);
 
% --- 7. Main comparison plot ---------------------------------------------
figure('Position', [50 50 1600 1000]);

ax1 = subplot(2,1,1);
imagesc(x_vec_pt, y_vec, U_pt);
set(gca,'YDir','normal'); axis image; colormap(gca, jet);
colorbar; clim([0 30]); hold on;
plot(x_vec_pt, floor_y_A,             'r.',  'MarkerSize', 4);
plot(x_vec_pt, floor_y_C,             'g.',  'MarkerSize', 4);
plot(x_vec_pt, polyval(p_A,x_vec_pt), 'r-',  'LineWidth', 2);
plot(x_vec_pt, polyval(p_C,x_vec_pt), 'g-',  'LineWidth', 2);
yline(y_lo, 'k:', 'LineWidth', 1);
yline(y_hi, 'k:', 'LineWidth', 1);
xlabel('X (mm)'); ylabel('Y (mm)');
title('Wall detections and linear fits');
legend('Method A raw','Method C raw','Method A fit','Method C fit', ...
    'Location','northeast');
ax2 = subplot(2,1,2);
hold on; 
resid_A =  polyval(p_A, x_vec_pt);%floor_y_A - polyval(p_A, x_vec_pt);
plot(x_vec_pt, resid_A, 'r.', 'MarkerSize', 4, 'DisplayName',sprintf('Method A residual  (rise=%.3f mm, %.5f°)', pitch_A_mm, pitch_A_deg)); hold on;
resid_C = polyval(p_C, x_vec_pt);%floor_y_C - polyval(p_C, x_vec_pt);
plot(x_vec_pt, resid_C, 'g.', 'MarkerSize', 4, 'DisplayName', sprintf('Method C residual  (rise=%.3f mm, %.5f°)', pitch_C_mm, pitch_C_deg)); hold on;
yline( 0,   'k--', 'LineWidth', 1.0);
yline( 0.5, 'k:',  'LineWidth', 0.8);
yline(-0.5, 'k:',  'LineWidth', 0.8);
xlabel('X (mm)'); ylabel('Residual (mm)');
grid on;
hold off; 
xlim([0 1400]); 
% ax3 = subplot(3,1,3);
% resid_C = polyval(p_C, x_vec_pt);%floor_y_C - polyval(p_C, x_vec_pt);
% plot(x_vec_pt, resid_C, 'g.', 'MarkerSize', 4); hold on;
% yline( 0,   'k--', 'LineWidth', 1.0);
% yline( 0.5, 'k:',  'LineWidth', 0.8);
% yline(-0.5, 'k:',  'LineWidth', 0.8);
% xlabel('X (mm)'); ylabel('Residual (mm)');
% title();
% grid on;

linkaxes([ax1, ax2], 'x');
sgtitle('Wall location and pitch: Method A vs Method C', 'FontSize', 14);

% --- 8. Cleanup ----------------------------------------------------------
clear U_band band_rows band_mask y_band;
clear floor_y_A floor_y_C valid_A valid_C resid_A resid_C;
clear p_A p_C pitch_A_mm pitch_A_deg pitch_C_mm pitch_C_deg;

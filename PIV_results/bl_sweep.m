% =========================================================================
% bl_sweep.m
% Ashley Kwong — March 2026
% Boundary layer integral parameters across the full PIV measurement window.
% Outputs: x, delta99 (hybrid/vel/vinuesa), delta_star, theta, H, U99, U_inf
%          + RMS uncertainty on each quantity from column-averaging spread.
% Vinuesa et al. (2016) method for delta99 detection.
% Hybrid: Vinuesa where diagnostic crosses 0.02, velocity-based fallback.
% =========================================================================
clear; clc; close all; 

%% ── 1. File paths ────────────────────────────────────────────────────────
mean_field = 'G:\Y235_AOAN04_AOAFN06_SmallerWindows\merge_instantaneousavg_20260429_083207\merged_meanUV_14loops_20260429_083207.mat';
rms_field  = 'G:\Y235_AOAN04_AOAFN06_SmallerWindows\merge_instantaneousavg_20260429_083207\merged_turbrms_hann_14loops_20260429_083207.mat';
savePath = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\16x96init\'; 
%% ── 1b. Smoothing parameters ─────────────────────────────────────────────
sg_order  = 3;    % Savitzky-Golay polynomial order
sg_window = 101;   % must be odd and > sg_order
x_search_start = 600;   % mm — don't look for crossover before this

%% ── 1c. Load data ────────────────────────────────────────────────────────
fprintf('[%s] Loading data...\n', datestr(now,'HH:MM:SS'));
tmp     = load(mean_field, 'U_hann_mean', 'worldX', 'worldY');
rms_tmp = load(rms_field,  'U_rms');
U_mean  = double(tmp.U_hann_mean);
U_rms   = double(rms_tmp.U_rms);
worldX  = double(tmp.worldX);
worldY  = double(tmp.worldY);
clear tmp rms_tmp

x_vec = worldX(1, :);      % [1 x Nx]  mm
y_vec = worldY(:, 1);      % [Ny x 1]  mm
[~, Nx] = size(U_mean);

fprintf('  U_rms NaN coverage: %.1f%%\n', 100*mean(isnan(U_rms(:))));

%% ── 2. Parameters ────────────────────────────────────────────────────────
mag_factor = (max(x_vec) - min(x_vec)) / Nx;
x_targets  = x_vec;
N_x        = length(x_targets);

%% ── 3. Pre-allocate outputs ──────────────────────────────────────────────
out_x                = nan(1, N_x);
out_Uinf             = nan(1, N_x);
out_U99_val          = nan(1, N_x);
out_delta99_hybrid   = nan(1, N_x);
out_delta99_vel      = nan(1, N_x);
out_delta99_vinuesa  = nan(1, N_x);
out_deltastar        = nan(1, N_x);
out_theta            = nan(1, N_x);
out_H                = nan(1, N_x);
out_vindiag          = nan(1, N_x);
out_delta99_rms      = nan(1, N_x);
out_deltastar_rms    = nan(1, N_x);
out_theta_rms        = nan(1, N_x);
out_H_rms            = nan(1, N_x);
out_Uinf_rms         = nan(1, N_x);

% pre-allocate
out_delta99_vinuesa_min = nan(1, N_x);
col_d99_vinuesa_min     = nan(1, N_x);

fprintf('[%s] Starting sweep over %d x-locations...\n', datestr(now,'HH:MM:SS'), N_x);
report_every = max(1, round(N_x / 20));

%% ── 4. Main loop ─────────────────────────────────────────────────────────
for i = 1:N_x

    x_target = x_targets(i);
    cols = find(x_vec >= x_target - mag_factor & ...
                x_vec <= x_target + mag_factor);
    if isempty(cols), continue; end

    out_x(i) = mean(x_vec(cols));
    n_cols    = length(cols);

    col_d99_hybrid  = nan(1, n_cols);
    col_d99_vel     = nan(1, n_cols);
    col_d99_vinuesa = nan(1, n_cols);
    col_ds          = nan(1, n_cols);
    col_th          = nan(1, n_cols);
    col_Uinf        = nan(1, n_cols);
    col_diagmin     = nan(1, n_cols);

    for c = 1:n_cols
        rms_col              = U_rms(:, cols(c));
        rms_col(isnan(rms_col)) = 0;
        urms_smooth          = sgolayfilt(rms_col, sg_order, sg_window);

        % Hybrid (Vinuesa + velocity fallback) — primary method
        [d99, ds, th, Ui, ~, dm] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'hybrid');
        col_d99_hybrid(c) = d99;
        col_ds(c)         = ds;
        col_th(c)         = th;
        col_Uinf(c)       = Ui;
        col_diagmin(c)    = dm;
        % vinuesa min method for all 
        [d99vm, ~, ~, ~, ~, ~, ~] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'vinuesa_min');
        col_d99_vinuesa_min(c) = d99vm;
        % Pure velocity-based
        [d99v, ~, ~, ~, ~, ~, ~] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'velocity');
        col_d99_vel(c) = d99v;

        % Pure Vinuesa (NaN where threshold never crossed) 
        [d99vin, ~, ~, ~, ~, ~, ~] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'vinuesa');
        col_d99_vinuesa(c) = d99vin;
    end

    % ── Column averages ───────────────────────────────────────────────────
    out_delta99_hybrid(i)  = mean(col_d99_hybrid(~isnan(col_d99_hybrid)));
    out_delta99_vel(i)     = mean(col_d99_vel(~isnan(col_d99_vel)));
    out_delta99_vinuesa(i) = mean(col_d99_vinuesa(~isnan(col_d99_vinuesa)));
    out_delta99_vinuesa_min(i) = mean(col_d99_vinuesa_min(~isnan(col_d99_vinuesa_min)));

    out_deltastar(i)       = mean(col_ds(~isnan(col_ds)));
    out_theta(i)           = mean(col_th(~isnan(col_th)));
    out_Uinf(i)            = mean(col_Uinf(~isnan(col_Uinf)));
    out_H(i)               = out_deltastar(i) / out_theta(i);
    out_U99_val(i)         =       0.99 * out_Uinf(i);
    out_vindiag(i)         = mean(col_diagmin(~isnan(col_diagmin)));

    % ── RMS uncertainty ───────────────────────────────────────────────────
    if n_cols > 1
        out_delta99_rms(i)   = std(col_d99_hybrid(~isnan(col_d99_hybrid)));
        out_deltastar_rms(i) = std(col_ds(~isnan(col_ds)));
        out_theta_rms(i)     = std(col_th(~isnan(col_th)));
        out_Uinf_rms(i)      = std(col_Uinf(~isnan(col_Uinf)));
        if out_deltastar(i) > 0 && out_theta(i) > 0
            out_H_rms(i) = out_H(i) * sqrt( ...
                (out_deltastar_rms(i)/out_deltastar(i))^2 + ...
                (out_theta_rms(i)    /out_theta(i)    )^2);
        end
    end

    if mod(i, report_every) == 0
        fprintf('  [%s] %d / %d (%.0f%%)\n', ...
            datestr(now,'HH:MM:SS'), i, N_x, 100*i/N_x);
    end
end

%% ── 5. Diagnostics ───────────────────────────────────────────────────────
n_valid = sum(~isnan(out_delta99_hybrid));
fprintf('  Valid delta99 points: %d / %d (%.1f%%)\n', n_valid, N_x, 100*n_valid/N_x);
if n_valid > 0
    fprintf('  First valid x: %.2f mm\n', min(out_x(~isnan(out_delta99_hybrid))));
    fprintf('  Last  valid x: %.2f mm\n', max(out_x(~isnan(out_delta99_hybrid))));
end
n_crossed  = sum(out_vindiag <= 0.02, 'omitnan');
n_fallback = sum(out_vindiag >  0.02, 'omitnan');
fprintf('  Threshold crossed (<=0.02): %d columns\n', n_crossed);
fprintf('  Fallback to velocity:       %d columns\n', n_fallback);

%% ── 5b. Smooth δ* and θ, compute H ──────────────────────────────────────
valid_hybrid  = ~isnan(out_x) & ~isnan(out_delta99_hybrid);
valid_vel     = ~isnan(out_x) & ~isnan(out_delta99_vel);
valid_vinuesa = ~isnan(out_x) & ~isnan(out_delta99_vinuesa);
valid_diag    = ~isnan(out_x) & ~isnan(out_vindiag);

ds_smooth            = nan(1, N_x);
th_smooth            = nan(1, N_x);
H_smooth             = nan(1, N_x);
ds_smooth(valid_hybrid) = sgolayfilt(out_deltastar(valid_hybrid), sg_order, sg_window);
th_smooth(valid_hybrid) = sgolayfilt(out_theta(valid_hybrid),     sg_order, sg_window);
H_smooth(valid_hybrid)  = ds_smooth(valid_hybrid) ./ th_smooth(valid_hybrid);
%% ── 5c. Find crossover point and update hybrid ───────────────────────────

% Find first x after x_search_start where d99_vinuesa and d99_vel cross
% i.e. Vinuesa transitions from below to above velocity estimate
search_mask  = out_x >= x_search_start & ...
               ~isnan(out_delta99_vinuesa) & ~isnan(out_delta99_vel);

diff_vin_vel = out_delta99_vinuesa - out_delta99_vel;  % +ve means Vin > vel

% Find first sign change from negative to positive (Vinuesa crosses above velocity)
sign_changes = find(search_mask & ...
    [false, diff_vin_vel(1:end-1) <= 0] & ...
    [false, diff_vin_vel(2:end)   >  0]);

if ~isempty(sign_changes)
    x_crossover = out_x(sign_changes(1));
    fprintf('  Vinuesa/velocity crossover at x = %.1f mm\n', x_crossover);
else
    % Fallback — use x_search_start if no crossover found
    x_crossover = x_search_start;
    fprintf('  No crossover found — switching at x = %.1f mm\n', x_crossover);
end

% Update hybrid: use velocity method from crossover point onward
switch_mask = out_x >= x_crossover & ~isnan(out_delta99_vel);
out_delta99_hybrid(switch_mask) = out_delta99_vel(switch_mask);

fprintf('  Hybrid switched to velocity for %d / %d columns\n', ...
    sum(switch_mask), N_x);

% Save metadata
blSweep.crossover.x_mm         = x_crossover;
blSweep.crossover.x_search_mm  = x_search_start;
blSweep.crossover.method       = 'velocity from crossover, Vinuesa upstream';
% %% ── Profile Check: all three methods at selected x stations ──────────────
% x_check = [800, 1000];   % <-- set your x locations in mm
% 
% for xi = 1:length(x_check)
%     [~, col] = min(abs(x_vec - x_check(xi)));
%     x_actual = x_vec(col);
% 
%     % ── Pull column data ──────────────────────────────────────────────────
%     U_col   = U_mean(:, col);
%     rms_col = U_rms(:, col);
%     rms_col(isnan(rms_col)) = 0;
%     urms_col = sgolayfilt(rms_col, sg_order, sg_window);
% 
%     % ── Run all four methods ──────────────────────────────────────────────
%     [d99_hyb,  ~, ~, Ue_hyb,  ~, dm_hyb,  ~] = compute_bl_params(U_col, urms_col, y_vec, 'hybrid');
%     [d99_vel,  ~, ~, Ue_vel,  ~, ~,        ~] = compute_bl_params(U_col, urms_col, y_vec, 'velocity');
%     [d99_vin,  ~, ~, Ue_vin,  ~, ~,        ~] = compute_bl_params(U_col, urms_col, y_vec, 'vinuesa');
%     [d99_vinm, ~, ~, Ue_vinm, ~, ~,        ~] = compute_bl_params(U_col, urms_col, y_vec, 'vinuesa_min');
% 
%     % ── Replicate internal profile prep ──────────────────────────────────
%     valid_above = (y_vec > 0) & (U_col ~= 0) & ~isnan(U_col) & ~isnan(urms_col);
%     if ~any(valid_above), continue; end
%     first_idx = find(valid_above, 1, 'first');
%     y_offset  = y_vec(first_idx);
% 
%     y_fit    = y_vec(valid_above) - y_offset;
%     U_fit    = U_col(valid_above);
%     urms_fit = urms_col(valid_above);
%     [y_fit, sidx] = sort(y_fit, 'ascend');
%     U_fit    = U_fit(sidx);
%     urms_fit = urms_fit(sidx);
%     y_fit    = reshape(y_fit,    [], 1);
%     U_fit    = reshape(U_fit,    [], 1);
%     urms_fit = reshape(urms_fit, [], 1);
% 
%     % Vinuesa diagnostic over full window using hybrid Ue and H
%     smooth_pts      = max(5, round(10 / mean(diff(y_fit))));
%     urms_fit_smooth = smoothdata(urms_fit, 'gaussian', smooth_pts);
%     H_ref           = 1.4;   % use nominal H for diagnostic display
%     vin_diag_full   = urms_fit_smooth ./ max(U_fit, 1e-6) / sqrt(H_ref);
% 
%     % Smoothed U profile for velocity method visualisation
%     U_fit_smooth = smoothdata(U_fit, 'gaussian', smooth_pts);
%     Ue_ref       = max(U_fit_smooth, [], 'omitnan');
% 
%     % ── Figure ────────────────────────────────────────────────────────────
%     figure('Position', [100, 100, 1400, 450], ...
%         'Name', sprintf('Profile check  x = %.1f mm', x_actual));
% 
%     tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
%     title(tl, sprintf(['x = %.1f mm  |  hybrid: %.2f mm  |  ' ...
%         'vel: %.2f mm  |  Vin (thresh): %.2f mm  |  Vin (min): %.2f mm'], ...
%         x_actual, d99_hyb, d99_vel, d99_vin, d99_vinm), 'FontSize', 10);
% 
% 
%     % ── Panel 1: velocity profile ─────────────────────────────────────────
%     nexttile;
%     hold on;
%     plot(y_fit,        U_fit ./ Ue_ref,       'b-',  'LineWidth', 1.5, 'DisplayName', 'U/U_e raw');
%     plot(y_fit,        U_fit_smooth ./ Ue_ref, 'b--', 'LineWidth', 1,   'DisplayName', 'U/U_e smoothed');
%     yline(      0.99, 'k:',  'LineWidth', 1, 'DisplayName', '      0.99');
%     yline(1.00, 'k--', 'LineWidth', 1, 'DisplayName', '1.00');
% 
%     % Mark where max(U_smooth) actually sits in y
%     [~, idx_Ue] = max(U_fit_smooth);
%     xline(y_fit(idx_Ue), 'g-', 'LineWidth', 2, 'DisplayName', 'y at max(U) = U_e');
% 
%     % Mark the actual U_e level on the y axis
%     yline(Ue_ref / Ue_ref, 'g--', 'LineWidth', 1, 'HandleVisibility', 'off');  % always 1.0 by definition
% 
%     if ~isnan(d99_hyb),  xline(d99_hyb,  'b-',  'LineWidth', 1.5, 'DisplayName', '\delta_{99} hybrid');     end
%     if ~isnan(d99_vel),  xline(d99_vel,  'r--', 'LineWidth', 1.5, 'DisplayName', '\delta_{99} velocity');   end
%     if ~isnan(d99_vin),  xline(d99_vin,  'c:',  'LineWidth', 1.5, 'DisplayName', '\delta_{99} Vin thresh'); end
%     if ~isnan(d99_vinm), xline(d99_vinm, 'm-.', 'LineWidth', 1.5, 'DisplayName', '\delta_{99} Vin min');    end
% 
%     xlabel('y – y_{wall} [mm]'); ylabel('U / U_e  [–]');
%     title(sprintf('Velocity profile  |  U_e = %.3f m/s  |  y(U_e) = %.1f mm', ...
%         Ue_ref, y_fit(idx_Ue)));
%     ylim([0 1.15]); grid on;
%     legend('Location', 'southeast', 'FontSize', 7);
%     % ── Panel 2: Vinuesa diagnostic ───────────────────────────────────────
%     nexttile;
%     hold on;
%     plot(y_fit, vin_diag_full, 'k-', 'LineWidth', 1.5);
%     yline(0.02, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold 0.02');
%     if ~isnan(d99_hyb),  xline(d99_hyb,  'b-',  'LineWidth', 1.5, 'DisplayName', '\delta_{99} hybrid');    end
%     if ~isnan(d99_vel),  xline(d99_vel,  'r--', 'LineWidth', 1.5, 'DisplayName', '\delta_{99} velocity');  end
%     if ~isnan(d99_vin),  xline(d99_vin,  'c:',  'LineWidth', 1.5, 'DisplayName', '\delta_{99} Vin thresh'); end
%     if ~isnan(d99_vinm), xline(d99_vinm, 'm-.', 'LineWidth', 1.5, 'DisplayName', '\delta_{99} Vin min');   end
%     xlabel('y – y_{wall} [mm]'); ylabel('u_{rms} / (U\surdH)  [–]');
%     title('Vinuesa diagnostic');
%     yl = ylim; ylim([0, min(yl(2), 0.5)]);
%     grid on; legend('Location', 'best', 'FontSize', 7);
% 
%     % ── Panel 3: RMS profile ──────────────────────────────────────────────
%     nexttile;
%     hold on;
%     plot(y_fit, urms_fit,        'b-',  'LineWidth', 1,   'DisplayName', 'u_{rms} sgolay');
%     plot(y_fit, urms_fit_smooth, 'r-',  'LineWidth', 1.5, 'DisplayName', 'u_{rms} gaussian');
%     if ~isnan(d99_hyb),  xline(d99_hyb,  'b-',  'LineWidth', 1.5, 'DisplayName', '\delta_{99} hybrid');    end
%     if ~isnan(d99_vel),  xline(d99_vel,  'r--', 'LineWidth', 1.5, 'DisplayName', '\delta_{99} velocity');  end
%     if ~isnan(d99_vin),  xline(d99_vin,  'c:',  'LineWidth', 1.5, 'DisplayName', '\delta_{99} Vin thresh'); end
%     if ~isnan(d99_vinm), xline(d99_vinm, 'm-.', 'LineWidth', 1.5, 'DisplayName', '\delta_{99} Vin min');   end
%     xlabel('y – y_{wall} [mm]'); ylabel('u_{rms}  [m/s]');
%     title('RMS profile');
%     grid on; legend('Location', 'best', 'FontSize', 7);
% end

%% ── 6. Package and save ──────────────────────────────────────────────────
blSweep.x_mm               = out_x;
blSweep.U_inf              = out_Uinf;
blSweep.U99                = out_U99_val;
blSweep.delta99_hybrid_mm  = out_delta99_hybrid;
blSweep.delta99_vel_mm     = out_delta99_vel;
blSweep.delta99_vinuesa_mm = out_delta99_vinuesa;
blSweep.deltastar_mm       = ds_smooth;
blSweep.theta_mm           = th_smooth;
blSweep.H                  = H_smooth;
blSweep.vinuesadiag        = out_vindiag;
blSweep.sg_order           = sg_order;
blSweep.sg_window          = sg_window;

blSweep.rms.U_inf          = out_Uinf_rms;
blSweep.rms.delta99_mm     = out_delta99_rms;
blSweep.rms.deltastar_mm   = out_deltastar_rms;
blSweep.rms.theta_mm       = out_theta_rms;
blSweep.rms.H              = out_H_rms;
blSweep.delta99_vinuesa_min_mm = out_delta99_vinuesa_min;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_file  = sprintf('blSweep_%s.mat', timestamp);
save(fullfile(savePath, out_file), 'blSweep', '-v7.3');
fprintf('[%s] Saved → %s\n', datestr(now,'HH:MM:SS'), out_file);

%% ── 7. Plots ─────────────────────────────────────────────────────────────
fprintf('[%s] Plotting...\n', datestr(now,'HH:MM:SS'));

% ── Figure 1: velocity field + all three delta99 lines ───────────────────
figure('Position', [100, 100, 1400, 500]);
imagesc(worldX(1,:), worldY(:,1), U_mean);
axis xy equal tight;
cb = colorbar; cb.Label.String = 'U [m/s]';
colormap(jet);
clim([0, max(U_mean(:), [], 'omitnan') * 1.05]);
hold on;
plot(out_x(valid_hybrid), sgolayfilt(out_delta99_hybrid(valid_hybrid),  sg_order, sg_window), ...
    'w-',  'LineWidth', 2, 'DisplayName', '\delta_{99} hybrid');
plot(out_x(valid_vel),    sgolayfilt(out_delta99_vel(valid_vel),        sg_order, sg_window), ...
    'y--', 'LineWidth', 2, 'DisplayName', '\delta_{99} velocity');
plot(out_x(valid_vinuesa), sgolayfilt(out_delta99_vinuesa_min(valid_vinuesa), sg_order, sg_window), ...
    'c:',  'LineWidth', 2, 'DisplayName', '\delta_{99} Vinuesa');



x_fill = [out_x(valid_hybrid), fliplr(out_x(valid_hybrid))];
y_fill = [out_delta99_hybrid(valid_hybrid) + out_delta99_rms(valid_hybrid), ...
          fliplr(out_delta99_hybrid(valid_hybrid) - out_delta99_rms(valid_hybrid))];
fill(x_fill, y_fill, 'w', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
yline(0, 'k-', 'LineWidth', 2, 'DisplayName', 'Wall');
xlabel('X [mm]'); ylabel('Y [mm]');
title('Mean U-Velocity Field — \delta_{99} method comparison');
legend('Location', 'northwest'); grid off;

% ── Figure 2: U_inf ───────────────────────────────────────────────────────
figure('Position', [100, 100, 1400, 500]);
plot(out_x(valid_hybrid), sgolayfilt(out_Uinf(valid_hybrid), sg_order, sg_window), ...
    'k-', 'LineWidth', 2);
xlabel('X [mm]'); ylabel('U_\infty [m/s]');
title('Freestream Velocity Development');
grid off;

% ── Figure 3: delta99 + H ─────────────────────────────────────────────────
figure('Position', [100, 650, 1400, 400]);
yyaxis left
plot(out_x(valid_hybrid), sgolayfilt(out_delta99_hybrid(valid_hybrid), sg_order, sg_window), ...
    'b-', 'LineWidth', 2);
ylabel('\delta_{99} [mm]');
yyaxis right
plot(out_x(valid_hybrid), H_smooth(valid_hybrid), 'r-', 'LineWidth', 2);
ylabel('H = \delta^* / \theta  [–]');
yline(1.4, 'r:', 'LineWidth', 1, 'DisplayName', 'H = 1.4 (ZPG ref)');
xlabel('X [mm]');
title('Boundary Layer Development');
legend('\delta_{99}', 'H', 'ZPG ref', 'Location', 'best');
grid on; hold off;

% ── Figure 4: δ*, θ, H ───────────────────────────────────────────────────
figure('Position', [100, 650, 1400, 400]);
yyaxis left
hold on;
plot(out_x(valid_hybrid), ds_smooth(valid_hybrid), 'b-', 'LineWidth', 2, 'DisplayName', '\delta^*');
plot(out_x(valid_hybrid), th_smooth(valid_hybrid), 'k-', 'LineWidth', 2, 'DisplayName', '\theta');
ylabel('\delta^*, \theta [mm]');
yyaxis right
plot(out_x(valid_hybrid), H_smooth(valid_hybrid),  'r-', 'LineWidth', 2, 'DisplayName', 'H');
ylabel('H = \delta^* / \theta  [–]');
yline(1.4, 'r:', 'LineWidth', 1, 'DisplayName', 'H = 1.4 (ZPG ref)');
xlabel('X [mm]');
title('Integral Parameters');
legend('Location', 'best'); grid on; hold off;

% ── Figure 5: delta99 method comparison ──────────────────────────────────
valid_vinuesa_min = ~isnan(out_x) & ~isnan(out_delta99_vinuesa_min);

figure('Position', [100, 550, 1400, 350]);
hold on;
plot(out_x(valid_hybrid),      sgolayfilt(out_delta99_hybrid(valid_hybrid),           sg_order, sg_window), ...
    'b-',  'LineWidth', 2, 'DisplayName', '\delta_{99} hybrid');
plot(out_x(valid_vel),         sgolayfilt(out_delta99_vel(valid_vel),                 sg_order, sg_window), ...
    'r--', 'LineWidth', 2, 'DisplayName', '\delta_{99} velocity');
plot(out_x(valid_vinuesa),     sgolayfilt(out_delta99_vinuesa(valid_vinuesa),         sg_order, sg_window), ...
    'c:',  'LineWidth', 2, 'DisplayName', '\delta_{99} Vinuesa (threshold only)');
plot(out_x(valid_vinuesa_min), sgolayfilt(out_delta99_vinuesa_min(valid_vinuesa_min), sg_order, sg_window), ...
    'm-.', 'LineWidth', 2, 'DisplayName', '\delta_{99} Vinuesa (min fallback)');
xlabel('X [mm]'); ylabel('\delta_{99} [mm]');
title('\delta_{99} Method Comparison');
legend('Location', 'best'); grid on; hold off;

% ── Figure 6: Vinuesa diagnostic vs x ────────────────────────────────────
figure();
crossed  = valid_diag & (out_vindiag <= 0.02);
fallback = valid_diag & (out_vindiag >  0.02);
ax1 = subplot(2,1,1);
scatter(ax1, out_x(crossed),  out_vindiag(crossed),  4, [0.2 0.6 0.2], 'filled', ...
    'DisplayName', 'Threshold crossed (\leq 0.02)');
hold on;
scatter(ax1, out_x(fallback), out_vindiag(fallback), 4, [0.8 0.2 0.2], 'filled', ...
    'DisplayName', 'Velocity fallback (> 0.02)');
yline(ax1, 0.02, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Vinuesa threshold');
xlabel(ax1, 'X [mm]');
ylabel(ax1, '$(u_{rms} / (U \sqrt{H}))_{\delta_{99}}$', 'Interpreter', 'latex');
title(ax1, 'Vinuesa Diagnostic Value at Detected \delta_{99}');
xlim(ax1, [min(worldX(1,:)) max(worldX(1,:))]);
legend(ax1, 'Location', 'best'); grid(ax1, 'on'); hold off;
ax2 = subplot(2,1,2);
imagesc(ax2, worldX(1,:), worldY(:,1), U_rms);
set(ax2, 'YDir', 'normal');
clim(ax2, [0 2]); colorbar();
axis(ax2, 'image');
xlim(ax2, [min(worldX(1,:)) max(worldX(1,:))]);
linkaxes([ax1, ax2], 'x');

% ── Figure 7: histogram of diagnostic values ─────────────────────────────
figure();
histogram(out_vindiag(valid_diag), 60, 'FaceColor', [0.3 0.5 0.8]);
xline(0.02, 'r--', 'LineWidth', 2, 'DisplayName', '0.02 threshold');
xlabel('$(u_{rms} / (U \sqrt{H}))_{\delta_{99}}$', 'Interpreter', 'latex');
ylabel('Count');
title('Distribution of Vinuesa Diagnostic at Detected Edge');
legend('Location', 'best'); grid on;

%% =========================================================================
function [delta99, deltastar, theta, U_e, U99, diag_min, used_fallback] = compute_bl_params(U_col, urms_col, y_col, method)
% method: 'hybrid' (default), 'vinuesa', 'velocity', 'vinuesa_min'
if nargin < 4, method = 'hybrid'; end

delta99       = NaN;  deltastar = NaN;
theta         = NaN;  U_e       = NaN;
U99           = NaN;  diag_min  = NaN;
used_fallback = false;

y_max_search   = 180;
diag_threshold = 0.02;
max_iter       = 10;
tol            = 1e-4;

% ── Wall offset ───────────────────────────────────────────────────────────
valid_above = (y_col > 0) & (U_col ~= 0) & ~isnan(U_col) & ~isnan(urms_col);
if ~any(valid_above), return; end

first_idx = find(valid_above, 1, 'first');
y_offset  = y_col(first_idx);

y_fit    = y_col(valid_above) - y_offset;
U_fit    = U_col(valid_above);
urms_fit = urms_col(valid_above);

[y_fit, sidx] = sort(y_fit, 'ascend');
U_fit    = U_fit(sidx);
urms_fit = urms_fit(sidx);
y_fit    = reshape(y_fit,    [], 1);
U_fit    = reshape(U_fit,    [], 1);
urms_fit = reshape(urms_fit, [], 1);

if length(y_fit) < 5, return; end

search_mask = y_fit <= y_max_search;
y_s    = reshape(y_fit(search_mask),    [], 1);
U_s    = reshape(U_fit(search_mask),    [], 1);
urms_s = reshape(urms_fit(search_mask), [], 1);

if length(y_s) < 5, return; end

smooth_pts  = max(5, round(10 / mean(diff(y_s))));
urms_smooth = smoothdata(urms_s, 'gaussian', smooth_pts);

H12_iter      = 1.4;
delta99_iter  = NaN;
U_e_iter      = NaN;
diag_min_iter = NaN;

% ── Iterative loop ────────────────────────────────────────────────────────
for iter = 1:max_iter

    safe_U       = max(U_s, 1e-6);
    vinuesa_diag = urms_smooth ./ (safe_U * sqrt(H12_iter));

    crosses = find(vinuesa_diag <= diag_threshold, 1, 'first');
    use_vin = ~isempty(crosses) && crosses >= 2;

    % ── Detection switch ──────────────────────────────────────────────────
    if strcmp(method, 'velocity') || (strcmp(method, 'hybrid') && ~use_vin)
        U_s_smooth = smoothdata(U_s, 'gaussian', smooth_pts);
        Ue_vel     = max(U_s_smooth, [], 'omitnan');
        cross_vel  = find(U_s_smooth >=       0.99 * Ue_vel, 1, 'first');

        if ~isempty(cross_vel) && cross_vel >= 2
            u1 = U_s_smooth(cross_vel-1);  y1 = y_s(cross_vel-1);
            u2 = U_s_smooth(cross_vel);    y2 = y_s(cross_vel);
            if u2 ~= u1
                delta99_new = interp1([u1,u2], [y1,y2],       0.99*Ue_vel);
            else
                delta99_new = y2;
            end
            U_e_new = Ue_vel;
        elseif ~isempty(cross_vel)
            delta99_new = y_s(cross_vel);
            U_e_new     = Ue_vel;
        else
            break
        end
        [diag_min_new, ~] = min(vinuesa_diag);
        used_fallback     = true;

    elseif strcmp(method, 'vinuesa') || (strcmp(method, 'hybrid') && use_vin)
        if strcmp(method, 'vinuesa') && ~use_vin
            break
        end
        d1 = vinuesa_diag(crosses-1);  d2 = vinuesa_diag(crosses);
        y1 = y_s(crosses-1);           y2 = y_s(crosses);
        if d1 ~= d2
            delta99_new = interp1([d1,d2], [y1,y2], diag_threshold);
            delta99_new = max(min(y1,y2), min(max(y1,y2), delta99_new));
            U_e_new     = interp1([y1,y2], [U_s(crosses-1), U_s(crosses)], delta99_new);
        else
            delta99_new = y2;
            U_e_new     = U_s(crosses);
        end
        diag_min_new  = diag_threshold;
        used_fallback = false;
    elseif strcmp(method, 'vinuesa_min') || (strcmp(method, 'hybrid') && use_vin)
        if strcmp(method, 'vinuesa_min') && ~use_vin
            % No crossing — take minimum as best estimate
            [diag_min_new, idx_min] = min(vinuesa_diag);
            if idx_min < 2, break; end
            delta99_new = y_s(idx_min);
            U_e_new     = U_s(idx_min);
        else
            % Genuine crossing
            d1 = vinuesa_diag(crosses-1);  d2 = vinuesa_diag(crosses);
            y1 = y_s(crosses-1);           y2 = y_s(crosses);
            if d1 ~= d2
                delta99_new = interp1([d1,d2], [y1,y2], diag_threshold);
                delta99_new = max(min(y1,y2), min(max(y1,y2), delta99_new));
                U_e_new     = interp1([y1,y2], [U_s(crosses-1), U_s(crosses)], delta99_new);
            else
                delta99_new = y2;
                U_e_new     = U_s(crosses);
            end
            diag_min_new = diag_threshold;
        end
        used_fallback = false;

    else
        break
    end

    % ── Integral thicknesses ──────────────────────────────────────────────
    in_bl = y_fit < delta99_new;
    if sum(in_bl) < 2, break; end

    y_int  = [y_fit(in_bl);  delta99_new];
    u_norm = [U_fit(in_bl);  U_e_new] / U_e_new;

    ds_new  = trapz(y_int, 1 - u_norm);
    th_new  = trapz(y_int, u_norm .* (1 - u_norm));
    H12_new = (th_new > 0) * (ds_new / th_new) + (th_new <= 0) * H12_iter;

    if ~isnan(delta99_iter) && abs(delta99_new - delta99_iter) < tol
        delta99_iter  = delta99_new;
        H12_iter      = H12_new;
        U_e_iter      = U_e_new;
        diag_min_iter = diag_min_new;
        break
    end

    delta99_iter  = delta99_new;
    H12_iter      = H12_new;
    U_e_iter      = U_e_new;
    diag_min_iter = diag_min_new;

end  % for iter

% ── Package outputs ───────────────────────────────────────────────────────
if ~isnan(delta99_iter) && ~isnan(U_e_iter)
    delta99  = delta99_iter;
    U_e      = U_e_iter;
    U99      =       0.99 * U_e;
    diag_min = diag_min_iter;

    in_bl = y_fit < delta99;
    if sum(in_bl) < 2, return; end
    y_int     = [y_fit(in_bl);  delta99];
    u_norm    = [U_fit(in_bl);  U_e] / U_e;
    deltastar = trapz(y_int, 1 - u_norm);
    theta     = trapz(y_int, u_norm .* (1 - u_norm));
end

end  % function
% =========================================================================
% bl_sweep_Tom.m
% Adapted from bl_sweep.m (Ashley Kwong, March 2026) for use with
% merged multi-position PIV fields produced by merge_positions_mean.m
%
% Input file contains: worldX, worldY, U_hann_mean, U_rms
% (variable names produced by merge_positions_mean.m)
% =========================================================================
clear; clc; close all;

%% ── 1. File paths ────────────────────────────────────────────────────────
% !! Point this at the file produced by merge_positions_mean.m !!
mean_field = 'G:\SW 500mm Mean Flow Fields\globalmerged_mean_field_20260424_215933.mat';

%% ── 1b. Smoothing parameters ─────────────────────────────────────────────
sg_order  = 3;    % Savitzky-Golay polynomial order
sg_window = 101;  % must be odd and > sg_order
% !! Review x_search_start once you know the x-extent of the merged field !!
% It sets the earliest x [mm] at which the Vinuesa/velocity crossover is
% searched for. Too small → spurious early crossover; too large → misses it.
x_search_start = 10000;   % mm — TUNE THIS after first run

%% ── 1c. Load data ────────────────────────────────────────────────────────
fprintf('[%s] Loading data...\n', datestr(now,'HH:MM:SS'));
tmp    = load(mean_field, 'U_hann_mean', 'U_rms', 'worldX', 'worldY');
U_mean = double(tmp.U_hann_mean);
U_rms  = double(tmp.U_rms);
worldX = double(tmp.worldX);
worldY = double(tmp.worldY);
clear tmp

x_vec = worldX(1, :);   % [1 x Nx]  mm
y_vec = worldY(:, 1);   % [Ny x 1]  mm
[~, Nx] = size(U_mean);

fprintf('  Grid: %d x %d   X:[%.1f, %.1f] mm   Y:[%.1f, %.1f] mm\n', ...
    size(U_mean,2), size(U_mean,1), min(x_vec), max(x_vec), min(y_vec), max(y_vec));
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

        [d99, ds, th, Ui, ~, dm] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'hybrid');
        col_d99_hybrid(c) = d99;
        col_ds(c)         = ds;
        col_th(c)         = th;
        col_Uinf(c)       = Ui;
        col_diagmin(c)    = dm;

        [d99vm, ~, ~, ~, ~, ~, ~] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'vinuesa_min');
        col_d99_vinuesa_min(c) = d99vm;

        [d99v, ~, ~, ~, ~, ~, ~] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'velocity');
        col_d99_vel(c) = d99v;

        [d99vin, ~, ~, ~, ~, ~, ~] = compute_bl_params( ...
            U_mean(:,cols(c)), urms_smooth, y_vec, 'vinuesa');
        col_d99_vinuesa(c) = d99vin;
    end

    out_delta99_hybrid(i)      = mean(col_d99_hybrid(~isnan(col_d99_hybrid)));
    out_delta99_vel(i)         = mean(col_d99_vel(~isnan(col_d99_vel)));
    out_delta99_vinuesa(i)     = mean(col_d99_vinuesa(~isnan(col_d99_vinuesa)));
    out_delta99_vinuesa_min(i) = mean(col_d99_vinuesa_min(~isnan(col_d99_vinuesa_min)));
    out_deltastar(i)           = mean(col_ds(~isnan(col_ds)));
    out_theta(i)               = mean(col_th(~isnan(col_th)));
    out_Uinf(i)                = mean(col_Uinf(~isnan(col_Uinf)));
    out_H(i)                   = out_deltastar(i) / out_theta(i);
    out_U99_val(i)             = 0.99 * out_Uinf(i);
    out_vindiag(i)             = mean(col_diagmin(~isnan(col_diagmin)));

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

%% ── 5b. Smooth δ* and θ, compute H ──────────────────────────────────────
valid_hybrid  = ~isnan(out_x) & ~isnan(out_delta99_hybrid);
valid_vel     = ~isnan(out_x) & ~isnan(out_delta99_vel);
valid_vinuesa = ~isnan(out_x) & ~isnan(out_delta99_vinuesa);
valid_diag    = ~isnan(out_x) & ~isnan(out_vindiag);

ds_smooth           = nan(1, N_x);
th_smooth           = nan(1, N_x);
H_smooth            = nan(1, N_x);
ds_smooth(valid_hybrid) = sgolayfilt(out_deltastar(valid_hybrid), sg_order, sg_window);
th_smooth(valid_hybrid) = sgolayfilt(out_theta(valid_hybrid),     sg_order, sg_window);
H_smooth(valid_hybrid)  = ds_smooth(valid_hybrid) ./ th_smooth(valid_hybrid);

%% ── 5c. Find crossover point and update hybrid ───────────────────────────
search_mask  = out_x >= x_search_start & ...
               ~isnan(out_delta99_vinuesa) & ~isnan(out_delta99_vel);
diff_vin_vel = out_delta99_vinuesa - out_delta99_vel;
sign_changes = find(search_mask & ...
    [false, diff_vin_vel(1:end-1) <= 0] & ...
    [false, diff_vin_vel(2:end)   >  0]);

if ~isempty(sign_changes)
    x_crossover = out_x(sign_changes(1));
    fprintf('  Vinuesa/velocity crossover at x = %.1f mm\n', x_crossover);
else
    x_crossover = x_search_start;
    fprintf('  No crossover found — switching at x = %.1f mm\n', x_crossover);
end

switch_mask = out_x >= x_crossover & ~isnan(out_delta99_vel);
out_delta99_hybrid(switch_mask) = out_delta99_vel(switch_mask);
fprintf('  Hybrid switched to velocity for %d / %d columns\n', sum(switch_mask), N_x);

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
blSweep.crossover.x_mm        = x_crossover;
blSweep.crossover.x_search_mm = x_search_start;
blSweep.crossover.method       = 'velocity from crossover, Vinuesa upstream';

out_file = sprintf('blSweep_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(out_file, 'blSweep', '-v7.3');
fprintf('[%s] Saved → %s\n', datestr(now,'HH:MM:SS'), out_file);

%% ── 7. Plots ─────────────────────────────────────────────────────────────
valid_vinuesa_min = ~isnan(out_x) & ~isnan(out_delta99_vinuesa_min);

figure('Position', [100, 100, 1400, 500]);
imagesc(worldX(1,:), worldY(:,1), U_mean);
axis xy equal tight; cb = colorbar; cb.Label.String = 'U [m/s]';
colormap(jet); clim([0, max(U_mean(:), [], 'omitnan') * 1.05]);
hold on;
plot(out_x(valid_hybrid),  sgolayfilt(out_delta99_hybrid(valid_hybrid),  sg_order, sg_window), ...
    'w-',  'LineWidth', 2, 'DisplayName', '\delta_{99} hybrid');
plot(out_x(valid_vel),     sgolayfilt(out_delta99_vel(valid_vel),        sg_order, sg_window), ...
    'y--', 'LineWidth', 2, 'DisplayName', '\delta_{99} velocity');
plot(out_x(valid_vinuesa), sgolayfilt(out_delta99_vinuesa_min(valid_vinuesa), sg_order, sg_window), ...
    'c:',  'LineWidth', 2, 'DisplayName', '\delta_{99} Vinuesa');
yline(0, 'k-', 'LineWidth', 2, 'DisplayName', 'Wall');
xlabel('X [mm]'); ylabel('Y [mm]');
title('Mean U-Velocity — \delta_{99} comparison');
legend('Location', 'northwest');

figure('Position', [100, 650, 1400, 400]);
yyaxis left
plot(out_x(valid_hybrid), sgolayfilt(out_delta99_hybrid(valid_hybrid), sg_order, sg_window), ...
    'b-', 'LineWidth', 2);
ylabel('\delta_{99} [mm]');
yyaxis right
plot(out_x(valid_hybrid), H_smooth(valid_hybrid), 'r-', 'LineWidth', 2);
ylabel('H = \delta^* / \theta  [–]');
yline(1.4, 'r:', 'LineWidth', 1);
xlabel('X [mm]'); title('Boundary Layer Development');
legend('\delta_{99}', 'H', 'ZPG ref', 'Location', 'best'); grid on;

%% =========================================================================
function [delta99, deltastar, theta, U_e, U99, diag_min, used_fallback] = compute_bl_params(U_col, urms_col, y_col, method)
if nargin < 4, method = 'hybrid'; end

delta99       = NaN;  deltastar = NaN;
theta         = NaN;  U_e       = NaN;
U99           = NaN;  diag_min  = NaN;
used_fallback = false;

y_max_search   = 180;
diag_threshold = 0.02;
max_iter       = 10;
tol            = 1e-4;

valid_above = (y_col > 0) & (U_col ~= 0) & ~isnan(U_col) & ~isnan(urms_col);
if ~any(valid_above), return; end

valid_idxs = find(valid_above);
[~, min_loc] = min(y_col(valid_idxs));
first_idx = valid_idxs(min_loc);

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

for iter = 1:max_iter
    safe_U       = max(U_s, 1e-6);
    vinuesa_diag = urms_smooth ./ (safe_U * sqrt(H12_iter));
    crosses      = find(vinuesa_diag <= diag_threshold, 1, 'first');
    use_vin      = ~isempty(crosses) && crosses >= 2;

    if strcmp(method, 'velocity') || (strcmp(method, 'hybrid') && ~use_vin)
        U_s_smooth = smoothdata(U_s, 'gaussian', smooth_pts);
        Ue_vel     = max(U_s_smooth, [], 'omitnan');
        cross_vel  = find(U_s_smooth >= 0.99 * Ue_vel, 1, 'first');
        if ~isempty(cross_vel) && cross_vel >= 2
            u1 = U_s_smooth(cross_vel-1); y1 = y_s(cross_vel-1);
            u2 = U_s_smooth(cross_vel);   y2 = y_s(cross_vel);
            if u2 ~= u1
                delta99_new = interp1([u1,u2],[y1,y2], 0.99*Ue_vel);
            else
                delta99_new = y2;
            end
            U_e_new = Ue_vel;
        elseif ~isempty(cross_vel)
            delta99_new = y_s(cross_vel); U_e_new = Ue_vel;
        else
            break
        end
        [diag_min_new, ~] = min(vinuesa_diag);
        used_fallback = true;

    elseif strcmp(method, 'vinuesa') || (strcmp(method, 'hybrid') && use_vin)
        if strcmp(method, 'vinuesa') && ~use_vin, break; end
        d1 = vinuesa_diag(crosses-1); d2 = vinuesa_diag(crosses);
        y1 = y_s(crosses-1);          y2 = y_s(crosses);
        if d1 ~= d2
            delta99_new = interp1([d1,d2],[y1,y2], diag_threshold);
            delta99_new = max(min(y1,y2), min(max(y1,y2), delta99_new));
            U_e_new     = interp1([y1,y2],[U_s(crosses-1),U_s(crosses)], delta99_new);
        else
            delta99_new = y2; U_e_new = U_s(crosses);
        end
        diag_min_new  = diag_threshold;
        used_fallback = false;

    elseif strcmp(method, 'vinuesa_min') || (strcmp(method, 'hybrid') && use_vin)
        if strcmp(method, 'vinuesa_min') && ~use_vin
            [diag_min_new, idx_min] = min(vinuesa_diag);
            if idx_min < 2, break; end
            delta99_new = y_s(idx_min); U_e_new = U_s(idx_min);
        else
            d1 = vinuesa_diag(crosses-1); d2 = vinuesa_diag(crosses);
            y1 = y_s(crosses-1);          y2 = y_s(crosses);
            if d1 ~= d2
                delta99_new = interp1([d1,d2],[y1,y2], diag_threshold);
                delta99_new = max(min(y1,y2), min(max(y1,y2), delta99_new));
                U_e_new     = interp1([y1,y2],[U_s(crosses-1),U_s(crosses)], delta99_new);
            else
                delta99_new = y2; U_e_new = U_s(crosses);
            end
            diag_min_new = diag_threshold;
        end
        used_fallback = false;
    else
        break
    end

    in_bl  = y_fit < delta99_new;
    if sum(in_bl) < 2, break; end
    y_int  = [y_fit(in_bl);  delta99_new];
    u_norm = [U_fit(in_bl);  U_e_new] / U_e_new;
    ds_new = trapz(y_int, 1 - u_norm);
    th_new = trapz(y_int, u_norm .* (1 - u_norm));
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
end

if ~isnan(delta99_iter) && ~isnan(U_e_iter)
    delta99  = delta99_iter;
    U_e      = U_e_iter;
    U99      = 0.99 * U_e;
    diag_min = diag_min_iter;
    in_bl    = y_fit < delta99;
    if sum(in_bl) < 2, return; end
    y_int     = [y_fit(in_bl);  delta99];
    u_norm    = [U_fit(in_bl);  U_e] / U_e;
    deltastar = trapz(y_int, 1 - u_norm);
    theta     = trapz(y_int, u_norm .* (1 - u_norm));
end
end
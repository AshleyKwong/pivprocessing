%% comparePG_PIVvsTaps_SW.m
clear; clc; close all;

%% -------------------------------------------------------------------------
% 0. USER SETTINGS
% -------------------------------------------------------------------------
sg_order    = 3;
c_chord     = 1.25;        % chord length [m]
x_LE        = 6.53;        % leading edge global x [m]
LE_offset_m = x_LE;        % x/delta origin shift to LE

% PIV x offset — set to 0 until known
piv_x_offset_m = 0;        % !! UPDATE when global PIV origin is known !!

% S-G windows
sg_win_piv_override = [];
sg_win_tap_override = 5;

piv_path = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\blSweep_20260424_220155_global.mat';   % !! UPDATE: path to blSweep .mat from bl_sweep_SW.m !!
tap_path = 'G:\Mean_Cp_Data_SW\Cp_SW_h_0.5m.csv';
out_path = 'PIV_vs_taps_results_SW.mat';

%% -------------------------------------------------------------------------
% 1. LOAD
% -------------------------------------------------------------------------
load(piv_path, 'blSweep');

tap_raw = readmatrix(tap_path);   % 21 x 4
xc_tap  = tap_raw(2:end, 2);          % x/c
Cp_tap_raw_all = tap_raw(2:end, 4);   % Cp

%% -------------------------------------------------------------------------
% 2. PIV: valid data + global x
% -------------------------------------------------------------------------
x_piv_mm = blSweep.x_mm;
U99      = blSweep.U99;
 
valid_piv = ~isnan(x_piv_mm) & ~isnan(U99);
x_piv_mm    = x_piv_mm(valid_piv);
U99         = U99(valid_piv);
delta99_mm  = blSweep.delta99_hybrid_mm(valid_piv);
deltastar_mm = blSweep.deltastar_mm(valid_piv);
theta_mm    = blSweep.theta_mm(valid_piv);
H           = blSweep.H(valid_piv);
 
x_piv_m = x_piv_mm/1000 + piv_x_offset_m;   % global x [m]
 
% Uinf from last 40 mm of PIV domain
xf_piv   = max(x_piv_m);
inf_idx  = x_piv_m >= (xf_piv - 0.04);
Uinf_piv = mean(U99(inf_idx), 'omitnan');
fprintf('PIV Uinf (outlet mean): %.4f m/s\n', Uinf_piv);
 
%% -------------------------------------------------------------------------
% 3. PIV: BL thickness at inlet (first 40 mm)
% -------------------------------------------------------------------------
x0_BL   = min(x_piv_m);
BL_idx  = x_piv_m <= (x0_BL + 0.04);
BL_0_mm = mean(delta99_mm(BL_idx), 'omitnan');
BL_0_m  = BL_0_mm / 1000;
fprintf('Reference BL thickness delta0: %.2f mm\n', BL_0_mm);
 
%% -------------------------------------------------------------------------
% 4. PIV: Cp and raw dCp/d(x/delta)
% -------------------------------------------------------------------------
Cp_piv     = 1 - (U99 / Uinf_piv).^2;
xdelta_piv = x_piv_m / BL_0_m;
dCp_piv_raw = gradient(Cp_piv, xdelta_piv);
 
%% -------------------------------------------------------------------------
% 5. TAPS: convert x/c to global x, then to x/delta
% -------------------------------------------------------------------------
x_tap_m = xc_tap * c_chord + x_LE;   % global x [m]
 
[x_tap_sorted, sort_idx] = sort(x_tap_m, 'ascend');
Cp_tap_sorted = Cp_tap_raw_all(sort_idx);
xdelta_tap_sorted = x_tap_sorted / BL_0_m;
 
valid_tap    = ~isnan(Cp_tap_sorted) & ~isnan(xdelta_tap_sorted) & ~isinf(xdelta_tap_sorted);
xdelta_tap_v = xdelta_tap_sorted(valid_tap);
Cp_tap_v     = Cp_tap_sorted(valid_tap);
x_tap_v      = x_tap_sorted(valid_tap);
 
% Deduplicate (average any repeated x/delta locations)
[xdelta_tap_u, ~, ic] = unique(xdelta_tap_v, 'sorted');
Cp_tap_u              = accumarray(ic, Cp_tap_v, [], @mean);
x_tap_u               = accumarray(ic, x_tap_v,  [], @mean);
 
dCp_tap_raw = gradient(Cp_tap_u, xdelta_tap_u);
 
fprintf('Taps loaded: %d valid points, x=[%.3f, %.3f] m\n', ...
    numel(x_tap_u), min(x_tap_u), max(x_tap_u));
 
%% -------------------------------------------------------------------------
% 6. S-G WINDOWS — matched physical width or user override
% -------------------------------------------------------------------------
tap_spacing_mean = mean(diff(xdelta_tap_u));
piv_spacing_mean = mean(diff(xdelta_piv));
 
default_physical_width = (sg_order + 2) * tap_spacing_mean;
 
default_win_piv = round(default_physical_width / piv_spacing_mean);
if mod(default_win_piv, 2) == 0, default_win_piv = default_win_piv + 1; end
default_win_piv = max(default_win_piv, sg_order + 2);
 
default_win_tap = sg_order + 2;
if mod(default_win_tap, 2) == 0, default_win_tap = default_win_tap + 1; end
 
if ~isempty(sg_win_piv_override)
    sg_win_piv = sg_win_piv_override;
    fprintf('PIV S-G window: OVERRIDE = %d points\n', sg_win_piv);
else
    sg_win_piv = default_win_piv;
    fprintf('PIV S-G window: AUTO = %d points (%.2f x/delta)\n', ...
            sg_win_piv, sg_win_piv * piv_spacing_mean);
end
 
if ~isempty(sg_win_tap_override)
    sg_win_tap = sg_win_tap_override;
    fprintf('Tap S-G window: OVERRIDE = %d points\n', sg_win_tap);
else
    sg_win_tap = default_win_tap;
    fprintf('Tap S-G window: AUTO = %d points (%.2f x/delta)\n', ...
            sg_win_tap, sg_win_tap * tap_spacing_mean);
end
 
% Guard windows
sg_win_piv = min(sg_win_piv, 2*floor(numel(dCp_piv_raw)/2)-1);
sg_win_piv = max(sg_win_piv, sg_order + 2);
if mod(sg_win_piv, 2) == 0, sg_win_piv = sg_win_piv + 1; end
 
sg_win_tap = min(sg_win_tap, 2*floor(numel(dCp_tap_raw)/2)-1);
sg_win_tap = max(sg_win_tap, sg_order + 2);
if mod(sg_win_tap, 2) == 0, sg_win_tap = sg_win_tap + 1; end
 
dCp_piv_sg = sgolayfilt(dCp_piv_raw, sg_order, sg_win_piv);
dCp_tap_sg = sgolayfilt(dCp_tap_raw, sg_order, sg_win_tap);
 
fprintf('Physical S-G width — PIV: %.2f x/delta | Taps: %.2f x/delta\n', ...
        sg_win_piv * piv_spacing_mean, sg_win_tap * tap_spacing_mean);
 
%% -------------------------------------------------------------------------
% 7. CLIP TO OVERLAP REGION
% -------------------------------------------------------------------------
xd_lo = max(min(xdelta_piv), min(xdelta_tap_u));
xd_hi = min(max(xdelta_piv), max(xdelta_tap_u));
fprintf('Overlap x/delta region: %.2f to %.2f\n', xd_lo, xd_hi);
 
piv_overlap = xdelta_piv   >= xd_lo & xdelta_piv   <= xd_hi;
tap_overlap = xdelta_tap_u >= xd_lo & xdelta_tap_u <= xd_hi;
 
xdelta_piv_cl  = xdelta_piv(piv_overlap);
Cp_piv_cl      = Cp_piv(piv_overlap);
dCp_piv_sg_cl  = dCp_piv_sg(piv_overlap);
 
xdelta_tap_cl  = xdelta_tap_u(tap_overlap);
Cp_tap_cl      = Cp_tap_u(tap_overlap);
dCp_tap_raw_cl = dCp_tap_raw(tap_overlap);
dCp_tap_sg_cl  = dCp_tap_sg(tap_overlap);
 
%% -------------------------------------------------------------------------
% 7b. SHIFT x/delta ORIGIN TO LE
% -------------------------------------------------------------------------
LE_xdelta = LE_offset_m / BL_0_m;
 
xdelta_piv_cl = xdelta_piv_cl - LE_xdelta;
xdelta_tap_cl = xdelta_tap_cl - LE_xdelta;
 
%% -------------------------------------------------------------------------
% 8. DIAGNOSTICS
% -------------------------------------------------------------------------
fprintf('\n=== DIAGNOSTICS ===\n');
fprintf('Cp_tap range (overlap):  %.4f to %.4f\n', min(Cp_tap_cl),  max(Cp_tap_cl));
fprintf('Cp_piv range (overlap):  %.4f to %.4f\n', min(Cp_piv_cl),  max(Cp_piv_cl));
fprintf('Tap points in overlap:   %d\n',            sum(tap_overlap));
fprintf('PIV points in overlap:   %d\n',            sum(piv_overlap));
 
if sum(piv_overlap) == 0
    warning(['No PIV-tap overlap found. ' ...
        'Check piv_x_offset_m — PIV global x range is [%.3f, %.3f] m, ' ...
        'tap global x range is [%.3f, %.3f] m.'], ...
        min(x_piv_m), max(x_piv_m), min(x_tap_u), max(x_tap_u));
end
 
%% -------------------------------------------------------------------------
% 9. PLOTS
% -------------------------------------------------------------------------
 
% --- Cp ---
figure('Name','Cp: PIV vs Taps SW','Position',[100 100 900 450]);
plot(xdelta_piv_cl, Cp_piv_cl, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'PIV (Bernoulli)'); hold on;
plot(xdelta_tap_cl, Cp_tap_cl, 'ro', 'MarkerSize', 6, ...
     'DisplayName', 'Taps (SW)');
xlabel('$\Delta x/\delta_0,\ x_0 = x_{LE}$', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('Pressure coefficient: PIV vs Taps (SW)');
legend('Location', 'best'); grid on; box on;
 
% --- dCp/dx ---
figure('Name','dCp/dx: PIV vs Taps SW','Position',[100 600 900 450]);
plot(xdelta_piv_cl, dCp_piv_sg_cl, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'PIV S-G'); hold on;
plot(xdelta_tap_cl, dCp_tap_raw_cl, 'rs', 'MarkerSize', 5, 'LineWidth', 1.0, ...
     'DisplayName', 'Taps raw');
plot(xdelta_tap_cl, dCp_tap_sg_cl, 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Taps S-G win=%d', sg_win_tap));
xlabel('$\Delta x/\delta_0,\ x_0 = x_{LE}$', 'Interpreter', 'latex');
ylabel('$dC_p/d(x/\delta_0)$', 'Interpreter', 'latex');
title('$dC_p/d(x/\delta_0)$: PIV vs Taps (SW)', 'Interpreter', 'latex');
legend('Location', 'best'); grid on; box on;
 
%% -------------------------------------------------------------------------
% 9b. ADDITIONAL PLOTS — physical x [mm]
% -------------------------------------------------------------------------
x_tap_mm = xc_tap * 1250 + 6530;   % x/c → mm
 
% Sort and clip taps to valid only (same as above but in mm)
[x_tap_mm_sorted, ~] = sort(x_tap_mm, 'ascend');
x_tap_mm_sorted = x_tap_mm_sorted(valid_tap);
[x_tap_mm_u, ~, ic2] = unique(x_tap_mm_sorted, 'sorted');
Cp_tap_u_mm   = accumarray(ic2, Cp_tap_v, [], @mean);
dCp_tap_raw_mm = gradient(Cp_tap_u_mm, x_tap_mm_u);
dCp_tap_sg_mm  = sgolayfilt(dCp_tap_raw_mm, sg_order, sg_win_tap);
 
% Full PIV range — no overlap clipping, x_piv_mm already in global coords
x_piv_mm_cl  = x_piv_mm;
Cp_piv_mm_cl = Cp_piv;
 
% dCp/dx in mm^-1 for PIV
dCp_piv_mm_raw = gradient(Cp_piv_mm_cl, x_piv_mm_cl);
dCp_piv_mm_sg  = sgolayfilt(dCp_piv_mm_raw, sg_order, sg_win_piv);
 
% Show all tap points — no overlap clipping
tap_overlap_mm = true(size(x_tap_mm_u));
 
% --- No-data band limits [mm] ---
gap_lo = 7200;
gap_hi = 7242;
 
% --- Cp physical ---
figure('Name','Cp: PIV vs Taps SW [mm]','Position',[100 100 900 450]);
plot(x_piv_mm_cl, Cp_piv_mm_cl, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'PIV (Bernoulli)'); hold on;
plot(x_tap_mm_u(tap_overlap_mm), Cp_tap_u_mm(tap_overlap_mm), 'ro', 'MarkerSize', 6, ...
     'DisplayName', 'Taps (SW)');
yl = ylim;
patch([gap_lo gap_hi gap_hi gap_lo], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'DisplayName', 'No data (gap)');
ylim(yl);
xlabel('$x$ [mm]', 'Interpreter', 'latex');
ylabel('$C_p$', 'Interpreter', 'latex');
title('Pressure coefficient: PIV vs Taps (SW) — physical x');
legend('Location', 'best'); grid on; box on;
 
% --- dCp/dx physical ---
figure('Name','dCp/dx: PIV vs Taps SW [mm]','Position',[100 600 900 450]);
plot(x_piv_mm_cl, dCp_piv_mm_sg.*1000.*0.08, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'PIV S-G'); hold on;
plot(x_tap_mm_u(tap_overlap_mm), dCp_tap_raw_mm(tap_overlap_mm).*1000*0.08, 'rs', ...
     'MarkerSize', 5, 'LineWidth', 1.0, 'DisplayName', 'Taps raw');
plot(x_tap_mm_u(tap_overlap_mm), dCp_tap_sg_mm(tap_overlap_mm).*1000*0.08, 'r-o', ...
     'MarkerFaceColor', 'r', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Taps S-G win=%d', sg_win_tap));
xline(BL_0_mm*3 + (x_LE+c_chord)*1000); 
yl = [-0.21 0.1];
xlim([5666 9000])
% patch([gap_lo gap_hi gap_hi gap_lo], [yl(1) yl(1) yl(2) yl(2)], ...
%     [1 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
%     'DisplayName', 'No data (gap)');
ylim(yl);
% xline(5866) % [Inlet]  → Pos_1
% xline(6480) %  [FPG max]  → Pos_2
% xline(6921) %  [Crossover]  → Pos_2
% xline(7320) %  [APG max]  → Pos_3
% xline(8000) %  [TE relative]  → Pos_3
% xline(8573) %  [ZPG recovery]  → Pos_4
xline([5866 6480 6921 7320 8000], "k", LineWidth=1.5, HandleVisibility='off') 

xlabel('$x$ [mm]', 'Interpreter', 'latex');
ylabel('$dC_p/dx$', 'Interpreter', 'latex');
title('$dC_p/dx$: PIV vs Taps (SW) — physical x', 'Interpreter', 'latex');
legend('Location', 'best'); grid on; box on;
 
%% -------------------------------------------------------------------------
% 10. SAVE
% -------------------------------------------------------------------------
piv_results.x_m         = x_piv_m(piv_overlap);
piv_results.xdelta      = xdelta_piv_cl;
piv_results.Cp          = Cp_piv_cl;
piv_results.dCp_dxd     = dCp_piv_sg_cl;
piv_results.Uinf        = Uinf_piv;
piv_results.BL0_mm      = BL_0_mm;
piv_results.sg_win      = sg_win_piv;
piv_results.x_offset_m  = piv_x_offset_m;
 
tap_results.x_m          = x_tap_u(tap_overlap);
tap_results.xdelta        = xdelta_tap_cl;
tap_results.Cp            = Cp_tap_cl;
tap_results.dCp_dxd_raw   = dCp_tap_raw_cl;
tap_results.dCp_dxd_sg    = dCp_tap_sg_cl;
tap_results.sg_win        = sg_win_tap;
tap_results.c_chord_m     = c_chord;
tap_results.x_LE_m        = x_LE;
 
% save(out_path, 'piv_results', 'tap_results');
fprintf('\nResults saved to %s\n', out_path);
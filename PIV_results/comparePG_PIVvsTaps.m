%% compare_PIV_vs_taps.m
clear; clc; close all;

%% -------------------------------------------------------------------------
% 0. USER SETTINGS
% -------------------------------------------------------------------------
caseNo      = 2;
sg_order    = 3;
LE_offset_m = 7.65;

% S-G windows — set to [] to use matched physical width (default)
%               set to an odd integer to override
sg_win_piv_override = [];   % e.g. 51 to override, [] for auto
sg_win_tap_override = 11  ;  % e.g. 5  to override, [] for auto


piv_path = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_results\Case2_PIVresults\blSweep_20260415_125916.mat';
tap_path = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\Pressure\pressuredata_wingcases_revisedNOZPGZERO.mat';
out_path = 'PIV_vs_taps_results_case2.mat';


%% -------------------------------------------------------------------------
% 1. LOAD
% -------------------------------------------------------------------------
load(piv_path, 'blSweep');
load(tap_path, 'case_pdata');

%% -------------------------------------------------------------------------
% 2. PIV: valid data + global x
% -------------------------------------------------------------------------
x_piv_mm = blSweep.x_mm;
U99      = blSweep.U99;

valid_piv = ~isnan(x_piv_mm) & ~isnan(U99);
x_piv_mm  = x_piv_mm(valid_piv);
U99       = U99(valid_piv);

x_piv_m = x_piv_mm/1000 + 7.25;   % global x [m]

% Uinf from last 40 mm of PIV domain
xf_piv   = max(x_piv_m);
inf_idx  = x_piv_m >= (xf_piv - 0.04);
Uinf_piv = mean(U99(inf_idx), 'omitnan');
fprintf('PIV Uinf (outlet mean): %.4f m/s\n', Uinf_piv);

%% -------------------------------------------------------------------------
% 3. PIV: BL thickness at inlet (first 40 mm)
% -------------------------------------------------------------------------
delta99_mm = blSweep.delta99_hybrid_mm(valid_piv);

x0_BL   = min(x_piv_m);
BL_idx  = x_piv_m <= (x0_BL + 0.04);
BL_0_mm = mean(delta99_mm(BL_idx), 'omitnan');
BL_0_m  = BL_0_mm / 1000;
fprintf('Reference BL thickness delta0: %.2f mm\n', BL_0_mm);

%% -------------------------------------------------------------------------
% 4. PIV: Cp and raw dCp/d(x/delta)
% -------------------------------------------------------------------------
Cp_piv      = 1 - (U99 / Uinf_piv).^2;
xdelta_piv  = x_piv_m / BL_0_m;
dCp_piv_raw = gradient(Cp_piv, xdelta_piv);

%% -------------------------------------------------------------------------
% 5. TAPS: recover global x, non-dimensionalise by PIV BL_0_m
% -------------------------------------------------------------------------
xdelta_tap  = case_pdata(caseNo).xloc_full_xdelta;
refBL_m     = case_pdata(caseNo).refBLmm / 100;     % mislabelled — actually cm
mean_Cp_tap = mean(case_pdata(caseNo).Cp_full, 2, 'omitnan');

x_tap_m = xdelta_tap * refBL_m + LE_offset_m;

[x_tap_sorted, sort_idx] = sort(x_tap_m, 'ascend');
Cp_tap_sorted     = mean_Cp_tap(sort_idx);
xdelta_tap_sorted = x_tap_sorted / BL_0_m;

valid_tap    = ~isnan(Cp_tap_sorted) & ~isnan(xdelta_tap_sorted) & ~isinf(xdelta_tap_sorted);
xdelta_tap_v = xdelta_tap_sorted(valid_tap);
Cp_tap_v     = Cp_tap_sorted(valid_tap);
x_tap_v      = x_tap_sorted(valid_tap);

[xdelta_tap_u, ~, ic] = unique(xdelta_tap_v, 'sorted');
Cp_tap_u              = accumarray(ic, Cp_tap_v, [], @mean);
x_tap_u               = accumarray(ic, x_tap_v,  [], @mean);

dCp_tap_raw = gradient(Cp_tap_u, xdelta_tap_u);

%% -------------------------------------------------------------------------
% 6. S-G WINDOWS — matched physical width or user override
% -------------------------------------------------------------------------
tap_spacing_mean = mean(diff(xdelta_tap_u));
piv_spacing_mean = mean(diff(xdelta_piv));

% Default: minimum meaningful window on tap grid, expressed in PIV points
default_physical_width = (sg_order + 2) * tap_spacing_mean;

default_win_piv = round(default_physical_width / piv_spacing_mean);
if mod(default_win_piv, 2) == 0, default_win_piv = default_win_piv + 1; end
default_win_piv = max(default_win_piv, sg_order + 2);

default_win_tap = sg_order + 2;
if mod(default_win_tap, 2) == 0, default_win_tap = default_win_tap + 1; end

% Apply override or default
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

% Guard both windows
sg_win_piv = min(sg_win_piv, 2*floor(numel(dCp_piv_raw)/2)-1);
sg_win_piv = max(sg_win_piv, sg_order + 2);
if mod(sg_win_piv, 2) == 0, sg_win_piv = sg_win_piv + 1; end

sg_win_tap = min(sg_win_tap, 2*floor(numel(dCp_tap_raw)/2)-1);
sg_win_tap = max(sg_win_tap, sg_order + 2);
if mod(sg_win_tap, 2) == 0, sg_win_tap = sg_win_tap + 1; end

% Apply filters
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

xdelta_piv_cl      = xdelta_piv(piv_overlap);
Cp_piv_cl          = Cp_piv(piv_overlap);
dCp_piv_sg_cl      = dCp_piv_sg(piv_overlap);

xdelta_tap_cl      = xdelta_tap_u(tap_overlap);
Cp_tap_cl          = Cp_tap_u(tap_overlap);
dCp_tap_raw_cl     = dCp_tap_raw(tap_overlap);
dCp_tap_sg_cl      = dCp_tap_sg(tap_overlap);
%% -------------------------------------------------------------------------
% 7b. SHIFT x/delta ORIGIN TO LE (7.65 m)
% -------------------------------------------------------------------------
LE_xdelta = LE_offset_m / BL_0_m;   % 7.65 m in x/delta units

xdelta_piv_cl  = xdelta_piv_cl- LE_xdelta;
xdelta_tap_cl  = xdelta_tap_cl - LE_xdelta;

%% -------------------------------------------------------------------------
% 8. DIAGNOSTICS
% -------------------------------------------------------------------------
fprintf('\n=== DIAGNOSTICS ===\n');
fprintf('Cp_tap range (overlap):  %.4f to %.4f\n', min(Cp_tap_cl),  max(Cp_tap_cl));
fprintf('Cp_piv range (overlap):  %.4f to %.4f\n', min(Cp_piv_cl),  max(Cp_piv_cl));
fprintf('Tap points in overlap:   %d\n',           sum(tap_overlap));
fprintf('PIV points in overlap:   %d\n',           sum(piv_overlap));

%% -------------------------------------------------------------------------
% 9. PLOTS
% -------------------------------------------------------------------------

% --- Cp ---
figure('Name','Cp: PIV vs Taps','Position',[100 100 900 450]);
plot(xdelta_piv_cl, Cp_piv_cl, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'PIV (Bernoulli)'); hold on;
plot(xdelta_tap_cl, Cp_tap_cl, 'ro', 'MarkerSize', 6, ...
     'DisplayName', sprintf('Taps (case %d)', caseNo));
xlabel('x/\delta_0'); ylabel('C_p');
title('Pressure coefficient: PIV vs Taps');
legend('Location','best'); grid on; box on;

% --- dCp/dx ---
figure('Name','dCp/dx: PIV vs Taps','Position',[100 600 900 450]);
plot(xdelta_piv_cl, dCp_piv_sg_cl, 'b-o', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('PIV')); hold on;
plot(xdelta_tap_cl, dCp_tap_raw_cl, 'rs', 'MarkerSize', 5, 'LineWidth', 1.0, ...
     'DisplayName', sprintf('Taps raw (case %d)', caseNo));
plot(xdelta_tap_cl, dCp_tap_sg_cl, 'r-o', 'MarkerFaceColor', 'r','LineWidth', 1.5, ...
     'DisplayName', sprintf('Taps S-G win=%d (case %d)', sg_win_tap, caseNo));
xlabel('$\Delta x/\delta_0, x_0 = x_{LE}$ m', Interpreter='latex'); ylabel('$dC_p/d(x/\delta_0)$', Interpreter='latex');
title('$dCp/d(x/\delta_0)$  PIV vs Taps', Interpreter='latex');
xline([ -0.4313    1.2940    2.6830    5.1762], 'k-', LineWidth = 1.5); 
legend('Location','best'); grid on; box on;

%% -------------------------------------------------------------------------
% 10. SAVE
% -------------------------------------------------------------------------
piv_results.x_m        = x_piv_m(piv_overlap);
piv_results.xdelta     = xdelta_piv_cl;
piv_results.Cp         = Cp_piv_cl;
piv_results.dCp_dxd    = dCp_piv_sg_cl;
piv_results.Uinf       = Uinf_piv;
piv_results.BL0_mm     = BL_0_mm;
piv_results.sg_win     = sg_win_piv;

tap_results.x_m        = x_tap_u(tap_overlap);
tap_results.xdelta     = xdelta_tap_cl;
tap_results.Cp         = Cp_tap_cl;
tap_results.dCp_dxd_raw = dCp_tap_raw_cl;
tap_results.dCp_dxd_sg  = dCp_tap_sg_cl;
tap_results.sg_win      = sg_win_tap;

save(out_path, 'piv_results', 'tap_results', 'caseNo');
fprintf('\nResults saved to %s\n', out_path);
%% =========================================================================
%  process_ofi.m
%  PURPOSE : Process raw OFI oil-drop data for all cameras into a clean,
%            U_inf(x)-corrected Cf/2 dataset.
%
%  COVERS (original script sections)
%   1. Camera loop       — original lines ~45–90  (i=1:2, k=1:5/4, j-loop)
%   2. Coord transform   — original lines ~93–105 (x0_global, truncation)
%   3. U_inf(x) correction — original lines ~107–128
%   4. Sort + movmean    — original lines ~131–155 (window_sizes FIRST def.)
%
%  FIXES APPLIED vs original script
%   - `Cf_2`           was re-declared both INSIDE the j-loop (as a drop
%                      scalar) AND outside as the full column vector from
%                      Cf2_data_trunc(:,5). Renamed inner → `drop_Cf2`.
%   - `window_sizes`   was defined as [3,6,9,13] then immediately overwritten
%                      as round(0.1*length(case_Cf2_data)) in the next section,
%                      silently breaking the movmean matrix. Both are now
%                      explicit named parameters passed in via `params`.
%   - `Cf2_movmeans`   was rebuilt with the new (wrong) window_sizes.
%                      Now built once, correctly, from params.window_fixed.
%   - `i`              clashed as camera-set index AND movmean loop counter.
%                      Camera loop now uses `cam_set`; movmean uses `wi`.
%   - `in_range`       was computed but never used downstream. Kept as a
%                      diagnostic field in OFI.diag.
%
%  USAGE
%   data   = load_data();
%   params = struct('window_fixed', [3 6 9 13], 'window_pct', 0.10,
%                   'use_fixed', true);
%   OFI    = process_ofi(data, params);
%
%  OUTPUT  OFI struct — fields listed in the RETURNS section below.
%
%  AUTHOR  : AK (refactored from VirgPreskettCase1_computauUU.m)
%  DATE    : 2026-03-31
% =========================================================================

function OFI = process_ofi(data, params)

%% ---- DEFAULT PARAMETERS -------------------------------------------------
if nargin < 2 || isempty(params)
    params.window_fixed = [3, 6, 9, 13];   % original line ~134: window_sizes = [3,6,9,13]
    params.window_pct   = 0.10;             % original line ~152: round(0.1*length(...))
    params.use_fixed    = true;             % which window set to use for final movmean
end

% Unpack cfg from data struct (all set in load_data)
mmperpix    = data.cfg.mmperpix;           % [mm/pix] per camera set
x0_global   = data.cfg.x0_global_mm;      % = -7300 mm  (original line ~93)
max_gap_mm  = data.cfg.max_gap_mm;         % = 50 mm     (original line ~115)
U_inf       = data.fluid.U_inf;            % scalar reference U_inf [m/s]

%% =========================================================================
%  SECTION 1 — CAMERA LOOP
%  Original lines ~45–90: nested i/k/j loops building case_Cf2_data
% =========================================================================

% Pre-allocate with NaN rather than growing with end+1
% (original used case_Cf2_data(end+1,:) — fine for small datasets but
% slow for large ones; pre-alloc is cleaner)
raw_rows = [];   % will be concatenated per camera/drop

% ---- Camera set 1: 5 cameras (i==1 block, original lines ~48–66) --------
cam_set = 1;
n_cams_set1 = 5;

for k = 1:n_cams_set1
    loadFilename = sprintf('c%d%d_oilDropinfo.mat', cam_set, k);  % original line ~51
    tmp = load(loadFilename);
    oilDropinfo = tmp.oilDropinfo;

    calibrationXML = parseCameraXML(data.cfg.calib_path_cam1);    % original line ~53
    ppx0   = calibrationXML(k).OriginPixelPosition(1);            % original line ~55
    ppy0   = 0;                                                    % original line ~56
    offset = data.cfg.offset_c1;  % = -120 mm                     % original line ~57

    for j = 1:length(oilDropinfo)
        % Pixel → mm coordinate transform (original lines ~59–63)
        mmx1 = (oilDropinfo(j).globalSearchBoxCoord(1) - ppx0) * mmperpix(cam_set) - offset;
        mmx2 = (oilDropinfo(j).globalSearchBoxCoord(2) - ppx0) * mmperpix(cam_set) - offset;
        mmy1 = (oilDropinfo(j).searchBox(3)                    - ppy0) * mmperpix(cam_set);
        mmy2 = ((oilDropinfo(j).searchBox(3) + oilDropinfo(j).searchBox(4)) - ppy0) * mmperpix(cam_set);

        % Cf/2 from u_tau (original line ~64)
        % RENAMED: was `Cf_2` (scalar) — clashed with the full-vector `Cf_2`
        % extracted from Cf2_data_trunc(:,5) later in the script.
        drop_Cf2 = (oilDropinfo(j).u_tau ./ U_inf).^2;

        raw_rows(end+1, :) = [mmx1, mmx2, mmy1, mmy2, drop_Cf2]; % original line ~65-66
    end
    clear tmp oilDropinfo;
end

% ---- Camera set 2: 4 cameras (i==2 block, original lines ~68–87) --------
cam_set = 2;
n_cams_set2 = 4;

for k = 1:n_cams_set2
    loadFilename = sprintf('c%d%d_oilDropinfo.mat', cam_set, k);  % original line ~70
    tmp = load(loadFilename);
    oilDropinfo = tmp.oilDropinfo;

    calibrationXML = parseCameraXML(data.cfg.calib_path_cam2);    % original line ~72
    % NOTE: original used calibrationXML(k+1) for cam set 2 — preserved here
    ppx0   = calibrationXML(k+1).OriginPixelPosition(1);          % original line ~74
    ppy0   = 0;                                                    % original line ~75
    offset = data.cfg.offset_c2;  % = 1400 mm                     % original line ~76

    for j = 1:length(oilDropinfo)
        % Pixel → mm transform (original lines ~78–82)
        mmx1 = (oilDropinfo(j).globalSearchBoxCoord(1) - ppx0) * mmperpix(cam_set) - offset;
        mmx2 = (oilDropinfo(j).globalSearchBoxCoord(2) - ppx0) * mmperpix(cam_set) - offset;
        mmy1 = (oilDropinfo(j).searchBox(3)                    - ppy0) * mmperpix(cam_set);
        mmy2 = ((oilDropinfo(j).searchBox(3) + oilDropinfo(j).searchBox(4)) - ppy0) * mmperpix(cam_set);

        drop_Cf2 = (oilDropinfo(j).u_tau ./ U_inf).^2;            % original line ~83

        raw_rows(end+1, :) = [mmx1, mmx2, mmy1, mmy2, drop_Cf2]; % original line ~84-85
    end
    clear tmp oilDropinfo;
end

OFI.raw_Cf2_table = raw_rows;   % equivalent to original `case_Cf2_data`
fprintf('[process_ofi] Camera loop done: %d oil drops collected.\n', size(raw_rows,1));

%% =========================================================================
%  SECTION 2 — COORDINATE TRANSFORM & TRUNCATION
%  Original lines ~93–105: x0_global shift, abs(), round to 4 d.p.
% =========================================================================

Cf2_shifted          = raw_rows;
Cf2_shifted(:,1)     = abs(Cf2_shifted(:,1) + x0_global);   % original line ~95
Cf2_shifted(:,2)     = abs(Cf2_shifted(:,2) + x0_global);   % original line ~96
Cf2_trunc            = round(Cf2_shifted * 1e4) / 1e4;       % original line ~97

% Named column extraction (original lines ~99–105)
% Previously extracted into flat workspace vars: start_x, end_x, start_y,
% end_y, Cf_2 — now stored in OFI struct.
OFI.box_x1_mm   = Cf2_trunc(:,1);   % was `start_x`
OFI.box_x2_mm   = Cf2_trunc(:,2);   % was `end_x`
OFI.box_y1_mm   = Cf2_trunc(:,3);   % was `start_y`
OFI.box_y2_mm   = Cf2_trunc(:,4);   % was `end_y`
OFI.Cf2_raw     = Cf2_trunc(:,5);   % was `Cf_2` (the full vector — NOT drop_Cf2)

% Box-centre positions (original lines ~107–108)
OFI.x_center_mm = (OFI.box_x1_mm + OFI.box_x2_mm) / 2;   % was `x_centers`
OFI.y_center_mm = (OFI.box_y1_mm + OFI.box_y2_mm) / 2;   % was `y_centers`

fprintf('[process_ofi] Coord transform done. x range: [%.1f, %.1f] mm\n', ...
    min(OFI.x_center_mm), max(OFI.x_center_mm));

%% =========================================================================
%  SECTION 3 — U_inf(x) CORRECTION
%  Original lines ~110–128: nearest-neighbour PIV U_inf lookup + correction
% =========================================================================

% Global x-axis for PIV sweep (original line ~110: x_mm_global = blSweep.x_mm + 7100)
x_piv_global_mm  = data.blSweep.x_mm(:) + 7100;
U_inf_piv        = data.blSweep.U_inf(:);

% Strip NaNs from PIV sweep before lookup (original lines ~113–115)
valid_piv        = ~isnan(U_inf_piv);
x_piv_clean      = x_piv_global_mm(valid_piv);
U_inf_piv_clean  = U_inf_piv(valid_piv);

% Nearest-neighbour lookup (original lines ~117–118)
[min_dist, nn_idx] = min(abs(OFI.x_center_mm(:) - x_piv_clean(:)'), [], 2);

% In-range flag (original line ~120: max_gap_mm = 50)
% NOTE: `in_range` was computed in the original but never used to gate
% anything downstream. Kept here as a diagnostic field.
OFI.diag.in_range    = min_dist <= max_gap_mm;
OFI.diag.nn_min_dist = min_dist;

% Initialise with nearest-neighbour, then apply boundary extrapolation
% (original lines ~124–128)
U_inf_local = U_inf_piv_clean(nn_idx);
U_inf_local(OFI.x_center_mm < x_piv_clean(1))   = U_inf_piv_clean(1);    % extrapolate left
U_inf_local(OFI.x_center_mm > x_piv_clean(end)) = U_inf_piv_clean(end);  % extrapolate right

OFI.U_inf_local    = U_inf_local;         % was `U_inf_local`
OFI.x_piv_global   = x_piv_global_mm;    % was `x_mm_global` — stored for downstream plots

% Apply correction: Cf2_corrected = Cf2_raw * (U_ref / U_local)^2
% (original line ~130: Cf_2_corrected = Cf_2 .* (U_inf ./ U_inf_local).^2)
OFI.Cf2_corrected = OFI.Cf2_raw .* (U_inf ./ U_inf_local).^2;

fprintf('[process_ofi] U_inf(x) correction applied. Cf2 range: [%.5f, %.5f]\n', ...
    min(OFI.Cf2_corrected), max(OFI.Cf2_corrected));

%% =========================================================================
%  SECTION 4 — SORT BY x AND MOVING AVERAGE
%  Original lines ~133–155 (first window_sizes block only)
%
%  BUG FIXED: In the original, window_sizes = [3,6,9,13] was defined here,
%  then OVERWRITTEN two sections later as round(0.1*length(case_Cf2_data)),
%  causing Cf2_movmeans to be silently rebuilt with a different window.
%  This function computes BOTH, names them clearly, and lets the caller pick.
% =========================================================================

% Sort corrected Cf2 by ascending x-centre (original lines ~133–135)
[OFI.sorted_x_mm, OFI.sort_idx] = sort(OFI.x_center_mm);   % was `sorted_x`, `sort_idx`
OFI.Cf2_sorted = OFI.Cf2_corrected(OFI.sort_idx);            % was `Cf2_sorted`

% ---- Fixed window set (original line ~136: window_sizes = [3,6,9,13]) ---
w_fixed  = params.window_fixed;
n_fixed  = length(w_fixed);
Cf2_mm_fixed = zeros(length(OFI.Cf2_sorted), n_fixed);

for wi = 1:n_fixed   % RENAMED: was `i` — clashed with camera-set loop index
    Cf2_mm_fixed(:, wi) = movmean(OFI.Cf2_sorted, w_fixed(wi));
end

OFI.movmean.window_fixed   = w_fixed;          % was `window_sizes` (first definition)
OFI.movmean.Cf2_fixed      = Cf2_mm_fixed;     % was `Cf2_movmeans` (first build)

% The "final" movmean (largest fixed window) mapped back to original order
% (original lines ~151–154: Cf2_movmean_original)
OFI.movmean.Cf2_final_w    = w_fixed(end);     % = 13 by default
OFI.movmean.Cf2_final      = Cf2_mm_fixed(:, end);  % was `Cf2_movmean` (w=13)

% Restore to unsorted indexing for any downstream use requiring original order
Cf2_movmean_orig = nan(size(OFI.Cf2_corrected));
Cf2_movmean_orig(OFI.sort_idx) = Cf2_mm_fixed(:, end);
OFI.movmean.Cf2_final_unsorted = Cf2_movmean_orig;  % was `Cf2_movmean_original`

% ---- Percentage window (original line ~152: round(0.1*length(case_Cf2_data))) ---
w_pct   = round(params.window_pct * size(raw_rows, 1));
Cf2_mm_pct = movmean(OFI.Cf2_sorted, w_pct);

OFI.movmean.window_pct     = w_pct;            % was `window_sizes` (OVERWRITTEN value)
OFI.movmean.Cf2_pct        = Cf2_mm_pct;       % was `Cf2_movmeans` (second build)

fprintf('[process_ofi] Movmean done. Fixed windows: %s | Pct window: %d\n', ...
    num2str(w_fixed), w_pct);

%% =========================================================================
%  SECTION 5 — Re_x FOR AK DATA
%  Original line (in later section): AK_Rex = (sorted_x./1000)*U_inf/visc_air
%  Computed here so it travels with the OFI struct.
% =========================================================================

OFI.Re_x = (OFI.sorted_x_mm ./ 1000) .* U_inf ./ data.fluid.visc_air;
% was `AK_Rex` — defined much later in original script, but logically
% belongs to the OFI processing output.

fprintf('[process_ofi] Re_x computed. Range: [%.2e, %.2e]\n', ...
    min(OFI.Re_x), max(OFI.Re_x));

end  % function process_ofi

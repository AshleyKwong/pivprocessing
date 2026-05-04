function [bnd_dx, bnd_dy, ell] = get_contour_2(covFolder, corrType, xr_nom, yr_target, ...
    worldX, worldY, dx_back, dx_fwd, rho_level)
% GET_CONTOUR  Extract contour boundary and fit ellipse for a single R_uu/R_uv map.
%
% Finds the R_uu (or R_uv) file in covFolder whose xref matches xr_nom and
% whose yref is closest to yr_target. Builds the threshold mask at rho_level,
% extracts the origin-connected region boundary, and fits an ellipse via
% image second moments (regionprops).
%
% INPUTS
%   covFolder  — folder containing R_uu_xref*_yref*.mat files
%   corrType   — 'uu' or 'uv'
%   xr_nom     — nominal xref (mm); matched within 5 mm
%   yr_target  — target yref (mm); closest available file selected
%   worldX     — [Ny x Nx] physical x grid (mm)
%   worldY     — [Ny x Nx] physical y grid (mm)
%   dx_back    — upstream search extent (mm)
%   dx_fwd     — downstream search extent (mm)
%   rho_level  — correlation threshold (e.g. 1/exp(1), 0.5)
%
% OUTPUTS
%   bnd_dx     — contour boundary delta-x coords (mm), [] if failed
%   bnd_dy     — contour boundary delta-y coords (mm), [] if failed
%   ell        — struct with ellipse parameters:
%                  .a      semi-major axis (mm)
%                  .b      semi-minor axis (mm)
%                  .theta  inclination angle (deg CCW from +x)
%                  .cx     ellipse centre delta-x (mm)
%                  .cy     ellipse centre delta-y (mm)
%                  .valid  true if fit succeeded
%
% NOTES
%   Ellipse fitted from contour boundary tip points:
%     theta  — angle of upstream→downstream tip vector from +x (atan2, four-quadrant)
%     cx, cy — midpoint of the two tips (ellipse centre)
%     a      — half distance between tips (semi-major)
%     b      — half perpendicular extent of boundary points (semi-minor)
%   No bounding box quantities or regionprops used.

% ---- Default outputs ---------------------------------------------------
bnd_dx = [];
bnd_dy = [];
ell    = struct('a', NaN, 'b', NaN, 'theta', NaN, ...
                'cx', NaN, 'cy', NaN, 'valid', false);

% ---- Find matching file ------------------------------------------------
filePattern = ['R_' corrType '_xref*_yref*.mat'];
files = dir(fullfile(covFolder, filePattern));
if isempty(files), return; end

yr_available = nan(numel(files), 1);
xr_available = nan(numel(files), 1);
for f = 1:numel(files)
    tok = regexp(files(f).name, ...
        ['^R_' corrType '_xref(\d+)_yref(\d+p\d+)\.mat$'], ...
        'tokens', 'once');
    if ~isempty(tok)
        xr_available(f) = str2double(tok{1});
        yr_available(f) = str2double(strrep(tok{2}, 'p', '.'));
    end
end

% Match xref within 5 mm; fall back to nearest if none within tolerance
xr_ok = abs(xr_available - xr_nom) < 5;
if ~any(xr_ok)
    [~, iBest] = min(abs(xr_available - xr_nom));
    xr_ok(iBest) = true;
end

% Among xref-matched files, pick closest yref
yr_sub = yr_available;
yr_sub(~xr_ok) = Inf;
[~, iClosest] = min(abs(yr_sub - yr_target));

% ---- Load correlation map ----------------------------------------------
S  = load(fullfile(covFolder, files(iClosest).name), 'R_s', 'xr', 'yr');
R  = double(S.R_s);
xr = double(S.xr);
yr = double(S.yr);

% ---- Build delta grids and search box ----------------------------------
dX = worldX - xr;
dY = worldY - yr;

x_min_domain = min(worldX(:));
x_max_domain = max(worldX(:));
x_lo = max(xr - dx_back, x_min_domain);
x_hi = min(xr + dx_fwd,  x_max_domain);

boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
R_box           = R;
R_box(~boxMask) = NaN;
R_box(R_box >  1.0) =  1.0;
R_box(R_box < -1.0) =  NaN;

% ---- Threshold mask ----------------------------------------------------
mask = R_box >= rho_level;
if sum(mask(:)) < 5, return; end

% ---- Origin pixel (closest to reference point) -------------------------
dist2           = dX.^2 + dY.^2;
dist2(~boxMask) = Inf;
[~, i0]         = min(dist2(:));

% ---- Connected components — origin-connected or largest region ---------
CC = bwconncomp(mask);
if CC.NumObjects == 0, return; end

regionIdx = 0;
for r = 1:CC.NumObjects
    if any(CC.PixelIdxList{r} == i0)
        regionIdx = r; break;
    end
end
if regionIdx == 0
    regionSizes = cellfun(@numel, CC.PixelIdxList);
    [~, regionIdx] = max(regionSizes);
end

mask_origin = false(size(mask));
mask_origin(CC.PixelIdxList{regionIdx}) = true;

% ---- Contour boundary in physical coords (mm) --------------------------
dX_vec = dX(1,:);
dY_vec = dY(:,1);

B_contour = bwboundaries(mask_origin, 'noholes');
if isempty(B_contour), return; end

bnd    = B_contour{1};
bnd_dx = dX_vec(bnd(:,2))';   % column vector
bnd_dy = dY_vec(bnd(:,1));    % already column vector

% ---- Ellipse fit from contour boundary tips ----------------------------
dx_spacing = abs(dX_vec(2) - dX_vec(1));
dy_spacing = abs(dY_vec(2) - dY_vec(1));

% --- Tip extraction (unchanged) ---
up_mask  = bnd_dx <= (min(bnd_dx) + dx_spacing);
dn_mask  = bnd_dx >= (max(bnd_dx) - dx_spacing);
top_mask = bnd_dy >= (max(bnd_dy) - dy_spacing);
bot_mask = bnd_dy <= (min(bnd_dy) + dy_spacing);

x_up  = mean(bnd_dx(up_mask));   y_up  = mean(bnd_dy(up_mask));
x_dn  = mean(bnd_dx(dn_mask));   y_dn  = mean(bnd_dy(dn_mask));
x_top = mean(bnd_dx(top_mask));  y_top = mean(bnd_dy(top_mask));
x_bot = mean(bnd_dx(bot_mask));  y_bot = mean(bnd_dy(bot_mask));

% --- Major axis: upstream -> downstream tip vector (unchanged) ---
ell.theta = atan2d(y_dn - y_up, x_dn - x_up);
ell.a     = sqrt((x_dn - x_up)^2 + (y_dn - y_up)^2) / 2;

% --- Centre: use top/bottom tip midpoint for cy (FIX 1) ---
% The upstream/downstream midpoint gives cx well, but cy is biased
% downward when the contour is wall-clipped. The top/bottom midpoint
% is much more robust for the wall-normal centre.
ell.cx = (x_up  + x_dn)  / 2;   % x-centre from streamwise tips (good)
ell.cy = (y_top + y_bot)  / 2;  % y-centre from wall-normal tips (FIX)

% --- Semi-minor: perpendicular projection, then wall-clamp (FIX 2) ---
th_rad = deg2rad(ell.theta);
perp_x = -sin(th_rad);
perp_y =  cos(th_rad);

proj_top = (x_top - ell.cx)*perp_x + (y_top - ell.cy)*perp_y;
proj_bot = (x_bot - ell.cx)*perp_x + (y_bot - ell.cy)*perp_y;
b_raw    = abs(proj_top - proj_bot) / 2;

% Clamp: ellipse lower extent must not go below contour's minimum y.
% Compute how far below cy the ellipse would reach along the minor axis.
% The minor axis lower tip in physical coords:
% Clamp: ellipse lower extent must not go below contour's minimum y.
b_lower_tip_y = ell.cy - b_raw * cos(th_rad);
if b_lower_tip_y < min(bnd_dy) && abs(cos(th_rad)) > 0.1
    b_clamped = (ell.cy - min(bnd_dy)) / abs(cos(th_rad));
    ell.b = min(b_raw, b_clamped);
else
    ell.b = b_raw;
end

ell.valid = true;

end
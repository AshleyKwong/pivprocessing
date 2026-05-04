function [bnd_dx, bnd_dy, ell] = get_contour(covFolder, corrType, xr_nom, yr_target, ...
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

B_contour = bwboundaries(mask_origin);
if isempty(B_contour), return; end

bnd    = B_contour{1};
bnd_dx = dX_vec(bnd(:,2))';   % column vector
bnd_dy = dY_vec(bnd(:,1));    % already column vector

% ---- Ellipse fit from contour boundary tips ----------------------------
% Method: find the leftmost and rightmost points on the boundary,
% fit the major axis through them, then compute the semi-minor axis
% as the half-width of the contour perpendicular to that axis.
%
% theta  = angle of left→right tip vector from +x axis (deg)
% cx, cy = midpoint of left and right tips (ellipse centre)
% a      = half the distance between left and right tips (semi-major)
% b      = half the perpendicular extent of boundary points (semi-minor)

% Upstream tip: mean position of boundary points within one grid cell of
% the leftmost x. Downstream tip: same for rightmost x.
dx_spacing = abs(dX_vec(2) - dX_vec(1));

up_mask = bnd_dx <= (min(bnd_dx) + dx_spacing);
dn_mask = bnd_dx >= (max(bnd_dx) - dx_spacing);

x_up = mean(bnd_dx(up_mask));   y_up = mean(bnd_dy(up_mask));
x_dn = mean(bnd_dx(dn_mask));   y_dn = mean(bnd_dy(dn_mask));

% Ellipse centre and major axis angle
ell.cx    = (x_up + x_dn) / 2;
ell.cy    = (y_up + y_dn) / 2;
ell.theta = atan2d(y_dn - y_up, x_dn - x_up);   % deg, four-quadrant

% Semi-major: half distance between tips
ell.a = sqrt((x_dn - x_up)^2 + (y_dn - y_up)^2) / 2;

% Semi-minor: project all boundary points onto the perpendicular axis,
% take half the total perpendicular extent
th_rad   = deg2rad(ell.theta);
perp_vec = [-sin(th_rad), cos(th_rad)];   % unit vector perpendicular to major axis
proj     = (bnd_dx - ell.cx) * perp_vec(1) + (bnd_dy - ell.cy) * perp_vec(2);
ell.b    = (max(proj) - min(proj)) / 2;

ell.valid = true;

end
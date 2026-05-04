function [Lx, Ly] = compute_length_scale_2D(R, x, y, ix_ref, iy_ref)
% COMPUTE_LENGTH_SCALE_2D
% Computes streamwise (Lx) and wall-normal (Ly) integral length scales
% from a pre-normalised 2D correlation map R(Ny x Nx).
%
% INPUTS:
%   R      - (Ny x Nx) normalised correlation field R(x_ref,y_ref; x,y)
%   x      - (1 x Nx) streamwise coordinate vector (mm or m, auto-detected)
%   y      - (Ny x 1) wall-normal coordinate vector (mm or m, auto-detected)
%   ix_ref - column index of reference point in R
%   iy_ref - row index of reference point in R
%
% OUTPUTS:
%   Lx - streamwise integral length scale (metres)
%   Ly - wall-normal integral length scale (metres)
%
% METHOD:
%   One-sided integration forward from the reference point to the first
%   zero-crossing of R in each direction.
%
% *** TODO: Review choices before using Lx/Ly in analysis: ***
% *** 1. One-sided (current) vs two-sided integration      ***
% *** 2. Exponential fit instead of direct trapz           ***
% *** 3. Behaviour when R never crosses zero               ***
% *** 4. Normalise output by delta, theta, or y+           ***

% --- Input checks ---
assert(ismatrix(R),            'compute_length_scale_2D: R must be a 2D matrix.');
assert(isvector(x),            'compute_length_scale_2D: x must be a vector.');
assert(isvector(y),            'compute_length_scale_2D: y must be a vector.');
assert(ix_ref >= 1 && ix_ref <= size(R,2), 'compute_length_scale_2D: ix_ref out of range.');
assert(iy_ref >= 1 && iy_ref <= size(R,1), 'compute_length_scale_2D: iy_ref out of range.');

% Ensure correct orientation
x = x(:)';   % force (1 x Nx) row
y = y(:);    % force (Ny x 1) column

% --- Auto unit detection and conversion ---
if max(abs(x)) > 10
    fprintf('compute_length_scale_2D: x appears to be in mm (max = %.1f mm) — converting to m.\n', max(abs(x)));
    x = x ./ 1000;
end
if max(abs(y)) > 10
    fprintf('compute_length_scale_2D: y appears to be in mm (max = %.1f mm) — converting to m.\n', max(abs(y)));
    y = y ./ 1000;
end

% ----------------------------------------------------------
% Lx: streamwise length scale
% Slice R along x at fixed y = y_ref → (1 x Nx) vector
% Integrate forward from ix_ref to first zero-crossing
% ----------------------------------------------------------
R_x     = R(iy_ref, :);           % (1 x Nx) row at y_ref
R_x_fwd = R_x(ix_ref:end);        % forward (downstream) from x_ref
x_fwd   = x(ix_ref:end);

zero_idx_x = find(R_x_fwd <= 0, 1, 'first');
if isempty(zero_idx_x)
    warning('compute_length_scale_2D: Lx — no zero-crossing found downstream of x_ref. Lx = NaN.');
    Lx = NaN;
else
    Lx = abs(trapz(x_fwd(1:zero_idx_x), R_x_fwd(1:zero_idx_x)));
end

% ----------------------------------------------------------
% Ly: wall-normal length scale
% Slice R along y at fixed x = x_ref → (Ny x 1) vector
% Integrate away from wall (increasing y) from iy_ref to first zero-crossing
% ----------------------------------------------------------
R_y     = R(:, ix_ref);           % (Ny x 1) column at x_ref
R_y_fwd = R_y(iy_ref:end);        % away from wall from y_ref
y_fwd   = y(iy_ref:end);

zero_idx_y = find(R_y_fwd <= 0, 1, 'first');
if isempty(zero_idx_y)
    warning('compute_length_scale_2D: Ly — no zero-crossing found above y_ref. Ly = NaN.');
    Ly = NaN;
else
    Ly = abs(trapz(y_fwd(1:zero_idx_y), R_y_fwd(1:zero_idx_y)));
end

fprintf('compute_length_scale_2D: Lx = %.4f m | Ly = %.4f m\n', Lx, Ly);
end

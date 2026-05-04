function [worldX, worldY, U_merged, V_merged] = merge_cameras_python_style_mean(...
    windowCenters, velocityU, velocityV, masks, window_type, taper_param, worldGrid)
% MERGE_CAMERAS_PYTHON_STYLE Merge multi-camera PIV using Python-style blending
%
% Inputs:
%   windowCenters - struct with .x1_mm{cam}, .x2_mm{cam} (coordinates per camera)
%   velocityU - cell array {cam} of U velocity fields
%   velocityV - cell array {cam} of V velocity fields
%   masks - cell array {cam} of boolean masks (true = invalid)
%   window_type - 'tukey', 'hann', 'cosine', 'linear', 'flat' (default: 'tukey')
%   taper_param - window parameter (default: 0.5 for tukey alpha)
%
% Outputs:
%   worldX, worldY - unified grid coordinates
%   U_merged, V_merged - blended velocity fields

if nargin < 5, window_type = 'tukey'; end
if nargin < 6, taper_param = 0.5; end

nCams = numel(velocityU);
%% STEP 1: Create unified grid from all camera bounds
% Add optional worldGrid parameter
if nargin < 7 || isempty(worldGrid)
    x_min = inf; x_max = -inf;
    y_min = inf; y_max = -inf;

    for cam = 1:nCams
        x_min = min(x_min, min(windowCenters.x1_mm{cam}(:)));
        x_max = max(x_max, max(windowCenters.x1_mm{cam}(:)));
        y_min = min(y_min, min(windowCenters.x2_mm{cam}(:)));
        y_max = max(y_max, max(windowCenters.x2_mm{cam}(:)));
    end

    % Get grid spacing from first camera
    % dx = abs(mean(diff(windowCenters.x1_mm{1}(1,:))));
    % dy = abs(mean(diff(windowCenters.x2_mm{1}(:,1))));
    dx_all = zeros(1, nCams);
    dy_all = zeros(1, nCams);
    for cam = 1:nCams
        dx_all(cam) = abs(mean(diff(windowCenters.x1_mm{cam}(1,:))));
        dy_all(cam) = abs(mean(diff(windowCenters.x2_mm{cam}(:,1))));
    end
    dx = mean(dx_all);
    dy = mean(dy_all);


    % Create unified grid with UNIFORM SPACING using linspace (matches Python)
    nx = round((x_max - x_min) / dx) + 1;
    ny = round((y_max - y_min) / dy) + 1;
    x_vec = linspace(x_min, x_max, nx);
    y_vec = linspace(y_min, y_max, ny);  % Ascending for cartesian coords for image convention

    [worldX, worldY] = meshgrid(x_vec, y_vec);

    fprintf('Unified grid: %d x %d, X:[%.1f, %.1f], Y:[%.1f, %.1f]\n', ...
        size(worldX,2), size(worldX,1), x_min, x_max, y_min, y_max);

else
    % Use provided grid
    worldX = worldGrid.x1;
    worldY = worldGrid.x2;
    fprintf('Using provided world grid: %d x %d\n', size(worldX,2), size(worldX,1));
end


%% STEP 2: Interpolate all cameras to unified grid
camera_interp = cell(1, nCams);

for cam = 1:nCams
    fprintf('Interpolating camera %d...\n', cam);
    
    % Get camera's native grid
    x_cam = windowCenters.x1_mm{cam};
    y_cam = windowCenters.x2_mm{cam};
    u_cam = velocityU{cam};
    v_cam = velocityV{cam};
    
    % IMPROVEMENT 2: Store original NaN mask for accurate tracking
    original_nan_mask = isnan(u_cam) | isnan(v_cam);
    
    % Handle mask
    if ~isempty(masks) && ~isempty(masks{cam})
        mask_cam = masks{cam};
        original_nan_mask = original_nan_mask | mask_cam;
    end
    
    % Ensure x and y are ascending
    x_vec_cam = x_cam(1,:);
    y_vec_cam = y_cam(:,1);
    
    if any(diff(x_vec_cam) < 0)
        x_vec_cam = fliplr(x_vec_cam);
        u_cam = fliplr(u_cam);
        v_cam = fliplr(v_cam);
        original_nan_mask = fliplr(original_nan_mask);
    end
    
    if any(diff(y_vec_cam) < 0)
        y_vec_cam = flipud(y_vec_cam);
        u_cam = flipud(u_cam);
        v_cam = flipud(v_cam);
        original_nan_mask = flipud(original_nan_mask);
    end
    
    % Force uniform spacing using linspace (handles slight variations)
    nx_cam = length(x_vec_cam);
    ny_cam = length(y_vec_cam);
    x_vec_uniform = linspace(x_vec_cam(1), x_vec_cam(end), nx_cam);
    y_vec_uniform = linspace(y_vec_cam(1), y_vec_cam(end), ny_cam);
    
    % IMPROVEMENT 1: Better NaN handling - fill with 0 for interpolation
    u_cam_filled = u_cam;
    v_cam_filled = v_cam;
    u_cam_filled(isnan(u_cam)) = 0;
    v_cam_filled(isnan(v_cam)) = 0;
    
    % IMPROVEMENT 1: Create interpolants with CUBIC method
    F_u = griddedInterpolant({y_vec_uniform, x_vec_uniform}, u_cam_filled, 'cubic', 'none');
    F_v = griddedInterpolant({y_vec_uniform, x_vec_uniform}, v_cam_filled, 'cubic', 'none');
    
    % IMPROVEMENT 3: Accurate valid region tracking with nearest-neighbor mask interpolation
    valid_mask_data = double(~original_nan_mask);
    F_mask = griddedInterpolant({y_vec_uniform, x_vec_uniform}, valid_mask_data, 'nearest', 'none');
    
    % Interpolate to unified grid
    u_interp = F_u(worldY, worldX);
    v_interp = F_v(worldY, worldX);
    mask_interp = F_mask(worldY, worldX);
    
    % IMPROVEMENT 3: Valid mask includes original data validity
    valid = ~isnan(u_interp) & ~isnan(v_interp) & (mask_interp > 0.5);
    
    camera_interp{cam}.u = u_interp;
    camera_interp{cam}.v = v_interp;
    camera_interp{cam}.valid = valid;
end
%% DIAGNOSTIC: Check actual x-extent of valid regions per camera
fprintf('\n=== Valid region x-extents on merged grid ===\n');
for cam = 1:nCams
    valid_cols = any(camera_interp{cam}.valid, 1);  % any valid row in each column
    x_valid = worldX(1, valid_cols);
    fprintf('  Cam%d: x=[%.2f, %.2f] mm, %d valid columns\n', ...
        cam, min(x_valid), max(x_valid), sum(valid_cols));
end

%% STEP 3: Compute weights using distance transform - FIXED TAPER
camera_weights = cell(1, nCams);

% Fixed taper width in grid points - set to ~half your overlap width
% Your overlap looks ~200-300 px wide, so 150 is a good starting point
taper_width_pts = 25;

for cam = 1:nCams
    valid_mask = camera_interp{cam}.valid;
    
    % Distance from edge of valid region (in grid points)
    edge_dist = bwdist(~valid_mask);
    
    % Normalise to [0,1] using FIXED taper width, cap at 1
    norm_dist = min(edge_dist / taper_width_pts, 1.0);
    
    % Apply window function
    weight = apply_window(norm_dist, window_type, taper_param);
    
    % Zero weight outside valid region
    weight(~valid_mask) = 0;
    
    camera_weights{cam} = weight;
    
    fprintf('Camera %d: taper_width=%d px, weight range=[%.3f, %.3f]\n', ...
        cam, taper_width_pts, min(weight(valid_mask)), max(weight(valid_mask)));
end

%% STEP 4: Normalize weights to sum to 1
total_weight = zeros(size(worldX));
for cam = 1:nCams
    total_weight = total_weight + camera_weights{cam};
end

for cam = 1:nCams
    camera_weights{cam} = camera_weights{cam} ./ max(total_weight, eps);
end
%% DIAGNOSTIC: Inspect total_weight before and after normalisation
figure;
subplot(2,1,1);
imagesc(total_weight);
colorbar;
clim([0 2]);
title('total\_weight before normalisation (should be >1 in overlaps)');
xlabel('x grid point'); ylabel('y grid point');

subplot(2,1,2);
% Show a horizontal slice through the middle of the domain
mid_row = round(size(total_weight, 1) / 2);
plot(total_weight(mid_row, :));
yline(1.0, 'r--', 'sum = 1');
title(sprintf('Horizontal slice at row %d', mid_row));
xlabel('x grid point'); ylabel('total weight');
ylim([0 2.5]);

% Also print per-camera max_dist values for comparison
fprintf('\n=== Weight diagnostic ===\n');
for cam = 1:nCams
    fprintf('  Cam%d: max_dist in grid points = %.1f\n', cam, max(bwdist(~camera_interp{cam}.valid), [], 'all'));
end
%% STEP 6: Blend
U_merged = zeros(size(worldX));
V_merged = zeros(size(worldX));

for cam = 1:nCams
    u_clean = camera_interp{cam}.u;
    v_clean = camera_interp{cam}.v;
    
    u_clean(isnan(u_clean)) = 0;
    v_clean(isnan(v_clean)) = 0;
    
    U_merged = U_merged + camera_weights{cam} .* u_clean;
    V_merged = V_merged + camera_weights{cam} .* v_clean;
end

% Set to NaN where no camera has data
no_data = total_weight == 0;
U_merged(no_data) = NaN;
V_merged(no_data) = NaN;

fprintf('Blending complete.\n');
end


function weight = apply_window(norm_dist, window_type, param)
% APPLY_WINDOW Apply windowing function to normalized distance
%
% norm_dist: [0,1] where 0=edge, 1=center of valid region
% window_type: 'tukey', 'hann', 'cosine', 'linear', 'flat'
% param: window parameter (alpha for tukey, etc.)

switch lower(window_type)
    case 'tukey'
        % Tukey: flat in center, cosine taper at edges
        % param = alpha (0.5 = half cosine taper on each side)
        alpha = param;
        weight = ones(size(norm_dist));
        taper_region = norm_dist < (alpha / 2);
        weight(taper_region) = 0.5 * (1 - cos(2*pi * norm_dist(taper_region) / alpha));
        
    case 'hann'
        % Hann: pure cosine, no flat region
        % Maps [0,1] to [0,1] via raised cosine
        weight = 0.5 * (1 - cos(pi * norm_dist));
        
    case 'cosine'
        % Cosine taper: linear-ish in middle, smooth at edges
        weight = sin(pi/2 * norm_dist);
        
    case 'linear'
        % Linear ramp from edge to center
        weight = norm_dist;
        
    case 'flat'
        % No blending - sharp transition (for debugging)
        weight = ones(size(norm_dist));
        
    otherwise
        error('Unknown window type: %s', window_type);
end
end

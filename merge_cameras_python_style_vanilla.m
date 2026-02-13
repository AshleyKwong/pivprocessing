function [worldX, worldY, U_merged, V_merged] = merge_cameras_python_style_vanilla(...
    windowCenters, velocityU, velocityV, masks, window_type, taper_param)

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
x_min = inf; x_max = -inf;
y_min = inf; y_max = -inf;

for cam = 1:nCams
    x_min = min(x_min, min(windowCenters.x1_mm{cam}(:)));
    x_max = max(x_max, max(windowCenters.x1_mm{cam}(:)));
    y_min = min(y_min, min(windowCenters.x2_mm{cam}(:)));
    y_max = max(y_max, max(windowCenters.x2_mm{cam}(:)));
end

% Get grid spacing from first camera
dx = abs(mean(diff(windowCenters.x1_mm{1}(1,:))));
dy = abs(me(diff(windowCenters.x2_mm{1}(:,1))));

% Create unified grid
x_vec = x_min:dx:x_max;
y_vec = y_max:-dy:y_min; % Descending for image convention

[worldX, worldY] = meshgrid(x_vec, y_vec);

fprintf('Unified grid: %d x %d, X:[%.1f, %.1f], Y:[%.1f, %.1f]\n', ...
    size(worldX,2), size(worldX,1), x_min, x_max, y_min, y_max);

%% STEP 2: Interpolate all cameras to unified grid (CUBIC + Better NaN handling)
camera_interp = cell(1, nCams);

for cam = 1:nCams
    fprintf('Interpolating camera %d...\n', cam);
    
    % Get camera's native grid
    x_cam = windowCenters.x1_mm{cam};
    y_cam = windowCenters.x2_mm{cam};
    u_cam = velocityU{cam};
    v_cam = velocityV{cam};
    
    % Store original NaN mask for tracking
    original_nan_mask = isnan(u_cam) | isnan(v_cam);
    
    % Handle mask
    if ~isempty(masks) && ~isempty(masks{cam})
        mask_cam = masks{cam};
        original_nan_mask = original_nan_mask | mask_cam;
    end
    
    % Ensure x and y are ascending
    x_vec = x_cam(1,:);
    y_vec = y_cam(:,1);
    
    if any(diff(x_vec) < 0)
        x_vec = fliplr(x_vec);
        u_cam = fliplr(u_cam);
        v_cam = fliplr(v_cam);
        original_nan_mask = fliplr(original_nan_mask);
    end
    
    if any(diff(y_vec) < 0)
        y_vec = flipud(y_vec);
        u_cam = flipud(u_cam);
        v_cam = flipud(v_cam);
        original_nan_mask = flipud(original_nan_mask);
    end
    
    % === IMPROVED NaN HANDLING (Python approach) ===
    % Fill NaN with 0 for interpolation, then track valid regions separately
    u_cam_filled = u_cam;
    v_cam_filled = v_cam;
    u_cam_filled(isnan(u_cam)) = 0;
    v_cam_filled(isnan(v_cam)) = 0;
    
    % Create interpolants with CUBIC method
    F_u = griddedInterpolant({y_vec, x_vec}, u_cam_filled, 'cubic', 'none');
    F_v = griddedInterpolant({y_vec, x_vec}, v_cam_filled, 'cubic', 'none');
    
    % === ACCURATE VALID REGION TRACKING (Nearest-neighbor mask interpolation) ===
    % Create a binary mask interpolant to track original valid regions
    valid_mask_data = double(~original_nan_mask);
    F_mask = griddedInterpolant({y_vec, x_vec}, valid_mask_data, 'nearest', 'none');
    
    % Interpolate to unified grid
    u_interp = F_u(worldY, worldX);
    v_interp = F_v(worldY, worldX);
    mask_interp = F_mask(worldY, worldX);
    
    % Valid mask: within bounds AND originally valid data
    valid = ~isnan(u_interp) & ~isnan(v_interp) & (mask_interp > 0.5);
    
    % Store results
    camera_interp{cam}.u = u_interp;
    camera_interp{cam}.v = v_interp;
    camera_interp{cam}.valid = valid;
    
    fprintf('  Valid points: %d/%d (%.1f%%)\n', ...
        sum(valid(:)), numel(valid), 100*sum(valid(:))/numel(valid));
end

%% STEP 3: Compute weights using distance transform with Tukey/Hann window
camera_weights = cell(1, nCams);

for cam = 1:nCams
    valid_mask = camera_interp{cam}.valid;
    
    % Distance from edge of valid region (in pixels)
    edge_dist = bwdist(~valid_mask);
    
    % Normalize to [0, 1]
    max_dist = max(edge_dist(:));
    
    if max_dist > 0
        norm_dist = edge_dist / max_dist;
    else
        norm_dist = zeros(size(edge_dist));
    end
    
    % Apply window function (Tukey or Hann)
    weight = apply_window(norm_dist, window_type, taper_param);
    
    % Zero weight outside valid region
    weight(~valid_mask) = 0;
    
    camera_weights{cam} = weight;
    
    if any(valid_mask(:))
        fprintf('Camera %d: max_dist=%.1f px, weight range=[%.3f, %.3f]\n', ...
            cam, max_dist, min(weight(valid_mask)), max(weight(valid_mask)));
    else
        fprintf('Camera %d: No valid data\n', cam);
    end
end

%% STEP 4: Normalize weights to sum to 1
total_weight = zeros(size(worldX));
for cam = 1:nCams
    total_weight = total_weight + camera_weights{cam};
end

for cam = 1:nCams
    valid_total = total_weight > 0;
    camera_weights{cam}(valid_total) = camera_weights{cam}(valid_total) ./ total_weight(valid_total);
end

%% STEP 5: Blend
U_merged = zeros(size(worldX));
V_merged = zeros(size(worldX));

for cam = 1:nCams
    u_clean = camera_interp{cam}.u;
    v_clean = camera_interp{cam}.v;
    
    % Replace NaN with 0 for blending
    u_clean(isnan(u_clean)) = 0;
    v_clean(isnan(v_clean)) = 0;
    
    U_merged = U_merged + camera_weights{cam} .* u_clean;
    V_merged = V_merged + camera_weights{cam} .* v_clean;
end

% Set to NaN where no camera has data
no_data = total_weight == 0;
U_merged(no_data) = NaN;
V_merged(no_data) = NaN;

fprintf('Blending complete. Using %s window.\n', window_type);
end


function weight = apply_window(norm_dist, window_type, param)
% APPLY_WINDOW Apply windowing function to normalized distance
%
% norm_dist: [0,1] where 0=edge, 1=center of valid region
% window_type: 'tukey', 'hann', 'cosine', 'linear', 'flat'
% param: window parameter (alpha for tukey)

switch lower(window_type)
    case 'tukey'
        % Tukey: flat in center, cosine taper at edges
        % param = alpha (0.5 = half cosine taper on each side)
        alpha = param;
        weight = ones(size(norm_dist));
        taper_region = norm_dist < (alpha / 2);
        weight(taper_region) = 0.5 * (1 - cos(2*pi * norm_dist(taper_region) / alpha));
        
    case 'hann'
        % Hann: pure cosine, no flat region (matches Python)
        % Maps [0,1] to [0,1] via raised cosine
        weight = 0.5 * (1 - cos(pi * norm_dist));
        
    case 'cosine'
        % Cosine taper: smooth sine-based ramp
        weight = sin(pi/2 * norm_dist);
        
    case 'linear'
        % Linear ramp from edge to center
        weight = norm_dist;
        
    case 'flat'
        % No blending - sharp transition (for debugging)
        weight = ones(size(norm_dist));
        
    otherwise
        error('Unknown window type: %s. Options: tukey, hann, cosine, linear, flat', window_type);
end
end

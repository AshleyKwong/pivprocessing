function profiles = extractBoundaryLayerProfiles(x_locations, U_field, worldX, worldY, varargin)
% EXTRACTBOUNDARYLAYERPROFILES Extract velocity profiles at specified streamwise locations
%
% Syntax:
%   profiles = extractBoundaryLayerProfiles(x_locations, U_field, worldX, worldY)
%   profiles = extractBoundaryLayerProfiles(..., 'Name', Value)
%
% Inputs:
%   x_locations - Array of streamwise positions to extract profiles [mm]
%   U_field     - 2D velocity field (typically streamwise component U) [m/s]
%   worldX      - 2D grid of X-coordinates [mm]
%   worldY      - 2D grid of Y-coordinates [mm]
%
% Optional Name-Value Pairs:
%   'WallPosition'        - Y-coordinate of wall [mm] (default: 0)
%   'AveragingWidth'      - Width to average around each x_location [mm] (default: auto)
%   'Verbose'             - Print progress messages (default: true)
%   'ComputeThickness'    - Compute boundary layer thickness δ₉₉ (default: true)
%   'ComputeDisplacement' - Compute displacement thickness δ* (default: false)
%   'ComputeMomentum'     - Compute momentum thickness θ (default: false)
%   'ComputeShapeFactor'  - Compute shape factor H = δ*/θ (default: false)
%
% Outputs:
%   profiles - Structure array with fields:
%     .x_target   - Requested x-position [mm]
%     .x_actual   - Actual x-position (center of averaged region) [mm]
%     .x_range    - [x_min, x_max] of averaged region [mm]
%     .y          - Wall-normal distance (y - y_wall) [mm]
%     .U          - Streamwise velocity profile [m/s]
%     .num_cols   - Number of columns averaged
%     .U_max      - Maximum velocity in profile [m/s]
%     .delta_99   - Boundary layer thickness (99% U_max), interpolated [mm]
%     .delta_star - Displacement thickness [mm] (if ComputeDisplacement=true)
%     .theta      - Momentum thickness [mm] (if ComputeMomentum=true)
%     .H          - Shape factor H = δ*/θ (if ComputeShapeFactor=true)
%
% Example:
%   % Extract profiles with all integral thicknesses
%   profiles = extractBoundaryLayerProfiles(x_locs, U_hann_mean, worldX, worldY, ...
%       'WallPosition', 0, ...
%       'ComputeDisplacement', true, ...
%       'ComputeMomentum', true, ...
%       'ComputeShapeFactor', true);
%
% Author: Ashley Kwong
% Date: 02/17/2026

%% Parse inputs
p = inputParser;
addRequired(p, 'x_locations', @isnumeric);
addRequired(p, 'U_field', @isnumeric);
addRequired(p, 'worldX', @isnumeric);
addRequired(p, 'worldY', @isnumeric);
addParameter(p, 'WallPosition', 0, @isnumeric);
addParameter(p, 'AveragingWidth', [], @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'ComputeThickness', true, @islogical);
addParameter(p, 'ComputeDisplacement', false, @islogical);
addParameter(p, 'ComputeMomentum', false, @islogical);
addParameter(p, 'ComputeShapeFactor', false, @islogical);

parse(p, x_locations, U_field, worldX, worldY, varargin{:});

y_wall = p.Results.WallPosition;
averaging_width = p.Results.AveragingWidth;
VERBOSE = p.Results.Verbose;
COMPUTE_THICKNESS = p.Results.ComputeThickness;
COMPUTE_DISPLACEMENT = p.Results.ComputeDisplacement;
COMPUTE_MOMENTUM = p.Results.ComputeMomentum;
COMPUTE_SHAPE = p.Results.ComputeShapeFactor;

% Auto-enable dependencies
if COMPUTE_SHAPE
    COMPUTE_DISPLACEMENT = true;
    COMPUTE_MOMENTUM = true;
end

%% Extract coordinate vectors
x_vec = worldX(1, :);  % X-coordinates
y_vec = worldY(:, 1);  % Y-coordinates

% Auto-detect averaging width if not provided
if isempty(averaging_width)
    dx = abs(mean(diff(x_vec)));
    averaging_width = dx;  % Use grid spacing
end

if VERBOSE
    fprintf('\n=== Boundary Layer Profile Extraction ===\n');
    fprintf('Wall position: y = %.2f mm\n', y_wall);
    fprintf('Averaging width: ±%.2f mm around each x-location\n', averaging_width);
    fprintf('Number of profiles: %d\n', length(x_locations));
    fprintf('Field dimensions: [%.1f, %.1f] × [%.1f, %.1f] mm\n', ...
        min(x_vec), max(x_vec), min(y_vec), max(y_vec));
    fprintf('Calculations: δ₉₉=%d, δ*=%d, θ=%d, H=%d\n', ...
        COMPUTE_THICKNESS, COMPUTE_DISPLACEMENT, COMPUTE_MOMENTUM, COMPUTE_SHAPE);
end

%% Check input dimensions
if ~isequal(size(U_field), size(worldX), size(worldY))
    error('U_field, worldX, and worldY must have the same dimensions');
end

%% Extract profiles at each x-location
profiles = struct();

for i = 1:length(x_locations)
    x_target = x_locations(i);
    
    % Find columns within ±averaging_width of target x
    x_min = x_target - averaging_width;
    x_max = x_target + averaging_width;
    
    % Find indices of columns to average
    cols_to_avg = find(x_vec >= x_min & x_vec <= x_max);
    
    if isempty(cols_to_avg)
        if VERBOSE
            fprintf('  ⚠ Warning: No data found near x = %.2f mm\n', x_target);
        end
        
        % Store empty profile
        profiles(i) = createEmptyProfile(x_target, COMPUTE_THICKNESS, ...
            COMPUTE_DISPLACEMENT, COMPUTE_MOMENTUM, COMPUTE_SHAPE);
        continue;
    end
    
    % Actual x-position (center of averaged columns)
    x_actual = mean(x_vec(cols_to_avg));
    
    % Extract and average the columns
    U_profile = mean(U_field(:, cols_to_avg), 2, 'omitnan');
    
    % Wall-normal distance (y - y_wall)
    y_wall_normal = y_vec - y_wall;
    
    % Store basic profile data
    profiles(i).x_target = x_target;
    profiles(i).x_actual = x_actual;
    profiles(i).x_range = [x_vec(cols_to_avg(1)), x_vec(cols_to_avg(end))];
    profiles(i).num_cols = length(cols_to_avg);
    profiles(i).y = y_wall_normal;
    profiles(i).U = U_profile;
    
    % Calculate boundary layer statistics
    valid = ~isnan(U_profile);
    
    if ~any(valid)
        profiles(i).U_max = NaN;
        profiles(i).delta_99 = NaN;
        if COMPUTE_DISPLACEMENT, profiles(i).delta_star = NaN; end
        if COMPUTE_MOMENTUM, profiles(i).theta = NaN; end
        if COMPUTE_SHAPE, profiles(i).H = NaN; end
        continue;
    end
    
    % Find U_max
    U_max = max(U_profile(valid));
    profiles(i).U_max = U_max;
    
    %% Compute boundary layer thickness δ₉₉ (interpolated between 98% and 99%)
    if COMPUTE_THICKNESS
        delta_99 = computeBoundaryLayerThickness(U_profile(valid), ...
            y_wall_normal(valid), U_max);
        profiles(i).delta_99 = delta_99;
    else
        delta_99 = NaN;
        profiles(i).delta_99 = NaN;
    end
    
    %% Compute integral thicknesses (if requested)
    if COMPUTE_DISPLACEMENT || COMPUTE_MOMENTUM
        % Need valid δ₉₉ for integration
        if isnan(delta_99)
            if COMPUTE_DISPLACEMENT, profiles(i).delta_star = NaN; end
            if COMPUTE_MOMENTUM, profiles(i).theta = NaN; end
            if COMPUTE_SHAPE, profiles(i).H = NaN; end
        else
            % Prepare data for integration (from wall to BL edge)
            [y_int, U_int, Ue] = prepareIntegrationData(y_wall_normal(valid), ...
                U_profile(valid), delta_99, U_max);
            
            % Displacement thickness: δ* = ∫[0 to δ] (1 - U/Ue) dy
            if COMPUTE_DISPLACEMENT
                integrand_disp = (1 - U_int / Ue);
                delta_star = trapz(y_int, integrand_disp);
                profiles(i).delta_star = delta_star;
            end
            
            % Momentum thickness: θ = ∫[0 to δ] (U/Ue)(1 - U/Ue) dy
            if COMPUTE_MOMENTUM
                integrand_mom = (U_int / Ue) .* (1 - U_int / Ue);
                theta = trapz(y_int, integrand_mom);
                profiles(i).theta = theta;
            end
            
            % Shape factor: H = δ*/θ
            if COMPUTE_SHAPE
                if delta_star > 0 && theta > 0
                    profiles(i).H = delta_star / theta;
                else
                    profiles(i).H = NaN;
                end
            end
        end
    end
    
    %% Verbose output
    if VERBOSE
        output_str = sprintf('  Profile %d: x=%.1f mm, %d cols, Umax=%.2f m/s', ...
            i, x_actual, profiles(i).num_cols, U_max);
        
        if COMPUTE_THICKNESS
            output_str = [output_str, sprintf(', δ₉₉=%.2f mm', delta_99)];
        end
        if COMPUTE_DISPLACEMENT
            output_str = [output_str, sprintf(', δ*=%.2f mm', profiles(i).delta_star)];
        end
        if COMPUTE_MOMENTUM
            output_str = [output_str, sprintf(', θ=%.2f mm', profiles(i).theta)];
        end
        if COMPUTE_SHAPE
            output_str = [output_str, sprintf(', H=%.3f', profiles(i).H)];
        end
        
        fprintf('%s\n', output_str);
    end
end

if VERBOSE
    fprintf('✓ Extraction complete\n\n');
end

end

%% Helper Functions

function profile = createEmptyProfile(x_target, comp_thick, comp_disp, comp_mom, comp_shape)
% Create an empty profile structure
    profile.x_target = x_target;
    profile.x_actual = NaN;
    profile.x_range = [NaN, NaN];
    profile.y = [];
    profile.U = [];
    profile.num_cols = 0;
    profile.U_max = NaN;
    
    if comp_thick, profile.delta_99 = NaN; end
    if comp_disp, profile.delta_star = NaN; end
    if comp_mom, profile.theta = NaN; end
    if comp_shape, profile.H = NaN; end
end

function delta_99 = computeBoundaryLayerThickness(U_profile, y_profile, U_max)
% Compute boundary layer thickness by interpolating between 98% and 99% of U_max
% Inputs are assumed to be valid (no NaNs)
    
    U_normalized = U_profile / U_max;
    
    % Find points bracketing 99%
    idx_99 = find(U_normalized >= 0.99, 1, 'first');
    
    if isempty(idx_99)
        % No point reaches 99% - use maximum extent
        delta_99 = max(y_profile);
        return;
    end
    
    if idx_99 == 1
        % First point already >= 99%
        delta_99 = y_profile(1);
        return;
    end
    
    % Interpolate between point below and above 99%
    idx_98 = idx_99 - 1;
    
    U_98 = U_normalized(idx_98);
    U_99_actual = U_normalized(idx_99);
    y_98 = y_profile(idx_98);
    y_99 = y_profile(idx_99);
    
    % Linear interpolation to find y where U/U_max = 0.99
    delta_99 = interp1([U_98, U_99_actual], [y_98, y_99], 0.99, 'linear');
    
    % Handle extrapolation edge case
    if isnan(delta_99)
        delta_99 = y_99;
    end
end

function [y_int, U_int, Ue] = prepareIntegrationData(y_profile, U_profile, delta_99, U_max)
% Prepare velocity profile data for integration from y=0 to y=δ₉₉
% Ensures boundary layer edge point is included
    
    % Truncate to boundary layer thickness
    idx_within_BL = y_profile <= delta_99;
    
    y_int = y_profile(idx_within_BL);
    U_int = U_profile(idx_within_BL);
    
    % Add boundary layer edge point if not already present
    if ~any(y_int == delta_99)
        % Interpolate velocity at exact δ₉₉
        U_at_delta = interp1(y_profile, U_profile, delta_99, 'linear', 'extrap');
        
        % Append to arrays (sorted ascending in y)
        y_int = [y_int; delta_99];
        U_int = [U_int; U_at_delta];
        
        % Sort in case insertion breaks ordering
        [y_int, sort_idx] = sort(y_int, 'ascend');
        U_int = U_int(sort_idx);
    end
    
    % Edge velocity (at boundary layer edge)
    Ue = U_max;  % Could also use U_int(end) if interpolated
end

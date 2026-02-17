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
%   'WallPosition'    - Y-coordinate of wall [mm] (default: 0)
%   'AveragingWidth'  - Width to average around each x_location [mm] (default: auto from grid spacing)
%   'Verbose'         - Print progress messages (default: true)
%
% Outputs:
%   profiles - Structure array with fields:
%       .x_target       - Requested x-position [mm]
%       .x_actual       - Actual x-position (center of averaged region) [mm]
%       .x_range        - [x_min, x_max] of averaged region [mm]
%       .y              - Wall-normal distance (y - y_wall) [mm]
%       .U              - Streamwise velocity profile [m/s]
%       .num_cols       - Number of columns averaged
%       .delta_99       - Boundary layer thickness (99% U_max) [mm]
%       .U_max          - Maximum velocity in profile [m/s]
%
% Example:
%   % Extract profiles at 5 streamwise locations
%   x_locs = [1000, 1100, 1200, 1300, 1400];
%   profiles = extractBoundaryLayerProfiles(x_locs, U_hann_mean, worldX, worldY, ...
%                                           'WallPosition', 0, ...
%                                           'AveragingWidth', 5);
%
%   % Plot first profile
%   figure;
%   plot(profiles(1).U, profiles(1).y, 'o-');
%   xlabel('U [m/s]'); ylabel('y [mm]');
%
% Author: Ashley Kwong
% Date: 02/16/2026

%% Parse inputs
p = inputParser;
addRequired(p, 'x_locations', @isnumeric);
addRequired(p, 'U_field', @isnumeric);
addRequired(p, 'worldX', @isnumeric);
addRequired(p, 'worldY', @isnumeric);
addParameter(p, 'WallPosition', 0, @isnumeric);
addParameter(p, 'AveragingWidth', [], @isnumeric);
addParameter(p, 'Verbose', true, @islogical);

parse(p, x_locations, U_field, worldX, worldY, varargin{:});

y_wall = p.Results.WallPosition;
averaging_width = p.Results.AveragingWidth;
VERBOSE = p.Results.Verbose;

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
        profiles(i).x_target = x_target;
        profiles(i).x_actual = NaN;
        profiles(i).x_range = [NaN, NaN];
        profiles(i).y = [];
        profiles(i).U = [];
        profiles(i).num_cols = 0;
        profiles(i).delta_99 = NaN;
        profiles(i).U_max = NaN;
        continue;
    end
    
    % Actual x-position (center of averaged columns)
    x_actual = mean(x_vec(cols_to_avg));
    
    % Extract and average the columns
    U_profile = mean(U_field(:, cols_to_avg), 2, 'omitnan');
    
    % Wall-normal distance (y - y_wall)
    y_wall_normal = y_vec - y_wall;
    
    % Calculate boundary layer statistics
    valid = ~isnan(U_profile);
    if any(valid)
        U_max = max(U_profile(valid));
        
        % Find boundary layer edge (99% of U_max)
        idx_99 = find(U_profile(valid) >= 0.99*U_max, 1, 'first');
        if ~isempty(idx_99)
            y_valid = y_wall_normal(valid);
            delta_99 = y_valid(idx_99);
        else
            delta_99 = NaN;
        end
    else
        U_max = NaN;
        delta_99 = NaN;
    end
    
    % Store in structure
    profiles(i).x_target = x_target;
    profiles(i).x_actual = x_actual;
    profiles(i).x_range = [x_vec(cols_to_avg(1)), x_vec(cols_to_avg(end))];
    profiles(i).num_cols = length(cols_to_avg);
    profiles(i).y = y_wall_normal;
    profiles(i).U = U_profile;
    profiles(i).delta_99 = delta_99;
    profiles(i).U_max = U_max;
    
    if VERBOSE
        fprintf('  Profile %d: x=%.1f mm (range [%.1f, %.1f]), %d cols, δ99=%.1f mm, Umax=%.2f m/s\n', ...
            i, x_actual, profiles(i).x_range(1), profiles(i).x_range(2), ...
            profiles(i).num_cols, delta_99, U_max);
    end
end

if VERBOSE
    fprintf('✓ Extraction complete\n\n');

end

end

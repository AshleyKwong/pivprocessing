% =========================================================================
% find_required_yrefs_Tom.m
%
% For each x_ref target and each y/delta target, computes the exact y_ref
% (mm) required: y_ref = y/delta_target * delta99(x_target)
%
% Outputs the required y_ref values so you know what sweeps to run.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================
clear; clc;

%% ===================== USER INPUTS =====================
x_targets      = [7754.8]; % , 6850 7754.8
ydelta_targets = [0.02, 0.04, 0.05, 0.1, 0.5, 0.8];
blSweepFile    = ['C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\TOM PIV Processing\h500mm\blSweep_20260424_215811_Pos3.mat'];

%% ===================== LOAD delta99 =====================
B      = load(blSweepFile, 'blSweep');
bl_x   = B.blSweep.x_mm;
bl_d99 = B.blSweep.delta99_hybrid_mm;

% Remove non-finite points
valid_bl     = isfinite(bl_x) & isfinite(bl_d99);
bl_x_clean   = bl_x(valid_bl);
bl_d99_clean = bl_d99(valid_bl);

% Deduplicate — average delta99 at any repeated x locations
% (duplicates arise from mag_factor window overlap in bl_sweep)
[bl_x_clean, ~, ic] = unique(bl_x_clean, 'sorted');
bl_d99_clean        = accumarray(ic, bl_d99_clean, [], @mean);

fprintf('bl_x unique points: %d\n', numel(bl_x_clean));

delta_at_xtarget = interp1(bl_x_clean, bl_d99_clean, x_targets, 'linear', NaN);

%% ===================== COMPUTE REQUIRED y_refs =====================
nX  = numel(x_targets);
nYd = numel(ydelta_targets);

required_yref = nan(nYd, nX);
for iYd = 1:nYd
    for iX = 1:nX
        if isfinite(delta_at_xtarget(iX))
            required_yref(iYd, iX) = ydelta_targets(iYd) * delta_at_xtarget(iX);
        end
    end
end

%% ===================== PRINT TABLE =====================
fprintf('\n=== delta99 at each x_target ===\n');
for iX = 1:nX
    fprintf('  x = %6.1f mm  ->  delta99 = %.2f mm\n', ...
        x_targets(iX), delta_at_xtarget(iX));
end

fprintf('\n=== Required y_ref (mm) = y/delta * delta99 ===\n\n');
fprintf('%-12s', 'y/delta');
for iX = 1:nX
    fprintf('  x=%6.0fmm', x_targets(iX));
end
fprintf('\n%s\n', repmat('-', 1, 12 + nX*12));
for iYd = 1:nYd
    fprintf('%-12.3f', ydelta_targets(iYd));
    for iX = 1:nX
        if isfinite(required_yref(iYd, iX))
            fprintf('  %9.2f', required_yref(iYd, iX));
        else
            fprintf('  %9s', 'NaN');
        end
    end
    fprintf('\n');
end

%% ===================== UNIQUE y_refs TO RUN =====================
all_yrefs = unique(round(required_yref(isfinite(required_yref)), 2));
fprintf('\n=== Unique y_ref values to run (mm) ===\n');
fprintf('%.2f  ', all_yrefs);
fprintf('\n\nTotal: %d sweeps\n', numel(all_yrefs));
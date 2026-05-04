% save_position_mean_fields.m
% Ashley Kwong
%
% Loads Mean_Flow_Feilds.mat from each Pos_n subfolder, applies below-wall
% NaN masking (Y <= 0), and resaves with variable names compatible with
% bl_sweep.m and compute_two_point_covariance_v5.m:
%
%   worldX, worldY      — coordinate grids [mm]
%   U_hann_mean         — streamwise mean velocity [m/s]
%   V_hann_mean         — wall-normal mean velocity [m/s]
%   U_rms               — sqrt(uu_mean), streamwise RMS [m/s]
%   V_rms               — sqrt(vv_mean), wall-normal RMS [m/s]
%
% Output per position:
%   baseDir/Pos_n/mean_fields_blsweep.mat

clear; clc;

%% ── USER SETTINGS ────────────────────────────────────────────────────────
baseDir = 'G:\SW 500mm Mean Flow Fields\';
nPos    = 4;
%% ─────────────────────────────────────────────────────────────────────────

fprintf('\n=== save_position_mean_fields.m ===\n');

for p = 1:nPos
    posDir   = fullfile(baseDir, sprintf('Pos_%d', p));
    meanFile = fullfile(posDir, 'Mean_Flow_Feilds.mat');   % typo intentional

    fprintf('\nPos_%d: loading...', p);
    mf = load(meanFile, 'X', 'Y', 'U_mean', 'V_mean', 'uu_mean', 'vv_mean');

    % Convert coordinates metres -> mm
    worldX = mf.X * 1e3;
    worldY = mf.Y * 1e3;

    % Compute RMS from variance, clamp negatives from any processing artefacts
    U_hann_mean = double(mf.U_mean);
    V_hann_mean = double(mf.V_mean);
    U_rms       = sqrt(max(double(mf.uu_mean), 0));
    V_rms       = sqrt(max(double(mf.vv_mean), 0));

    % NaN-mask below-wall points (Y <= 0)
    below_wall = worldY <= 0;
    U_hann_mean(below_wall) = NaN;
    V_hann_mean(below_wall) = NaN;
    U_rms(below_wall)       = NaN;
    V_rms(below_wall)       = NaN;

    fprintf(' done.\n');
    fprintf('  Grid:        %d x %d\n',   size(worldX, 1), size(worldX, 2));
    fprintf('  X:           [%.2f, %.2f] mm\n', min(worldX(:)), max(worldX(:)));
    fprintf('  Y:           [%.2f, %.2f] mm\n', min(worldY(:)), max(worldY(:)));
    fprintf('  Below-wall:  %d / %d points NaN-masked (%.1f%%)\n', ...
        sum(below_wall(:)), numel(below_wall), 100*mean(below_wall(:)));
    fprintf('  U_hann_mean: [%.3f, %.3f] m/s\n', min(U_hann_mean(:), [], 'omitnan'), max(U_hann_mean(:), [], 'omitnan'));
    fprintf('  U_rms:       [%.3f, %.3f] m/s\n', min(U_rms(:), [], 'omitnan'),       max(U_rms(:), [], 'omitnan'));

    outFile = fullfile(posDir, 'mean_fields_blsweep.mat');
    save(outFile, 'worldX', 'worldY', ...
        'U_hann_mean', 'V_hann_mean', ...
        'U_rms',       'V_rms',       ...
        '-v7.3');
    fprintf('  Saved -> %s\n', outFile);
end

fprintf('\n=== Done ===\n');
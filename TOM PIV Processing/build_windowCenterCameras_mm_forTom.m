% build_windowCenterCameras_mm_forTom.m
% Ashley Kwong
%
% Builds windowCenterCameras_mm.mat from the X, Y coordinate grids
% in each Pos_n/Mean_Flow_Feilds.mat, converting metres -> mm.
%
% Output: baseDir/windowCenterCameras_mm.mat
%   windowCenterCameras_mm.x1_mm{p} — streamwise coordinate grid [mm]  (Ny x Nx)
%   windowCenterCameras_mm.x2_mm{p} — wall-normal coordinate grid [mm] (Ny x Nx)

clear; clc;

%% ── USER SETTINGS ────────────────────────────────────────────────────────
baseDir = 'G:\SW 500mm Mean Flow Fields\';
nPos    = 4;
%% ─────────────────────────────────────────────────────────────────────────
 
fprintf('\n=== build_windowCenterCameras_mm.m ===\n');
 
for p = 1:nPos
    posDir   = fullfile(baseDir, sprintf('Pos_%d', p));
    meanFile = fullfile(posDir, 'Mean_Flow_Feilds.mat');   % typo intentional
 
    mf = load(meanFile, 'X', 'Y');
  
    % {1 x 1} cell — single camera, matches {nFrames x 1} fluctuations
    windowCenterCameras_mm.x1_mm = {mf.X * 1e3};   % metres -> mm
    windowCenterCameras_mm.x2_mm = {mf.Y * 1e3};   % metres -> mm
 
    fprintf('Pos_%d  X:[%.2f, %.2f] mm   Y:[%.2f, %.2f] mm   grid:[%d x %d]\n', p, ...
        min(mf.X(:))*1e3, max(mf.X(:))*1e3, ...
        min(mf.Y(:))*1e3, max(mf.Y(:))*1e3, ...
        size(mf.X, 1), size(mf.X, 2));
 
    outFile = fullfile(posDir, 'windowCenterCameras_mm.mat');
    save(outFile, 'windowCenterCameras_mm', '-v7.3');
    fprintf('  Saved -> %s\n', outFile);
end
 
fprintf('\n=== Done ===\n');
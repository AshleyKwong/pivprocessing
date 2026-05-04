%% =========================================================
%  Two-point correlation figure around a central HWA profile
%  Ashley - template script (updated layout & styling)
% ==========================================================
% clear; clc; close all;

%% ---------------- USER INPUTS ----------------
folder = 'C:\Users\ak1u24\Downloads\xref_HWA1_1150mm';

% y_ref in mm (same order as files after sorting & manualFileOrder)
yref_mm = [ 3.5
            3.9
            5.0
            6.0
           18.1
           25.7
           89.1 ];

% Convert y_ref to y+ using your calibration
yref_yplus = yref_mm./1000 .* utau_ofi_matchedUinflocal_piv ./ nu_air;

% x-window for the correlation panels
xwin = [1000 1200];

% contour levels / colour scaling
nLevels = 21;
useSymmetricCaxis = true;
manualCLim = [0 1]    ;         % set e.g. [-0.3 0.3] to override automatic

% HWA profile data (already in y+, U+)
yplus_prof = HWA_y * HWA_utau / HWA_nu;
Uplus_prof = HWA_U / HWA_utau;

% Optional: manually specify order of files if alphabetical order is wrong
manualFileOrder = [2 3 5 1 4 6];

%% ---------------- LOAD FILES ----------------
files = dir(fullfile(folder, '*.mat'));
assert(~isempty(files), 'No .mat files found in folder.');

[~, idx] = sort({files.name});
files = files(idx);

if ~isempty(manualFileOrder)
    files = files(manualFileOrder);
end

nFiles = numel(files);

assert(numel(yref_mm)    == nFiles, 'Length of yref_mm must match number of .mat files.');
assert(numel(yref_yplus) == nFiles, 'Length of yref_yplus must match number of .mat files.');

disp('Files being used:')
for i = 1:nFiles
    fprintf('%d : %s\n', i, files(i).name);
end

%% ---------------- INSPECT ONE FILE ----------------
S0 = load(fullfile(folder, files(1).name));
disp('Variables in first .mat file:')
disp(fieldnames(S0))

gridData = load('C:\Users\ak1u24\Downloads\grid.mat');

% Assumed structure:
%   S.R_s              -> correlation field (ny x nx)
%   gridData.worldX_merged -> x grid (either [1 x nx] or ny x nx)
%   gridData.worldY_merged -> y grid (either [ny x 1] or ny x nx, in mm)

getCorrField = @(S) S.R_s;

%% ---------------- PRELOAD ALL PANELS ----------------
Rall = cell(nFiles,1);
Xall = cell(nFiles,1);
Yall = cell(nFiles,1);

globalMin = +inf;
globalMax = -inf;

for i = 1:nFiles
    S = load(fullfile(folder, files(i).name));
    R = getCorrField(S);

    X = gridData.worldX_merged;
    Y = gridData.worldY_merged;

    Rall{i} = R;
    Xall{i} = X;
    Yall{i} = Y;

    globalMin = min(globalMin, min(R(:), [], 'omitnan'));
    globalMax = max(globalMax, max(R(:), [], 'omitnan'));
end

if isempty(manualCLim)
    if useSymmetricCaxis
        cmax = max(abs([globalMin globalMax]));
        cLimVals = [-cmax cmax];
    else
        cLimVals = [globalMin globalMax];
    end
else
    cLimVals = manualCLim;
end

%% ---------------- CREATE FIGURE ----------------
fig = figure('Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [50 50 1700 950]);

% Central main axis
axMain = axes('Parent', fig, ...
    'Position', [0.25 0.25 0.40 0.45]);
hold(axMain, 'on')
box(axMain, 'on')
grid(axMain, 'on')

set(axMain, 'XScale', 'log', 'YScale', 'linear');  % enforce semilogx

assert(~isempty(yplus_prof) && ~isempty(Uplus_prof), ...
    'Please supply yplus_prof and Uplus_prof before plotting.');

semilogx(axMain, yplus_prof, Uplus_prof, 'k.-', ...
    'LineWidth', 1.2, 'MarkerSize', 12);

xlabel(axMain, 'y^+')
ylabel(axMain, 'U^+')
title(axMain, 'HWA profile with linked two-point correlation panels')

% Mark the y+ reference locations on the profile
Uref = interp1(yplus_prof, Uplus_prof, yref_yplus, 'linear', 'extrap');

plot(axMain, yref_yplus, Uref, 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 6)

for i = 1:nFiles
    text(axMain, yref_yplus(i), Uref(i), sprintf('  %g', yref_yplus(i)), ...
        'Color', 'r', 'FontSize', 9, 'Clipping', 'on');
end

%% ---------------- PANEL POSITIONS ----------------
% New layout: 3 panels top row, 3 panels bottom row
% More vertical space from main axes (top row high, bottom row low)
panelPos = [
    0.08 0.78 0.22 0.17   % top-left
    0.38 0.78 0.22 0.17   % top-middle
    0.68 0.78 0.22 0.17   % top-right
    0.08 0.001 0.22 0.17   % bottom-left
    0.38 0.001 0.22 0.17   % bottom-middle
    0.68 0.001 0.22 0.17   % bottom-right
];

assert(size(panelPos,1) >= nFiles, ...
    'Not enough panel positions defined for number of files.');

axCorr = gobjects(nFiles,1);

%% ---------------- PLOT CORRELATION PANELS ----------------
for i = 1:nFiles
    axCorr(i) = axes('Parent', fig, 'Position', panelPos(i,:));
    hold(axCorr(i), 'on')
    box(axCorr(i), 'on')

    X = Xall{i};
    Y = Yall{i};
    R = Rall{i};
    % Filled field with smooth colormap
    contourf(axCorr(i), X, Y, R, nLevels, 'LineColor', 'none');
    hold(axCorr(i), 'on')

    % Overlay thin black contour lines (fewer levels so it's not busy)
    contour(axCorr(i), X, Y, R, 10, 'LineColor', 'k', 'LineWidth', 0.5);

    
    % Contour plot
    % imagesc(axCorr(i), X(1,:), Y(:,1), R);  % assuming X,Y are regular grids
    % set(axCorr(i),'YDir','normal');        % fix inversion from imagesc

    % caxis(axCorr(i), cLimVals);
    % colormap(axCorr(i), redblue(10));

    % hold(axCorr(i),'on')
    contour(axCorr(i), X, Y, R, 10, 'k', 'LineWidth', 0.5);
    % Shared colormap & colour limits
    % Ensure redblue(10) exists; otherwise use parula(10) or similar
    try
        colormap(axCorr(i), redblue(10));
    catch
        warning('redblue(10) not found on path; using parula(10) instead.');
        colormap(axCorr(i), parula(10));
    end
    caxis(axCorr(i), cLimVals);

    % Axes limits
    xlim(axCorr(i), xwin);
    ylim(axCorr(i), [0 120]);

    % Title with y_ref and y+
    title(axCorr(i), sprintf('y_{ref}=%.1f mm, y^+=%.0f', ...
        yref_mm(i), yref_yplus(i)), 'FontSize', 9);

    set(axCorr(i), 'FontSize', 8);

    % Remove all x tick labels & x label
    axCorr(i).XTickLabel = [];
    axCorr(i).XLabel.String = '';

    % y labels: show ticks in mm only (use existing Y tick locations)
    ylabel(axCorr(i), 'y [mm]');
end

linkaxes(axCorr, 'x');

% Shared colourbar (right side)
cb = colorbar(axCorr(end), 'Position', [0.92 0.22 0.012 0.30]);
cb.Label.String = 'R_s';

%% ---------------- CONNECT MAIN PROFILE TO PANELS ----------------
for i = 1:nFiles
    [xf, yf] = data2norm(axMain, yref_yplus(i), Uref(i));

    p = panelPos(i,:);
    xPanel = p(1) + 0.5*p(3);
    yPanel = p(2) + 0.5*p(4);

    annotation(fig, 'line', [xf xPanel], [yf yPanel], ...
        'Color', [0 0.45 0.45], 'LineWidth', 1.2);
end

%% ---------------- EXPORT ----------------
% exportgraphics(fig, 'two_point_correlation_layout.png', 'Resolution', 300);

%% =========================================================
% Helper function: convert data coords -> normalized figure coords
% Works with semilogx main axis
%% =========================================================
function [xfig, yfig] = data2norm(ax, xdata, ydata)

    axUnits = ax.Units;
    ax.Units = 'normalized';
    axPos = ax.Position;
    ax.Units = axUnits;

    xlim_ = xlim(ax);
    ylim_ = ylim(ax);

    if strcmp(ax.XScale, 'log')
        xnorm = (log10(xdata) - log10(xlim_(1))) / (log10(xlim_(2)) - log10(xlim_(1)));
    else
        xnorm = (xdata - xlim_(1)) / (xlim_(2) - xlim_(1));
    end

    if strcmp(ax.YScale, 'log')
        ynorm = (log10(ydata) - log10(ylim_(1))) / (log10(ylim_(2)) - log10(ylim_(1)));
    else
        ynorm = (ydata - ylim_(1)) / (ylim_(2) - ylim_(1));
    end

    xfig = axPos(1) + axPos(3) * xnorm;
    yfig = axPos(2) + axPos(4) * ynorm;
end
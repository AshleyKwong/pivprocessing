% inspect_position_extents.m
% Ashley Kwong
%
% Loads each Pos_n mean field, plots U_mean and U_rms side by side
% WITHOUT merging, and identifies which Pos_n subfolder each user-defined
% x location belongs to. xline markers are added to all plots.

clear; clc; close all;

%% ── USER SETTINGS ────────────────────────────────────────────────────────
baseDir = 'G:\SW 500mm Mean Flow Fields';
nPos    = 4;

% x locations of interest [mm] — will be identified to a position
x_query_mm = [5853, 6850 7754.8];
% x_labels   = {'Inlet', 'FPG max', 'Crossover', 'APG max', 'TE relative', 'ZPG recovery'};
x_labels   = {'Inlet', 'Crossover','TE relative', };

%% ─────────────────────────────────────────────────────────────────────────

% Colour per position for xlines
pos_colours = lines(nPos);

%% STEP 1: Load all positions
X_mm     = cell(1, nPos);
Y_mm     = cell(1, nPos);
U_mean   = cell(1, nPos);
U_rms    = cell(1, nPos);
x_extents = nan(nPos, 2);   % [x_min, x_max] per position

for p = 1:nPos
    posDir   = fullfile(baseDir, sprintf('Pos_%d', p));
    meanFile = fullfile(posDir, 'Mean_Flow_Feilds.mat');

    mf = load(meanFile, 'X', 'Y', 'U_mean', 'uu_mean');

    X_mm{p}   = mf.X * 1e3;
    Y_mm{p}   = mf.Y * 1e3;
    U_mean{p} = mf.U_mean;
    U_rms{p}  = sqrt(max(mf.uu_mean, 0));

    x_extents(p, :) = [min(X_mm{p}(:)), max(X_mm{p}(:))];

    fprintf('Pos_%d  X:[%.2f, %.2f] mm   Y:[%.2f, %.2f] mm\n', p, ...
        x_extents(p,1), x_extents(p,2), ...
        min(Y_mm{p}(:)), max(Y_mm{p}(:)));
end

%% STEP 2: Identify which position each query x belongs to
fprintf('\n=== Query x locations ===\n');
x_owner = nan(size(x_query_mm));   % which position owns each query x

for q = 1:numel(x_query_mm)
    xq = x_query_mm(q);
    candidates = find(xq >= x_extents(:,1) & xq <= x_extents(:,2));
    if isempty(candidates)
        fprintf('  x = %g mm  [%s]  → NOT in any position domain\n', xq, x_labels{q});
    else
        x_owner(q) = candidates(1);   % take first if overlap
        fprintf('  x = %g mm  [%s]  → Pos_%d\n', xq, x_labels{q}, candidates(1));
        if numel(candidates) > 1
            fprintf('    (also in overlap with: %s)\n', ...
                strjoin(arrayfun(@(c) sprintf('Pos_%d',c), candidates(2:end), ...
                'UniformOutput', false), ', '));
        end
    end
end

%% STEP 3: Per-position xline colours — assign query points to their position colour
% xline colour matches the position the point belongs to
xline_colours = zeros(numel(x_query_mm), 3);
for q = 1:numel(x_query_mm)
    if ~isnan(x_owner(q))
        xline_colours(q, :) = pos_colours(x_owner(q), :);
    else
        xline_colours(q, :) = [0.5 0.5 0.5];   % grey if unowned
    end
end

%% STEP 4: Plot U_mean per position
figure('Name', 'U_mean per position', 'Position', [50 50 1400 900]);
t = tiledlayout(nPos, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'U_{mean} per position — unmerged', 'FontSize', 13);

for p = 1:nPos
    nexttile;
    pcolor(X_mm{p}, Y_mm{p}, U_mean{p}); shading interp;
    cb = colorbar; cb.Label.String = 'U [m/s]';
    clim([0 40]);
    hold on;
    % Position boundary lines
    xline(x_extents(p,1), '--', 'Color', pos_colours(p,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    xline(x_extents(p,2), '--', 'Color', pos_colours(p,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    % Query x locations
    for q = 1:numel(x_query_mm)
        xq = x_query_mm(q);
        if xq >= x_extents(p,1) && xq <= x_extents(p,2)
            xline(xq, '-', 'Color', 'k', 'LineWidth', 1.5, ...
                'Label', x_labels{q}, 'LabelVerticalAlignment', 'bottom', ...
                'FontSize', 8);
        end
    end
    ylabel('Y [mm]');
    title(sprintf('Pos_%d  |  X: [%.1f, %.1f] mm', p, x_extents(p,1), x_extents(p,2)));
    axis tight; box on;
end
xlabel('X [mm]');

%% STEP 5: Plot U_rms per position
figure('Name', 'U_rms per position', 'Position', [100 100 1400 900]);
t2 = tiledlayout(nPos, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t2, 'U_{rms} per position — unmerged', 'FontSize', 13);

for p = 1:nPos
    nexttile;
    pcolor(X_mm{p}, Y_mm{p}, U_rms{p}); shading interp;
    cb = colorbar; cb.Label.String = 'U_{rms} [m/s]';
    clim([0 3]);
    hold on;
    xline(x_extents(p,1), '--', 'Color', pos_colours(p,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    xline(x_extents(p,2), '--', 'Color', pos_colours(p,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    for q = 1:numel(x_query_mm)
        xq = x_query_mm(q);
        if xq >= x_extents(p,1) && xq <= x_extents(p,2)
            xline(xq, '-', 'Color', 'k', 'LineWidth', 1.5, ...
                'Label', x_labels{q}, 'LabelVerticalAlignment', 'bottom', ...
                'FontSize', 8);
        end
    end
    ylabel('Y [mm]');
    title(sprintf('Pos_%d  |  X: [%.1f, %.1f] mm', p, x_extents(p,1), x_extents(p,2)));
    axis tight; box on;
end
xlabel('X [mm]');
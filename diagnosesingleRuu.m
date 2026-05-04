% =========================================================================
% diagnose_single_Ruu.m
%
% Loads a single R_uu correlation map and visualises the threshold masks,
% bounding boxes and fitted ellipses at all three rho levels side by side.
%
% Inclination angle method follows Marusic & Heuer (2007) / Volino et al.
% (2007): least-squares line through the upstream-most and downstream-most
% contour tip points, pooled across all rho levels.
%
% Semi-axes come directly from the bounding box:
%   a = (Lxu + Lxd) / 2   (semi-major, streamwise)
%   b = (Lyt + Lyb) / 2   (semi-minor, wall-normal)
%
% No external helper functions required.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc;

%% ===================== USER INPUTS ======================================

matFile  = 'G:\SW 400mm Processed Data Single Snapshots\-8\Pos_3\two_point_covariance_20260423_113156\R_uu_xref7320\R_uu_xref7320_yref003p20.mat';
gridFile = 'G:\SW 400mm Processed Data Single Snapshots\-8\Pos_3\two_point_covariance_20260423_113156\grid.mat';

rho_levels     = [1/exp(1), 0.6, 0.7];
dx_back        = 300;
dx_fwd         = 1000;
clip_warn_frac = 0.1;

%% ===================== LOAD =============================================

G      = load(gridFile, 'worldX_merged', 'worldY_merged');
worldX = double(G.worldX_merged);
worldY = double(G.worldY_merged);

S  = load(matFile, 'R_s', 'xr', 'yr');
R  = double(S.R_s);
xr = double(S.xr);
yr = double(S.yr);

fprintf('Loaded: %s\n', matFile);
fprintf('Reference point: xr = %.2f mm, yr = %.2f mm\n', xr, yr);

%% ===================== BUILD DELTA GRIDS + SEARCH BOX ==================

dX = worldX - xr;
dY = worldY - yr;

x_min_domain = min(worldX(:));   x_max_domain = max(worldX(:));
y_min_domain = min(worldY(:));   y_max_domain = max(worldY(:));

x_lo = max(xr - dx_back, x_min_domain);
x_hi = min(xr + dx_fwd,  x_max_domain);

clipped = (xr - x_lo) < (1 - clip_warn_frac)*dx_back || ...
          (x_hi - xr) < (1 - clip_warn_frac)*dx_fwd;
if clipped
    fprintf('WARNING: search box clipped by domain boundary\n');
    fprintf('  Requested: [-%.0f, +%.0f] mm\n', dx_back, dx_fwd);
    fprintf('  Actual:    [-%.0f, +%.0f] mm\n', xr-x_lo, x_hi-xr);
end

boxMask         = (dX >= -(xr - x_lo)) & (dX <= (x_hi - xr));
R_box           = R;
R_box(~boxMask) = NaN;
R_box(R_box >  1.0) =  1.0;
R_box(R_box < -1.0) =  NaN;

dX_box = dX;  dX_box(~boxMask) = NaN;
dY_box = dY;  dY_box(~boxMask) = NaN;

dist2           = dX.^2 + dY.^2;
dist2(~boxMask) = Inf;
[~, i0]         = min(dist2(:));

%% ===================== FIGURE 1: RAW CORRELATION MAP ===================

figure('Color','w','Position',[50 50 900 400],'Name','Raw R_uu map');
imagesc(dX(1,:), dY(:,1), R_box);
axis xy; axis tight; clim([-0.3 1]);
colormap(redblue(10)); colorbar; hold on;
plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
xlo_rel = -(xr - x_lo);  xhi_rel = (x_hi - xr);
yext    = [min(dY(:)) max(dY(:))];
plot([xlo_rel xlo_rel], yext, 'w--', 'LineWidth', 1.2);
plot([xhi_rel xhi_rel], yext, 'w--', 'LineWidth', 1.2);
xlabel('\Deltax  (mm)');  ylabel('\Deltay  (mm)');
title(sprintf('R_{uu}  |  x_{ref} = %.1f mm,  y_{ref} = %.1f mm', xr, yr));
set(gca, 'YDir', 'normal');

%% ===================== PER-LEVEL PROCESSING ============================

nLevels = numel(rho_levels);

cmapRho = [0.85 0.33 0.10;   % rho = 0.2  (orange-red)
           0.47 0.67 0.19;   % rho = 1/e  (green)
           0.13 0.47 0.71];  % rho = 0.5  (blue)

rhoLabels = {sprintf('\\rho = %.2f',    rho_levels(1)), ...
             sprintf('\\rho = 1/e = %.2f', rho_levels(2)), ...
             sprintf('\\rho = %.2f',    rho_levels(3))};

dx_vec     = dX(1,:);
dy_vec     = dY(:,1);
dx_spacing = abs(dx_vec(2) - dx_vec(1));
dy_spacing = abs(dy_vec(2) - dy_vec(1));

bbox     = nan(nLevels, 4);   % [Lxu, Lxd, Lyt, Lyb]
tip_pts  = nan(nLevels, 4);   % [x_up, y_up, x_dn, y_dn]  (upstream/downstream tips)
nReg     = nan(nLevels, 1);
bndStore = cell(nLevels, 1);
touchFlags(nLevels) = struct('up',false,'dn',false,'top',false,'bot',false);

for iL = 1:nLevels

    mask = R_box >= rho_levels(iL);
    if sum(mask(:)) < 5, continue; end

    CC = bwconncomp(mask);
    nReg(iL) = CC.NumObjects;

    % Identify origin-connected (or largest) region
    regionIdx = 0;
    for r = 1:CC.NumObjects
        if any(CC.PixelIdxList{r} == i0), regionIdx = r; break; end
    end
    if regionIdx == 0
        [~, regionIdx] = max(cellfun(@numel, CC.PixelIdxList));
    end

    mask_origin = false(size(mask));
    mask_origin(CC.PixelIdxList{regionIdx}) = true;
    bndStore{iL} = bwboundaries(mask_origin);

    dX_pts = dX_box(mask_origin);
    dY_pts = dY_box(mask_origin);
    if numel(dX_pts) < 5, continue; end

    % Bounding box extents
    Lxu = abs(min(dX_pts));
    Lxd = max(dX_pts);
    Lyt = max(dY_pts);
    Lyb = abs(min(dY_pts));
    bbox(iL,:) = [Lxu, Lxd, Lyt, Lyb];

    % Upstream tip = mean of all points within one grid cell of the leftmost x
    % Downstream tip = mean of all points within one grid cell of the rightmost x
    up_mask   = dX_pts <= (min(dX_pts) + dx_spacing);
    dn_mask   = dX_pts >= (max(dX_pts) - dx_spacing);
    tip_pts(iL,:) = [mean(dX_pts(up_mask)),  mean(dY_pts(up_mask)), ...
                     mean(dX_pts(dn_mask)),   mean(dY_pts(dn_mask))];

    % Edge-touch flags
    touchFlags(iL).up  = Lxu >= (xr - x_lo)       - dx_spacing;
    touchFlags(iL).dn  = Lxd >= (x_hi - xr)        - dx_spacing;
    touchFlags(iL).top = Lyt >= (y_max_domain - yr) - dy_spacing;
    touchFlags(iL).bot = Lyb >= (yr - y_min_domain) - dy_spacing;
end

%% ===================== GLOBAL INCLINATION ANGLE (LS FIT) ===============
% Pool upstream + downstream tips across all valid rho levels.
% Fit y = p(1)*x + p(2) → theta = atand(p(1)).
% This follows Marusic & Heuer (2007) / Volino et al. (2007).

valid = find(~any(isnan(tip_pts), 2));
x_tips = [tip_pts(valid, 1); tip_pts(valid, 3)];
y_tips = [tip_pts(valid, 2); tip_pts(valid, 4)];

if numel(x_tips) >= 2
    p_fit     = polyfit(x_tips, y_tips, 1);
    theta_deg = atand(p_fit(1));
else
    p_fit     = [0 0];
    theta_deg = NaN;
    warning('Not enough valid tip points to fit inclination angle.');
end

fprintf('\nGlobal inclination angle (LS, pooled tips): %.2f deg\n', theta_deg);

%% ===================== FIGURE 2: MASKS + BBOX + ELLIPSE ================

figure('Color','w','Position',[50 500 1400 420], ...
    'Name','Threshold masks, bounding boxes and ellipses');

for iL = 1:nLevels

    ax = subplot(1, nLevels, iL);
    hold on; box on;

    imagesc(ax, dx_vec, dy_vec, R_box);
    colormap(ax, gray); clim(ax, [0 1]);
    axis(ax,'xy'); axis(ax,'tight');

    if all(isnan(bbox(iL,:)))
        title(ax, sprintf('%s\n(insufficient points)', rhoLabels{iL}));
        set(ax,'YDir','normal'); continue
    end

    % All regions – grey outlines
    mask_all = R_box >= rho_levels(iL);
    CC_all   = bwconncomp(mask_all);
    for r = 1:CC_all.NumObjects
        tmp = false(size(mask_all));
        tmp(CC_all.PixelIdxList{r}) = true;
        Braw = bwboundaries(tmp);
        for b = 1:numel(Braw)
            bnd = Braw{b};
            plot(ax, dx_vec(bnd(:,2)), dy_vec(bnd(:,1)), '-', ...
                'Color',[0.7 0.7 0.7],'LineWidth',0.8);
        end
    end

    % Origin-region boundary
    for b = 1:numel(bndStore{iL})
        bnd = bndStore{iL}{b};
        plot(ax, dx_vec(bnd(:,2)), dy_vec(bnd(:,1)), '-', ...
            'Color', cmapRho(iL,:), 'LineWidth', 2);
    end

    % Bounding box
    Lxu = bbox(iL,1);  Lxd = bbox(iL,2);
    Lyt = bbox(iL,3);  Lyb = bbox(iL,4);
    Lx  = Lxu + Lxd;  Ly  = Lyt + Lyb;

    plot(ax, [-Lxu, Lxd, Lxd, -Lxu, -Lxu], [-Lyb,-Lyb, Lyt, Lyt,-Lyb], ...
        '--', 'Color', cmapRho(iL,:), 'LineWidth', 1.4);
    plot(ax, [-Lxu 0],   [0 0],    '-', 'Color',cmapRho(iL,:),'LineWidth',1.0);
    plot(ax, [0  Lxd],   [0 0],    '-', 'Color',cmapRho(iL,:),'LineWidth',1.0);
    plot(ax, [0    0],   [0 Lyt],  '-', 'Color',cmapRho(iL,:),'LineWidth',1.0);
    plot(ax, [0    0],   [-Lyb 0], '-', 'Color',cmapRho(iL,:),'LineWidth',1.0);

    % Ellipse: semi-axes from bbox, inclination from global LS fit
    a_ell = Lx / 2;
    b_ell = Ly / 2;
    cx    = (Lxd - Lxu) / 2;   % centre offset due to streamwise asymmetry
    cy    = (Lyt - Lyb) / 2;   % centre offset due to wall-normal asymmetry

    t         = linspace(0, 2*pi, 361);
    R_rot     = [cosd(theta_deg), -sind(theta_deg);
                 sind(theta_deg),  cosd(theta_deg)];
    pts_world = R_rot * [a_ell.*cos(t); b_ell.*sin(t)];

    xe = cx + pts_world(1,:);
    ye = cy + pts_world(2,:);

    fill(ax, xe, ye, cmapRho(iL,:), 'FaceAlpha', 0.12, 'EdgeColor','none');
    plot(ax, xe, ye, '-', 'Color', cmapRho(iL,:), 'LineWidth', 2.2);

    % Major axis through ellipse centre
    plot(ax, cx + a_ell*[-cosd(theta_deg), cosd(theta_deg)], ...
             cy + a_ell*[-sind(theta_deg), sind(theta_deg)], ...
        ':', 'Color', cmapRho(iL,:), 'LineWidth', 1.5);

    % Tip markers used in LS fit
    plot(ax, tip_pts(iL,1), tip_pts(iL,2), 'v', ...
        'Color',cmapRho(iL,:),'MarkerFaceColor',cmapRho(iL,:),'MarkerSize',7);
    plot(ax, tip_pts(iL,3), tip_pts(iL,4), '^', ...
        'Color',cmapRho(iL,:),'MarkerFaceColor',cmapRho(iL,:),'MarkerSize',7);

    % Annotation
    flagStr = '';
    if touchFlags(iL).up,  flagStr = [flagStr, '⚠ upstream '];   end
    if touchFlags(iL).dn,  flagStr = [flagStr, '⚠ downstream ']; end
    if touchFlags(iL).top, flagStr = [flagStr, '⚠ top '];        end
    if touchFlags(iL).bot, flagStr = [flagStr, '⚠ bottom '];     end
    if isempty(flagStr), flagStr = 'no clipping'; end

    text(ax, 0.05, 0.97, ...
        sprintf(['L_x^u=%.1f  L_x^d=%.1f mm\n', ...
                 'L_y^t=%.1f  L_y^b=%.1f mm\n', ...
                 'a=%.1f mm   b=%.1f mm\n',      ...
                 'L_x/L_y=%.2f  |  N_{reg}=%d\n',...
                 '%s'],                           ...
            Lxu, Lxd, Lyt, Lyb, a_ell, b_ell, Lx/Ly, nReg(iL), flagStr), ...
        'Units','normalized','VerticalAlignment','top', ...
        'FontSize',7.5,'Color','w','BackgroundColor',[0 0 0 0.5]);

    plot(ax, 0, 0, 'w+', 'MarkerSize',10,'LineWidth',2);
    plot(ax, [-(xr-x_lo) -(xr-x_lo)],[min(dy_vec) max(dy_vec)],'w:','LineWidth',1);
    plot(ax, [ (x_hi-xr)  (x_hi-xr)],[min(dy_vec) max(dy_vec)],'w:','LineWidth',1);

    xlabel(ax,'\Deltax  (mm)');
    if iL == 1, ylabel(ax,'\Deltay  (mm)'); end
    title(ax, rhoLabels{iL});
    set(ax,'YDir','normal');
end

sgtitle(sprintf('Masks, boxes & ellipses  |  \\theta_{LS} = %.1f°  |  x_{ref}=%.1f, y_{ref}=%.1f mm', ...
    theta_deg, xr, yr), 'FontSize',12);

%% ===================== FIGURE 3: SUMMARY CHARTS ========================

a_vals = (bbox(:,1) + bbox(:,2)) / 2;
b_vals = (bbox(:,3) + bbox(:,4)) / 2;

figure('Color','w','Position',[50 100 1000 400],'Name','Summary');

% Left: bounding box extents
ax1 = subplot(1,3,1);
hold(ax1,'on'); box(ax1,'on');
b_bbox = bar(ax1, bbox','grouped');
for iL = 1:nLevels, b_bbox(iL).FaceColor = cmapRho(iL,:); end
set(ax1,'XTickLabel',{'L_x^u','L_x^d','L_y^t','L_y^b'});
ylabel(ax1,'Extent  (mm)'); title(ax1,'Bounding box extents');
legend(ax1, arrayfun(@(r) sprintf('\\rho=%.2f',r),rho_levels,'UniformOutput',false), ...
    'Location','best');
grid(ax1,'on');

% Middle: ellipse semi-axes
ax2 = subplot(1,3,2);
hold(ax2,'on'); box(ax2,'on');
b_ab = bar(ax2, [a_vals, b_vals],'grouped');
b_ab(1).FaceColor = [0.2 0.2 0.2];
b_ab(2).FaceColor = [0.7 0.7 0.7];
set(ax2,'XTick',1:nLevels,'XTickLabel',rhoLabels);
ylabel(ax2,'Semi-axis  (mm)'); title(ax2,'Ellipse semi-axes');
legend(ax2,{'a  (semi-major)','b  (semi-minor)'},'Location','best');
grid(ax2,'on');

% Right: LS inclination angle — scatter of tip pairs + fitted line
ax3 = subplot(1,3,3);
hold(ax3,'on'); box(ax3,'on');

x_line = linspace(min(x_tips)*1.3, max(x_tips)*1.3, 100);
plot(ax3, x_line, polyval(p_fit,x_line), 'k--', 'LineWidth',1.5);

for iL = 1:nLevels
    if any(isnan(tip_pts(iL,:))), continue; end
    plot(ax3, [tip_pts(iL,1), tip_pts(iL,3)], ...
              [tip_pts(iL,2), tip_pts(iL,4)], 'o-', ...
        'Color',cmapRho(iL,:),'MarkerFaceColor',cmapRho(iL,:), ...
        'LineWidth',1.5,'MarkerSize',7);
end
plot(ax3, 0, 0, 'k+','MarkerSize',10,'LineWidth',2);

xlabel(ax3,'\Deltax  (mm)'); ylabel(ax3,'\Deltay  (mm)');
title(ax3, sprintf('Inclination angle: \\theta = %.1f°', theta_deg));
legend(ax3, ['LS fit  (\theta='+string(sprintf('%.1f°',theta_deg))+')'; ...
    arrayfun(@(r) sprintf('\\rho=%.2f tips',r),rho_levels,'UniformOutput',false)'], ...
    'Location','best');
axis(ax3,'equal'); grid(ax3,'on');

sgtitle('Correlation structure summary','FontSize',12);

%% ===================== CONSOLE SUMMARY =================================

fprintf('\n%s\n', repmat('=',1,66));
fprintf('  Ellipse summary  (semi-axes from bbox, angle from LS tip fit)\n');
fprintf('%s\n', repmat('=',1,66));
fprintf('  Inclination angle  theta = %.2f deg\n\n', theta_deg);
fprintf('%-8s  %8s  %8s  %8s  %8s  %8s  %8s  %6s\n', ...
    'rho','Lxu(mm)','Lxd(mm)','Lyt(mm)','Lyb(mm)','a(mm)','b(mm)','a/b');
fprintf('%s\n', repmat('-',1,66));
for iL = 1:nLevels
    if any(isnan(bbox(iL,:))), continue; end
    fprintf('%-8.3f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %6.2f\n', ...
        rho_levels(iL), bbox(iL,1), bbox(iL,2), bbox(iL,3), bbox(iL,4), ...
        a_vals(iL), b_vals(iL), a_vals(iL)/b_vals(iL));
end
fprintf('%s\n\n', repmat('=',1,66));
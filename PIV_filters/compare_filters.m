% =========================================================
% USER CONFIG — edit these variables only
% =========================================================
clc; clear; close all; 
base_dir = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\PIV\PIV_filters\';
conditions = {
'timepodminmax_fulldataset_c45timeonly',...
'fulltimeandpodonly_allcams', ...
'fulltimeonly', ...
};
cmap = 'gray';

% --- Column 5 config ---
col5_mode = 'merged'; % 'piv' | 'merged'
col5_piv_file = '00075.mat';
col5_merged_pattern = 'merged_turbvariance_hann*.mat';
% =========================================================
% END OF USER CONFIG
% =========================================================

% Helper: extract filter type label from filename stem
%   e.g. '03_after_time_A' -> 'after time'
%        '04_after_pod_B'  -> 'after pod'
get_stage_label = @(fname) strrep( ...
    regexprep(fname, '^[0-9]+_(.+)_[AB]$', '$1'), '_', ' ');

n_rows = numel(conditions);
n_cols = 7; % col1=Before(A), col2=Before(B), col3=last(b0000 A),
            % col4=last(b0000 B), col5=last(b0001 A), col6=last(b0001 B), col7=merged

fig = figure('Units','normalized','Position',[0.05 0.05 0.98 0.85]);
tl = tiledlayout(n_rows, n_cols, 'TileSpacing','compact','Padding','compact');

fprintf('\n=== NaN percentages (piv_result(end).nan_mask) ===\n');

for r = 1:n_rows

    % ---- Paths for the two batches ----
    path_b0 = fullfile(base_dir, conditions{r}, 'batch_0000');
    path_b1 = fullfile(base_dir, conditions{r}, 'batch_0001');

    % ---- Find last filter stage for each batch ----
    % Lists all after_*.mat files, sorted, picks the last by prefix number
    find_last_stage = @(bpath, ab) find_last_filter(bpath, ab);

    [file_b0A, label_b0A] = find_last_filter(path_b0, 'A');
    [file_b0B, label_b0B] = find_last_filter(path_b0, 'B');
    [file_b1A, label_b1A] = find_last_filter(path_b1, 'A');
    [file_b1B, label_b1B] = find_last_filter(path_b1, 'B');

    % ---- Load all images ----
    imgs = cell(1,6);
    % Cols 1-2: before filtering
    for c = 1:2
        ab = 'A'; if c == 2, ab = 'B'; end
        fname = fullfile(path_b0, sprintf('00_before_filtering_%s.mat', ab));
        if isfile(fname)
            S = load(fname); flds = fieldnames(S);
            imgs{c} = double(S.(flds{1}));
        end
    end
    % Cols 3-6: last filter for each batch/channel
    src_files = {file_b0A, file_b0B, file_b1A, file_b1B};
    for c = 1:4
        if ~isempty(src_files{c}) && isfile(src_files{c})
            S = load(src_files{c}); flds = fieldnames(S);
            imgs{c+2} = double(S.(flds{1}));
        end
    end

    % ---- Column titles (only on first row) ----
    col_titles = {
        'Before (A)', 'Before (B)', ...
        sprintf('b0000 A\n[%s]', label_b0A), ...
        sprintf('b0000 B\n[%s]', label_b0B), ...
        sprintf('b0001 A\n[%s]', label_b1A), ...
        sprintf('b0001 B\n[%s]', label_b1B)
    };

    % ---- Clims ----
    clim_raw  = compute_clim(imgs(1:2));
    clim_filt = compute_clim(imgs(3:6));

    % Symmetric clim if last stage is 'pod'
    if contains(label_b0A, 'pod')
        abs_max = max(abs(clim_filt));
        clim_filt = [-abs_max, abs_max];
        cmap_filt = 'gray'; % diverging-friendly
    else
        cmap_filt = cmap;
    end

    % ---- Plot cols 1-6 ----
    for c = 1:6
        nexttile((r-1)*n_cols + c);
        if isempty(imgs{c})
            text(0.5,0.5,'File not found','HorizontalAlignment','center', ...
                'Units','normalized','FontSize',8,'Color','r');
            axis off;
        else
            imagesc(imgs{c});
            if c <= 2
                colormap(gca, cmap);   clim(clim_raw);
            else
                colormap(gca, cmap_filt); clim(clim_filt);
            end
            colorbar; axis image off;
        end
        if r == 1
            title(col_titles{c}, 'FontSize', 8, 'FontWeight', 'bold');
        end
        if c == 1
            ylabel(strrep(conditions{r},'_',' '), ...
                'FontSize',7,'FontWeight','bold','Interpreter','none');
        end
    end

    % ---- Col 7: merged field ----
    nexttile((r-1)*n_cols + 7);
    if strcmp(col5_mode, 'piv')
        col5_fpath = fullfile(base_dir, conditions{r}, col5_piv_file);
        if ~isfile(col5_fpath)
            text(0.5,0.5,'File not found','HorizontalAlignment','center', ...
                'Units','normalized','FontSize',8,'Color','r'); axis off;
        else
            N = load(col5_fpath);
            nan_mask = double(N.piv_result(3).nan_mask);
            nan_pct  = 100 * sum(nan_mask(:)) / numel(nan_mask);
            imagesc(nan_mask,[0 1]); colormap(gca,'hot');
            axis image off; colorbar;
            title_str = sprintf('NaN mask\n%.2f%% NaN', nan_pct);
            if r==1, title(title_str,'FontSize',9,'FontWeight','bold');
            else,    title(title_str,'FontSize',8); end
            fprintf(' [%-40s] NaN = %.2f%%\n', conditions{r}, nan_pct);
        end

    elseif strcmp(col5_mode, 'merged')
        col5_search = dir(fullfile(base_dir, conditions{r}, col5_merged_pattern));
        if isempty(col5_search)
            text(0.5,0.5,'File not found','HorizontalAlignment','center', ...
                'Units','normalized','FontSize',8,'Color','r'); axis off;
        else
            % col5_fpath = fullfile(col5_search(1).folder, col5_search(1).name);
            % N = load(col5_fpath);
            % quant_u = double(N.U_variance);
            % imagesc(N.worldX(1,:), N.worldY(:,1), quant_u);
            % set(gca,'YDir','normal'); colormap(gca,'gray');
            % axis image; colorbar;
            % clim([prctile(quant_u(:),1), prctile(quant_u(:),99)]);
            % if r==1, title('U_{variance}','FontSize',9,'FontWeight','bold'); end
        end
    end
end

fprintf('=================================================\n\n');

% Link axes within each column
for c = 1:n_cols
    col_axes = gobjects(n_rows,1);
    for r = 1:n_rows
        col_axes(r) = nexttile(tl, (r-1)*n_cols + c);
    end
    linkaxes(col_axes,'xy');
end

sgtitle(sprintf('Filtering comparison — last stage | %s', strjoin(conditions,' / ')), ...
    'FontSize',11,'Interpreter','none');


% =========================================================
% Separate figure: merged U_variance fields
% =========================================================
fig2 = figure('Units','normalized','Position',[0.05 0.05 0.95 0.5]);
t2 = tiledlayout(n_rows,1,'TileSpacing','compact','Padding','compact');
title(t2,'Merged U_{variance} fields','FontSize',12,'FontWeight','bold');

for r = 1:n_rows
    col5_search = dir(fullfile(base_dir, conditions{r}, col5_merged_pattern));
    nexttile;
    if isempty(col5_search)
        text(0.5,0.5,'File not found','HorizontalAlignment','center', ...
            'Units','normalized','FontSize',9,'Color','r'); axis off;
    else
        col5_fpath = fullfile(col5_search(1).folder, col5_search(1).name);
        N = load(col5_fpath);
        quant_u = double(N.U_variance);
        imagesc(N.worldX(1,:), N.worldY(:,1), quant_u);
        set(gca,'YDir','normal'); colormap(gca,'turbo');
        axis image; colorbar;
        % clim([prctile(quant_u(:),1), prctile(quant_u(:),99)]);
        clim([0 6]);
        xlabel('x (world)','FontSize',8); ylabel('y (world)','FontSize',8);
    end
    title(strrep(conditions{r},'_',' '),'FontSize',9,'FontWeight','bold','Interpreter','none');
end
% =========================================================
% PIV results figure: ux and peak_mag for camera4 & camera5
% Only runs when col5_mode = 'piv'
% =========================================================
if strcmp(col5_mode, 'piv')

    cameras    = {'camera4', 'camera5'};
    n_cams     = numel(cameras);
    n_cols_piv = n_cams * 2; % ux + peak_mag per camera = 4 cols

    fig3 = figure('Units','normalized','Position',[0.05 0.05 0.98 0.85]);
    tl3  = tiledlayout(n_rows, n_cols_piv, ...
                       'TileSpacing','compact','Padding','compact');
    title(tl3, sprintf('PIV results — %s (last window pass)', col5_piv_file), ...
          'FontSize',11,'FontWeight','bold');

    % Column headers
    col_headers = {};
    for k = 1:n_cams
        col_headers{end+1} = sprintf('%s — u_x',      cameras{k});
        col_headers{end+1} = sprintf('%s — peak mag',  cameras{k});
    end

    for r = 1:n_rows
        for k = 1:n_cams

            piv_fpath = fullfile(base_dir, conditions{r}, cameras{k}, col5_piv_file);

            tile_ux = (r-1)*n_cols_piv + (k-1)*2 + 1;
            tile_pk = (r-1)*n_cols_piv + (k-1)*2 + 2;

            % --- ux ---
            nexttile(tile_ux);
            if ~isfile(piv_fpath)
                text(0.5,0.5,'File not found','HorizontalAlignment','center', ...
                    'Units','normalized','FontSize',8,'Color','r');
                axis off;
            else
                N  = load(piv_fpath);
                ux = double(N.piv_result(end).ux);
                imagesc(ux);
                colormap(gca,'turbo');
                clim([prctile(ux(:),1), prctile(ux(:),99)]);
                colorbar; 
            end
            if r == 1
                title(col_headers{(k-1)*2+1}, 'FontSize',8,'FontWeight','bold');
            end
            if k == 1
                ylabel(strrep(conditions{r},'_',' '), ...
                       'FontSize',7,'FontWeight','bold','Interpreter','none');
            end

            % --- peak_mag ---
            nexttile(tile_pk);
            if ~isfile(piv_fpath)
                text(0.5,0.5,'File not found','HorizontalAlignment','center', ...
                    'Units','normalized','FontSize',8,'Color','r');
                axis off;
            else
                % peak_mag is 3×W×H — take primary peak (index 1)
                pk = squeeze(double(N.piv_result(end).peak_mag(1,:,:)));
                imagesc(pk);
                colormap(gca,'hot');
                clim([prctile(pk(:),1), prctile(pk(:),99)]);
                colorbar;
            end
            if r == 1
                title(col_headers{(k-1)*2+2}, 'FontSize',8,'FontWeight','bold');
            end

        end
    end

    % Link axes across rows per column
    for col = 1:n_cols_piv
        ax_col = gobjects(n_rows,1);
        for r = 1:n_rows
            ax_col(r) = nexttile(tl3, (r-1)*n_cols_piv + col);
        end
        linkaxes(ax_col,'xy');
    end

end

% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function [fpath, label] = find_last_filter(batch_path, ab)
    fpath = '';
    label = 'N/A';

    % Use plain wildcard — no bracket character classes
    hits = dir(fullfile(batch_path, sprintf('*_after_*_%s.mat', ab)));

    if isempty(hits)
        return;
    end

    % Filter to only files that start with a numeric prefix (e.g. 01_, 03_)
    valid = ~cellfun(@isempty, regexp({hits.name}, '^\d{2}_after_', 'match'));
    hits  = hits(valid);

    if isempty(hits)
        return;
    end

    % Sort by filename (numeric prefix guarantees correct order)
    [~, idx] = sort({hits.name});
    last  = hits(idx(end));
    fpath = fullfile(last.folder, last.name);

    % Extract label: '03_after_time_A.mat' -> 'after time'
    stem  = regexprep(last.name, '\.mat$', '');
    label = regexprep(stem, '^\d+_(.+)_[AB]$', '$1');
    label = strrep(label, '_', ' ');
end

function cl = compute_clim(img_cell)
    vals = [];
    for k = 1:numel(img_cell)
        if ~isempty(img_cell{k})
            v = img_cell{k}(:);
            vals = [vals; v(isfinite(v) & v ~= 0)]; %#ok<AGROW>
        end
    end
    if isempty(vals)
        cl = [0 1];
    elseif numel(vals) > 100
        cl = [prctile(vals,1), prctile(vals,99)];
    else
        cl = [min(vals), max(vals)];
    end
    if cl(1) >= cl(2), cl(2) = cl(1) + 1; end
end

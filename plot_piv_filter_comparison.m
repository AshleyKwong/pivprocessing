% Script to compare before/after PIV filter images
% Loads .mat files from batch folders and plots them in 2x2 subplots

function plot_piv_filter_comparison(basePath)
    % plot_piv_filter_comparison - Compare PIV images before and after filtering
    %
    % Input:
    %   basePath - Path to the basic_filters folder
    %              Example: 'E:\Case 6 PIV pivtools Processed\loop=00\basic_filters\150'

    % Get all batch folders
    batchFolders = dir(fullfile(basePath, 'batch_*'));
    batchFolders = batchFolders([batchFolders.isdir]);

    if isempty(batchFolders)
        error('No batch folders found in: %s', basePath);
    end

    fprintf('Found %d batch folders\n', length(batchFolders));

    % Process each batch folder
    for i = 1:length(batchFolders)
        batchPath = fullfile(basePath, batchFolders(i).name);
        fprintf('Processing %s...\n', batchFolders(i).name);

        % Get all .mat files in this batch
        matFiles = dir(fullfile(batchPath, '*.mat'));

        if isempty(matFiles)
            fprintf('  No .mat files found in %s\n', batchFolders(i).name);
            continue;
        end

        fprintf('  Found %d .mat files\n', length(matFiles));

        % Create figure for this batch
        figHandle = figure('Name', sprintf('%s - Filter Comparison', batchFolders(i).name), ...
                           'NumberTitle', 'off', 'Position', [100, 100, 1400, 1000]);

        % Process each .mat file
        for j = 1:min(4, length(matFiles))  % Plot up to 4 files in 2x2 grid
            matFile = fullfile(batchPath, matFiles(j).name);
            [~, fileName, ~] = fileparts(matFiles(j).name);

            fprintf('    Loading %s...\n', matFiles(j).name);

            % Load the data
            data = load(matFile);
            fieldNames = fieldnames(data);

            % Get the image array (assume first field contains the data)
            imgData = data.(fieldNames{1});

            % Determine subplot position based on file name
            if contains(lower(fileName), 'framea') && contains(lower(fileName), 'before')
                subplotIdx = 1;
                titleStr = 'Frame A - Before Filter';
            elseif contains(lower(fileName), 'framea') && contains(lower(fileName), 'after')
                subplotIdx = 2;
                titleStr = 'Frame A - After Filter';
            elseif contains(lower(fileName), 'frameb') && contains(lower(fileName), 'before')
                subplotIdx = 3;
                titleStr = 'Frame B - Before Filter';
            elseif contains(lower(fileName), 'frameb') && contains(lower(fileName), 'after')
                subplotIdx = 4;
                titleStr = 'Frame B - After Filter';
            else
                % If naming convention is different, plot sequentially
                subplotIdx = j;
                titleStr = strrep(fileName, '_', ' ');
            end

            % Plot the data
            subplot(2, 2, subplotIdx);
            imagesc(imgData);
            colormap gray;

            % Adjust color limits for better visibility
            % Use percentile-based scaling to enhance contrast
            imgMin = prctile(imgData(:), 1);   % 1st percentile
            imgMax = prctile(imgData(:), 99);  % 99th percentile
            caxis([imgMin, imgMax]);

            colorbar;
            axis image;

            title(titleStr, 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel('X (pixels)');
            ylabel('Y (pixels)');

            % Display intensity range in subtitle
            subtitle(sprintf('Range: [%.1f, %.1f]', imgMin, imgMax), 'FontSize', 9);
        end

        % Add main title
        sgtitle(sprintf('Batch %s - PIV Filter Comparison', batchFolders(i).name), ...
                'FontSize', 14, 'FontWeight', 'bold');

        % Save the figure
        savePath = fullfile(batchPath, sprintf('%s_comparison.png', batchFolders(i).name));
        saveas(figHandle, savePath);
        fprintf('  Saved comparison figure to: %s\n', savePath);
    end

    fprintf('\nProcessing complete!\n');
end

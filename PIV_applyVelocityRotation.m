function PIV_applyVelocityRotation(camInstDir, p_C, opts)

    arguments
        camInstDir          (1,1) string
        p_C                 (1,2) double
        opts.BackupOld      (1,1) logical = false
        opts.SaveUnrotated  (1,1) logical = false
    end

    camInstDir     = char(camInstDir);
    doBackup       = opts.BackupOld;
    doUnrotated    = opts.SaveUnrotated;

    if ~isfolder(camInstDir)
        error('PIV_applyVelocityRotation:folderNotFound', ...
              'Folder not found:\n  %s', camInstDir);
    end

    % ------------------------------------------------------------------
    %  1.  Rotation constants
    % ------------------------------------------------------------------
    theta = atan(p_C(1));
    cosT  = single(cos(theta));
    sinT  = single(sin(theta));

    fprintf('Camera folder : %s\n', camInstDir);
    fprintf('Floor slope   : %.6f  →  theta = %.6f rad (%.5f deg)\n', ...
            p_C(1), theta, rad2deg(theta));

    % ------------------------------------------------------------------
    %  2.  Discover frame files
    % ------------------------------------------------------------------
    allFiles   = dir(fullfile(camInstDir, '*.mat'));
    frameFiles = allFiles(~contains({allFiles.name}, 'coordinates'));

    if isempty(frameFiles)
        warning('PIV_applyVelocityRotation:noFiles', ...
                'No frame .mat files found in:\n  %s', camInstDir);
        return
    end

    nFiles    = numel(frameFiles);
    filePaths = fullfile({frameFiles.folder}, {frameFiles.name})';
    fprintf('Found %d frame files. Rotating (parfor)...\n', nFiles);

    % ------------------------------------------------------------------
    %  3.  Copy originals to 'unrotated/' BEFORE touching anything
    % ------------------------------------------------------------------
    unrotatedDir = fullfile(camInstDir, 'unrotated');

    if doUnrotated
        if ~isfolder(unrotatedDir)
            mkdir(unrotatedDir);
            fprintf('  Created backup folder: %s\n', unrotatedDir);
        else
            fprintf('  Backup folder already exists, skipping copy: %s\n', unrotatedDir);
        end
        fprintf('  Copying %d original files to unrotated/ ...', nFiles);
        for k = 1:nFiles
            [~, fname, ext] = fileparts(filePaths{k});
            destPath = fullfile(unrotatedDir, [fname ext]);
            if ~isfile(destPath)
                copyfile(filePaths{k}, destPath);
            end
        end
        fprintf(' done.\n');
    end

    % ------------------------------------------------------------------
    %  4.  Parallel rotation loop
    %      save() with a string variable name causes a parfor transparency
    %      violation — solved by delegating to helper functions below.
    % ------------------------------------------------------------------
    parfor k = 1:nFiles

        fPath      = filePaths{k};
        loaded     = load(fPath, 'piv_result');
        piv_result = loaded.piv_result;  %#ok<PFBNS>

        % Optional per-file backup (also uses helper to avoid violation)
        if doBackup
            [fdr, fname, ext] = fileparts(fPath);
            piv_result_bak    = piv_result;
            savePivResultBak(fullfile(fdr, [fname '_old' ext]), piv_result_bak);
        end

        % Rotate in-place, keeping single precision
        U = piv_result{end}.ux;
        V = piv_result{end}.uy;

        piv_result{end}.ux =  cosT .* U + sinT .* V;
        piv_result{end}.uy = -sinT .* U + cosT .* V;

        % Use helper — avoids transparency violation
        savePivResult(fPath, piv_result);

    end

    fprintf('✓ Done. All %d frames rotated by %.5f deg.\n', nFiles, rad2deg(theta));
    if doUnrotated
        fprintf('  Unrotated originals preserved in: %s\n', unrotatedDir);
    end

end

% ======================================================================
%  Private helper functions — save() is transparent here because MATLAB
%  can see exactly which variable is being written in each function scope.
% ======================================================================

function savePivResult(fPath, piv_result)
    save(fPath, 'piv_result', '-v7');
end

function savePivResultBak(fPath, piv_result_bak)
    save(fPath, 'piv_result_bak', '-v7');
end

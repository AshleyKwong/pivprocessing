function PIV_applyVelocityRotation(camInstDir, p_C, opts)
% PIV_applyVelocityRotation  Rotate PIV velocity vectors into the floor-aligned frame.
%   Key changes vs v1:
%     - parfor over files (uses all allocated SLURM cores)
%     - '-v7' instead of '-v7.3' (no HDF5 overhead for small arrays)
%     - arithmetic kept in single (no double upcast needed)

    arguments
        camInstDir     (1,1) string
        p_C            (1,2) double
        opts.BackupOld (1,1) logical = false
    end

    camInstDir = char(camInstDir);

    if ~isfolder(camInstDir)
        error('PIV_applyVelocityRotation:folderNotFound', ...
              'Folder not found:\n  %s', camInstDir);
    end

    % ------------------------------------------------------------------
    %  1.  Rotation constants  (broadcast into parfor workers)
    % ------------------------------------------------------------------
    theta = atan(p_C(1));
    cosT  = single(cos(theta));     % single: matches ux/uy dtype, avoids upcast
    sinT  = single(sin(theta));
    doBackup = opts.BackupOld;      % extract struct field — parfor requires plain vars

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
    filePaths = fullfile({frameFiles.folder}, {frameFiles.name})';  % cell col-vec
    fprintf('Found %d frame files. Rotating (parfor)...\n', nFiles);

    % ------------------------------------------------------------------
    %  3.  Parallel loop — each file is fully independent
    % ------------------------------------------------------------------
    parfor k = 1:nFiles

        fPath      = filePaths{k};
        loaded     = load(fPath, 'piv_result');
        piv_result = loaded.piv_result;   %#ok<PFBNS>

        % Optional backup — done before any modification
        if doBackup
            [fdr, fname, ext] = fileparts(fPath);
            piv_result_bak    = piv_result;
            save(fullfile(fdr, [fname '_old' ext]), 'piv_result_bak', '-v7');
        end

        % Rotate in-place, keeping single precision throughout
        U = piv_result{end}.ux;
        V = piv_result{end}.uy;

        piv_result{end}.ux =  cosT .* U + sinT .* V;
        piv_result{end}.uy = -sinT .* U + cosT .* V;

        % '-v7' is ~3-5x faster than '-v7.3' for small arrays (no HDF5 overhead)
        save(fPath, 'piv_result', '-v7');

    end

    fprintf('✓ Done. All %d frames rotated by %.5f deg.\n', nFiles, rad2deg(theta));

end

function PIV_applyVelocityRotation(camInstDir, p_C, opts)
% PIV_applyVelocityRotation  Rotate PIV velocity vectors into the floor-aligned frame.
%
%   PIV_applyVelocityRotation(camInstDir, p_C)
%   PIV_applyVelocityRotation(camInstDir, p_C, BackupOld=true)
%
%   Inputs
%   ------
%   camInstDir : char | string
%       Full path to the 'instantaneous' subfolder for one camera, e.g.:
%       'E:\ProcessedPIV_case2fullpipe\loop=05\calibrated_piv\150\Cam1\instantaneous'
%
%   p_C        : 1x2 double
%       Linear floor-fit coefficients [slope, intercept] from PIV_detectFloor.
%       The rotation angle is derived as theta = atan(p_C(1)).
%
%   Options (name-value)
%   --------------------
%   BackupOld  : logical, default false
%       If true, saves each original file as 00001_old.mat etc. before
%       overwriting. WARNING: this doubles disk usage — only use on a
%       small subset to verify the rotation is correct first.
%
%   What it does
%   ------------
%   For every 000XX.mat file in camInstDir:
%     1. Loads piv_result  (1x3 cell array)
%     2. Rotates piv_result{end}.ux and piv_result{end}.uy by -theta
%        (transforms from tilted lab frame into floor-parallel / wall-normal frame)
%     3. Leaves piv_result{1}, piv_result{2}, and b_mask untouched
%     4. Overwrites the original .mat file with the modified piv_result
%
%   Rotation applied (same theta for all frames — it is a property of the rig):
%       U_rot =  U*cos(theta) + V*sin(theta)
%       V_rot = -U*sin(theta) + V*cos(theta)

    arguments
        camInstDir  (1,1) string
        p_C         (1,2) double
        opts.BackupOld (1,1) logical = false
    end

    camInstDir = char(camInstDir);

    if ~isfolder(camInstDir)
        error('PIV_applyVelocityRotation:folderNotFound', ...
              'Folder not found:\n  %s', camInstDir);
    end

    % ------------------------------------------------------------------
    %  1.  Derive rotation angle from floor fit slope
    % ------------------------------------------------------------------
    theta = atan(p_C(1));   % radians; same for every frame in this camera
    fprintf('Camera folder : %s\n', camInstDir);
    fprintf('Floor slope   : %.6f  →  theta = %.6f rad (%.5f deg)\n', ...
            p_C(1), theta, rad2deg(theta));

    cosT =  cos(theta);
    sinT =  sin(theta);

    % ------------------------------------------------------------------
    %  2.  Discover all frame files (exclude coordinates*.mat)
    % ------------------------------------------------------------------
    allFiles  = dir(fullfile(camInstDir, '*.mat'));
    frameFiles = allFiles(~contains({allFiles.name}, 'coordinates'));

    if isempty(frameFiles)
        warning('PIV_applyVelocityRotation:noFiles', ...
                'No frame .mat files found in:\n  %s', camInstDir);
        return
    end

    fprintf('Found %d frame files. Rotating...\n', numel(frameFiles));

    % ------------------------------------------------------------------
    %  3.  Loop over frames
    % ------------------------------------------------------------------
    for k = 1:numel(frameFiles)

        fPath = fullfile(frameFiles(k).folder, frameFiles(k).name);

        % Load
        loaded     = load(fPath, 'piv_result');
        piv_result = loaded.piv_result;   %#ok<NASGU> — will be modified below

        % Extract U and V from the last cell (index 3 = end)
        U = double(piv_result{end}.ux);
        V = double(piv_result{end}.uy);

        % Rotate
        U_rot =  cosT .* U + sinT .* V;
        V_rot = -sinT .* U + cosT .* V;

        % Write back (preserve original numeric class)
        piv_result{end}.ux = cast(U_rot, class(piv_result{end}.ux));
        piv_result{end}.uy = cast(V_rot, class(piv_result{end}.uy));

        % Optional backup
        if opts.BackupOld
            [~, fname, ext] = fileparts(fPath);
            backupPath = fullfile(frameFiles(k).folder, [fname '_old' ext]);
            loaded_orig = load(fPath, 'piv_result'); %#ok — re-read before overwrite
            piv_result_orig = loaded_orig.piv_result;
            save(backupPath, 'piv_result_orig', '-v7.3');
        end

        % Overwrite
        save(fPath, 'piv_result', '-v7.3');

        if mod(k, 50) == 0 || k == numel(frameFiles)
            fprintf('  Processed %d / %d files\n', k, numel(frameFiles));
        end

    end

    fprintf('✓ Done. All frames rotated by %.5f deg.\n', rad2deg(theta));

end

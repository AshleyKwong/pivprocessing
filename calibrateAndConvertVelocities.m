function [allCameras, windowCenterCameras_mm] = calibrateAndConvertVelocities( ...
        base_dir, calib, windowCenterCameras_mm, camIdx, ...
        final_pass_index, skip_frames, dt, skipCalibration)
% CALIBRATEANDCONVERTVELOCITIES
%   Applies polynomial calibration to uncalibrated PIV pixel displacements
%   and converts to physical velocities in m/s.
%
%   INPUTS:
%     base_dir               - dir() struct of .mat files for this camera
%     calib                  - calibration struct (calib_tools, poly_2d2c)
%     windowCenterCameras_mm - struct with .x1_mm / .x2_mm cell arrays
%                              (pre-populated if skipCalibration = true)
%     camIdx                 - camera index (integer, 1–5)
%     final_pass_index       - PIVtools multi-pass index to read (e.g. 5)
%     skip_frames            - vector of frame indices to skip ([] if none)
%     dt                     - inter-frame time [s]
%     skipCalibration        - logical; true = load centres from struct,
%                              false = compute from first frame
%
%   OUTPUTS:
%     allCameras             - struct with fields:
%                                .u  {nFrames × 1} single, streamwise [m/s]
%                                .v  {nFrames × 1} single, wall-normal [m/s]
%     windowCenterCameras_mm - input struct, updated with x1_mm / x2_mm
%                              for camIdx (only modified if ~skipCalibration)

% ---- 1. Physical window centres (pixel → mm) ----------------------------
if ~skipCalibration
    tmp   = load(fullfile(base_dir(1).folder, base_dir(1).name));
    win_y = single(tmp.piv_result(final_pass_index).win_ctrs_y);
    win_x = single(tmp.piv_result(final_pass_index).win_ctrs_x);
    clear tmp;

    [px1, px2] = meshgrid(win_x, win_y);

    cx1_mm = reshape( ...
        polyvaln(calib.fit_cx1{camIdx}.p, [px1(:), px2(:)]), size(px1));
    cx2_mm = reshape( ...
        polyvaln(calib.fit_cx2{camIdx}.p, [px1(:), px2(:)]), size(px2));

    windowCenterCameras_mm.x1_mm{camIdx} = cx1_mm;
    windowCenterCameras_mm.x2_mm{camIdx} = cx2_mm;
    fprintf('    ✓ Computed window centres [%d×%d]\n', size(cx1_mm,1), size(cx1_mm,2));
else
    cx1_mm = windowCenterCameras_mm.x1_mm{camIdx};
    cx2_mm = windowCenterCameras_mm.x2_mm{camIdx};
    [px1, px2] = meshgrid(1:size(cx1_mm,2), 1:size(cx1_mm,1));
    fprintf('    ✓ Loaded window centres [%d×%d]\n', size(cx1_mm,1), size(cx1_mm,2));
end

% ---- 2. Preallocate outputs ---------------------------------------------
nTotalFrames = length(base_dir) - 1;   % last file is often a summary
nValidFrames = nTotalFrames - length(skip_frames);

allCameras.u = cell(nValidFrames, 1);
allCameras.v = cell(nValidFrames, 1);

% ---- 3. Frame loop: displace → calibrate → velocity --------------------
counter = 0;
for b = 1:nTotalFrames
    if ismember(b, skip_frames), continue; end
    counter = counter + 1;

    fd = load(fullfile(base_dir(b).folder, base_dir(b).name), 'piv_result');
    uy = single(fd.piv_result(final_pass_index).uy);
    ux = single(fd.piv_result(final_pass_index).ux);
    clear fd;

    % Displaced pixel positions (centre + displacement)
    px1_disp = px1 + uy;
    px2_disp = px2 + ux;

    % Physical positions at displaced pixel grid
    vx1_mm = reshape( ...
        polyvaln(calib.fit_cx1{camIdx}.p, [px1_disp(:), px2_disp(:)]), size(uy));
    vx2_mm = reshape( ...
        polyvaln(calib.fit_cx2{camIdx}.p, [px1_disp(:), px2_disp(:)]), size(uy));

    % Velocity = physical displacement / dt  [m/s]
    allCameras.u{counter} = (vx2_mm - cx2_mm) * 1e-3 / dt;
    allCameras.v{counter} = (vx1_mm - cx1_mm) * 1e-3 / dt;
end
fprintf('    ✓ %d / %d frames processed\n', counter, nTotalFrames);

end

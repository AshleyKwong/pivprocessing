function PIV_applyFloorShift(coordFilePath, p_C)
% PIV_applyFloorShift  Shift PIV window-centre Y coordinates to a floor-referenced frame.
%
%   PIV_applyFloorShift(coordFilePath, p_C)
%
%   Inputs
%   ------
%   coordFilePath : char | string
%       Full path to the coordinates.mat file for one camera subfolder,
%       e.g. 'E:\...\Cam1\instantaneous\coordinates.mat'
%
%   p_C           : 1x2 double
%       Linear polynomial coefficients [slope, intercept] of the floor fit
%       as returned by PIV_detectFloor  (i.e. polyval(p_C, x) gives the
%       floor height in mm at each x position).
%       Determined once from the ensemble-averaged field and reused here.
%
%   Behaviour
%   ---------
%   - Saves the original file as coordinates_old.mat (sibling of input file).
%   - Overwrites coordinates.mat with floor-shifted y values.
%   - coordinates{end}.x is left unchanged.
%   - coordinates{end}.y  <-- y_original - polyval(p_C, x)
%
%   Example
%   -------
%   p_C = [-0.00123, 12.456];   % from PIV_detectFloor on your mean field
%   PIV_applyFloorShift('E:\ProcessedPIV_case2fullpipe\loop = 1\calibrated_piv\150\Cam1\instantaneous\coordinates.mat', p_C);

    % ------------------------------------------------------------------ %
    %  0.  Validate inputs                                                 %
    % ------------------------------------------------------------------ %
    arguments
        coordFilePath  (1,1) string
        p_C            (1,2) double
    end

    coordFilePath = char(coordFilePath);

    if ~isfile(coordFilePath)
        error('PIV_applyFloorShift:fileNotFound', ...
              'Cannot find coordinates file:\n  %s', coordFilePath);
    end

    % ------------------------------------------------------------------ %
    %  1.  Load                                                            %
    % ------------------------------------------------------------------ %
    loaded     = load(coordFilePath);        % loads variable 'coordinates'
    coordinates = loaded.coordinates;

    % Accept both cell-array and plain struct storage
    if iscell(coordinates)
        coords_end = coordinates{end};
    else
        coords_end = coordinates(end);
    end

    % ------------------------------------------------------------------ %
    %  2.  Compute the floor fit on the x grid of this camera             %
    % ------------------------------------------------------------------ %
    x_vec  = coords_end.x(1, :);            % [1 x nX]  x positions (mm)
    C_fit  = polyval(p_C, x_vec);           % [1 x nX]  floor height at each x

    nY          = size(coords_end.y, 1);
    C_fit_grid  = repmat(C_fit, nY, 1);     % [nY x nX]  broadcast to full grid

    y_shifted = coords_end.y - C_fit_grid;  % floor-referenced y

    % ------------------------------------------------------------------ %
    %  3.  Save old file, then overwrite with shifted coordinates         %
    % ------------------------------------------------------------------ %
    [folder, ~, ~]  = fileparts(coordFilePath);
    oldPath         = fullfile(folder, 'coordinates_old.mat');
    newPath         = coordFilePath;          % same path — overwrite in place

    % Back-up original
    save(oldPath, 'coordinates');
    fprintf('  Saved original  →  %s\n', oldPath);

    % Write shifted y back into the same struct/cell structure
    if iscell(coordinates)
        coordinates{end}.y = y_shifted;
    else
        coordinates(end).y = y_shifted;
    end

    save(newPath, 'coordinates');
    fprintf('  Saved shifted   →  %s\n', newPath);
    fprintf('  Floor pitch: slope = %.6f  intercept = %.4f mm\n', p_C(1), p_C(2));
    fprintf('  Y range after shift: [%.3f, %.3f] mm\n', min(y_shifted(:)), max(y_shifted(:)));

end

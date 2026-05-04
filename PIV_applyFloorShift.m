function PIV_applyFloorShift(coordFilePath, p_C)
    arguments
        coordFilePath  (1,1) string
        p_C            (1,2) double
    end

    coordFilePath = char(coordFilePath);

    if ~isfile(coordFilePath)
        error('PIV_applyFloorShift:fileNotFound', ...
              'Cannot find coordinates file:\n  %s', coordFilePath);
    end

    % ------------------------------------------------------------------
    %  1.  Load
    % ------------------------------------------------------------------
    loaded      = load(coordFilePath);
    coordinates = loaded.coordinates;

    if iscell(coordinates)
        coords_end = coordinates{end};
    else
        coords_end = coordinates(end);
    end

    % ------------------------------------------------------------------
    %  2.  Compute floor fit and shift y
    % ------------------------------------------------------------------
    x_vec      = coords_end.x(1, :);
    C_fit      = polyval(p_C, x_vec);
    nY         = size(coords_end.y, 1);
    C_fit_grid = repmat(C_fit, nY, 1);
    y_shifted  = coords_end.y - C_fit_grid;

    % ------------------------------------------------------------------
    %  3.  Save old, then overwrite — both explicitly '-v7'
    % ------------------------------------------------------------------
    [folder, ~, ~] = fileparts(coordFilePath);
    oldPath        = fullfile(folder, 'coordinates_old.mat');

    save(oldPath, 'coordinates', '-v7');                % ← was: no flag
    fprintf('  Saved original  →  %s\n', oldPath);

    if iscell(coordinates)
        coordinates{end}.y = y_shifted;
    else
        coordinates(end).y = y_shifted;
    end

    save(coordFilePath, 'coordinates', '-v7');           % ← was: no flag
    fprintf('  Saved shifted   →  %s\n', coordFilePath);
    fprintf('  Floor pitch: slope = %.6f  intercept = %.4f mm\n', p_C(1), p_C(2));
    fprintf('  Y range after shift: [%.3f, %.3f] mm\n', min(y_shifted(:)), max(y_shifted(:)));

end

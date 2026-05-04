%% Load coordinates.mat for each camera
% Root path for the specific loop and resolution\
clear ;clc ; close all;
%%
savePath   = 'E:\ProcessedPIV_case2fullpipe';   % root folder
cameraList = ["Cam1","Cam2","Cam3","Cam4","Cam5"];
%% DISCOVER LOOP FOLDERS
% -------------------------------------------------------------------------
d = dir(savePath);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));

loopPattern = '^loop\s*=\s*\d+$';
validLoops  = false(size(d));
for i = 1:length(d)
    validLoops(i) = ~isempty(regexpi(d(i).name, loopPattern));
end
totalLoops = d(validLoops);

if isempty(totalLoops)
    error('No "loop=XX" folders found in %s', savePath);
else
    fprintf('Found %d loop folders:\n', length(totalLoops));
    for i = 1:length(totalLoops)
        fprintf('  %s\n', totalLoops(i).name);
    end
end
mean_U_per_cam = cell(1, numel(cameraList));   % [nY x nX] single, one per camera
mean_V_per_cam = cell(1, numel(cameraList));
coords_per_cam = cell(1, numel(cameraList));   % coordinates struct, one per camera

%% LOAD PIVTOOLS CALIBRATED/MERGED FIELDS
% -------------------------------------------------------------------------
% pivtools_inst_all = {}; % instantaneous: {frame, loopNo} cell of structs

for cameraNo = 1:numel(cameraList) % per each camera we are going to average within each loop, and then average over all the loops.
    pivtools_mean_U     = [];
    pivtools_mean_V     = [];
    pivtools_mean_count = 0;
 
    for loopNo = 1%:length(totalLoops)
        loopU = [];
        loopV = [];
        fprintf('\n%s\n', fullfile(savePath, totalLoops(loopNo).name, ...
            'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous'));


        camDir = fullfile(savePath, totalLoops(loopNo).name, ...
            'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous/');

        if ~isfolder(camDir)
            fprintf('  ⚠  No %s folder in %s – skipping\n', ...
                cameraList(cameraNo), totalLoops(loopNo).name);
            continue;
        end

        mergedFiles = dir(fullfile(camDir, '*.mat'));
        mergedFiles = mergedFiles(~strcmp({mergedFiles.name}, 'coordinates.mat'));

        if isempty(mergedFiles)
            fprintf('  ⚠  %s folder empty in %s – skipping\n', ...
                cameraList(cameraNo), totalLoops(loopNo).name);
            continue;
        end

        coordinates_file = fullfile(savePath, totalLoops(loopNo).name, ...
            'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous/coordinates.mat');
        
        if ~isfile(coordinates_file)
            fprintf('   ⚠  No coordinates.mat for %s in %s\n', ...
                cameraList(cameraNo), totalLoops(loopNo).name);
            continue;
        else
            tmp = load(coordinates_file).coordinates;
            if iscell(tmp)
                coords_per_cam{cameraNo} = tmp{end};
            else
                coords_per_cam{cameraNo} = tmp;   % already a struct, use directly
            end

        end


        for f = 1:length(mergedFiles)
            pivData = load(fullfile(mergedFiles(f).folder, mergedFiles(f).name));

            % --- Flexible field detection ------------------------------------
            % PIVtools may store fields under different variable names.
            % Common candidates: piv_result, U, V, u, v, Ux, Uy, vel
            if isfield(pivData, 'piv_result')
                if iscell(pivData.piv_result)
                    U_frame = single(pivData.piv_result{end}.ux);   % cell array → brace index
                    V_frame = single(pivData.piv_result{end}.uy);
                else
                    U_frame = single(pivData.piv_result(end).ux);   % struct array → parent index
                    V_frame = single(pivData.piv_result(end).uy);
                end
            elseif isfield(pivData, 'U') && isfield(pivData, 'V')
                U_frame = single(pivData.U);
                V_frame = single(pivData.V);
            elseif isfield(pivData, 'u') && isfield(pivData, 'v')
                U_frame = single(pivData.u);
                V_frame = single(pivData.v);
            elseif isfield(pivData, 'Ux') && isfield(pivData, 'Uy')
                U_frame = single(pivData.Ux);
                V_frame = single(pivData.Uy);
            else
                fprintf('  ⚠  Frame %d loop %s: unrecognised field names. Fields: %s\n', ...
                    f, totalLoops(loopNo).name, strjoin(fieldnames(pivData)', ', '));
                continue;
            end

            % Store instantaneous frame
            % pivtools_inst_all{end+1} = struct('U', U_frame, 'V', V_frame, ...
            %     'loop', loopNo, 'frame', f); %#ok<SAGROW>

            % Accumulate for mean
            if isempty(loopU)
                loopU = double(U_frame);
                loopV = double(V_frame);
            else
                loopU = loopU + double(U_frame);
                loopV = loopV + double(V_frame);
            end
        end

        % Running mean accumulation across loops
        if ~isempty(loopU)
            loopU = loopU / length(mergedFiles);
            loopV = loopV / length(mergedFiles);
            if isempty(pivtools_mean_U)
                pivtools_mean_U = loopU;
                pivtools_mean_V = loopV;
            else
                pivtools_mean_U = pivtools_mean_U + loopU;
                pivtools_mean_V = pivtools_mean_V + loopV;
            end
            pivtools_mean_count = pivtools_mean_count + 1;
            fprintf('  ✓ Loop %s: loaded %d merged frames\n', ...
                totalLoops(loopNo).name, length(mergedFiles));
        end
    end % loopNo

    if pivtools_mean_count > 0
        mean_U_per_cam{cameraNo} = pivtools_mean_U / pivtools_mean_count;
        mean_V_per_cam{cameraNo} = pivtools_mean_V / pivtools_mean_count;
        fprintf('✓ PIVtools mean computed for %s, %d total frames\n', ...
            cameraList(cameraNo), length(mergedFiles));
    else
        warning('No PIVtools Merged data found. Check folder structure.');
    end
end % end of cam

%-------------------------------------------------------------------------
%% loading the floor
close all;
% cameraList = ["Cam1","Cam2","Cam3","Cam4","Cam5"];

% % --- Robust load: handles both flat and nested .mat structures ----------
% matPath  = 'C:\Users\ak1u24\Downloads\pivtools_mean_UV_allcams.mat';
% raw      = load(matPath);
% 
% % Show what's actually in the file so you can see the structure
% disp('Variables found in .mat file:');
% disp(fieldnames(raw));
% 
% if isfield(raw, 'mean_pivtools')
%     % Nested case: everything lives inside mean_pivtools struct
%     mean_pivtools  = raw.mean_pivtools;
%     coords_per_cam = mean_pivtools.coords_per_cam;
%     mean_U_per_cam = mean_pivtools.mean_U_per_cam;
%     mean_V_per_cam = mean_pivtools.mean_V_per_cam;
% else
%     % Flat case: variables saved directly at top level
%     coords_per_cam = raw.coords_per_cam;
%     mean_U_per_cam = raw.mean_U_per_cam;
%     mean_V_per_cam = raw.mean_V_per_cam;
%     % Reconstruct mean_pivtools as a container for new outputs
%     mean_pivtools  = raw;
% end
% % -----------------------------------------------------------------------
% 
cam_y_bands = [ ...
   -14,   10;   % Cam1
    -7,   10;   % Cam2
   -15,   10;   % Cam3
   -21,   10;   % Cam4
    -0.5, 10 ]; % Cam5

for cameraNo = 1:numel(cameraList)

    % --- Coordinates (struct inside cell — note {}, not .) ----------
    x_vec_pt = coords_per_cam{cameraNo}.x(1, :);   % [1 x nX]
    y_vec    = coords_per_cam{cameraNo}.y(:, 1);    % [nY x 1]

    % --- Floor detection --------------------------------------------
    [~, p_C, ~, floor_y_C, ~, pitch_C_mm] = PIV_detectFloor( ...
        x_vec_pt, y_vec, mean_U_per_cam{cameraNo}, cam_y_bands(cameraNo,1), cam_y_bands(cameraNo,2)); %cam_y_bands(cameraNo,1), cam_y_bands(cameraNo,2)

    p_C_per_cam(cameraNo, :)       = p_C;           % ← was cut off before
    mean_pivtools.pC{cameraNo}     = p_C;
    pitch_C_deg                    = atand(p_C(1));
    theta_rad                      = atan(p_C(1));
    theta_per_cam(cameraNo)        = theta_rad;

    % --- Floor shift on coordinates ---------------------------------
    C_fit      = polyval(p_C, x_vec_pt);            % [1 x nX]
    nY         = size(coords_per_cam{cameraNo}.y, 1);
    C_fit_grid = repmat(C_fit, nY, 1);              % [nY x nX]

    coords_per_cam_corrected.x1_mm{cameraNo} = coords_per_cam{cameraNo}.x;
    coords_per_cam_corrected.x2_mm{cameraNo} = coords_per_cam{cameraNo}.y - C_fit_grid;

    % --- Velocity rotation ------------------------------------------
    U = mean_U_per_cam{cameraNo};
    V = mean_V_per_cam{cameraNo};
    mean_pivtools.mean_U_per_cam_rotated{cameraNo} =  U * cos(theta_rad) + V * sin(theta_rad);
    mean_pivtools.mean_V_per_cam_rotated{cameraNo} = -U * sin(theta_rad) + V * cos(theta_rad);

    % ================================================================
    %  PLOTS
    % ================================================================
    figure('Name', sprintf('%s – Floor correction summary', cameraList(cameraNo)), ...
           'Units','normalized', 'Position',[0.1 0.05 0.8 0.9]);

    % --- Subplot 1: Raw <U>, original coordinates, unrotated --------
    ax1 = subplot(3,1,1);
    imagesc(x_vec_pt, y_vec, mean_U_per_cam{cameraNo});
    set(gca,'YDir','normal'); axis image; colormap(gca, jet);
    colorbar; clim([0 30]); hold on;
    plot(x_vec_pt, floor_y_C, 'g.', 'MarkerSize', 4, ...
        'DisplayName', 'Detected floor points');
    plot(x_vec_pt, C_fit,     'g-', 'LineWidth', 2,  ...
        'DisplayName', sprintf('Linear fit  (%.5f°)', pitch_C_deg));
    xlabel('X (mm)'); ylabel('Y (mm)');
    title(sprintf('%s  |  Raw \\langle U \\rangle — original lab coordinates, unrotated vectors', ...
        cameraList(cameraNo)));
    legend('Location','northeast');

    % --- Subplot 2: Corrected <U>, floor-shifted coords, rotated ----
    ax2 = subplot(3,1,2);
    imagesc(coords_per_cam_corrected.x1_mm{cameraNo}(1,:), ...
            coords_per_cam_corrected.x2_mm{cameraNo}(:,1), ...
            mean_pivtools.mean_U_per_cam_rotated{cameraNo});
    set(gca,'YDir','normal'); axis image; colormap(gca, jet);
    colorbar; clim([0 30]); hold on;
    yline(0, 'w--', 'LineWidth', 1.5, 'DisplayName', 'Wall  (y = 0)');
    xlabel('X (mm)'); ylabel('y – y_{wall}  (mm)');
    title(sprintf('%s  |  Corrected \\langle U \\rangle — floor-referenced coordinates, rotated vectors', ...
        cameraList(cameraNo)));
    legend('Location','northeast');

    % --- Subplot 3: Floor detection diagnostic ----------------------
    ax3 = subplot(3,1,3);
    hold on;
    plot(x_vec_pt, floor_y_C, 'b.', 'MarkerSize', 6, ...
        'DisplayName', 'Detected floor points');
    plot(x_vec_pt, C_fit,     'g-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Linear fit:  slope=%.5f,  intercept=%.4f mm,  rise=%.3f mm,  %.5f°', ...
        p_C(1), p_C(2), pitch_C_mm, pitch_C_deg));
    xlabel('X (mm)'); ylabel('Y (mm)');
    title(sprintf('%s  |  Floor detection — detected points vs linear fit', cameraList(cameraNo)));
    legend('Location','best'); grid on; hold off;

    linkaxes([ax1, ax2, ax3], 'x');

end
% Add this just before the save line
mean_pivtools.coords_per_cam_corrected = coords_per_cam_corrected;


save('pivtools_mean_UV_allcams_corrected.mat', 'mean_pivtools');

%% now need to interpolate onto a grid

masks = {}; 
[worldX, worldY, U_hann, ~] = merge_cameras_python_style_mean(...
    coords_per_cam_corrected, mean_pivtools.mean_U_per_cam_rotated, mean_pivtools.mean_V_per_cam_rotated, masks, 'hann', []);


figure();
imagesc(worldX(1,:), (worldY(:, 1)), U_hann); % no need to flip - just set ydir as normal.
title('Hann'); colorbar;
axis image;
set(gca, 'YDir', 'normal');
colormap(jet);
% now perform a final floor check
[~, p_C, ~, floor_y_C,  ~, pitch_C_mm ] = PIV_detectFloor(worldX(1,:), worldY(:,1),U_hann);

figure();
imagesc(worldX(1,:), (worldY(:, 1)), U_hann); % no need to flip - just set ydir as normal.
set(gca,'YDir','normal'); axis image; colormap(gca, jet);
colorbar; clim([0 30]); hold on;
plot(worldX(1,:), floor_y_C,'r.',  'MarkerSize', 4);
plot(worldX(1,:), polyval(p_C,worldX(1,:)), 'r-',  'LineWidth', 2);
xlabel('X (mm)'); ylabel('Y (mm)');
title('Wall detections and linear fits');


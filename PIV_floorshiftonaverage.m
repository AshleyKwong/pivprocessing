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
 
    for loopNo = 2%:length(totalLoops)
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
mean_pivtools = load('C:\Users\ak1u24\Downloads\pivtools_mean_UV_allcams.mat'); 
coords_per_cam = mean_pivtools.coords_per_cam; 
mean_U_per_cam = mean_pivtools.mean_U_per_cam; 

% *** NEW: pre-allocate a row per camera to store [slope, intercept]
p_C_per_cam = nan(numel(cameraList), 2);
% now we will load the floor correction per camera
coords_per_cam_corrected = struct();
coords_per_cam_corrected.x1_mm = cell(1, numel(cameraList));
coords_per_cam_corrected.x2_mm = cell(1, numel(cameraList));
cam_y_bands = [ ...
   -12,  -10;   % Cam1
   -6,  -3;   % Cam2
   -14,  -10;   % Cam3
   -21,  -16;   % Cam4
   -0.5,  10 ]; % Cam5
for cameraNo = 1:numel(cameraList)
    x_vec_pt = coords_per_cam{cameraNo}.x(1, :);
    y_vec    = coords_per_cam{cameraNo}.y(:, 1);
    [~, p_C, ~, floor_y_C, ~, pitch_C_mm] = PIV_detectFloor( ...
        x_vec_pt, y_vec, mean_U_per_cam{cameraNo}, ...
        cam_y_bands(cameraNo, 1), cam_y_bands(cameraNo, 2));   % ← per-camera     p_C_per_cam(cameraNo, :) = p_C;
    
    
    mean_pivtools.pC{cameraNo} = p_C; % save a copy.
    pitch_C_deg = atand(p_C(1));
    theta_rad         = atan(p_C(1)); 
    theta_per_cam(cameraNo) = theta_rad;

    C_fit = polyval(p_C, x_vec_pt);
    nY          = size(coords_per_cam{cameraNo}.y, 1);
    C_fit_grid  = repmat(C_fit, nY, 1);     % [nY x nX]  broadcast to full grid

    coords_per_cam_corrected.x2_mm{cameraNo} = coords_per_cam{cameraNo}.y - C_fit_grid;  % floor-referenced y
    coords_per_cam_corrected.x1_mm{cameraNo} =coords_per_cam{cameraNo}.x ;  % floor-referenced y
    
    % --- 2. ROTATION: align velocity vectors with floor frame ------------
    U = mean_U_per_cam{cameraNo};   % [nY x nX]
    V = mean_V_per_cam{cameraNo};   % [nY x nX]

    % theta > 0 means floor rises in +X  →  rotate frame by -theta
    mean_pivtools.mean_U_per_cam_rotated{cameraNo} =  U * cos(theta_rad) + V * sin(theta_rad);
    mean_pivtools.mean_V_per_cam_rotated{cameraNo} = -U * sin(theta_rad) + V * cos(theta_rad);


    % --- plots (unchanged) ---
    figure();
    ax1 = subplot(3,1,1);
    imagesc(coords_per_cam_corrected.x1_mm{cameraNo}(1, :), coords_per_cam_corrected.x2_mm{cameraNo}(:,1), mean_U_per_cam{cameraNo});
    set(gca,'YDir','normal'); axis image; colormap(gca, jet);
    colorbar; clim([0 30]); hold on;
    plot(x_vec_pt, floor_y_C,             'g.',  'MarkerSize', 4);
    plot(x_vec_pt, polyval(p_C,x_vec_pt), 'g-',  'LineWidth', 2);
    xlabel('X (mm)'); ylabel('Y (mm)');
    title('Wall detections and linear fits');
    legend('Method C raw','Method C fit', 'Location','northeast');
    
    ax2 = subplot(3,1,2);
    imagesc(coords_per_cam_corrected.x1_mm{cameraNo}(1, :), coords_per_cam_corrected.x2_mm{cameraNo}(:,1),  mean_pivtools.mean_U_per_cam_rotated{cameraNo});
    set(gca,'YDir','normal'); axis image; colormap(gca, jet);
    colorbar; clim([0 30]); hold on;
    plot(x_vec_pt, floor_y_C,             'g.',  'MarkerSize', 4);
    plot(x_vec_pt, polyval(p_C,x_vec_pt), 'g-',  'LineWidth', 2);
    xlabel('X (mm)'); ylabel('Y (mm)');
    title('Rotated plane');
    legend('Method C raw','Method C fit', 'Location','northeast');
    ax3 = subplot(3, 1, 3);
    hold on;
    plot(x_vec_pt, C_fit, 'g.', 'MarkerSize', 4, ...
        'DisplayName', sprintf('Method C residual  (rise=%.3f mm, %.5f°)', pitch_C_mm, pitch_C_deg));
    plot(x_vec_pt, floor_y_C, 'bo', 'MarkerFaceColor','b', 'MarkerEdgeColor','b', ...
        'MarkerSize', 10, ...
        'DisplayName', sprintf('Method C Raw y  (rise=%.3f mm, %.5f°)', pitch_C_mm, pitch_C_deg));
    xlabel('X (mm)'); ylabel('Residual (mm)');
    grid on; hold off;
    linkaxes([ax1, ax2,  ax3], 'x');
    
    % --- apply floor shift via function ---
    % camDir = fullfile(savePath, totalLoops(loopNo).name, ...
    %     'calibrated_piv', '150', cameraList(cameraNo), 'instantaneous/');
    % PIV_applyFloorShift(fullfile(camDir, 'coordinates.mat'), p_C_per_cam(cameraNo, :)); 
end
save('pivtools_mean_UV_allcams.mat', 'mean_pivtools')

% % *** NEW: save all p_C coefficients alongside the loop data
% floorFitPath = fullfile(savePath, totalLoops(loopNo).name, 'floor_fits.mat');
% save(floorFitPath, 'p_C_per_cam', 'cameraList');
% fprintf('Saved floor fit coefficients → %s\n', floorFitPath);

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


function results = compute_integral_length_scales_hpc( resultFolder, params)
if nargin < 2 || isempty(params)
    params = default_integral_length_params();
else
    params = apply_default_params(params, default_integral_length_params());
end

assert(isfolder(resultFolder), 'Folder not found: %s', resultFolder);
assert(ismember(lower(params.marchDirection), {'dx','dy'}), 'params.marchDirection must be ''dx'' or ''dy''.');

gridFile = fullfile(resultFolder, params.gridFileName);
assert(exist(gridFile, 'file') == 2, 'grid.mat not found in %s', resultFolder);

Sg = load(gridFile, 'worldX_merged', 'worldY_merged');
assert(isfield(Sg, 'worldX_merged') && isfield(Sg, 'worldY_merged'), 'grid.mat must contain worldX_merged and worldY_merged');
worldX = double(Sg.worldX_merged);
worldY = double(Sg.worldY_merged);

files = dir(fullfile(resultFolder, params.filePattern));
assert(~isempty(files), 'No correlation files found matching %s', params.filePattern);
N = numel(files);

% if params.useParallel
%     pool = gcp('nocreate');
%     if isempty(pool)
%         if isempty(params.numWorkers)
%             parpool;
%         else
%             parpool(params.numWorkers);
%         end
%     end
% end

file_name      = strings(N,1);
corr_type      = strings(N,1);
xref_from_name = nan(N,1);
yref_from_name = nan(N,1);
xr_exact       = nan(N,1);
yr_exact       = nan(N,1);
L_int          = nan(N,1);
zero_crossing  = nan(N,1);
n_points_used  = nan(N,1);
line_index     = nan(N,1);
origin_index   = nan(N,1);
threshold_scales = nan(N, numel(params.thresholds));
status_code      = zeros(N,1);
status_message   = strings(N,1);

if params.saveDecayCurves
    curve_coord_raw = cell(N,1);
    curve_rho_raw   = cell(N,1);
    curve_coord_smooth = cell(N,1);
    curve_rho_smooth   = cell(N,1);
else
    curve_coord_raw = {};
    curve_rho_raw = {};
    curve_coord_smooth = {};
    curve_rho_smooth = {};
end

if params.useParallel
    parfor k = 1:N
        [file_name(k), corr_type(k), xref_from_name(k), yref_from_name(k), xr_exact(k), yr_exact(k), ...
            L_int(k), zero_crossing(k), n_points_used(k), line_index(k), origin_index(k), threshold_scales(k,:), ...
            status_code(k), status_message(k), ccr, rrr, ccs, rrs] = process_one_file(files(k), resultFolder, worldX, worldY, params);
        if params.saveDecayCurves
            curve_coord_raw{k} = ccr;
            curve_rho_raw{k} = rrr;
            curve_coord_smooth{k} = ccs;
            curve_rho_smooth{k} = rrs;
        end
    end
else
    for k = 1:N
        [file_name(k), corr_type(k), xref_from_name(k), yref_from_name(k), xr_exact(k), yr_exact(k), ...
            L_int(k), zero_crossing(k), n_points_used(k), line_index(k), origin_index(k), threshold_scales(k,:), ...
            status_code(k), status_message(k), ccr, rrr, ccs, rrs] = process_one_file(files(k), resultFolder, worldX, worldY, params);
        if params.saveDecayCurves
            curve_coord_raw{k} = ccr;
            curve_rho_raw{k} = rrr;
            curve_coord_smooth{k} = ccs;
            curve_rho_smooth{k} = rrs;
        end
    end
end

T = table(file_name, corr_type, xref_from_name, yref_from_name, xr_exact, yr_exact, L_int, zero_crossing, n_points_used, line_index, origin_index, status_code, status_message);
for j = 1:numel(params.thresholds)
    varName = matlab.lang.makeValidName(sprintf('threshold_%g', params.thresholds(j)));
    T.(varName) = threshold_scales(:,j);
end
[~, orderOut] = sortrows([yref_from_name, xref_from_name], [1 2]);
T = T(orderOut,:);

results = struct();
results.table = T;
results.params = params;
results.resultFolder = resultFolder;
results.created = datetime('now');
if params.saveDecayCurves
    results.curve_file_name = file_name;
    results.curve_coord_raw = curve_coord_raw;
    results.curve_rho_raw = curve_rho_raw;
    results.curve_coord_smooth = curve_coord_smooth;
    results.curve_rho_smooth = curve_rho_smooth;
end

summaryMat = fullfile(resultFolder, sprintf('integral_length_summary_%s.mat', lower(params.marchDirection)));
summaryCsv = fullfile(resultFolder, sprintf('integral_length_summary_%s.csv', lower(params.marchDirection)));
save(summaryMat, 'results', '-v7.3');
writetable(T, summaryCsv);

if params.saveDecayCurves
    decayMat = fullfile(resultFolder, sprintf('decay_curves_%s.mat', lower(params.marchDirection)));
    save(decayMat, 'file_name', 'curve_coord_raw', 'curve_rho_raw', 'curve_coord_smooth', 'curve_rho_smooth', 'params', '-v7.3');
end

end

function params = default_integral_length_params()
params = struct();
params.gridFileName = 'grid.mat';
params.filePattern = 'R_*_xref*_yref*.mat';
params.marchDirection = 'dx';
params.useParallel = false;
params.numWorkers = [];
params.saveDecayCurves = true;
params.useAbsSymmetry = true;
params.integratePositiveOnly = true;
params.thresholds = [exp(-1), 0.5, 0.1];
params.allowNoZeroCrossing = true;
params.verbose = true;
params.smoothCurve = false;
params.smoothMethod = 'sgolay';
params.sgolayOrder = 3;
params.sgolayFrameLength = 9;
end

function params = apply_default_params(params, defaults)
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(params, fn{i}) || isempty(params.(fn{i}))
        params.(fn{i}) = defaults.(fn{i});
    end
end
end

function [file_name, corr_type, xref_from_name, yref_from_name, xr_exact, yr_exact, L_int, zero_crossing, n_points_used, line_index, origin_index, threshold_row, status_code, status_message, curve_coord_raw, curve_rho_raw, curve_coord_smooth, curve_rho_smooth] = process_one_file(fileInfo, resultFolder, worldX, worldY, params)
file_name = string(fileInfo.name);
corr_type = ""; xref_from_name = nan; yref_from_name = nan; xr_exact = nan; yr_exact = nan; L_int = nan; zero_crossing = nan; n_points_used = nan; line_index = nan; origin_index = nan; threshold_row = nan(1, numel(params.thresholds)); status_code = -999; status_message = "";
curve_coord_raw = []; curve_rho_raw = []; curve_coord_smooth = []; curve_rho_smooth = [];
try
    tok = regexp(fileInfo.name, '^R_([A-Za-z0-9]+)_xref([-+]?\d*\.?\d+)_yref([-+]?\d*\.?\d+)\.mat$', 'tokens', 'once');
    if isempty(tok)
        status_code = -1; status_message = "Bad filename pattern"; return
    end
    corr_type = string(tok{1}); xref_from_name = str2double(tok{2}); yref_from_name = str2double(tok{3});
    S = load(fullfile(resultFolder, fileInfo.name), 'xr', 'yr', 'R_s');
    if ~isfield(S, 'xr') || ~isfield(S, 'yr') || ~isfield(S, 'R_s')
        status_code = -2; status_message = "Missing xr, yr, or R_s"; return
    end
    xr_exact = double(S.xr); yr_exact = double(S.yr); R = double(S.R_s);
    if ~isequal(size(R), size(worldX)) || ~isequal(size(R), size(worldY))
        status_code = -3; status_message = "Size mismatch with grid"; return
    end
    dX = worldX - xr_exact; dY = worldY - yr_exact;
    switch lower(params.marchDirection)
        case 'dx'
            [~, iy0] = min(abs(dY(:,1))); line_index = iy0; coord_line = dX(iy0, :); rho_line = R(iy0, :);
        case 'dy'
            [~, ix0] = min(abs(dX(1,:))); line_index = ix0; coord_line = dY(:, ix0); rho_line = R(:, ix0);
    end
    [coord_pos_raw, rho_pos_raw, idx0_local] = build_1d_curve(coord_line, rho_line, params.useAbsSymmetry);
    origin_index = idx0_local;
    if isempty(coord_pos_raw) || isempty(rho_pos_raw) || numel(coord_pos_raw) < 2
        status_code = -4; status_message = "Insufficient 1D curve"; return
    end
    [coord_pos_raw, order] = sort(coord_pos_raw(:)); rho_pos_raw = rho_pos_raw(order);
    valid = isfinite(coord_pos_raw) & isfinite(rho_pos_raw); coord_pos_raw = coord_pos_raw(valid); rho_pos_raw = rho_pos_raw(valid);
    if numel(coord_pos_raw) < 2
        status_code = -5; status_message = "Too few valid points"; return
    end
    [~, i0] = min(abs(coord_pos_raw));
    % minOriginRho = 0.01; 
    % if abs(rho_pos_raw(i0)) < minOriginRho
    %     status_code = -6;
    %     status_message = sprintf("Origin rho too small (%.4f < %.4f) — field likely uncorrelated", ...
    %         abs(rho_pos_raw(i0)), minOriginRho);
    %     return
    % end
    % rho_pos_raw = rho_pos_raw ./ rho_pos_raw(i0); rho_pos_raw(i0) = 1;
    [~, i0] = min(abs(coord_pos_raw));

% local fitting window around zero
fitHalfWidth = 5;   % points on each side; tune as needed
i1 = max(1, i0 - fitHalfWidth);
i2 = min(numel(coord_pos_raw), i0 + fitHalfWidth);

x_fit = coord_pos_raw(i1:i2);
y_fit = rho_pos_raw(i1:i2);

% remove any bad values
goodFit = isfinite(x_fit) & isfinite(y_fit);
x_fit = x_fit(goodFit);
y_fit = y_fit(goodFit);

% defaults in case fit fails
rho_peak = max(y_fit);
x_peak   = x_fit(y_fit == rho_peak);
x_peak   = x_peak(1);

% try local Gaussian peak fit
try
    % start point: [amplitude, center, width]
    a0 = max(y_fit);
    [~, imax] = max(y_fit);
    b0 = x_fit(imax);

    % rough width guess from local window
    c0 = max(range(x_fit)/4, eps);

    ft = fittype('a1*exp(-((x-b1)/c1)^2)', ...
        'independent', 'x', 'coefficients', {'a1','b1','c1'});

    opts = fitoptions(ft);
    opts.StartPoint = [a0, b0, c0];
    opts.Lower      = [0, min(x_fit), eps];
    opts.Upper      = [2*max(y_fit), max(x_fit), range(x_fit)];

    fobj = fit(x_fit, y_fit, ft, opts);

    rho_peak = fobj.a1;
    x_peak   = fobj.b1;

catch
    % fallback: local discrete maximum
end

% reject weak peaks
minOriginRho = 0.05;
if rho_peak < minOriginRho
    status_code = -6;
    status_message = sprintf("Peak rho too small (%.4f < %.4f) — field likely uncorrelated", ...
        rho_peak, minOriginRho);
    return
end

% normalize entire curve by fitted local peak amplitude
rho_pos_raw = rho_pos_raw ./ rho_peak;

% optional sanity clamp
rhoClamp = 1.5;
if max(rho_pos_raw) > rhoClamp
    status_code = -11;
    status_message = sprintf("Normalized rho exceeds %.1f (max = %.2f) — bad slice", ...
        rhoClamp, max(rho_pos_raw));
    return
end
    coord_pos_s = coord_pos_raw; rho_pos_s = rho_pos_raw;
    if params.smoothCurve
        switch lower(params.smoothMethod)
            case 'sgolay'
                frameLen = params.sgolayFrameLength; polyOrd = params.sgolayOrder;
                if mod(frameLen,2) == 0, frameLen = frameLen + 1; end
                if frameLen <= polyOrd, frameLen = polyOrd + 2 + mod(polyOrd+2,2); end
                if frameLen > numel(rho_pos_s)
                    frameLen = numel(rho_pos_s);
                    if mod(frameLen,2) == 0, frameLen = frameLen - 1; end
                end
                if frameLen >= 3 && polyOrd < frameLen
                    rho_pos_s = sgolayfilt(rho_pos_s, polyOrd, frameLen);
                    rho_pos_s(i0) = 1;
                end
            otherwise
                error('Unknown smoothMethod: %s', params.smoothMethod)
        end
    end
    curve_coord_raw = coord_pos_raw; curve_rho_raw = rho_pos_raw;
    curve_coord_smooth = coord_pos_s; curve_rho_smooth = rho_pos_s;
    coord_pos = coord_pos_s; rho_pos = rho_pos_s;
    maxReasonableL = max(abs(coord_pos));   % can't exceed the domain itself

    posMask = coord_pos >= 0; coord_pos = coord_pos(posMask); rho_pos = rho_pos(posMask);
    if numel(coord_pos) < 2
        status_code = -7; status_message = "Too few positive-separation points"; return
    end
    idxCross = find(rho_pos(1:end-1) >= 0 & rho_pos(2:end) < 0, 1, 'first');
    if ~isempty(idxCross)
        x1 = coord_pos(idxCross); x2 = coord_pos(idxCross+1); y1 = rho_pos(idxCross); y2 = rho_pos(idxCross+1);
        x0 = x1 - y1*(x2-x1)/(y2-y1);
        coord_int = [coord_pos(1:idxCross); x0]; rho_int = [rho_pos(1:idxCross); 0];
        zero_crossing = x0; L_int = trapz(coord_int, rho_int); n_points_used = numel(coord_int);
    else
        if params.allowNoZeroCrossing
            if params.integratePositiveOnly
                posLobe = rho_pos >= 0; coord_int = coord_pos(posLobe); rho_int = rho_pos(posLobe);
            else
                coord_int = coord_pos; rho_int = rho_pos;
            end
            if numel(coord_int) >= 2
                L_int = trapz(coord_int, rho_int); zero_crossing = coord_int(end); n_points_used = numel(coord_int);
                status_code    = 2;   % valid but no zero crossing found in domain
                status_message = "OK — no zero crossing in domain";
            else
                status_code = -8; status_message = "No zero crossing and too few points"; return
            end
        else
            status_code = -9; status_message = "No zero crossing"; return
        end
    end
    % ---- SANITY CHECK: L_int must not exceed the domain half-extent ----
    if L_int > maxReasonableL
        status_code = -10;
        status_message = sprintf("L_int (%.2f) exceeds domain extent (%.2f) — outlier rejected", ...
            L_int, maxReasonableL);
        L_int = NaN;
        return
    end
    for j = 1:numel(params.thresholds)
        thr = params.thresholds(j);
        idxThr = find(rho_pos(1:end-1) >= thr & rho_pos(2:end) < thr, 1, 'first');
        if ~isempty(idxThr)
            x1 = coord_pos(idxThr); x2 = coord_pos(idxThr+1); y1 = rho_pos(idxThr); y2 = rho_pos(idxThr+1);
            threshold_row(j) = x1 + (thr-y1)*(x2-x1)/(y2-y1);
        end
    end
    status_code = 1; status_message = "OK";
catch ME
    status_code = -100; status_message = string(ME.message);
end
end

function [coord_pos, rho_pos, idx0] = build_1d_curve(coord_line, rho_line, useAbsSymmetry)
coord_line = double(coord_line(:)); rho_line = double(rho_line(:)); [coord_line, ord] = sort(coord_line); rho_line = rho_line(ord); [~, idx0] = min(abs(coord_line));
if useAbsSymmetry
    posCoord = coord_line(coord_line >= 0); posRho = rho_line(coord_line >= 0);
    negCoord = -flipud(coord_line(coord_line <= 0)); negRho = flipud(rho_line(coord_line <= 0));
    if isempty(posCoord), coord_pos = []; rho_pos = []; return; end
    if isempty(negCoord), coord_pos = posCoord; rho_pos = posRho; return; end
    maxCommon = min(max(posCoord), max(negCoord)); commonMask = posCoord <= maxCommon; coord_pos = posCoord(commonMask); posRho = posRho(commonMask); negRhoInterp = interp1(negCoord, negRho, coord_pos, 'linear', 'extrap'); rho_pos = 0.5*(posRho + negRhoInterp);
else
    mask = coord_line >= 0; coord_pos = coord_line(mask); rho_pos = rho_line(mask);
end
end
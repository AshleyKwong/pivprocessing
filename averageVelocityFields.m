%% LOCAL FUNCTION: averageVelocityFields
% =========================================================================
function meanCameras = averageVelocityFields(savePath, totalLoops)
tic;
fprintf('\n=== Averaging Velocity Fields ===\n');

firstPath = fullfile(savePath, totalLoops(1).name, 'processedvelocityfields.mat');
firstData = load(firstPath, 'allCameras');
[nFrames, nCameras] = size(firstData.allCameras.u);
sumU = firstData.allCameras.u;
sumV = firstData.allCameras.v;
clear firstData;
fprintf('  Initialised from loop 1/%d\n', length(totalLoops));

for k = 2:length(totalLoops)
    nextPath = fullfile(savePath, totalLoops(k).name, 'processedvelocityfields.mat');
    nextData = load(nextPath, 'allCameras');
    for fr = 1:nFrames
        for cam = 1:nCameras
            sumU{fr,cam} = sumU{fr,cam} + nextData.allCameras.u{fr,cam};
            sumV{fr,cam} = sumV{fr,cam} + nextData.allCameras.v{fr,cam};
        end
    end
    clear nextData;
    fprintf('  Added loop %d/%d\n', k, length(totalLoops));
end

nL = length(totalLoops);
meanU = cell(1, nCameras);
meanV = cell(1, nCameras);
for cam = 1:nCameras
    stackU = cat(3, sumU{:,cam}) / nL;
    stackV = cat(3, sumV{:,cam}) / nL;
    meanU{cam} = mean(stackU, 3);
    meanV{cam} = mean(stackV, 3);
end

meanCameras.u = meanU;
meanCameras.v = meanV;

S = whos('meanCameras');
fprintf('✓ Averaged %d loops × %d frames  (%.2f GB) in %.1f s\n', ...
    nL, nFrames, S.bytes/1e9, toc);

tag      = sprintf('averagedvelfields_uv_%d', nL*150);
filename = fullfile(savePath, [tag '.mat']);
save(filename, 'meanCameras', '-v7.3');
fprintf('✓ Saved %s\n', tag);
end

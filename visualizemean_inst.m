clc; close all;
if ~exist('rms_field', 'var')
    % rms_field = load('D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\merge_instantaneousavg_20260403_011751\merged_turbrms_hann_14loops_20260403_011751.mat');
    % rms_field = load('C:\Users\ak1u24\Downloads\PG_fixedcal\loop=06\merge_instantaneousavg_20260410_001222\merged_turbrms_hann_1loops_20260410_001222.mat');
    % rms_field = load('C:\Users\ak1u24\Downloads\ZPG_fixedcal\merge_instantaneousavg_20260410_013235\merged_turbrms_hann_11loops_20260410_013235.mat');
    rms_field = load('C:\Users\ak1u24\Downloads\PG_fixedmergandcal\loop=06\merge_instantaneousavg_20260410_133714\merged_turbrms_hann_1loops_20260410_133714.mat');

end
if ~exist('mean_field', 'var')
    % mean_field = load('D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\merge_instantaneousavg_20260403_011751\merged_meanUV_14loops_20260403_011751.mat');
    % mean_field = load('C:\Users\ak1u24\Downloads\PG_fixedcal\loop=06\merge_instantaneousavg_20260410_001222\merged_meanUV_1loops_20260410_001222.mat');
    % mean_field = load('C:\Users\ak1u24\Downloads\ZPG_fixedcal\merge_instantaneousavg_20260410_013235\merged_meanUV_11loops_20260410_013235.mat');
    mean_field = load('C:\Users\ak1u24\Downloads\PG_fixedmergandcal\loop=06\merge_instantaneousavg_20260410_133714\merged_meanUV_1loops_20260410_133714.mat');


end
if ~exist('inst_frames', 'var')
    % inst_frames = load('D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\loop=06\mergedvelocityfields_20260403_011751.mat');
    % inst_frames = load('C:\Users\ak1u24\Downloads\PG_fixedcal\loop=06\mergedvelocityfields_20260410_001222.mat');
    % inst_frames = load('C:\Users\ak1u24\Downloads\ZPG_fixedcal\mergedvelocityfields_20260410_013235.mat');
    inst_frames = load('C:\Users\ak1u24\Downloads\PG_fixedmergandcal\loop=06\mergedvelocityfields_20260410_133714.mat');


end


nFrame = 130 ;

figure(1)
hold on; 
ax1 = subplot(2,1,1); 
imagesc(ax1 , mean_field.worldX(1,:), mean_field.worldY(:,1), mean_field.U_hann_mean)
ylabel( '$ \bar{U} [m/s] $', Interpreter='latex'); 
set(ax1, 'YDir', 'normal'); 
axis image; 

colormap(jet); clim([ 0 30]); colorbar(); 
ax2 = subplot(2,1,2); 
imagesc(ax2 , mean_field.worldX(1,:), mean_field.worldY(:,1), inst_frames.mergedFrames.U{nFrame}); 
set(ax2, 'YDir', 'normal'); 
axis image; 

colormap(jet); clim([ 0 30]); colorbar(); 
xlabel( '$ x [mm] $', Interpreter='latex'); 
ylabel( '$ u [m/s] $', Interpreter='latex'); 
hold off

figure(2)
hold on; 
ax1 = subplot(2,1,1); 
imagesc(ax1 , mean_field.worldX(1,:), mean_field.worldY(:,1), mean_field.V_hann_mean)
ylabel( '$ \bar{V} [m/s] $', Interpreter='latex'); 
set(ax1, 'YDir', 'normal'); 
axis image; 
colormap(jet); clim([ -2 5]); colorbar(); 
ax2 = subplot(2,1,2); 
imagesc(ax2 , mean_field.worldX(1,:), mean_field.worldY(:,1),inst_frames.mergedFrames.V{nFrame}); 
set(ax2, 'YDir', 'normal'); 
colormap(jet); clim([ -2 5]); colorbar(); 
xlabel( '$ x [mm] $', Interpreter='latex'); 
ylabel( '$ v [m/s] $', Interpreter='latex'); 
axis image; 
hold off


figure(3)
hold on;

% Auto clim for U_rms
u_rms_clean = rms_field.U_rms(~isnan(rms_field.U_rms));
clim_u = [prctile(u_rms_clean, 2), prctile(u_rms_clean, 98)];

% Auto clim for V_rms
v_rms_clean = rms_field.V_rms(~isnan(rms_field.V_rms));
clim_v = [prctile(v_rms_clean, 2), prctile(v_rms_clean, 98)];

ax1 = subplot(2,1,1);
imagesc(ax1, rms_field.worldX(1,:), rms_field.worldY(:,1), rms_field.U_rms);
clim(ax1, clim_u);
ylabel('$u_{rms}$ [m/s]', 'Interpreter', 'latex');
set(ax1, 'YDir', 'normal');
axis image;
colormap(ax1, jet); colorbar(ax1);

ax2 = subplot(2,1,2);
imagesc(ax2, rms_field.worldX(1,:), rms_field.worldY(:,1), rms_field.V_rms);
clim(ax2, clim_v);
set(ax2, 'YDir', 'normal');
colormap(ax2, jet); colorbar(ax2);
xlabel('$x$ [mm]', 'Interpreter', 'latex');
ylabel('$v_{rms}$ [m/s]', 'Interpreter', 'latex');
axis image;

hold off

%%
if ~exist('individual_frames', 'var')
    individual_frames = load('D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\loop=06\vel_fluctuations\fluctuations_all_frames.mat'); 
    % individual_frames = load('C:\Users\ak1u24\Downloads\ZPG_fixedcal\fluctuations_all_frames.mat'); 
end 
if ~exist('windowCenters', 'var')
    windowCenters = load('D:\FULLYPROCESSEDY235AOAN04AOAFN06PIVDATA\loop=06\windowCenterCameras_mm.mat'); 
    % windowCenters = load('C:\Users\ak1u24\Downloads\ZPG_fixedcal\windowCenterCameras_mm.mat');

end 
nFrame = 75; 
u_fluc = individual_frames.fluctuations.u_prime(nFrame, :); 
v_fluc = individual_frames.fluctuations.v_prime(nFrame, :); 

% [worldX, worldY, u_fluc_merged, v_fluc_merged] = merge_cameras_python_style_mean(...
%     windowCenters.windowCenterCameras_mm, u_fluc, v_fluc, {}, 'hann', []);
[worldX, worldY, u_fluc_merged, v_fluc_merged] = merge_cameras_python_style_mean(...
    windowCenters.windowCenterCameras_mm, u_fluc, v_fluc, {}, 'tukey', 0.5);
figure(1001)
hold on; 
imagesc(worldX(1,:), worldY(:,1), v_fluc_merged)
ylabel( '$ v'' [m/s] $', Interpreter='latex'); 
set(gca, 'YDir', 'normal'); 
axis image; 
colormap(jet); clim([ -2 5]); colorbar(); axis image; 
hold off
%% 

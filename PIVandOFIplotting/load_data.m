%% =========================================================================
%  load_data.m
%  PURPOSE : Centralised data loading module for Case 1 Virgilio/Preskett
%            comparison analysis.
%
%  USAGE   : data = load_data();
%            Returns a struct `data` containing all raw loaded datasets.
%
%  NOTES   :
%   - All addpath / load calls from VirgPreskettCase1_computauUU.m are
%     reproduced here with their original-script line references.
%   - Hardcoded paths are isolated to this ONE file. To run on a different
%     machine, only edit the cfg struct at the top.
%   - Loaded variables that were previously dumped into base workspace are
%     now namespaced into sub-structs (data.blSweep, data.VP, data.tunnel,
%     data.pressure).
%   - Variables requiring fluid property computation are also handled here
%     and returned in data.fluid.
%
%  AUTHOR  : AK (refactored from VirgPreskettCase1_computauUU.m)
%  DATE    : 2026-03-29
% =========================================================================

function data = load_data()

%% ---- PATH CONFIGURATION -------------------------------------------------
% Original lines 14-17: addpath calls
% Edit these paths for your machine / HPC environment.

cfg.paths.ofi_cf_dev  = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\Case 1 Streamwise Cf development\';
cfg.paths.matlab_root = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB';
cfg.paths.ofi_codes   = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\OFI\OFI Codes Tak';

addpath(cfg.paths.ofi_cf_dev);    % Original line 14
addpath(cfg.paths.matlab_root);   % Original line 15
addpath(cfg.paths.ofi_codes);     % Original line 16 — fluid_prop() lives here

%% ---- AK PIV SWEEP -------------------------------------------------------
% Original line 18: load('...blSweep_20260318_163411.mat')
% Was: dumped into base workspace as a flat struct `blSweep`
% Now: data.blSweep

blSweep_path = fullfile(cfg.paths.matlab_root, ...
    'Experimental Campaign 1\PIV\PIV_results\Case1_PIVresults\blSweep_20260318_163411.mat');
data.blSweep = load(blSweep_path);
% Access fields as: data.blSweep.U99, data.blSweep.x_mm, etc.
% (blSweep was itself a struct, so fields are at data.blSweep.blSweep.xxx)
% Unwrap one level for convenience:
if isfield(data.blSweep, 'blSweep')
    data.blSweep = data.blSweep.blSweep;
end

%% ---- VIRGILIO / PRESKETT 2025 DATA --------------------------------------
% Original lines 20-22: VirgilioPreskett_cp, Cf_minus8 (Cf_vec1 dumped raw)

vp_base = 'C:\Users\ak1u24\OneDrive - University of Southampton\Desktop\Preskett_Virgilio_Data\Data_Virgilio25\Dataset_for_Pressure_gradient_history_effects\Data OpenShare\';

% Cp mean (named load — original line 20)
data.VP.cp = load(fullfile(vp_base, 'Cp_mean.mat'));
% fields: data.VP.cp.tap_m, data.VP.cp.meight_500mm

% Cf at alpha=-8 (raw dump of Cf_vec1 — original line 21)
% Was: load(...Cf_minus8.mat) → Cf_vec1 dumped into workspace
% Then: VirgilioPreskett_cfmatch = mean(Cf_vec1, 2)  [original line 22]
tmp_cf = load(fullfile(vp_base, 'Cf_minus8.mat'));
data.VP.Cf_match = mean(tmp_cf.Cf_vec1, 2);   % replaces VirgilioPreskett_cfmatch
clear tmp_cf;

% Re_x at alpha=-8 — original line in %% Create combined figure section
data.VP.Re = load(fullfile(vp_base, 'Re_minus8.mat'));
% field: data.VP.Re.Re_vec1

% SW PIV summary — original load inside %% Create combined figure section
data.VP.SW_PIV = load(fullfile(vp_base, 'SW_PIV_Summary_-8.mat'));
% fields: data.VP.SW_PIV.X_m, .U99, .theta_m, .deltastar_m, .H, .delta_m

%% ---- PRESKETT PRESSURE DATA (THOMAS OLD WORK) ---------------------------
% Original line 24: load('...\case_pressuredata.mat')
% Was: dumped several variables into base workspace
% Now: data.tomPressure

tom_path = 'C:\Users\ak1u24\OneDrive - University of Southampton\Thomas Preskett Old Work\TOM MODEL SHI\Pressure Data - Upload\case_pressuredata.mat';
data.tomPressure = load(tom_path);

%% ---- TUNNEL CONDITIONS (OFI NIDAQ) --------------------------------------
% Original lines 25-27: tunnelConditions = load(...), Pdyn, T_atm, P_atm

tunnel_path = 'F:\LAB7 COMPUTER\OFI NIDAQ PRESSURES\AK_OFI\06072025\y250_aoan04_aoafn06\y250_aoan04_aoafn06_U20_secondplate_02.mat';
data.tunnel = load(tunnel_path);
% Access as: data.tunnel.qu, data.tunnel.T, data.tunnel.P0

%% ---- AK PRESSURE SWEEP DATA ---------------------------------------------
% Original load inside %% Create combined figure section:
% load('...pressuredata_wingcases_revisedNOZPGZERO.mat') → case_pdata dumped
pres_path = 'C:\Users\ak1u24\OneDrive - University of Southampton\MATLAB\Experimental Campaign 1\Pressure\pressuredata_wingcases_revisedNOZPGZERO.mat';
data.pressure = load(pres_path);
% Access as: data.pressure.case_pdata(2).xloc_testsection_xdelta, etc.

%% ---- FLUID PROPERTIES ---------------------------------------------------
% Original lines 28-34: P_atm, T_atm, fluid_prop(), visc_air override
% Now computed once and stored in data.fluid

Pdyn   = data.tunnel.qu;
T_atm  = mean(data.tunnel.T);
P_atm  = data.tunnel.P0 / 1013;   % convert to atm (original line 29: P0 in mbar)

[rho_air, ~, visc_oil, rho_w, g, ~] = fluid_prop(T_atm, P_atm);

data.fluid.rho_air  = rho_air;
data.fluid.visc_air = 1.5e-5;        % Original line 32: manual override
data.fluid.visc_oil = visc_oil;
data.fluid.rho_w    = rho_w;
data.fluid.g        = g;
data.fluid.T_atm    = T_atm;
data.fluid.P_atm    = P_atm;
data.fluid.Pdyn     = Pdyn;
data.fluid.U_inf    = sqrt(Pdyn * 2 / rho_air);   % original line 34

%% ---- VIRGILIO REFERENCE CONSTANTS ---------------------------------------
% Original lines 35-36: Virg_Uinf0, Virgilio_nuair

data.VP.Uinf0   = 19.7;       % [m/s]  — original line 35
data.VP.nu_air  = 1.51e-5;    % [m²/s] — original line 36

%% ---- GEOMETRY & CALIBRATION CONSTANTS -----------------------------------
% Original lines 12-13
data.cfg.mmperpix  = 1 ./ [15.95, 14.95];   % [mm/pix] for camera set 1 and 2
data.cfg.panel     = 0:1200:12000;

% OFI calibration file paths (used inside camera loops — original lines ~55, 70)
data.cfg.calib_path_cam1 = 'D:\ALK_OFI_250530_105600\Properties\Calibration\Calibration.xml';
data.cfg.calib_path_cam2 = 'D:\ALK_OFI_06072025\Properties\Calibration\Calibration.xml';

% Pixel offsets per camera set (original lines ~58, 73)
data.cfg.ppy0      = 0;
data.cfg.offset_c1 = -120;    % mm offset for camera set 1
data.cfg.offset_c2 =  1400;   % mm offset for camera set 2

% Global x origin for coordinate transform (original line after %% block 2)
data.cfg.x0_global_mm = -7300;

% Wing leading/trailing edge positions [mm] (used in multiple plot sections)
data.cfg.wing_LE_mm = 7650;
data.cfg.wing_TE_mm = 7950;

% Chord lengths [m]
data.cfg.c_ak_m  = 0.30;    % AK wing chord
data.cfg.c_VP_m  = 1.25;    % Virgilio wing chord

% Virgilio wing LE position [m]
data.cfg.VP_LE_m = 6.53;

% % Reference BL thickness for Virgilio [m] (used in dCpdx normalisation)
% data.cfg.VP_BL0_m = 0.08;

% Maximum spatial gap for nearest-neighbour U_inf matching [mm]
data.cfg.max_gap_mm = 50;

fprintf('[load_data] All data loaded successfully.\n');
end

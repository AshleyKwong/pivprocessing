% =========================================================================
% compute_uncertainty_forcompute_quad_spectra_skewkurt.m
%
% Post-processes turbulence statistics to compute uncertainty estimates.
%
% Quantities:
%   - u'rms uncertainty          from var(u'^2) and N_eff
%   - v'rms uncertainty          from var(v'^2) and N_eff
%   - Reynolds shear stress      from var(u'v') and N_eff
%   - Quadrant fraction          from binomial statistics
%   - Spectral confidence bounds from N_frames
%
% N_eff accounts for temporal correlation between frames using the
% integral time scale approach. If frames are statistically independent
% (e.g. low repetition rate PIV) then N_eff = N_frames.
%
% Author:  Ashley Kwong
% Date:    April 2026
% =========================================================================

clear; clc;

%% ===================== USER INPUTS =====================

turbStatsFile = '/iridisfs/scratch/ak1u24/case1_fullmask/turbulence_statistics_PIV/turbstats_case2.mat';
outDir        = '/iridisfs/scratch/ak1u24/case1_fullmask/turbulence_statistics_PIV/';
caseLabel     = 'case2';

% Confidence level
confLevel = 0.95;   % 95%

% Effective sample size factor
% If your PIV acquisition rate is low enough that frames are independent:
%   N_eff = totalFrames
% If frames are correlated (high rep rate), set this < 1 to reduce N_eff
% A conservative estimate: N_eff = totalFrames * dt / T_integral
% where T_integral is the integral time scale of the flow.
% Set to 1.0 if unsure — gives optimistic (lower) uncertainty bounds.
N_eff_factor = 1.0;

%% ===================== LOAD =====================

S      = load(turbStatsFile, 'results', 'totalFrames', 'xref_targets');
res    = S.results;
N      = S.totalFrames;
N_eff  = round(N * N_eff_factor);

fprintf('Total frames: %d\n', N);
fprintf('N_eff:        %d  (factor = %.2f)\n', N_eff, N_eff_factor);
fprintf('Confidence:   %.0f%%\n\n', confLevel*100);

% z-score for confidence level
z = norminv(1 - (1-confLevel)/2);   % 1.96 for 95%

nXref = numel(res);

%% ===================== COMPUTE UNCERTAINTY =====================

unc = struct();

for iX = 1:nXref

    if ~res(iX).valid
        unc(iX).xr    = res(iX).xr;
        unc(iX).valid = false;
        continue;
    end

    xr   = res(iX).xr;
    N_e  = N_eff;

    fprintf('xref = %.1f mm\n', xr);

    % --- Convenience ---
    urms = res(iX).urms;
    vrms = res(iX).vrms;
    u2   = urms.^2;          % <u'^2>
    v2   = vrms.^2;          % <v'^2>
    uv   = -res(iX).RS;      % <u'v'> (note RS = -<u'v'>)
    uv2  = res(iX).uv2;      % <(u'v')^2>
    Q1   = res(iX).Q1;
    Q2   = res(iX).Q2;
    Q3   = res(iX).Q3;
    Q4   = res(iX).Q4;

    % -------------------------------------------------------
    % u'rms uncertainty
    % var(u'^2) = <u'^4> - <u'^2>^2
    % We have <u'^4> from kurtosis: kurt = <u'^4>/<u'^2>^2
    % So <u'^4> = kurt * <u'^2>^2
    % -------------------------------------------------------
    u4       = res(iX).kurt .* u2.^2;          % <u'^4>
    var_u2   = max(u4 - u2.^2, 0);             % var(u'^2), floor at 0
    se_u2    = sqrt(var_u2 / N_e);             % std error of <u'^2>

    % Propagate to urms = sqrt(<u'^2>):
    % d(urms)/d(<u'^2>) = 1/(2*urms)
    safe_urms      = max(urms, eps);
    se_urms        = se_u2 ./ (2 * safe_urms);
    ci_urms        = z * se_urms;

    % -------------------------------------------------------
    % v'rms uncertainty — same approach
    % -------------------------------------------------------
    v4       = res(iX).kurt .* v2.^2;          % approximate using u kurt
    % Better: if you have v kurtosis separately use that
    var_v2   = max(v4 - v2.^2, 0);
    se_v2    = sqrt(var_v2 / N_e);
    safe_vrms      = max(vrms, eps);
    se_vrms        = se_v2 ./ (2 * safe_vrms);
    ci_vrms        = z * se_vrms;

    % -------------------------------------------------------
    % Reynolds shear stress uncertainty
    % var(u'v') = <(u'v')^2> - <u'v'>^2
    % -------------------------------------------------------
    var_uv   = max(uv2 - uv.^2, 0);
    se_RS    = sqrt(var_uv / N_e);
    ci_RS    = z * se_RS;

    % -------------------------------------------------------
    % Quadrant fraction uncertainty — binomial
    % var(Q_i) = p*(1-p)/N  where p is the fraction
    % -------------------------------------------------------
    ci_Q1 = z * sqrt(Q1.*(1-Q1) / N_e);
    ci_Q2 = z * sqrt(Q2.*(1-Q2) / N_e);
    ci_Q3 = z * sqrt(Q3.*(1-Q3) / N_e);
    ci_Q4 = z * sqrt(Q4.*(1-Q4) / N_e);

    % -------------------------------------------------------
    % Spectral uncertainty
    % For a periodogram averaged over N_e realisations,
    % the spectral estimate follows a chi-squared distribution
    % with 2*N_e degrees of freedom.
    % 95% CI: [Phi * 2*N_e/chi2inv(0.975, 2*N_e),
    %          Phi * 2*N_e/chi2inv(0.025, 2*N_e)]
    % Approximation for large N_e: CI ~ Phi * (1 +/- z/sqrt(N_e))
    % -------------------------------------------------------
    dof          = 2 * N_e;
    chi2_lo      = chi2inv((1-confLevel)/2,     dof);
    chi2_hi      = chi2inv(1-(1-confLevel)/2,   dof);
    Puu          = res(iX).kx_Puu;
    spec_ci_lo   = Puu * dof ./ chi2_hi;   % lower bound
    spec_ci_hi   = Puu * dof ./ chi2_lo;   % upper bound

    % Premultiplied bounds (skip DC)
    kx              = res(iX).kx;
    spec_ci_lo_plot = spec_ci_lo(:, 2:end) .* kx(2:end);
    spec_ci_hi_plot = spec_ci_hi(:, 2:end) .* kx(2:end);

    % -------------------------------------------------------
    % Store
    % -------------------------------------------------------
    unc(iX).xr           = xr;
    unc(iX).valid        = true;
    unc(iX).y_mm         = res(iX).y_mm;

    % Half-widths of confidence intervals (+ and -)
    unc(iX).ci_urms      = ci_urms;
    unc(iX).ci_vrms      = ci_vrms;
    unc(iX).ci_RS        = ci_RS;
    unc(iX).ci_Q1        = ci_Q1;
    unc(iX).ci_Q2        = ci_Q2;
    unc(iX).ci_Q3        = ci_Q3;
    unc(iX).ci_Q4        = ci_Q4;

    % Spectral bounds (full arrays)
    unc(iX).spec_ci_lo      = spec_ci_lo;
    unc(iX).spec_ci_hi      = spec_ci_hi;
    unc(iX).spec_ci_lo_plot = spec_ci_lo_plot;
    unc(iX).spec_ci_hi_plot = spec_ci_hi_plot;

    % Also store relative uncertainties (%) for quick inspection
    unc(iX).rel_urms_pct = 100 * ci_urms ./ max(urms, eps);
    unc(iX).rel_RS_pct   = 100 * ci_RS   ./ max(abs(res(iX).RS), eps);

    fprintf('  urms CI:  mean ±%.2f%%\n', mean(unc(iX).rel_urms_pct, 'omitnan'));
    fprintf('  RS    CI: mean ±%.2f%%\n', mean(unc(iX).rel_RS_pct,   'omitnan'));

end

%% ===================== SAVE =====================

saveFile = fullfile(outDir, sprintf('uncertainty_%s.mat', caseLabel));
save(saveFile, 'unc', 'xref_targets', 'N_eff', 'N_eff_factor', ...
    'confLevel', 'z', 'caseLabel', '-v7.3');

fprintf('\nSaved -> %s\n', saveFile);
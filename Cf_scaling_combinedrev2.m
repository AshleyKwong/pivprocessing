% =========================================================================
% Cf_scaling_combined.m
%
% Generalised multi-case Cf scaling pipeline.
%
% Pipeline:
%   1. Define N cases in the `cases` struct array (USER CONFIG)
%   2. ΔCp matching  — case 1 is the reference; all other cases have their
%      PG end matched to case 1's ΔCp
%   3. Compute Cf normalisations and Δx scalings for every case
%   4. Plot:
%        Figure A — Cp verification (one panel per case pair)
%        Figure B — Cf/Cf0 vs Δx/L_PG  (left)  and  Δx/δ0  (right)
%
% Adding a new case
% -----------------
%   Copy the last `cases(N)` block in USER CONFIG, increment the index,
%   and fill in the fields.  Nothing else needs changing.
%
% Workspace prerequisites — must exist before running
% ---------------------------------------------------
%   Each case struct references workspace variables by name (string fields).
%   The helper `ws` retrieves them via evalin('caller',...).
%   See field descriptions in USER CONFIG.
% =========================================================================

close all;

% =========================================================================
%  USER CONFIG
%  Each cases(i) describes one experimental dataset.
%
%  Required fields
%  ---------------
%  label        string   display name in legends
%  color        [R G B]  line/marker colour
%
%  Cp data (workspace variable names as strings)
%  ------
%  cp_x_expr    expression that evaluates to x-positions of Cp taps (m)
%  cp_val_expr  expression that evaluates to Cp values
%  cp_rm_idx    indices to remove (bad taps); [] for none
%
%  PG bounds (metres)
%  --
%  pg_x_start   x of PG start (used as Cf0 / δ0 reference)
%  pg_x_end     x of PG end   (only used for the REFERENCE case, i=1;
%               all other cases get their end from ΔCp matching)
%
%  Cf data (workspace variable names as strings)
%  ---
%  cf_expr          Cf array (local Uinf scaling, already movmean'd etc.)
%  uinf_local_expr  local Uinf array (same length as cf_expr)
%  uinf_outlet_expr scalar or expression for outlet Uinf
%  uinf_inlet_expr  scalar or expression for inlet Uinf
%  x_m_expr         x-positions of Cf data (m)
%
%  PIV boundary-layer data
%  -----------------------
%  piv_x_expr       x-positions of PIV sweep (m)
%  piv_d99_expr     δ99 at those positions (m)
%
%  Geometry
%  --------
%  chord            chord length c (m)
%  wing_le_x        wing LE x-position (m)
%
%  Smoothing
%  ---------
%  sg_win           Savitzky-Golay window (odd int ≥ 5); [] = no smoothing
% =========================================================================

% ---- Case 1 : AK (reference — ΔCp is taken from this case) -------------
cases(1).label          = 'AK';
cases(1).color          = [0.85 0.15 0.15];

cases(1).cp_x_expr      = '(case_pdata(2).xloc_testsection_xdelta .* case_pdata(2).refBLmm) ./ 100 + 7.65';
cases(1).cp_val_expr    = 'mean(case_pdata(2).Cp_testsection, 2)';
cases(1).cp_rm_idx      = 22;

cases(1).pg_x_start     = 7.44;
cases(1).pg_x_end       = 8.28;

cases(1).cf_expr        = '2 .* Cf2_movmean';
cases(1).uinf_local_expr  = 'U_inf_local(sort_idx)';
cases(1).uinf_outlet_expr = 'mean(blSweep.U_inf(end-10:end), ''omitnan'')';
cases(1).uinf_inlet_expr  = 'sqrt(Pdyn*2/rho_air)';
cases(1).x_m_expr       = 'x_AK_m';

cases(1).piv_x_expr     = 'blSweep.x_mm ./ 1000 + 7.1';
cases(1).piv_d99_expr   = 'blSweep.delta99_mm(~isnan(blSweep.delta99_mm)) ./ 1000';

cases(1).chord          = 0.3;
cases(1).wing_le_x      = 7.65;
cases(1).sg_win         = 5;

% ---- Case 2 : VP (end matched via ΔCp) ----------------------------------
cases(2).label          = 'Virgilio';
cases(2).color          = [0.15 0.35 0.85];

cases(2).cp_x_expr      = 'Virg_Preskett_xm';
cases(2).cp_val_expr    = 'Virg_Preskett_Cp';
cases(2).cp_rm_idx      = [];

cases(2).pg_x_start     = 6.12;
cases(2).pg_x_end       = NaN;   % filled by ΔCp matching

cases(2).cf_expr        = '2 .* filterVirgilioCf_localUinf';
cases(2).uinf_local_expr  = 'Virg_Uinflocal_PIV';
cases(2).uinf_outlet_expr = 'mean(SW_PIVData.U99(end-10:end)./0.99, ''omitnan'')';
cases(2).uinf_inlet_expr  = '19.7';
cases(2).x_m_expr       = 'Re_x_Virg_abs .* Virgilio_nuair ./ Virg_Uinf0';

cases(2).piv_x_expr     = 'SW_PIVData.X_m';
cases(2).piv_d99_expr   = 'SW_PIVData.delta_m';

cases(2).chord          = 1.25;
cases(2).wing_le_x      = 6.53;
cases(2).sg_win         = 5;



% ---- Case 3 : add more cases here ---------------------------------------
% cases(3).label = 'NewCase'; ...

% =========================================================================
%  END USER CONFIG
% =========================================================================

N = numel(cases);

%% ---- Step 1: load Cp and Cf arrays from workspace ----------------------
for i = 1:N
    c = cases(i);

    % Cp
    cp_x   = evalin('base', c.cp_x_expr);
    cp_val = evalin('base', c.cp_val_expr);
    if ~isempty(c.cp_rm_idx)
        cp_x(c.cp_rm_idx)   = [];
        cp_val(c.cp_rm_idx) = [];
    end
    cases(i).cp_x   = cp_x(:);
    cases(i).cp_val = cp_val(:);

    % Cf and velocities
    cases(i).Cf_local   = evalin('base', c.cf_expr);
    cases(i).Uinf_local = evalin('base', c.uinf_local_expr);
    cases(i).Uinf_out   = evalin('base', c.uinf_outlet_expr);
    cases(i).Uinf_in    = evalin('base', c.uinf_inlet_expr);
    cases(i).x_m        = evalin('base', c.x_m_expr);

    % PIV δ99
    % piv_x will never contain NaNs — derive the valid mask from δ99 only,
    % then index piv_x with the same mask (handles the case where piv_x is
    % longer than piv_d99 because δ99 was computed on a NaN-stripped subset)
    piv_d99_raw = evalin('base', c.piv_d99_expr);
    piv_x_raw   = evalin('base', c.piv_x_expr);
    valid        = ~isnan(piv_d99_raw);
    cases(i).piv_d99 = piv_d99_raw(valid);
    cases(i).piv_x   = piv_x_raw(valid);
end

%% ---- Step 2: ΔCp matching (case 1 = reference) -------------------------

% Reference ΔCp from case 1
ref  = cases(1);
[~, i0_ref] = min(abs(ref.cp_x - ref.pg_x_start));
[~, i1_ref] = min(abs(ref.cp_x - ref.pg_x_end));
delta_Cp_ref = ref.cp_val(i1_ref) - ref.cp_val(i0_ref);
cases(1).pg_x_end_matched = ref.pg_x_end;
cases(1).idx_cp_start     = i0_ref;
cases(1).idx_cp_end       = i1_ref;

fprintf('\n--- Pressure-gradient length summary ---\n');
fprintf('[Ref] %s:  x_start=%.4f  x_end=%.4f  L_PG=%.4f m  ΔCp=%.4f\n', ...
        ref.label, ref.pg_x_start, ref.pg_x_end, diff([ref.pg_x_start ref.pg_x_end]), delta_Cp_ref);

for i = 2:N
    cp_x   = cases(i).cp_x;
    cp_val = cases(i).cp_val;
    [~, i0] = min(abs(cp_x - cases(i).pg_x_start));
    Cp0     = cp_val(i0);
    Cp_tgt  = Cp0 + delta_Cp_ref;

    % Search only after start index
    [~, rel] = min(abs(cp_val(i0:end) - Cp_tgt));
    i1 = i0 + rel - 1;

    cases(i).pg_x_end_matched = cp_x(i1);
    cases(i).idx_cp_start     = i0;
    cases(i).idx_cp_end       = i1;
    cases(i).Cp_target        = Cp_tgt;

    L_PG_i = cp_x(i1) - cases(i).pg_x_start;
    fprintf('[Mtch] %s:  x_start=%.4f  x_end=%.4f  L_PG=%.4f m  Cp_tgt=%.4f\n', ...
            cases(i).label, cases(i).pg_x_start, cp_x(i1), L_PG_i, Cp_tgt);
end

%% ---- Step 3: compute δ0, L_PG, Cf normalisations, and Δx axes ---------

for i = 1:N
    c = cases(i);
    L_pgxm = [c.pg_x_start, c.pg_x_end_matched];
    L_PG   = diff(L_pgxm);

    % Incoming δ0 from PIV sweep at PG start
    [~, id0] = min(abs(c.piv_x - c.pg_x_start));
    d0 = c.piv_d99(id0);   % metres

    % Cf normalisations — normalise by value at PG start index
    [~, icf0] = min(abs(c.x_m - c.pg_x_start));

    Cf_out = c.Cf_local .* (c.Uinf_local ./ c.Uinf_out).^2;
    Cf_in  = c.Cf_local .* (c.Uinf_local ./ c.Uinf_in ).^2;
    Cf_loc = c.Cf_local;

    cases(i).Cf_outlet_norm = Cf_out ./ Cf_out(icf0);
    cases(i).Cf_inlet_norm  = Cf_in  ./ Cf_in(icf0);
    cases(i).Cf_local_norm  = Cf_loc ./ Cf_loc(icf0);

    % Δx axes
    dx = c.x_m - c.pg_x_start;
    cases(i).dx_lpg  = dx ./ L_PG;
    cases(i).dx_d0   = dx ./ d0;
    cases(i).dx_chord= dx ./ c.chord;
    cases(i).lpg_coords = [0 1]; 
    cases(i).lpg_coords_d0 = [0 L_PG]./ d0; 

    % Wing patch extents
    dw = c.wing_le_x - c.pg_x_start;
    cases(i).wing_lpg = [dw/L_PG,  (dw+c.chord)/L_PG ];
    cases(i).wing_d0  = [dw/d0,    (dw+c.chord)/d0   ];


    % Store for annotations
    cases(i).L_PG = L_PG;
    cases(i).d0   = d0;
    cases(i).ratio_lpg_to_d0 = L_PG / d0;   % how many δ0 per L_PG
end

%% ====================================================================
%  FIGURE A — Cp distributions + matched end-point verification
%  ====================================================================
figA = figure('Name','Cp — ΔCp matching verification', ...
              'Units','normalized','Position',[0.05 0.55 0.90 0.38]);

for i = 1:N
    ax = subplot(1, N, i, 'Parent', figA);
    hold(ax,'on');

    % Full Cp curve
    plot(ax, cases(i).cp_x, cases(i).cp_val, '-o', ...
        'Color', cases(i).color, 'MarkerFaceColor', cases(i).color, ...
        'MarkerSize', 4, 'LineWidth', 1.4, 'DisplayName', [cases(i).label, ' $C_p$']);

    % Start marker
    i0 = cases(i).idx_cp_start;
    plot(ax, cases(i).cp_x(i0), cases(i).cp_val(i0), 'v', ...
        'Color', cases(i).color, 'MarkerFaceColor', cases(i).color, ...
        'MarkerSize', 9, 'DisplayName', 'PG start');

    % End marker
    i1 = cases(i).idx_cp_end;
    plot(ax, cases(i).cp_x(i1), cases(i).cp_val(i1), 's', ...
        'Color', cases(i).color, 'MarkerFaceColor', cases(i).color, ...
        'MarkerSize', 9, 'DisplayName', 'PG end');

    if i > 1
        yline(ax, cases(i).Cp_target, '--', 'Color', cases(i).color, ...
            'LineWidth', 1.0, 'Label', '$C_{p,\mathrm{target}}$', ...
            'Interpreter','latex','HandleVisibility','off');
    end

    xlabel(ax, '$x$ (m)', 'Interpreter','latex');
    ylabel(ax, '$C_p$',   'Interpreter','latex');
    title(ax, cases(i).label, 'Interpreter','latex');
    legend(ax, 'Location','best','Interpreter','latex','FontSize',8);
    grid(ax,'on'); box(ax,'on');
    hold(ax,'off');
end
sgtitle(figA, '$C_p$ - $\Delta C_p$ end-point matching verification', ...
        'Interpreter','latex','FontSize',11);

%% ====================================================================
%  FIGURE B — Cf/Cf0:  left = Δx/L_PG,  right = Δx/δ0
%  ====================================================================
figB = figure('Name','Cf scaling', ...
              'Units','normalized','Position',[0.05 0.05 0.85 0.40]);

ax_lpg = subplot(1,2,1,'Parent',figB); hold(ax_lpg,'on');
ax_d0  = subplot(1,2,2,'Parent',figB); hold(ax_d0, 'on');

for i = 1:N
    c   = cases(i);
    sw  = c.sg_win;
    use_sg = ~isempty(sw);
    if use_sg
        if mod(sw,2)==0, sw = sw+1; end
        Cf_plt = sgolayfilt(c.Cf_outlet_norm(:), 3, sw);
    else
        Cf_plt = c.Cf_outlet_norm(:);
    end

    pale = c.color * 0.45 + 0.55;

    % Raw scatter underneath (if smoothed)
    if use_sg
        plot(ax_lpg, c.dx_lpg, c.Cf_outlet_norm, 'o', ...
            'Color',pale,'MarkerFaceColor',pale,'MarkerSize',3, ...
            'HandleVisibility','off');
        plot(ax_d0,  c.dx_d0,  c.Cf_outlet_norm, 'o', ...
            'Color',pale,'MarkerFaceColor',pale,'MarkerSize',3, ...
            'HandleVisibility','off');
    end

    % Ratio annotation for legend
    ratio_str = sprintf('%s: $L_{\\mathrm{PG}} = %.1f\\,\\delta_0$', ...
                        c.label, c.ratio_lpg_to_d0);

    plot(ax_lpg, c.dx_lpg, Cf_plt, '-o', ...
        'Color',c.color,'MarkerFaceColor',c.color,'MarkerSize',4, ...
        'LineWidth',1.5,'DisplayName', ratio_str);
    plot(ax_d0,  c.dx_d0,  Cf_plt, '-o', ...
        'Color',c.color,'MarkerFaceColor',c.color,'MarkerSize',4, ...
        'LineWidth',1.5,'DisplayName', ratio_str);

    % xline(ax_d0 ,(L_pgvp_xm(2)- L_pgvp_xm(1))/BL_VP_incoming, 'Color',c.color, 'DisplayName' ,'Case 1 PG end');
    % xline(ax_d0 , (L_pgak_xm(2)-L_pgak_xm(1))/BL_ak_incoming, 'Color',c.color, 'DisplayName' , 'VP PG end');
end

% Wing patches
for i = 1:N
    % add_wing_patch(ax_lpg, cases(i).wing_lpg, cases(i).color);
    % add_wing_patch(ax_d0,  cases(i).wing_d0,  cases(i).color);

    add_wing_patch(ax_lpg, cases(i).lpg_coords, cases(i).color);
    add_wing_patch(ax_d0,  cases(i).lpg_coords_d0,  cases(i).color);

end

xlabel(ax_lpg, '$\Delta x / L_{\mathrm{PG}}$', 'Interpreter','latex');
ylabel(ax_lpg, '$C_f / C_{f,0}$',              'Interpreter','latex');
grid(ax_lpg,'on'); box(ax_lpg,'on');
legend(ax_lpg,'Location','best','Interpreter','latex','FontSize',9);
hold(ax_lpg,'off');

% xline(ax_d0 ,0, "k--", 'DisplayName' ,'PG start ='); 
 

xlabel(ax_d0, '$\Delta x / \delta_0$',  'Interpreter','latex');
ylabel(ax_d0, '$C_f / C_{f,0}$',        'Interpreter','latex');
grid(ax_d0,'on'); box(ax_d0,'on');
legend(ax_d0,'Location','best','Interpreter','latex','FontSize',9);
hold(ax_d0,'off');

%% ====================================================================
%  LOCAL FUNCTION — wing shading
%  ====================================================================
function add_wing_patch(ax, wing_ext, col)
    yl = ylim(ax);
    patch(ax, [wing_ext(1) wing_ext(2) wing_ext(2) wing_ext(1)], ...
              [yl(1)       yl(1)       yl(2)       yl(2)      ], col, ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    uistack(findobj(ax,'Type','patch'),'bottom');
end

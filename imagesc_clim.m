% ── Helper: normalised 2-D autocorrelation ──────────────────────────────────

% After imagesc(imgA) calls, replace with this helper:
function h = imagesc_clim(ax, img, lo_pct, hi_pct)
% Displays img on axes ax with clim set by percentiles lo_pct / hi_pct
% Typical values: lo_pct = 1, hi_pct = 99
flat  = double(img(:));
clo   = prctile(flat, lo_pct);
chi   = prctile(flat, hi_pct);
h     = imagesc(ax, img);
clim(ax, [clo, chi]);
end


"""
report.py -- figure export (PDF / PNG) and the one-page per-cluster report (Phase 2.5).

TERMINOLOGY
  Figure helper : an existing function that RETURNS a matplotlib Figure without touching pyplot / show() (gmm.make_gmm_figures, imf.plot_imf,
                  extinction_isochrone.plot_cmd). This module only adds one new helper, msr_figure (no MSR plot existed).
  Bundle        : plain dict with the objects a report needs (see collect_figures); filled by cli.py from the finished session.
  One-page report: A4 PDF = header + key numbers + six figure slots + warnings. The six panels are the helper figures rasterised at REPORT_DPI
                  (the helpers draw whole Figures, not Axes, so they cannot be vector-composed); the stand-alone PDFs stay vector.

WORKFLOW
  1) collect_figures(bundle): call each existing helper that has its inputs, skip (with a reason) the ones that do not;
  2) export_figure(): save every figure as PDF + PNG;
  3) one_page_report(): rasterise the figures into a 3 x 2 grid, add the numbers and every warning (isochrone source, lower-limit mass, ...), save PDF.
"""
from __future__ import annotations

import io
import os

import matplotlib.image as mpimg
import numpy as np
from matplotlib.figure import Figure

from . import extinction_isochrone as EXT
from . import gmm as GMM
from . import imf as IMF

# ===================== GLOBAL CONFIGURATION =====================
FIGURE_FORMATS = ('pdf', 'png')           # formats written for every figure
FIGURE_DPI = 200                          # png resolution of stand-alone figures
REPORT_DPI = 130                          # resolution of the rasterised panels inside the one-page report
PAGE_SIZE_IN = (8.27, 11.69)              # A4 portrait
PANEL_ORDER = ('pm', 'plx', 'cmd', 'imf', 'msr', 'cmd_member')   # 3 rows x 2 columns, left to right
PANEL_TITLES = {'pm': 'Proper motion', 'plx': 'Parallax', 'cmd': 'CMD + isochrone', 'imf': 'Mass function', 'msr': 'Mass segregation', 'cmd_member': 'Member CMD (observed)'}
SUMMARY_FIELDS = (('N members', 'N', '{:d}'), ('d (pc)', 'd', '{:.1f}'), ('E(B-V)', 'E(B-V)', '{:.3f}'), ('log age', 'logAge', '{:.2f}'),
                  ('Gamma', 'Gamma', '{:+.2f}'), ('sigma_Gamma', 'sigma_Gamma', '{:.2f}'), ('MSR', 'MSR', '{:.2f}'), ('v_out (km/s)', 'v_out', '{:+.2f}'),
                  ('sigma_1D (km/s)', 'sigma_1D', '{:.2f}'))
WARNING_COLOUR = '#b00020'
# =================================================================


def msr_figure(msr_records: list, title: str = '') -> Figure:
    """Mass segregation ratio Lambda_MSR vs number of heaviest stars (Allison et al. 2009). Lambda = 1: no segregation. Error bars: 68 % of random draws."""
    f = Figure(figsize=(7, 5)); ax = f.subplots()
    r = [x for x in msr_records if np.isfinite(x.get('Lambda', np.nan))]
    n = np.array([x['N'] for x in r], float); lam = np.array([x['Lambda'] for x in r])
    lo = np.array([x['Lambda_lo'] for x in r]); hi = np.array([x['Lambda_hi'] for x in r])
    ax.errorbar(n, lam, yerr=[lam - lo, hi - lam], fmt='o-', color='tab:blue', capsize=4)
    ax.axhline(1.0, color='k', ls=':', lw=1)
    ax.set(xscale='log', xlabel='N heaviest stars', ylabel=r'$\Lambda_{\rm MSR}$', title=title or 'Mass segregation ratio')
    return f


def collect_figures(bundle: dict) -> tuple[dict, dict]:
    """Run every existing figure helper whose inputs are present. Returns (figures, skipped) with skipped = {panel: reason}."""
    figs, skipped = {}, {}
    name = bundle.get('name', '')
    res = bundle.get('gmm_res')
    if res is not None:
        g = GMM.make_gmm_figures(res, name)
        figs['pm'], figs['plx'], figs['cmd_member'] = g['pm'], g['plx'], g['cmd_member']
    else:
        skipped.update(pm='GMM did not run', plx='GMM did not run', cmd_member='GMM did not run')
    df, iso = bundle.get('analysis_df'), bundle.get('iso')
    if df is not None and iso is not None:
        figs['cmd'] = EXT.plot_cmd(df, iso, f'{name} CMD', reddened=bool(bundle.get('reddened')))
    else:
        skipped['cmd'] = 'extinction / isochrone step did not run'
    ir = bundle.get('imf_res')
    if ir is not None:
        figs['imf'] = IMF.plot_imf(ir, name, bundle.get('mf_kind', IMF.MF_KIND))
    else:
        skipped['imf'] = 'IMF step did not run'
    if bundle.get('msr'):
        figs['msr'] = msr_figure(bundle['msr'], f'{name} mass segregation')
    else:
        skipped['msr'] = 'dynamics step did not run'
    return figs, skipped


def export_figure(fig: Figure, base_path: str, formats=FIGURE_FORMATS, dpi: int = FIGURE_DPI) -> list:
    """Save one figure as base_path.<fmt> for each format; returns the written paths."""
    os.makedirs(os.path.dirname(os.path.abspath(base_path)), exist_ok=True)
    out = []
    for fmt in formats:
        p = f'{base_path}.{fmt}'
        fig.savefig(p, format=fmt, dpi=dpi, bbox_inches='tight')
        out.append(p)
    return out


def export_all(figs: dict, out_dir: str, prefix: str, formats=FIGURE_FORMATS) -> list:
    return [p for k, f in figs.items() for p in export_figure(f, os.path.join(out_dir, f'{prefix}_{k}'), formats)]


def _raster(fig: Figure, dpi: int) -> np.ndarray:
    buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight'); buf.seek(0)
    return mpimg.imread(buf)


def one_page_report(path: str, name: str, summary: dict, figs: dict, skipped: dict, warnings: list, provenance: str = '') -> str:
    """A4 one-page PDF for one cluster. `summary` uses the summary.csv column names; missing / NaN values print as '-'."""
    page = Figure(figsize=PAGE_SIZE_IN)
    gs = page.add_gridspec(5, 2, height_ratios=[0.55, 0.55, 1.4, 1.4, 1.4], hspace=0.22, wspace=0.08, left=0.04, right=0.96, top=0.97, bottom=0.12)
    ax = page.add_subplot(gs[0, :]); ax.axis('off')
    ax.text(0, 0.8, f'StarNest report: {name}', fontsize=16, fontweight='bold', va='top')
    ax.text(0, 0.18, provenance, fontsize=7.5, color='#444', va='top', wrap=True)
    ax = page.add_subplot(gs[1, :]); ax.axis('off')
    cells = []
    for label, key, fmt in SUMMARY_FIELDS:
        v = summary.get(key)
        try:
            cells.append(f'{label}: ' + (fmt.format(v) if v is not None and np.isfinite(float(v)) else '-'))
        except (TypeError, ValueError):
            cells.append(f'{label}: -')
    for i, c in enumerate(cells):
        ax.text((i % 3) / 3.0, 0.85 - 0.33 * (i // 3), c, fontsize=9, va='top', family='monospace')
    for i, key in enumerate(PANEL_ORDER):
        a = page.add_subplot(gs[2 + i // 2, i % 2]); a.axis('off')
        if key in figs:
            a.imshow(_raster(figs[key], REPORT_DPI))
        else:
            a.text(0.5, 0.5, f'{PANEL_TITLES[key]}\nnot available:\n{skipped.get(key, "no data")}', ha='center', va='center', fontsize=9, color='#777', transform=a.transAxes)
    page.text(0.04, 0.105, 'Warnings' if warnings else 'No warnings raised', fontsize=9, fontweight='bold', va='top', color=WARNING_COLOUR if warnings else '#2e7d32')
    y = 0.088
    for w in warnings[:9]:
        page.text(0.04, y, '- ' + w[:150], fontsize=7, va='top', color=WARNING_COLOUR); y -= 0.0095
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    page.savefig(path, format='pdf')
    return path

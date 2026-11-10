/* footer.js -- two-mode footer, citation buttons and citation modal.
   Self-contained: it only READS the page state (it never edits app.js). */

// ======================= GLOBAL CONFIGURATION =======================
const GITHUB_REPO_URL   = 'https://github.com/NIKHILKAUSHIK9/demo_trial.github.io';
const GITHUB_ISSUES_URL = GITHUB_REPO_URL + '/issues';
const PAPER_DOI_URL     = '#';              // set to the preprint / DOI link when available ('#' = button shows "coming soon")
const COPIED_MS         = 2000;             // how long the "Copied!" state stays
const WORKSPACE_PAGES   = { membership: 'Membership', isochrone: 'Isochrone', kinematics: 'Kinematics', imf: 'IMF' };

// TODO(NK): the supplied text ended at "author = {Panwar, Neel" - author order, venue and DOI below are placeholders; edit before release.
const BIBTEX_REPO = `@software{starnest_github_2026,
  author       = {Kaushik, Nikhil and Basur, Pryanshu and {Monu} and Prasad, Shruti and Panwar, Neelam and Kundu, Richa},
  title        = {{StarNest: Open-Source Web Pipeline for Open Star Cluster Photometry and Kinematics}},
  year         = {2026},
  version      = {1.0.0},
  publisher    = {GitHub},
  url          = {${GITHUB_REPO_URL}}
}`;
const BIBTEX_PAPER = `@article{starnest_clusters_2026,
  author       = {Panwar, Neelam and Kundu, Richa and Kaushik, Nikhil and Basur, Pryanshu and {Monu} and Prasad, Shruti},
  title        = {Photometric and Kinematic Study of Open Star Clusters},
  journal      = {In preparation},
  year         = {2026},
  note         = {Haffner 14, Melotte 111, NGC 2910, Ruprecht 147, UBC 52, UPK 545; Gaia DR3}
}`;
const BIBITEM_PAPER = `\\bibitem{starnest_clusters_2026} Panwar N., Kundu R., Kaushik N., Basur P., Monu, Prasad S., 2026, Photometric and Kinematic Study of Open Star Clusters, in preparation`;
// =====================================================================

(function () {
  const $  = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
  const CITES = { repo: BIBTEX_REPO, paper: BIBTEX_PAPER, bibitem: BIBITEM_PAPER };

  // ---------- 1. mode switching: Home -> full footer, workspace tabs -> slim status bar ----------
  function currentPage() {
    // is-entering = page sliding in (switch early); is-active = settled page
    const el = $('.page-section.is-entering') || $('.page-section.is-active');
    return el ? el.id.replace(/^page-/, '') : 'home';
  }
  function applyMode(page) {
    page = page || currentPage();
    const ws = page in WORKSPACE_PAGES;
    const home = $('#footer-home'), bar = $('#footer-workspace');
    if (home) home.style.display = ws ? 'none' : '';
    if (bar) { bar.classList.toggle('hidden', !ws); bar.classList.toggle('flex', ws); }
    document.body.classList.toggle('ws-footer', ws);              // adds bottom padding so the bar never covers content
    const st = $('#footer-stage'); if (st && ws) st.textContent = WORKSPACE_PAGES[page];
  }
  const sections = $$('.page-section');
  const mo = new MutationObserver(() => applyMode());
  sections.forEach(s => mo.observe(s, { attributes: true, attributeFilter: ['class'] }));
  $$('.nav-btn').forEach(b => b.addEventListener('click', () => applyMode(b.dataset.page)));
  document.addEventListener('click', (e) => {                      // any [data-navigate] link on the page
    const t = e.target.closest('[data-navigate]'); if (t) applyMode(t.dataset.navigate);
  });
  applyMode();

  // ---------- 2. links from constants ----------
  $$('[data-repo-link],[data-repo-cite]').forEach(a => { a.href = GITHUB_REPO_URL; });
  $$('[data-issues-link]').forEach(a => { a.href = GITHUB_ISSUES_URL; });
  $$('[data-paper-link]').forEach(a => {
    a.href = PAPER_DOI_URL;
    if (PAPER_DOI_URL === '#') { a.removeAttribute('target'); a.title = 'Paper link coming soon'; a.addEventListener('click', e => e.preventDefault()); }
  });

  // ---------- 3. footer jump links -> real nav buttons (reuses app.js navigation) ----------
  $$('[data-footer-jump]').forEach(a => a.addEventListener('click', (e) => {
    e.preventDefault();
    const btn = $(`.nav-btn[data-page="${a.dataset.footerJump}"]`); if (btn) btn.click();
  }));

  // ---------- 4. copy buttons with 2 s "Copied!" state ----------
  async function copyText(txt) {
    try { await navigator.clipboard.writeText(txt); return true; } catch (e) { /* fall through (iframes / http) */ }
    const ta = document.createElement('textarea'); ta.value = txt; ta.style.cssText = 'position:fixed;opacity:0';
    document.body.appendChild(ta); ta.select(); let ok = false;
    try { ok = document.execCommand('copy'); } catch (e) {} ta.remove(); return ok;
  }
  $$('[data-copy]').forEach(btn => btn.addEventListener('click', async () => {
    const lab = $('[data-label]', btn) || btn, orig = btn.dataset.orig || lab.textContent; btn.dataset.orig = orig;
    const ok = await copyText(CITES[btn.dataset.copy]);
    lab.textContent = ok ? 'Copied!' : 'Press Ctrl+C';
    btn.classList.toggle('is-copied', ok);
    clearTimeout(btn._t); btn._t = setTimeout(() => { lab.textContent = orig; btn.classList.remove('is-copied'); }, COPIED_MS);
  }));

  // ---------- 5. citation modal ----------
  const modal = $('#cite-modal');
  function openModal() {
    $('#cite-pre-repo').textContent = BIBTEX_REPO; $('#cite-pre-paper').textContent = BIBTEX_PAPER; $('#cite-pre-bibitem').textContent = BIBITEM_PAPER;
    modal.classList.remove('hidden'); modal.classList.add('flex');
  }
  function closeModal() { modal.classList.add('hidden'); modal.classList.remove('flex'); }
  $$('[data-open-cite]').forEach(b => b.addEventListener('click', openModal));
  $$('[data-close-cite]').forEach(b => b.addEventListener('click', closeModal));
  modal.addEventListener('click', (e) => { if (e.target === modal) closeModal(); });
  document.addEventListener('keydown', (e) => { if (e.key === 'Escape') closeModal(); });
})();

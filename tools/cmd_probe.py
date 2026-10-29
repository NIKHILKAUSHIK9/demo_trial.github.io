# ===== GLOBAL CONFIGURATION =====
PAGES = ["https://stev.oapd.inaf.it/cgi-bin/cmd_3.8", "https://stev.oapd.inaf.it/cgi-bin/cmd"]
TIMEOUT_S = 30
OUT_DIR = "cmd_probe_out"        # raw HTML of each page is saved here (send me these files if the output is unclear)
# ================================
# Workflow: GET each CMD form page once -> print title, form action/method and every field with its default
# -> save the raw HTML. Purpose: see exactly which form fields CMD 3.8 / 3.9 expect today.
import os, re, sys, requests
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))   # import the project modules from the repo root
os.makedirs(OUT_DIR, exist_ok=True)
for i, u in enumerate(PAGES):
    print("=" * 70, "\n", u)
    try:
        r = requests.get(u, timeout=(30, TIMEOUT_S)); html = r.text
    except Exception as e:
        print("  FAILED:", type(e).__name__, e); continue
    path = os.path.join(OUT_DIR, f"page{i}.html"); open(path, "w", encoding="utf-8").write(html)
    print("  HTTP", r.status_code, "| saved", path, "|", len(html), "chars")
    print("  title:", (re.findall(r"<title>(.*?)</title>", html, re.S | re.I) or ["?"])[0].strip())
    for f in re.findall(r"<form[^>]*>", html, re.I): print("  form:", re.sub(r"\s+", " ", f))
    for m in re.finditer(r"<(input|select|textarea)\b([^>]*)>", html, re.I | re.S):
        a = m.group(2); nm = re.search(r'name\s*=\s*["\']?([^"\'\s>]+)', a, re.I)
        if not nm: continue
        val = re.search(r'value\s*=\s*["\']([^"\']*)', a, re.I)
        print(f"  {m.group(1).lower():8s} {nm.group(1):22s} = {val.group(1) if val else ''}")

# ---- test POST: exactly what the app sends (logAge 7.5, [M/H] 0, A_V 0), as form body and as URL query ----
import ocmass_extinction_isochrone as E
U = "https://stev.oapd.inaf.it/cgi-bin/cmd_3.8"
params = E.build_cmd_params(7.5, 0.0, 0.0)
for variant in ("body", "query"):
    print("=" * 70, "\n POST", variant, U)
    try:
        r = requests.post(U, timeout=(30, 120), **(dict(data=params) if variant == "body" else dict(params=params)))
    except Exception as e:
        print("  FAILED:", type(e).__name__, e); continue
    path = os.path.join(OUT_DIR, f"post_{variant}.html"); open(path, "w", encoding="utf-8").write(r.text)
    print("  HTTP", r.status_code, "| saved", path, "| output link:", re.findall(r"output\d+", r.text)[:3] or "NONE")
    print("  summary:", E._reply_summary(r.text, 600))
    print("  links:", re.findall(r'href\s*=\s*["\']([^"\']+)', r.text)[:15])
    if re.findall(r"output\d+", r.text): break

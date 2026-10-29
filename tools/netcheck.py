# ===== GLOBAL CONFIGURATION =====
HOSTS = ["irsa.ipac.caltech.edu", "stev.oapd.inaf.it"]
URLS = ["https://irsa.ipac.caltech.edu/applications/DUST/", "https://stev.oapd.inaf.it/cgi-bin/cmd_3.8"]
PORT = 443
TIMEOUT_S = 15
# ================================
# Workflow: (1) raw TCP connect per address family (ignores proxies) -> (2) HTTPS GET with system proxy settings
# -> (3) HTTPS GET ignoring all proxies. Compare the three to see whether the network, IPv6 or a proxy is the problem.
import socket, urllib.request, requests
print("proxies seen by Python:", urllib.request.getproxies() or "none")
for h in HOSTS:
    print("== TCP", h)
    try: infos = socket.getaddrinfo(h, PORT, proto=socket.IPPROTO_TCP)
    except Exception as e: print("  DNS FAILED:", e); continue
    for fam, _, _, _, addr in sorted(set(infos), key=lambda x: x[0].name):
        s = socket.socket(fam, socket.SOCK_STREAM); s.settimeout(TIMEOUT_S)
        try: s.connect(addr); print(f"  {fam.name:8s} {addr[0]} connect OK")
        except Exception as e: print(f"  {fam.name:8s} {addr[0]} FAILED: {e}")
        finally: s.close()
for u in URLS:
    for label, trust in (("with system proxy", True), ("direct, no proxy", False)):
        s = requests.Session(); s.trust_env = trust
        try: print(f"{label:18s} {u} -> HTTP {s.get(u, timeout=TIMEOUT_S).status_code}")
        except Exception as e: print(f"{label:18s} {u} -> FAILED {type(e).__name__}")

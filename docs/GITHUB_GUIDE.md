# Publishing the project on GitHub: step by step

## 0. Before you publish (checklist)

| Item | Why | Done |
|---|---|---|
| Team members agreed to their photos and names being public | Photos of people in a public repo are hard to take back | ☐ |
| Supervisors / department agree to the license (MIT by default) | Research code can fall under institute rules | ☐ |
| Replace `<your-user>/<your-repo>` in `README.md` and `CITATION.cff` | Links and citation | ☐ |
| No private data in the repo (your Gaia CSVs are in `.gitignore` under `data/`) | Unpublished results | ☐ |
| Add `images/team/<name>.jpg` photos (lower-case names) | Shown on the Home page | ☐ |

## 1. Install tools (Windows)

1. Install **Git for Windows**: https://git-scm.com/download/win (defaults are fine).
2. Create a free account at https://github.com.
3. In a terminal (PowerShell or Git Bash) set your identity once:
   ```
   git config --global user.name  "Your Name"
   git config --global user.email "you@example.com"
   ```

## 2. Create the repository on GitHub

1. GitHub → **New repository**.
2. Name: e.g. `starnest-ocmass`. Visibility: **Public** (or Private first, then switch later).
3. **Do not** tick "Add a README / .gitignore / license" (the project already has them).
4. Click **Create repository** and keep the page open (it shows the repository URL).

## 3. Push the project

Run inside the unzipped project folder:

```
git init
git branch -M main
git add .
git status                    # check: no data files, no __pycache__, photos present
git commit -m "Initial release: StarNest / OCMass pipeline"
git remote add origin https://github.com/<your-user>/<your-repo>.git
git push -u origin main
```
When asked to sign in, use the browser prompt; if a password is requested, create a **Personal Access Token**
(GitHub → Settings → Developer settings → Tokens) and paste it instead of the password.

## 4. Publish the web page (GitHub Pages)

1. Repository → **Settings → Pages**.
2. Source: **Deploy from a branch**, Branch: `main`, folder `/ (root)`, **Save**.
3. After about a minute the page is at `https://<your-user>.github.io/<your-repo>/`.

| Runs on Pages | Does not run on Pages |
|---|---|
| GMM, kinematics, IMF, plots, team page | IRSA extinction lookup and the PARSEC file (need the Python backend, which Pages cannot host) |

State this in the README (already done) so visitors know to run `uvicorn main:app` locally for the full chain.

## 5. Release a version

```
git tag -a v1.0.0 -m "First public release"
git push origin v1.0.0
```
GitHub → **Releases → Draft a new release** → choose tag `v1.0.0` → add notes → **Publish**. For a citable DOI, link the repository to **Zenodo** (zenodo.org → GitHub) and publish the release again.

## 6. Day-to-day workflow

| Task | Commands |
|---|---|
| Save a change | `git add -A` → `git commit -m "what changed"` → `git push` |
| See what changed | `git status`, `git diff` |
| Try something risky | `git switch -c experiment` … `git switch main` |
| Update from GitHub | `git pull` |
| Undo an uncommitted edit of one file | `git restore path/to/file` |

## 7. Good practice

- One logical change per commit; write messages in the present tense ("Add IMF completeness limit").
- Keep outputs (`ocmass_output/`, CSV results) out of git; share them as release attachments if needed.
- Add team members as **collaborators** (Settings → Collaborators) or use pull requests for review.
- Photos: keep each under about 300 KB (resize before committing) so the repository stays small.
- Never commit tokens or passwords. If you do, revoke the token first, then clean the history.
- Pin library versions (`pip freeze > requirements.txt`) at release time so the results can be reproduced.

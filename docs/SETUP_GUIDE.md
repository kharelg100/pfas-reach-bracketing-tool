# Hosting Setup Guide

Step-by-step instructions to deploy the PFAS Reach Bracketing Tool on **GitHub Pages**, **Zenodo** (DOI), and **Netlify**.

---

## 1. GitHub Pages (Primary Hosting)

### Initial Setup (one time)

**Step 1: Create the repository**
1. Go to [github.com/new](https://github.com/new)
2. Repository name: `pfas-reach-bracketing-tool`
3. Description: `Interactive browser-based PFAS reach bracketing analysis tool with multiple imputation`
4. Set to **Public**
5. Check **"Add a README file"** (we'll replace it)
6. License: **MIT**
7. Click **Create repository**

**Step 2: Upload all files**
1. On the repository page, click **"Add file" → "Upload files"**
2. Drag the entire contents of this repository folder into the upload area:
   - `index.html` (this is the browser tool)
   - `README.md`
   - `LICENSE`
   - `CITATION.cff`
   - `.zenodo.json`
   - `netlify.toml`
   - `python/` folder (with the Python script)
   - `templates/` folder (with CSV templates)
   - `docs/` folder (with this guide)
3. Commit message: `Initial release v1.0.0`
4. Click **Commit changes**

**Step 3: Enable GitHub Pages**
1. Go to **Settings** → **Pages** (left sidebar)
2. Source: **Deploy from a branch**
3. Branch: **main** / **/ (root)**
4. Click **Save**
5. Wait 1–2 minutes. Your site will be live at:
   ```
   https://YOUR_USERNAME.github.io/pfas-reach-bracketing-tool/
   ```
   For example: `https://gkharel.github.io/pfas-reach-bracketing-tool/`

**Step 4: Verify**
- Visit the URL — you should see the tool's purple header and upload boxes
- Click **"▶ Run with demo data"** to confirm it works
- Test on mobile to verify responsive layout

### Updating the Tool

To push updates after the initial setup:
1. Edit files on GitHub directly (click any file → pencil icon → edit → commit), or
2. Use Git from your computer:
   ```bash
   git clone https://github.com/YOUR_USERNAME/pfas-reach-bracketing-tool.git
   cd pfas-reach-bracketing-tool
   # Make edits...
   git add .
   git commit -m "Description of changes"
   git push
   ```
   GitHub Pages redeploys automatically within ~1 minute.

---

## 2. Zenodo (DOI for Citation)

Zenodo assigns a permanent DOI that you can cite in the ES&T Letters manuscript.

**Step 1: Connect GitHub to Zenodo**
1. Go to [zenodo.org](https://zenodo.org) and log in with your GitHub account
2. Click your name (top right) → **GitHub**
3. Find `pfas-reach-bracketing-tool` in the list
4. Flip the toggle to **ON**

**Step 2: Create a GitHub Release**
1. On your GitHub repository page, click **Releases** (right sidebar)
2. Click **"Create a new release"**
3. Tag version: `v1.0.0`
4. Release title: `PFAS Reach Bracketing Tool v1.0.0`
5. Description:
   ```
   Initial release accompanying:
   Kharel, G. Machine-readable reach bracketing reveals persistent 
   tributary-driven PFAS loading in surface water networks. 
   ES&T Letters (submitted 2026).
   
   Features:
   - Browser-based MI analysis (250–10,000 replicates)
   - Machine-readable reach experiment contrasts
   - RL/2 sensitivity comparison
   - SVG figure exports and Excel workbook
   - Python companion script
   ```
6. Click **Publish release**

**Step 3: Get your DOI**
1. Wait ~5 minutes, then go to [zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/)
2. Your release should appear with a DOI badge
3. Click the DOI — it will look like: `10.5281/zenodo.1234567`

**Step 4: Update your files**
1. Replace `XXXXXXX` in `README.md` with your actual Zenodo DOI number
2. Replace `XXXXXXX` in `CITATION.cff` if you wish
3. Add your ORCID to `CITATION.cff` and `.zenodo.json`
4. Commit and push the changes

**Step 5: Add to manuscript**
In your ES&T Letters manuscript, cite as:
```
The interactive tool is available at 
https://gkharel.github.io/pfas-reach-bracketing-tool/ 
(archived at DOI: 10.5281/zenodo.XXXXXXX).
```

### Important Notes
- Each new GitHub Release creates a **new DOI version** on Zenodo
- Zenodo also provides a **concept DOI** that always resolves to the latest version
- Use the concept DOI in publications for longevity

---

## 3. Netlify (Alternative/Backup Hosting)

Netlify provides a second hosting option with instant deployment and optional custom domains.

### Option A: Drag-and-Drop (fastest, no account needed at first)

1. Go to [app.netlify.com/drop](https://app.netlify.com/drop)
2. Drag your entire repository folder onto the page
3. Within seconds you'll get a URL like: `https://random-name-abc123.netlify.app`
4. Create a free account to claim the site and customize the URL

### Option B: Connect to GitHub (auto-deploys on push)

**Step 1: Connect**
1. Go to [app.netlify.com](https://app.netlify.com) → Sign up with GitHub
2. Click **"Add new site" → "Import an existing project"**
3. Select **GitHub** → Authorize → Choose `pfas-reach-bracketing-tool`
4. Build settings (auto-detected from `netlify.toml`):
   - Build command: *(leave blank)*
   - Publish directory: `.`
5. Click **Deploy site**

**Step 2: Customize URL**
1. Go to **Site configuration → Site details → Change site name**
2. Set to: `pfas-reach-bracketing` 
3. Your URL becomes: `https://pfas-reach-bracketing.netlify.app`

**Step 3: Custom Domain (optional)**
If TCU provides a subdomain:
1. Go to **Domain management → Add custom domain**
2. Enter: `pfas-tool.tcu.edu` (or whatever TCU provides)
3. Follow DNS instructions (add a CNAME record pointing to Netlify)
4. Netlify auto-provisions HTTPS

### Netlify vs GitHub Pages

| Feature | GitHub Pages | Netlify |
|---------|-------------|---------|
| Cost | Free | Free (basic) |
| Custom domain | Yes (manual DNS) | Yes (guided setup) |
| HTTPS | Automatic | Automatic |
| Deploy speed | ~1 min | ~10 sec |
| Auth/password protect | No | Yes (paid plan) |
| Analytics | No | Yes (paid plan) |
| Best for | Academic + DOI workflow | Custom domain + polish |

---

## Recommended Workflow for ES&T Letters

1. **Now:** Set up GitHub repository + GitHub Pages (Steps 1.1–1.4)
2. **Before submission:** Create GitHub Release + Zenodo DOI (Steps 2.1–2.5)
3. **In manuscript:** Include both the live URL and DOI
4. **Optional:** Set up Netlify as backup/custom-domain host

### What to Include in the Manuscript

**Supporting Information section:**
```
An interactive browser-based tool implementing the reach bracketing 
framework is available at https://gkharel.github.io/pfas-reach-bracketing-tool/ 
(archived: DOI 10.5281/zenodo.XXXXXXX). The tool and companion Python 
script are provided as Supporting Information files S1 (HTML) and S2 (Python).
```

**Data Availability Statement:**
```
The PFAS Reach Bracketing Tool (browser-based and Python versions), 
CSV templates, and demonstration data are openly available at 
https://github.com/gkharel/pfas-reach-bracketing-tool under the MIT license 
(DOI: 10.5281/zenodo.XXXXXXX).
```

---

© 2025–2026 Dr. Gehendra Kharel, Texas Christian University.

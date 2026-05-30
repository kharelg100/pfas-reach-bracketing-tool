# Deploy PFAS Reach Bracketing Tool v1.1.1 via PowerShell

This package uses `index.html` as a full GitHub Pages entry page and also includes `PFAS_Reach_Bracketing_Tool.html` as an intuitive local filename. Both files are fully styled tool pages.

## 1. Extract the package

```powershell
$ZipPath = "$env:USERPROFILE\Downloads\PFAS_Reach_Bracketing_Browser_Tool_v1_1_1_Final.zip"
$ReleaseDir = "$env:USERPROFILE\Desktop\PFAS_Reach_Bracketing_Tool_v1_1_1"
Remove-Item -Recurse -Force $ReleaseDir -ErrorAction SilentlyContinue
Expand-Archive -Path $ZipPath -DestinationPath $ReleaseDir -Force
```

## 2. Test locally

```powershell
Start-Process "$ReleaseDir\pfas_reach_bracketing_browser_v1_1_1_final\PFAS_Reach_Bracketing_Tool.html"
```

Click **Run with manuscript demo data** and confirm that TW91 is approximately 15,667 ng L-1 and WF1 `f_trib` is approximately 0.047.

## 3. Clone or open the GitHub Pages repository

```powershell
$RepoUrl = "https://github.com/kharelg100/pfas-reach-bracketing-tool.git"
$RepoDir = "$env:USERPROFILE\Desktop\pfas-reach-bracketing-tool"
if (-not (Test-Path $RepoDir)) { git clone $RepoUrl $RepoDir }
```

## 4. Copy files to the repository root

```powershell
$Source = "$ReleaseDir\pfas_reach_bracketing_browser_v1_1_1_final"
robocopy $Source $RepoDir /E /XD .git /XF .DS_Store
Set-Location $RepoDir
git status
git add .
git commit -m "Deploy PFAS Reach Bracketing Tool v1.1.1 final"
git push origin main
```

## 5. GitHub Pages settings

In GitHub, set **Settings → Pages → Deploy from a branch → main → /(root)**, unless your repository uses `/docs` as the Pages source.

## 6. Zenodo

If archiving this exact deployment package, create a **new version** of the existing Zenodo record so the concept DOI remains `10.5281/zenodo.18937739`. Update the version DOI in the tool only after Zenodo assigns it.


Validation note: When running the manuscript demo data, the summary metric cards should display TW91 ΣPFAS40, WF1 ΔΣPFAS40, and WF1 ftrib. WF1 ftrib appears in the summary metrics and in the ReachExperiments.csv output.

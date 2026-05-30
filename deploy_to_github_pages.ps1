# PFAS Reach Bracketing Tool v1.1.1 deployment helper
$ErrorActionPreference = "Stop"
$RepoUrl = "https://github.com/kharelg100/pfas-reach-bracketing-tool.git"
$RepoDir = "$env:USERPROFILE\Desktop\pfas-reach-bracketing-tool"
$Source = Split-Path -Parent $MyInvocation.MyCommand.Path
if (-not (Test-Path $RepoDir)) { git clone $RepoUrl $RepoDir }
robocopy $Source $RepoDir /E /XD .git /XF .DS_Store
Set-Location $RepoDir
git status
Write-Host "Review changes. Then run: git add .; git commit -m 'Deploy PFAS Reach Bracketing Tool v1.1.1 final'; git push origin main"

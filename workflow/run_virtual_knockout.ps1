[CmdletBinding()]
param(
    [string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }),
    [string]$Rscript = $(if ($env:ONFH_RSCRIPT) { $env:ONFH_RSCRIPT } else { "Rscript" }),
    [ValidateRange(1, 128)] [int]$Cores = 1
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$env:ONFH_ROOT = $RepoRoot
$env:ONFH_PROJECT_ROOT = $RepoRoot
$env:PYTHONPATH = Join-Path $RepoRoot "analysis"

function Invoke-Checked {
    param([string]$Executable, [string[]]$Arguments)
    Write-Host "`n> $Executable $($Arguments -join ' ')" -ForegroundColor Cyan
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Stage failed with exit code $LASTEXITCODE" }
}

Push-Location $RepoRoot
try {
    Invoke-Checked $Rscript @("virtual_knockout/run_official_vko.R", "--profile=manuscript", "--cores=$Cores")
    Invoke-Checked $Rscript @("virtual_knockout/run_official_vko.R", "--profile=official-default", "--cores=$Cores")
    Invoke-Checked $Python @("virtual_knockout/postprocess_official_vko.py")
    Invoke-Checked $Rscript @("virtual_knockout/export_figure_data.R")
    Invoke-Checked $Python @("plotting/make_virtual_knockout_figure.py")
} finally {
    Pop-Location
}

Write-Host "Official-R virtual-knockout workflow completed." -ForegroundColor Green

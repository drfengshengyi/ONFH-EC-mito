[CmdletBinding()]
param([string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }))

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$env:ONFH_ROOT = $RepoRoot
$env:PYTHONPATH = Join-Path $RepoRoot "analysis"

function Invoke-Python {
    param([string]$Script)
    Write-Host "`n> $Python $Script" -ForegroundColor Cyan
    & $Python $Script
    if ($LASTEXITCODE -ne 0) { throw "$Script failed with exit code $LASTEXITCODE" }
}

Push-Location $RepoRoot
try {
    Invoke-Python "plotting/make_evidence_model.py"
    Invoke-Python "plotting/make_virtual_knockout_figure.py"
    Invoke-Python "plotting/assemble_manuscript_figures.py"
    Invoke-Python "plotting/make_reviewed_figures.py"
    Invoke-Python "plotting/correct_figure4_fdr.py"
} finally {
    Pop-Location
}

Write-Host "Final figures completed." -ForegroundColor Green

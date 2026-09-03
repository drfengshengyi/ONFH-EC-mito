[CmdletBinding()]
param(
    [string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }),
    [ValidateRange(0, 1000000)] [int]$MinimumUmi = 100
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$env:ONFH_ROOT = $RepoRoot
$env:PYTHONPATH = Join-Path $RepoRoot "analysis"

Push-Location $RepoRoot
try {
    Write-Host "`n> $Python analysis/spatial_contextualization.py --minimum-umi $MinimumUmi" -ForegroundColor Cyan
    & $Python "analysis/spatial_contextualization.py" --minimum-umi $MinimumUmi
    if ($LASTEXITCODE -ne 0) { throw "Spatial contextualization failed with exit code $LASTEXITCODE" }
} finally {
    Pop-Location
}

Write-Host "External femoral-head spatial contextualization completed." -ForegroundColor Green

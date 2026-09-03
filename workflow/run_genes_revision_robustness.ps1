[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)] [string]$SourceAnalysis,
    [Parameter(Mandatory = $true)] [string]$SourceData,
    [string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }),
    [string]$Rscript = $(if ($env:ONFH_RSCRIPT) { $env:ONFH_RSCRIPT } else { "Rscript" }),
    [ValidateRange(1, 128)] [int]$Jobs = 4
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$FigureInput = Join-Path $RepoRoot "results\figure_inputs"
$FgseaOutput = Join-Path $RepoRoot "results\participant_fgsea_stability"
$env:PYTHONPATH = Join-Path $RepoRoot "analysis"

function Invoke-Checked {
    param([string]$Exe, [string[]]$Arguments)
    Write-Host "> $Exe $($Arguments -join ' ')" -ForegroundColor Cyan
    & $Exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Exe failed with exit code $LASTEXITCODE"
    }
}

Push-Location $RepoRoot
try {
    $env:ONFH_DATA_ROOT = (Resolve-Path $SourceData).Path
    $env:ONFH_OUTPUT_DIR = $FigureInput
    Invoke-Checked $Python @(
        "analysis/serum_classifier.py", "--permutations", "1000",
        "--jobs", "$Jobs", "--no-plot"
    )
    Invoke-Checked $Rscript @("analysis/serum_paired_comparison.R")

    Invoke-Checked $Python @(
        "analysis/participant_fgsea_stability.py",
        "--source-analysis", (Resolve-Path $SourceAnalysis).Path,
        "--output-dir", $FgseaOutput
    )
    $env:ONFH_FGSEA_STABILITY_DIR = $FgseaOutput
    Invoke-Checked $Rscript @("analysis/participant_fgsea_stability.R")
} finally {
    Pop-Location
}

Write-Host "Genes revision robustness analyses completed." -ForegroundColor Green

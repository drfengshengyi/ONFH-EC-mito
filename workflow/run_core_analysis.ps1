[CmdletBinding()]
param(
    [string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }),
    [string]$Rscript = $(if ($env:ONFH_RSCRIPT) { $env:ONFH_RSCRIPT } else { "Rscript" }),
    [ValidateRange(0, 100000)] [int]$Permutations = 1000,
    [ValidateRange(1, 128)] [int]$Jobs = 4
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$env:ONFH_ROOT = $RepoRoot
$env:ONFH_PROJECT_ROOT = $RepoRoot
$analysisPath = Join-Path $RepoRoot "analysis"
$env:PYTHONPATH = if ($env:PYTHONPATH) { "$analysisPath;$env:PYTHONPATH" } else { $analysisPath }

function Invoke-Checked {
    param([string]$Executable, [string[]]$Arguments)
    Write-Host "`n> $Executable $($Arguments -join ' ')" -ForegroundColor Cyan
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Stage failed with exit code $LASTEXITCODE" }
}

Push-Location $RepoRoot
try {
    Invoke-Checked $Rscript @("analysis/prepare_matrices.R")
    foreach ($script in @(
        "atlas_normalize.py",
        "atlas_batch_correct.py",
        "atlas_cluster.py",
        "atlas_annotate.py",
        "endothelial_subset.py",
        "build_genesets.py",
        "endothelial_panels.py",
        "annotation_audit.py",
        "composition_inference.py",
        "pseudobulk.py"
    )) { Invoke-Checked $Python @("analysis/$script") }
    Invoke-Checked $Rscript @("analysis/pathway_enrichment.R")
    foreach ($script in @(
        "mitochondrial_effects.py",
        "cell_communication.py",
        "regulon_activity.py"
    )) { Invoke-Checked $Python @("analysis/$script") }
    Invoke-Checked $Python @("analysis/serum_classifier.py", "--permutations", "$Permutations", "--jobs", "$Jobs")
    foreach ($script in @(
        "serum_comparator.py",
        "supplementary_audits.py",
        "write_provenance.py"
    )) { Invoke-Checked $Python @("analysis/$script") }
} finally {
    Pop-Location
}

Write-Host "Core analysis completed." -ForegroundColor Green

param(
    [Parameter(Mandatory = $true)]
    [string]$SourceAnalysis,
    [string]$Python = "python",
    [string]$ProjectRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$OutputDir = ""
)

$ErrorActionPreference = "Stop"
$script = Join-Path $ProjectRoot "analysis\ec_preprocessing_sensitivity.py"
if (-not (Test-Path -LiteralPath $script)) {
    throw "Missing analysis script: $script"
}
if (-not (Test-Path -LiteralPath (Join-Path $SourceAnalysis "ec_final.h5ad"))) {
    throw "Missing ec_final.h5ad in SourceAnalysis: $SourceAnalysis"
}
if (-not (Test-Path -LiteralPath (Join-Path $SourceAnalysis "atlas_annotated.h5ad"))) {
    throw "Missing atlas_annotated.h5ad in SourceAnalysis: $SourceAnalysis"
}

$arguments = @(
    $script,
    "--source-analysis", $SourceAnalysis,
    "--project-root", $ProjectRoot
)
if ($OutputDir) {
    $arguments += @("--output-dir", $OutputDir)
}

& $Python @arguments
if ($LASTEXITCODE -ne 0) {
    throw "EC preprocessing sensitivity analysis failed with exit code $LASTEXITCODE"
}

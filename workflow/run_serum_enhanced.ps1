param(
    [ValidateSet("preflight", "smoke", "formal", "plot", "finalize", "all")]
    [string]$Phase = "all",
    [switch]$Resume,
    [string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }),
    [string]$Rscript = $(if ($env:ONFH_RSCRIPT) { $env:ONFH_RSCRIPT } else { "Rscript" }),
    [string]$ProjectRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$DataDir = "",
    [string]$OutputDir = "",
    [ValidateRange(1, 128)] [int]$Jobs = 4
)

$ErrorActionPreference = "Stop"
if (-not $DataDir) { $DataDir = Join-Path $ProjectRoot "data" }
if (-not $OutputDir) { $OutputDir = Join-Path $ProjectRoot "results\serum_enhanced_20260910" }
$Analysis = Join-Path $ProjectRoot "analysis\serum_classifier_enhanced.py"
$QA = Join-Path $ProjectRoot "qa\check_serum_enhanced.py"
$Plot = Join-Path $ProjectRoot "plotting\make_figure7_enhanced.R"

if (-not (Test-Path -LiteralPath $Python)) { throw "Missing isolated Python environment: $Python" }
if (-not (Test-Path -LiteralPath $Rscript)) { throw "Missing required Rscript: $Rscript" }

Set-Location -LiteralPath $ProjectRoot
$ResumeArg = @()
if ($Resume) { $ResumeArg = @("--resume") }
$AnalysisPaths = @("--project-root", $ProjectRoot, "--data-dir", $DataDir, "--output-dir", $OutputDir)
$QAPaths = @("--project-root", $ProjectRoot, "--output-dir", $OutputDir)

if ($Phase -in @("preflight", "all")) {
    & $Python $Analysis --phase preflight --jobs $Jobs @AnalysisPaths
    if ($LASTEXITCODE -ne 0) { throw "Preflight analysis failed" }
    & $Python $QA --stage preflight @QAPaths
    if ($LASTEXITCODE -ne 0) { throw "Preflight QC failed" }
}
if ($Phase -in @("smoke", "all")) {
    & $Python $Analysis --phase smoke --jobs $Jobs @ResumeArg @AnalysisPaths
    if ($LASTEXITCODE -ne 0) { throw "Smoke analysis failed" }
    & $Python $QA --stage smoke @QAPaths
    if ($LASTEXITCODE -ne 0) { throw "Smoke QC failed; formal analysis is prohibited" }
}
if ($Phase -in @("formal", "all")) {
    & $Python $Analysis --phase formal --jobs $Jobs @ResumeArg @AnalysisPaths
    if ($LASTEXITCODE -ne 0) { throw "Formal analysis failed" }
    & $Python $QA --stage formal @QAPaths
    if ($LASTEXITCODE -ne 0) { throw "Formal analysis QC failed" }
}
if ($Phase -in @("plot", "all")) {
    & $Rscript $Plot --project-root $ProjectRoot --output-dir $OutputDir
    if ($LASTEXITCODE -ne 0) { throw "Enhanced Figure 7 plotting failed" }
}
if ($Phase -in @("finalize", "all")) {
    & $Python $Analysis --phase finalize --jobs $Jobs @AnalysisPaths
    if ($LASTEXITCODE -ne 0) { throw "Report finalization failed" }
    & $Python $QA --stage final @QAPaths
    if ($LASTEXITCODE -ne 0) { throw "Final QA failed" }
}

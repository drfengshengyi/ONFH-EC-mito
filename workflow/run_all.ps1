[CmdletBinding()]
param(
    [string]$Python = $(if ($env:ONFH_PYTHON) { $env:ONFH_PYTHON } else { "python" }),
    [string]$Rscript = $(if ($env:ONFH_RSCRIPT) { $env:ONFH_RSCRIPT } else { "Rscript" }),
    [ValidateRange(0, 100000)] [int]$Permutations = 1000,
    [ValidateRange(1, 128)] [int]$Jobs = 4,
    [ValidateRange(1, 128)] [int]$Cores = 1
)

$ErrorActionPreference = "Stop"
& (Join-Path $PSScriptRoot "run_core_analysis.ps1") -Python $Python -Rscript $Rscript -Permutations $Permutations -Jobs $Jobs
if ($LASTEXITCODE -ne 0) { throw "Core analysis failed." }
& (Join-Path $PSScriptRoot "run_virtual_knockout.ps1") -Python $Python -Rscript $Rscript -Cores $Cores
if ($LASTEXITCODE -ne 0) { throw "Virtual knockout failed." }
& (Join-Path $PSScriptRoot "run_figures.ps1") -Python $Python -Rscript $Rscript
if ($LASTEXITCODE -ne 0) { throw "Figure generation failed." }

$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
$taskOut = Join-Path $taskRoot 'outputs/MFSD10_表达分布摸底_2026-10-26'
$env:LC_ALL = 'English_United States.utf8'
$env:LC_CTYPE = 'English_United States.utf8'
$env:OMP_NUM_THREADS = '2'
$env:OPENBLAS_NUM_THREADS = '2'
$env:TMPDIR = Join-Path $taskOut 'cache/r_tmp'
$env:TMP = $env:TMPDIR
$env:TEMP = $env:TMPDIR
foreach ($taskScript in @('00_inventory.R','01_extract_and_summarize.R','02_figures_and_checks.R','03_independent_check_and_digest.R','04_write_report.R')) {
    & 'D:/R4.5.3/R-4.5.3/bin/Rscript.exe' --vanilla --encoding=UTF-8 (Join-Path $PSScriptRoot $taskScript) $taskRoot
    if ($LASTEXITCODE -ne 0) { throw "Analysis failed: $taskScript" }
}
& (Join-Path $PSScriptRoot '05_finalize_audit.ps1')

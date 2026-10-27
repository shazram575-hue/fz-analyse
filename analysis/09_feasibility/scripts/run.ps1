$ErrorActionPreference='Stop'
$taskRoot=(Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
$env:LC_ALL='English_United States.utf8'
$env:LC_CTYPE='English_United States.utf8'
$env:OMP_NUM_THREADS='2'
$env:OPENBLAS_NUM_THREADS='2'
foreach($taskScript in @('01_patient_audit.R','02_assign_feasibility.R')) {
    & 'D:/R4.5.3/R-4.5.3/bin/Rscript.exe' --vanilla --encoding=UTF-8 (Join-Path $PSScriptRoot $taskScript) $taskRoot
    if($LASTEXITCODE -ne 0){throw "Feasibility audit failed: $taskScript"}
}
# This runner contains only audits; no NMF, cNMF, projection, or gene-program analysis.

$ErrorActionPreference='Stop'
$taskRoot=(Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
$taskOut=Join-Path $taskRoot 'outputs/MFSD10_稀疏与聚合探索_2026-10-26'
$taskPy=Join-Path $taskRoot 'outputs/MFSD10_程序发现与独立验证_2026-10-26/env/Scripts/python.exe'
$env:PYTHONUTF8='1';$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1'
if(Test-Path -LiteralPath (Join-Path $taskOut 'audit/COMPLETED.json')){Write-Output 'Complete. Stop rule applies; no automatic rerun.';exit 0}
foreach($taskScript in @('01_audit_reuse.py','02_build_units.py','03_estimate_effects.py','04_compare_report.py','05_final_audit.py')){
 & $taskPy -u (Join-Path $PSScriptRoot $taskScript) $taskRoot
 if($LASTEXITCODE -ne 0){throw "Failed: $taskScript"}
}

$ErrorActionPreference='Stop'
$taskRoot=(Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
$taskOut=Join-Path $taskRoot 'outputs/MFSD10_程序发现与独立验证_2026-10-26'
$taskPython=Join-Path $taskOut 'env/Scripts/python.exe'
$taskBundledPython='C:/Users/fangz/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$env:PYTHONUTF8='1';$env:LC_ALL='English_United States.utf8';$env:LC_CTYPE='English_United States.utf8'
$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
$env:NUMBA_CACHE_DIR=Join-Path $taskOut 'cache/numba';$env:TEMP=Join-Path $taskOut 'cache/tmp';$env:TMP=$env:TEMP
function Invoke-TaskPython([string]$taskScript){
    & $taskPython -u (Join-Path $PSScriptRoot $taskScript) $taskRoot
    if($LASTEXITCODE -ne 0){throw "Stage failed: $taskScript; inspect logs and resume successful checkpoints"}
}
function Invoke-TaskR([string]$taskScript){
    & 'D:/R4.5.3/R-4.5.3/bin/Rscript.exe' --vanilla --encoding=UTF-8 (Join-Path $PSScriptRoot $taskScript) $taskRoot
    if($LASTEXITCODE -ne 0){throw "R stage failed: $taskScript"}
}
if(Test-Path -LiteralPath (Join-Path $taskOut 'audit/PIPELINE_COMPLETED.json')) {
    Write-Output 'Pipeline already complete. See report and final audit; no automatic rerun.'
    exit 0
}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'metadata/Supplement_workbook_inventory.json'))) {
    & $taskBundledPython (Join-Path $PSScriptRoot '00_read_supporting_sources.py') $taskRoot
    if($LASTEXITCODE -ne 0){throw 'Source extraction failed'}
}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'cache/GSE161529_singlet_counts_gene_by_cell.mtx'))){Invoke-TaskR '01_export_discovery.R'}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'audit/GSE176078_import_audit.json'))){Invoke-TaskPython '02_import_validation.py'}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'objects/GSE176078_malignant_candidates_author_counts_clinical_v2.h5ad'))){Invoke-TaskPython '02b_reconcile_validation_clinical.py'}
Invoke-TaskPython '03_prepare_and_run_cnmf.py'
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'frozen/Programs_v1.json'))){Invoke-TaskPython '04_aggregate_programs.py'}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'cache/GSE161529_three_patients_all_confirmed.mtx'))){Invoke-TaskR '01b_export_doublet_sensitivity.R'}
Invoke-TaskPython '07_parameter_sensitivity.py'
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'discovery/Convergence_limit_followup.tsv'))){Invoke-TaskPython '07b_convergence_followup.py'}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'frozen/Analysis_freeze_v1.json'))){Invoke-TaskPython '05_discovery_associations_and_freeze.py'}
if(!(Test-Path -LiteralPath (Join-Path $taskOut 'audit/External_validation_summary.json'))){Invoke-TaskPython '06_external_validation.py'}
Invoke-TaskPython '08_figures_and_report.py'
Invoke-TaskPython '09_final_audit.py'

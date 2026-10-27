$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path (Join-Path $PSScriptRoot '../../..')).Path
$taskOut = Join-Path $taskRoot 'outputs/MFSD10_表达分布摸底_2026-10-26'
$taskValidation = Import-Csv -LiteralPath (Join-Path $taskOut 'audit/Analysis_validation.tsv') -Delimiter "`t"
if ($taskValidation.Count -ne 17 -or @($taskValidation | Where-Object pass -ne 'TRUE').Count -gt 0) { throw 'Analysis validation failed' }
$taskSpot = Import-Csv -LiteralPath (Join-Path $taskOut 'audit/Direct_Seurat_RNA_counts_spotcheck.tsv') -Delimiter "`t"
if ($taskSpot.Count -ne 4 -or @($taskSpot | Where-Object { $_.target_counts_identical -ne 'TRUE' -or $_.RNA_totals_identical -ne 'TRUE' -or $_.nFeature_identical -ne 'TRUE' }).Count -gt 0) { throw 'Direct RNA spot-check failed' }
$taskBaseline = Import-Csv -LiteralPath (Join-Path $taskRoot 'outputs/第一阶段_2026-10-23/audit/Final_artifact_SHA256_manifest.tsv') -Delimiter "`t"
$taskOld = foreach ($taskItem in $taskBaseline) {
    if ($taskItem.relative_path -like 'outputs*') {
        $taskActual = (Get-FileHash -LiteralPath (Join-Path $taskRoot $taskItem.relative_path) -Algorithm SHA256).Hash
        [pscustomobject]@{relative_path=$taskItem.relative_path; baseline_sha256=$taskItem.sha256; current_sha256=$taskActual; unchanged=($taskActual -eq $taskItem.sha256)}
    }
}
$taskOld | Export-Csv -LiteralPath (Join-Path $taskOut 'audit/Original_stage1_SHA256_preservation.tsv') -Delimiter "`t" -NoTypeInformation
if ($taskOld.Count -ne 86 -or @($taskOld | Where-Object unchanged -ne $true).Count -gt 0) { throw 'Frozen results changed' }
$taskSources = Import-Csv -LiteralPath (Join-Path $taskOut 'audit/Read_only_input_manifest.tsv') -Delimiter "`t"
$taskHashes = foreach ($taskItem in $taskSources) {
    $taskInfo = Get-Item -LiteralPath $taskItem.path
    if ($taskInfo.Length -ne [long]$taskItem.bytes) { throw "Input size changed: $($taskItem.path)" }
    [pscustomobject]@{path=$taskItem.path; bytes=$taskInfo.Length; sha256=(Get-FileHash -LiteralPath $taskItem.path -Algorithm SHA256).Hash}
}
$taskHashes | Export-Csv -LiteralPath (Join-Path $taskOut 'audit/Read_only_input_SHA256.tsv') -Delimiter "`t" -NoTypeInformation
$taskPDFs = @(Get-ChildItem -LiteralPath (Join-Path $taskOut 'figures') -Filter '*.pdf')
$taskPNGs = @(Get-ChildItem -LiteralPath (Join-Path $taskOut 'figures') -Filter '*.png')
if ($taskPDFs.Count -ne 6 -or $taskPNGs.Count -ne 6) { throw 'Expected six figure pairs' }
$taskPDFChecks = foreach ($taskPDF in $taskPDFs) {
    $taskBytes = [System.IO.File]::ReadAllBytes($taskPDF.FullName)
    $taskHead = [System.Text.Encoding]::ASCII.GetString($taskBytes,0,5)
    $taskTail = [System.Text.Encoding]::ASCII.GetString($taskBytes,[Math]::Max(0,$taskBytes.Length-1024),[Math]::Min(1024,$taskBytes.Length))
    if ($taskHead -ne '%PDF-' -or $taskTail -notmatch '%%EOF') { throw "Incomplete PDF: $($taskPDF.Name)" }
    [pscustomobject]@{file=$taskPDF.Name; bytes=$taskPDF.Length; header_ok=$true; eof_ok=$true}
}
$taskPDFChecks | Export-Csv -LiteralPath (Join-Path $taskOut 'audit/PDF_integrity.tsv') -Delimiter "`t" -NoTypeInformation
$taskDeliverables = foreach ($taskArtifact in Get-ChildItem -LiteralPath $taskOut -Recurse -File) {
    if ($taskArtifact.FullName -match '\\cache\\' -or $taskArtifact.Name -in @('Deliverable_SHA256.tsv','EXPRESSION_SURVEY_COMPLETED.json')) { continue }
    [pscustomobject]@{relative_path=$taskArtifact.FullName.Substring($taskOut.Length+1); bytes=$taskArtifact.Length; sha256=(Get-FileHash -LiteralPath $taskArtifact.FullName -Algorithm SHA256).Hash}
}
$taskDeliverables | Export-Csv -LiteralPath (Join-Path $taskOut 'audit/Deliverable_SHA256.tsv') -Delimiter "`t" -NoTypeInformation
$taskComplete = [ordered]@{status='complete'; completed_date='2026-10-26'; QC_cells=340979; libraries=69; patient_keys=54; malignant_cells=9486; malignant_patients=32; malignant_patients_min20=27; analysis_checks_passed=17; direct_Seurat_checks=4; protected_stage1_files_unchanged=86; figure_pairs=6; input_files_hashed=$taskHashes.Count; new_DE_or_mechanism_analysis=$false; background_task_remaining=$false}
$taskComplete | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $taskOut 'audit/EXPRESSION_SURVEY_COMPLETED.json') -Encoding utf8
$taskComplete | ConvertTo-Json

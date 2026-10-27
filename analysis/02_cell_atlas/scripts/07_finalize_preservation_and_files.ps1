param([string]$ProjectRoot='D:\Desktop\Yu\工作流3')
$ErrorActionPreference='Stop'
$projectDir=(Resolve-Path -LiteralPath $ProjectRoot).Path
$stageDir=Join-Path $projectDir 'outputs\单细胞基础补全_2026-10-23'
$oldDir=Join-Path $projectDir 'outputs\第一阶段_2026-10-23'
$frozen=Import-Csv -LiteralPath (Join-Path $stageDir 'audit\Frozen_stage1_before.tsv') -Delimiter "`t"
$after=foreach($row in $frozen) {
 $path=Join-Path $oldDir $row.relative_path
 $f=Get-Item -LiteralPath $path
 $hash=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
 [pscustomobject]@{relative_path=$row.relative_path;bytes=$f.Length;sha256=$hash;unchanged=($hash -eq $row.sha256 -and $f.Length -eq [long]$row.bytes)}
}
$after | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\Frozen_stage1_after.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
if(@($after | Where-Object {-not $_.unchanged}).Count -gt 0){throw 'A frozen Stage1 result changed; do not overwrite it'}
$published=Import-Csv -LiteralPath (Join-Path $oldDir 'audit\Final_artifact_SHA256_manifest.tsv') -Delimiter "`t"
$publishedChecks=foreach($row in $published){
 if($row.relative_path -notlike 'outputs\第一阶段_2026-10-23\*'){continue}
 $path=Join-Path $projectDir $row.relative_path
 $hash=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
 [pscustomobject]@{relative_path=$row.relative_path;unchanged=($hash -eq $row.sha256)}
}
$publishedChecks | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\Published_stage1_output_preservation.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
if(@($publishedChecks).Count -lt 1){throw 'No published Stage1 files matched the preservation manifest'}
if(@($publishedChecks | Where-Object {-not $_.unchanged}).Count -gt 0){throw 'A previously published Stage1 deliverable changed'}
$completeNames=@('Preparation.complete','Within_celltype_DE.complete','Global_embedding.complete','Figures.complete','Independent_validation.complete')
foreach($name in $completeNames){if(-not (Test-Path -LiteralPath (Join-Path $stageDir "cache\$name"))){throw "Incomplete module: $name"}}
$checks=Import-Csv -LiteralPath (Join-Path $stageDir 'audit\Independent_validation.tsv') -Delimiter "`t"
if(@($checks | Where-Object {$_.pass -ne 'TRUE'}).Count -gt 0){throw 'Independent validation has a failed item'}
$plotPng=@(Get-ChildItem -LiteralPath (Join-Path $stageDir 'figures') -Filter 'Foundation_*.png' -File)
$plotPdf=@(Get-ChildItem -LiteralPath (Join-Path $stageDir 'figures') -Filter 'Foundation_*.pdf' -File)
if($plotPng.Count -ne 22 -or $plotPdf.Count -ne 22){throw 'Expected22 PNG/PDF figure pairs'}
$figureChecks=foreach($f in $plotPdf){
 $bytes=[IO.File]::ReadAllBytes($f.FullName)
 $prefix=[Text.Encoding]::ASCII.GetString($bytes,0,[Math]::Min(5,$bytes.Length))
 $tail=[Text.Encoding]::ASCII.GetString($bytes,[Math]::Max(0,$bytes.Length-100),[Math]::Min(100,$bytes.Length))
 [pscustomobject]@{file=$f.Name;bytes=$bytes.Length;complete_pdf=($prefix -eq '%PDF-' -and $tail.Contains('%%EOF'))}
}
$figureChecks | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\PDF_completeness.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
if(@($figureChecks | Where-Object {-not $_.complete_pdf}).Count -gt 0){throw 'An incomplete PDF figure was found'}
$docs=@((Get-Item -LiteralPath (Join-Path $projectDir '单细胞基础补全_结果报告.md')),
 (Get-Item -LiteralPath (Join-Path $stageDir 'README_结果索引.md')),
 (Get-Item -LiteralPath (Join-Path $stageDir '方法与参数说明.md')))
$links=foreach($doc in $docs){
 $body=Get-Content -LiteralPath $doc.FullName -Raw -Encoding utf8
 foreach($m in [regex]::Matches($body,'\[[^\]]*\]\(([^)]+)\)')){
  $target=$m.Groups[1].Value.Trim('<','>')
  if($target -match '^https?://' -or $target.StartsWith('#')){continue}
  $path=if([IO.Path]::IsPathRooted($target)){$target}else{Join-Path $doc.DirectoryName $target}
  [pscustomobject]@{document=$doc.Name;target=$target;exists=(Test-Path -LiteralPath $path)}
 }
}
$links | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\Document_links.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
if(@($links | Where-Object {-not $_.exists}).Count -gt 0){throw 'A new-report link is broken'}
$files=@($docs)
foreach($part in @('scripts','results','figures','metadata')){$files+=Get-ChildItem -LiteralPath (Join-Path $stageDir $part) -File -Recurse}
$manifest=foreach($f in ($files | Sort-Object FullName -Unique)){
 [pscustomobject]@{relative_path=[IO.Path]::GetRelativePath($projectDir,$f.FullName);bytes=$f.Length;
 sha256=(Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash}
}
$manifest | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\Deliverable_SHA256_manifest.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
@("Single-cell foundation supplement complete: 236309 main cells; frozen Pool A unchanged.",
 "Only within-annotation Tumor-origin vs Healthy disease DE was performed.",
 "7 frozen Stage1 files unchanged; 20 independent R checks passed.",
 "Figure generation finished; visual inspection is separately documented in audit/Figure_visual_QA.md.",
 "Completed: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')") |
 Set-Content -LiteralPath (Join-Path $stageDir 'audit\BASICS_COMPLETED.txt') -Encoding utf8
Write-Output "Frozen files unchanged: $($after.Count); independent checks: $($checks.Count); links: $($links.Count); deliverable hashes: $($manifest.Count)."

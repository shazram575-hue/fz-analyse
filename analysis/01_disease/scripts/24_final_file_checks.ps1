param([string]$ProjectRoot='D:\Desktop\Yu\工作流3')
$ErrorActionPreference='Stop'
$projectDir=(Resolve-Path -LiteralPath $ProjectRoot).Path
$stageDir=Join-Path $projectDir 'outputs\第一阶段_2026-10-23'
$docFiles=@(Get-ChildItem -LiteralPath $projectDir -Filter '*.md' -File)
$docFiles+=Get-Item -LiteralPath (Join-Path $stageDir 'README_结果索引.md')
$linkChecks=foreach($doc in $docFiles) {
 $body=Get-Content -LiteralPath $doc.FullName -Raw -Encoding utf8
 foreach($m in [regex]::Matches($body,'\[[^\]]*\]\(([^)]+)\)')) {
  $target=$m.Groups[1].Value.Trim('<','>')
  if($target -match '^https?://' -or $target.StartsWith('#')) {continue}
  $resolved=if([IO.Path]::IsPathRooted($target)) {$target} else {Join-Path $doc.DirectoryName $target}
  [pscustomobject]@{document=$doc.Name;target=$target;exists=(Test-Path -LiteralPath $resolved)}
 }
}
$linkChecks | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\Final_document_links.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
if(@($linkChecks | Where-Object {-not $_.exists}).Count -gt 0) {throw 'A local document link is broken'}
$artifacts=@($docFiles)
$artifacts+=Get-ChildItem -LiteralPath (Join-Path $stageDir 'scripts') -File
$artifacts+=Get-ChildItem -LiteralPath (Join-Path $stageDir 'results') -File
$artifacts+=Get-ChildItem -LiteralPath (Join-Path $stageDir 'figures') -File -Filter 'Final_*'
$artifacts+=Get-ChildItem -LiteralPath (Join-Path $stageDir 'audit') -File -Filter 'Final_*'
$manifest=foreach($f in ($artifacts | Sort-Object FullName -Unique)) {
 [pscustomobject]@{relative_path=[IO.Path]::GetRelativePath($projectDir,$f.FullName);bytes=$f.Length;
  sha256=(Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash}
}
$manifest | Export-Csv -LiteralPath (Join-Path $stageDir 'audit\Final_artifact_SHA256_manifest.tsv') -Delimiter "`t" -NoTypeInformation -Encoding utf8
@("Phase 1 complete: Pool A 813 protein-coding genes.",
 "Independent R checks and document links passed.",
 "Final figures PNG visually inspected; PDF counterparts saved.",
 "Original input files unchanged; no drug prediction, pathway selection, spatial analysis or MD performed.",
 "Completed: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')") |
 Set-Content -LiteralPath (Join-Path $stageDir 'audit\PHASE1_COMPLETED.txt') -Encoding utf8
Write-Output "Passed $($linkChecks.Count) local links; hashed $($manifest.Count) deliverable files."

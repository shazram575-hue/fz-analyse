$ErrorActionPreference = 'Stop'
$stageRoot = Split-Path -Parent (Split-Path -Parent $PSCommandPath)
$metadataRoot = Join-Path $stageRoot 'metadata'
$localFiles = Import-Csv -LiteralPath (Join-Path $metadataRoot 'TCGA_local_file_ids.tsv') -Delimiter "`t"
$gdcFields = 'file_id,file_name,file_size,md5sum,state,analysis.workflow_type,cases.case_id,cases.submitter_id,cases.project.project_id,cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.sample_type_id,cases.samples.tissue_type,cases.samples.portions.analytes.aliquots.submitter_id,cases.demographic.gender,cases.demographic.race,cases.diagnoses.age_at_diagnosis,cases.diagnoses.primary_diagnosis,cases.diagnoses.ajcc_pathologic_stage'
$allHits = [System.Collections.Generic.List[object]]::new()
for ($chunkStart = 0; $chunkStart -lt $localFiles.Count; $chunkStart += 100) {
    $chunkEnd = [Math]::Min($chunkStart + 99, $localFiles.Count - 1)
    $fileIds = @($localFiles[$chunkStart..$chunkEnd] | ForEach-Object { $_.file_id })
    $filters = @{op='in'; content=@{field='files.file_id'; value=$fileIds}} | ConvertTo-Json -Depth 8 -Compress
    $requestUri = 'https://api.gdc.cancer.gov/files?filters=' + [Uri]::EscapeDataString($filters) + '&fields=' + [Uri]::EscapeDataString($gdcFields) + '&size=100&format=JSON'
    $pageNumber = [int]($chunkStart / 100) + 1
    $responsePath = Join-Path $metadataRoot ('GDC_files_response_{0:D2}.json' -f $pageNumber)
    if (Test-Path -LiteralPath $responsePath) {
        $responseText = Get-Content -LiteralPath $responsePath -Raw
    } else {
        $responseText = (Invoke-WebRequest -Uri $requestUri -TimeoutSec 90).Content
        $responseText | Set-Content -LiteralPath $responsePath -Encoding utf8
    }
    $responseObject = $responseText | ConvertFrom-Json
    if ($responseObject.data.pagination.total -ne $fileIds.Count) { throw "Incomplete metadata batch $pageNumber" }
    foreach ($hit in $responseObject.data.hits) { $allHits.Add($hit) }
    Write-Output "GDC batch $pageNumber retrieved $($responseObject.data.hits.Count) files"
}
$allHits | ConvertTo-Json -Depth 35 | Set-Content -LiteralPath (Join-Path $metadataRoot 'GDC_files_metadata_all.json') -Encoding utf8
Write-Output "TOTAL $($allHits.Count)"

$ErrorActionPreference = 'Stop'
$taskRoot = Split-Path -Parent $PSScriptRoot
$taskPython = Join-Path $taskRoot 'env\Scripts\python.exe'
$taskCache = Join-Path $taskRoot 'package_cache'
$taskLogs = Join-Path $taskRoot 'logs'
New-Item -ItemType Directory -Path $taskLogs -Force | Out-Null
$env:PYTHONUTF8 = '1'
$env:PIP_DISABLE_PIP_VERSION_CHECK = '1'
$env:OMP_NUM_THREADS = '2'
$env:MKL_NUM_THREADS = '2'
& $taskPython -m pip install --cache-dir $taskCache --log (Join-Path $taskLogs 'pip_torch.log') --index-url https://download.pytorch.org/whl/cpu 'torch==2.5.1+cpu'
if ($LASTEXITCODE -ne 0) { throw 'PyTorch installation failed; inspect pip_torch.log' }
& $taskPython -m pip install --cache-dir $taskCache --log (Join-Path $taskLogs 'pip_dependencies.log') 'numpy==1.26.4' 'pandas==2.2.3' 'scikit-learn==1.5.2' 'scipy==1.14.1' 'rdkit==2024.9.6' 'descriptastorus==2.8.0' 'subword-nmt==0.3.8' 'wget==3.2' 'lifelines==0.30.0' 'prettytable==3.12.0' 'tensorboard==2.18.0' 'matplotlib==3.9.4' 'tqdm==4.67.1' requests
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed; inspect pip_dependencies.log' }
& $taskPython -m pip install --no-deps --cache-dir $taskCache --log (Join-Path $taskLogs 'pip_deeppurpose.log') 'DeepPurpose==0.1.5'
if ($LASTEXITCODE -ne 0) { throw 'DeepPurpose installation failed; inspect pip_deeppurpose.log' }
& $taskPython -m pip freeze | Set-Content -LiteralPath (Join-Path $taskLogs 'pip_freeze.txt') -Encoding utf8
& $taskPython -c "import torch,rdkit; from DeepPurpose import DTI; print('IMPORT_OK',torch.__version__,rdkit.__version__, 'CUDA_AVAILABLE',torch.cuda.is_available())"
if ($LASTEXITCODE -ne 0) { throw 'Import check failed; environment not marked ready' }
'IMPORT_OK_ONLY_NOT_MODEL_VALIDATION' | Set-Content -LiteralPath (Join-Path $taskLogs 'environment_import.complete') -Encoding utf8

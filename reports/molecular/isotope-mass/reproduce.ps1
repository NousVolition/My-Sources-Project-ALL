param([int]$Workers = 4)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
python -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed' }
python test_science.py
if ($LASTEXITCODE -ne 0) { throw 'Scientific checks failed' }
python simulate.py --workers $Workers
if ($LASTEXITCODE -ne 0) { throw 'Simulation failed' }
python switch_assay.py --workers $Workers
if ($LASTEXITCODE -ne 0) { throw 'Mass-switch assay failed' }
python refine_response.py
if ($LASTEXITCODE -ne 0) { throw 'Response refinement failed' }
python analyze.py --workers $Workers
if ($LASTEXITCODE -ne 0) { throw 'Analysis failed' }
python make_report.py
if ($LASTEXITCODE -ne 0) { throw 'Report generation failed' }
python verify_results.py
if ($LASTEXITCODE -ne 0) { throw 'Artifact verification failed' }

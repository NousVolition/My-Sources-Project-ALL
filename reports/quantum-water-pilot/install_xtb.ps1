$ErrorActionPreference = 'Stop'
$taskTools = Join-Path $PSScriptRoot 'tools'
if (Test-Path -LiteralPath $taskTools) { throw 'tools already exists; use a clean directory.' }
New-Item -ItemType Directory -Path $taskTools | Out-Null
$taskArchive = Join-Path $taskTools 'xtb-windows.zip'
Invoke-WebRequest -Uri 'https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1pre-windows-x86_64.zip' -OutFile $taskArchive
$taskExpected = '043e578da4a7e114a4d584972959a875e3ffb9f2767a86723b95aa6719d28d9c'
if ((Get-FileHash -LiteralPath $taskArchive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskExpected) { throw 'Archive checksum mismatch; extraction stopped.' }
Expand-Archive -LiteralPath $taskArchive -DestinationPath $taskTools
Write-Output ('Verified xTB installed under ' + $taskTools)

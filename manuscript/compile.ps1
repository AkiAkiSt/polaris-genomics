param([string]$tex = "polaris.tex")
$tec = "C:\Users\Aki\OneDrive\Aki\Vs code\Carnegie Mellon\fp_test\data\external\tec_msvc\tectonic.exe"
$cache = "C:\Users\Aki\AppData\Local\Temp\claude\tectonic_cache"
New-Item -ItemType Directory -Force -Path $cache | Out-Null
$env:TECTONIC_CACHE_DIR = $cache
Set-Location "C:\Users\Aki\OneDrive\Aki\Vs code\Carnegie Mellon\fp_test\manuscript"
& $tec -X compile $tex > "tec_log.txt" 2>&1
$code = $LASTEXITCODE
$pdf = [System.IO.Path]::ChangeExtension($tex, ".pdf")
"exit=$code"
if (Test-Path $pdf) { "PDF OK: $pdf  " + [math]::Round((Get-Item $pdf).Length/1kb) + " KB" }
else { "NO PDF"; Get-Content "tec_log.txt" -Tail 16 }
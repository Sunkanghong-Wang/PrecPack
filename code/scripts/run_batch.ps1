$ErrorActionPreference = "Stop"

if (Get-Command py -ErrorAction SilentlyContinue) {
    & py -3 "$PSScriptRoot/run_batch.py" @args
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    & python "$PSScriptRoot/run_batch.py" @args
} else {
    throw "Python 3 was not found."
}
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

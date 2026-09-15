[CmdletBinding()]
param(
    [ValidateSet('archived', 'plot', 'check', 'raw', 'full', 'dqmc', 'ed', 'analyze', 'smoke')]
    [string]$Mode = 'full',
    [ValidateSet('both', 'number_conserving', 'pairing')]
    [string]$Model = 'both',
    [ValidateSet('all', 'main', 'supplement')]
    [string]$Scope = 'all',
    [string]$Output = '',
    [string]$PythonED = 'python3',
    [int]$Threads = 1,
    [switch]$Plan,
    [switch]$Resume
)

$ErrorActionPreference = 'Stop'
# Quote each argument for Bash; values are WSL paths, not native Windows paths.
function Quote-Bash([string]$Value) {
    $sq = [string][char]39
    $dq = [string][char]34
    return $sq + $Value.Replace($sq, $sq + $dq + $sq + $dq + $sq) + $sq
}
$arguments = @('python3', 'reproduce.py', '--mode', $Mode, '--model', $Model,
    '--scope', $Scope, '--python-ed', $PythonED, '--threads', [string]$Threads)
if ($Output) { $arguments += @('--output', $Output) }
if ($Plan) { $arguments += '--plan' }
if ($Resume) { $arguments += '--resume' }
$command = ($arguments | ForEach-Object { Quote-Bash $_ }) -join ' '
& (Join-Path $PSScriptRoot 'run-in-wsl.ps1') $command
exit $LASTEXITCODE

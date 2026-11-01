# Export the current editable F1 using Microsoft PowerPoint (no reconstruction).
param(
 [string]$Source = (Join-Path $PSScriptRoot 'figures/F1_overview.pptx'),
 [string]$OutputDir = (Join-Path $PSScriptRoot 'figures')
)
$ErrorActionPreference='Stop'
$sourcePath=(Resolve-Path -LiteralPath $Source).Path
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
$outputPath=(Resolve-Path -LiteralPath $OutputDir).Path
$app=New-Object -ComObject PowerPoint.Application
$deck=$app.Presentations.Open($sourcePath,$true,$false,$false)
try {
 if($deck.Slides.Count -ne 1){throw 'F1 must contain exactly one slide.'}
 $deck.SaveAs((Join-Path $outputPath 'F1_overview.pdf'),32)
 $previewHeight=[int][Math]::Round(1920*$deck.PageSetup.SlideHeight/$deck.PageSetup.SlideWidth)
 $deck.Slides.Item(1).Export((Join-Path $outputPath 'F1_overview.png'),'PNG',1920,$previewHeight)
 Write-Output "Exported F1 from $sourcePath"
} finally {
 $deck.Close()
 # Do not quit PowerPoint: the author may have other presentations open.
}

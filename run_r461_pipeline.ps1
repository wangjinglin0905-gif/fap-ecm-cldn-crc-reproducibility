$ErrorActionPreference = "Stop"

$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Rscript = if ($env:FAP_R461) {
    $env:FAP_R461
} else {
    "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
}

if (-not (Test-Path -LiteralPath $Rscript)) {
    throw "R 4.6.1 Rscript was not found. Set FAP_R461 to its full path."
}

& $Rscript --vanilla -e "if (as.character(getRversion()) != '4.6.1') stop('R 4.6.1 is required')"
if ($LASTEXITCODE -ne 0) {
    throw "The configured Rscript is not R 4.6.1."
}

$Stages = @(
    @{ Label = "01_tcga_l0"; Script = "work/reproducibility/scripts/01_tcga_l0.R" },
    @{ Label = "03_spatial_l2"; Script = "work/reproducibility/scripts/03_spatial_l2.R" },
    @{ Label = "cellchat_conventional"; Script = "work/reproducibility/scripts/cellchat_legacy/run_cellchat_reanalysis.R" },
    @{ Label = "cellchat_smc20_sensitivity"; Script = "work/reproducibility/scripts/cellchat_legacy/run_cellchat_smc20_sensitivity.R" },
    @{ Label = "04_nichenet_l3"; Script = "work/reproducibility/scripts/04_nichenet_l3.R" },
    @{ Label = "05_immune_l4"; Script = "work/reproducibility/scripts/05_immune_l4.R" },
    @{ Label = "07a_select_qi_hvg"; Script = "work/reproducibility/scripts/07a_select_qi_hvg.R" },
    @{ Label = "07c_harmony_integration"; Script = "work/reproducibility/scripts/07c_harmony_integration.R" },
    @{ Label = "08_rppa_validation"; Script = "work/reproducibility/scripts/08_rppa_validation.R" },
    @{ Label = "09_pathway_tf_l3"; Script = "work/reproducibility/scripts/09_pathway_tf_l3.R" },
    @{ Label = "10_cptac_protein_validation"; Script = "work/reproducibility/scripts/10_cptac_protein_validation.R" },
    @{ Label = "12b_gse166555_statistics"; Script = "work/reproducibility/scripts/12b_gse166555_statistics.R" },
    @{ Label = "13_valdeolivas_spatial"; Script = "work/reproducibility/scripts/13_valdeolivas_spatial_validation.R" },
    @{ Label = "14_spatial_cellchat"; Script = "work/reproducibility/scripts/14_spatial_cellchat_validation.R" },
    @{ Label = "14b_spatial_cellchat_postprocess"; Script = "work/reproducibility/scripts/14b_spatial_cellchat_postprocess.R" },
    @{ Label = "15_ualcan_cptac_parser"; Script = "work/reproducibility/scripts/15_ualcan_cptac_parser.R" },
    @{ Label = "17_tcga_primary_sensitivity"; Script = "work/reproducibility/scripts/17_tcga_primary_tumour_sensitivity.R" },
    @{ Label = "20_regenerate_all_figures"; Script = "work/reproducibility/scripts/20_regenerate_all_figures.R" }
)

$LogDirectory = Join-Path $ProjectRoot "work/reproducibility/run_logs"
New-Item -ItemType Directory -Force -Path $LogDirectory | Out-Null
$ManifestRows = @()

Push-Location $ProjectRoot
try {
    foreach ($Stage in $Stages) {
        $Started = Get-Date
        $LogPath = Join-Path $LogDirectory ($Stage.Label + ".log")
        & $Rscript --vanilla $Stage.Script 2>&1 | Tee-Object -FilePath $LogPath
        $ExitCode = $LASTEXITCODE
        $Ended = Get-Date
        $ManifestRows += [pscustomobject]@{
            label = $Stage.Label
            script = $Stage.Script
            started = $Started.ToString("o")
            ended = $Ended.ToString("o")
            duration_seconds = [math]::Round(($Ended - $Started).TotalSeconds, 3)
            exit_code = $ExitCode
            log = (Resolve-Path $LogPath).Path.Replace($ProjectRoot + "\", "")
        }
        if ($ExitCode -ne 0) {
            throw "Stage $($Stage.Label) failed with exit code $ExitCode."
        }
    }
} finally {
    Pop-Location
}

$ManifestRows | Export-Csv -NoTypeInformation -Encoding utf8 -Path (Join-Path $LogDirectory "run_manifest.csv")
Write-Host "Completed $($Stages.Count) R 4.6.1 stages."

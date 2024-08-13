# docs/make.jl
push!(LOAD_PATH, "../src/")
using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    modules = [MetaCommunityMetrics],
    format = Documenter.HTML(),
    pages = ["Home" => "index.md",
                "Beta Diverisity Functions" => "BetaDiversity.md",
                "Dispersal-Niche Continuum Index (DNCI) Functions" => "DNCI.md",
                "Niche Overlap Index Function" => "NicheOverlapIndex.md",
                "Occupied Patches Proportion Function" => "OccupiedPatchesProportion.md",
                "Variability Metrics Functions" => "VariabilityMetrics.md"]
)
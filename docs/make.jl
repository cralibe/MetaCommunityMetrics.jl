# docs/make.jl
push!(LOAD_PATH, "../src/")
using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    repo = "https://github.com/cralibe/MetaCommunityMetrics.jl", 
    modules = [MetaCommunityMetrics],
    clean = true,
    doctest = true,  
    format = Documenter.HTML(),
    logo = "assets/MetaCommunityMetrics_logo.png",  
    pages = ["Home" => "index.md",
                "Beta Diverisity Functions" => "BetaDiversity.md",
                "Dispersal-Niche Continuum Index (DNCI) Functions" => "DNCI.md",
                "Niche Overlap Index Function" => "NicheOverlapIndex.md",
                "Occupied Patches Proportion Function" => "OccupiedPatchesProportion.md",
                "Variability Metrics Functions" => "VariabilityMetrics.md"]
)
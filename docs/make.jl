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
    format = Documenter.HTML(edit_link = "main"),
    pages = ["Home" => "index.md",
                "Beta Diverisity Functions" => "BetaDiversity.md",
                "Dispersal-Niche Continuum Index (DNCI) Functions" => "DNCI.md",
                "Niche Overlap Index Function" => "NicheOverlapIndex.md",
                "Occupied Patches Proportion Function" => "OccupiedPatchesProportion.md",
                "Variability Metrics Functions" => "VariabilityMetrics.md",
                "Hypervolume Functions" => "Hypervolume.md",
                "Benchmarking Julia vs R" => "Benchmarking.md",
                "Acknowledgment" => "Acknowledgment.md"]
)
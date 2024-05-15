using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    format = Documenter.HTML(),
    modules = [MetaCommunityMetrics],
    pages = [
        "Home" => "index.md",
        "Beta Diversity" => "BetaDiversity.md"
    ]
)

deploydocs(
    repo = "github.com/cralibe/MetaCommunityMetrics.jl.git",
)

using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics",
    format = Documenter.HTML(),
    modules = [MetaCommunityMetrics],
    pages = [
        "Home" => "index.md",
        "Beta Diversity" => "BetaDiversity.md",
        "Dispersal–niche continuum index" => "DNCI.md"
    ],
    root = "docs/src"
)

if haskey(ENV, "CI")
    deploydocs(
        repo = "github.com/cralibe/MetaCommunityMetrics.jl.git",
        branch = "gh-pages",
        target = ".",  # Deploy to the root of the gh-pages branch
        deploy_config = Dict("GITHUB_TOKEN" => ENV["CRALIBE_TOKEN_1"])
    )
end

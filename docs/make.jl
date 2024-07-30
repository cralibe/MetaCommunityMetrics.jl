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
    ]
)

if haskey(ENV, "CI")
    deploydocs(
        repo = "github.com/cralibe/MetaCommunityMetrics.jl.git",  # Replace with your repo
        branch = "gh-pages",
        target = "site",  # Deploy to the root of the gh-pages branch
        deploy_config = Dict("GITHUB_TOKEN" => ENV["CRALIBE_TOKEN_1"])
    )
end
using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    format = Documenter.HTML(),
    modules = [MetaCommunityMetrics],
    pages = [
        "Home" => "index.md",
        "Beta Diversity" => "BetaDiversity.md"
        "Dispersal–niche continuum index" => "DNCI.md"
    ]
)

# Conditionally deploy only if in CI environment
if haskey(ENV, "CI")
    deploydocs(
        repo = "github.com/username/MetaCommunityMetrics.jl.git",
        target = "site",
        branch = "gh-pages",
        deploy_config = Dict("GITHUB_TOKEN" => ENV["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"]),
    )
end
using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    format = Documenter.HTML(),
    modules = [MetaCommunityMetrics],
    pages = [
        "Home" => "index.md"
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

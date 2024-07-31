push!(LOAD_PATH,"../src/")

# docs/make.jl

using Documenter, MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    modules = [MetaCommunityMetrics],
    format = Documenter.HTML()
)
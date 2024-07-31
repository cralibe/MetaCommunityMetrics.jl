# docs/make.jl
push!(LOAD_PATH, "../src/")
using Documenter
using MetaCommunityMetrics

makedocs(
    sitename = "MetaCommunityMetrics Documentation",
    modules = [MetaCommunityMetrics],
    format = Documenter.HTML()
)
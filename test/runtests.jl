using MetaCommunityMetrics
using Test

@testset "MetaCommunityMetrics.jl" begin
    @test MetaCommunityMetrics.greet_your_package_name() == "Hello MetaCommunityMetrics!"
    @test MetaCommunityMetrics.greet_your_package_name() != "Hello world!"
end

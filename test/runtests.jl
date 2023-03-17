using testProject
using Test
using TestItemRunner

@run_package_tests
@testset "testProject.jl" begin
    # Write your tests here.
    @test 1 + 1 < 4
end

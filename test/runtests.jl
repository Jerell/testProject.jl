using testProject
using Unitful
using Test

@testset "testProject.jl" begin
    # pressure = 16u"bar"
    # temperature = 16u"°C"
    # flowrate = 120u"kg/s"

    # @show a₀ = 4.0531u"bar*L^2*mol^-2"
    @show typeof(compoundparams[1].b)
    # @show f = Fluid(pressure, temperature, flowrate)
    # @show compoundparams
    @test 1 == 1
end

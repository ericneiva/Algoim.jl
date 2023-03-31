using Algoim
using Test

@testset "Algoim.jl" begin

  @time @testset "Algoim" begin include("AlgoimTests.jl") end

end

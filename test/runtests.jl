using Algoim
using Test

@testset "Algoim.jl" begin

  @time @testset "Quadrature" begin include("QuadratureTests.jl") end
  @time @testset "Closest-Point-Projections" begin include("ClosestPointTests.jl") end

end

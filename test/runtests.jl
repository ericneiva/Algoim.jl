using Algoim
using Test

@testset "Algoim.jl" begin

  @time @testset "Quadrature Single Polynomial" begin include("QuadratureSinglePolynomialTests.jl") end

  @time @testset "Quadrature Dual Polynomial" begin include("QuadratureDualPolynomialTests.jl") end
  
  @time @testset "Closest-Point-Projections" begin include("ClosestPointTests.jl") end

end

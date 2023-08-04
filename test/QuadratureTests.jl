module QuadratureTests

using Test
using CxxWrap
using StaticArrays
using Algoim

const IN = -1
const OUT = 1
const CUT = 0

function run_case(T,degree,phase)

  u(x) = (x[1]*x[1]/(1.0*1.0)+x[2]*x[2]/(0.5*0.5)+x[3]*x[3]/(0.33*0.33)) - 1.0
  function gradu(x::V) where {V}
    V([2.0*x[1]/(1.0*1.0),2.0*x[2]/(0.5*0.5),2.0*x[3]/(0.33*0.33)])
  end

  phi = AlgoimCallLevelSetFunction(u,gradu)

  xmin = T[-1.1,-1.1,-1.1]
  xmax = T[1.1,1.1,1.1]

  _, weights = fill_quad_data(phi,xmin,xmax,phase,degree)

  s = sum(weights)
  if phase == IN
    @test s ≈ 0.7494442830903447
  elseif phase == CUT
    @test s ≈ 4.7770866154803056
  else
    error()
  end

end

f(x) = (x[1]*x[1]/(1.0*1.0)+x[2]*x[2]/(0.5*0.5)+x[3]*x[3]/(0.33*0.33)) - 1.0
function gradf(x::V) where {V}
  V([2.0*x[1]/(1.0*1.0),2.0*x[2]/(0.5*0.5),2.0*x[3]/(0.33*0.33)])
end

phi = AlgoimCallLevelSetFunction(f,gradf)
@test phi([1.0,0.5,0.33]) == 2.0
@test phi(SA[1.0,0.5,0.33]) == 2.0
@test all(phi.∇φ([1.0,0.5,0.33]) .≈ [2.0, 4.0, 6.0606060606060606])
@test all(phi.∇φ(SA[1.0,0.5,0.33]) .≈ [2.0, 4.0, 6.0606060606060606])

run_case(Float64,3,IN)
run_case(Float64,3,CUT)

run_case(SA,3,IN)
run_case(SA,3,CUT)

end # module

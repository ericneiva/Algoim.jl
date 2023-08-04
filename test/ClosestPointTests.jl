module ClosestPointTests

using Test
using CxxWrap
using StaticArrays
using Algoim

const IN = -1
const OUT = 1
const CUT = 0

function run_case(T,degree)

  u(x) = (x[1]*x[1]/(0.5*0.5)+x[2]*x[2]/(0.5*0.5)+x[3]*x[3]/(0.5*0.5)) - 1.0
  function gradu(x::V) where {V}
    V([2.0*x[1]/(0.5*0.5),2.0*x[2]/(0.5*0.5),2.0*x[3]/(0.5*0.5)])
  end

  phi = AlgoimCallLevelSetFunction(u,gradu)

  xmin = T[-1.1,-1.1,-1.1]
  xmax = T[1.1,1.1,1.1]
  partition = Int32[16,16,16]

  coords = fill_cpp_data(phi,partition,xmin,xmax,degree)
  @test maximum(u.(coords)) < 1.0e-4

end

run_case(Float64,2)
run_case(Float64,3)
run_case(Float64,4)
run_case(Float64,5)
run_case(Float64,-1)
run_case(SA,2)

end # module

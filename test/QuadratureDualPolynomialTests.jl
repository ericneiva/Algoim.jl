module TMP

using Test
using CxxWrap
using StaticArrays
using Algoim

const IN  = -1
const OUT =  1
const CUT =  0

u¹(x) = ( (x[1]-0.25) * (x[1]-0.25) +
          (x[2]-0.25) * (x[2]-0.25) + 
          (x[3]-0.25) * (x[3]-0.25) ) / ( 0.20 * 0.20 ) - 1.0

function ∇u¹(x::V) where {V}
  V([ 2.0*(x[1]-0.25)/(0.20*0.20), 2.0*(x[2]-0.25)/(0.20*0.20), 2.0*(x[3]-0.25)/(0.20*0.20) ])
end

u²(x) = ( (x[1]-0.75) * (x[1]-0.75) +
          (x[2]-0.75) * (x[2]-0.75) + 
          (x[3]-0.75) * (x[3]-0.75) ) / ( 0.22 * 0.22 ) - 1.0

function ∇u²(x::V) where {V}
  V([ 2.0*(x[1]-0.75)/(0.22*0.22), 2.0*(x[2]-0.75)/(0.22*0.22), 2.0*(x[3]-0.75)/(0.22*0.22) ])
end

u³(x) = ( (x[1]-0.50) * (x[1]-0.50) +
          (x[2]-0.50) * (x[2]-0.50) + 
          (x[3]-0.50) * (x[3]-0.50) ) / ( 0.40 * 0.40 ) - 1.0

function ∇u³(x::V) where {V}
  V([ 2.0*(x[1]-0.50)/(0.40*0.40), 2.0*(x[2]-0.50)/(0.40*0.40), 2.0*(x[3]-0.50)/(0.40*0.40) ])
end

u⁴(x) = ( (x[1]-0.50) * (x[1]-0.50) +
          (x[2]-0.50) * (x[2]-0.50) + 
          (x[3]-0.50) * (x[3]-0.50) ) / ( 0.20 * 0.20 ) - 1.0

function ∇u⁴(x::V) where {V}
  V([ 2.0*(x[1]-0.50)/(0.20*0.20), 2.0*(x[2]-0.50)/(0.20*0.20), 2.0*(x[3]-0.50)/(0.20*0.20) ])
end

function run_case(T,degree,phase¹,phase²,u,∇u,v,∇v)

  φ¹ = AlgoimCallLevelSetFunction(u,∇u)
  φ² = AlgoimCallLevelSetFunction(v,∇v)

  xmin = T[0.0,0.0,0.0]
  xmax = T[1.1,1.1,1.1]

  coords, weights = fill_quad_data(φ¹,φ²,xmin,xmax,phase¹,phase²,degree)

  s = sum(weights)
  if ( phase¹ == IN ) && ( phase² == OUT )
    @test s ≈ 0.03633669251347109
  elseif ( phase¹ == OUT ) && ( phase² == IN )
    @test s ≈ 0.048364137735430045
  elseif ( phase¹ == IN ) && ( phase² == CUT )
    @test s ≈ 0.11481813064130597
  elseif ( phase¹ == OUT ) && ( phase² == CUT )
    @test s ≈ 1.9536966722642424
  end

end

run_case(SA,3,IN,OUT,u¹,∇u¹,u²,∇u²)
run_case(SA,3,OUT,IN,u¹,∇u¹,u²,∇u²)

run_case(SA,3,IN,CUT,u¹,∇u¹,u³,∇u³)
run_case(SA,3,OUT,CUT,u¹,∇u¹,u³,∇u³)

end # module

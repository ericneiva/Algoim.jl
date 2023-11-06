module Algoim

# Load algoimWrapper module

using CxxWrap

using algoimWrapper_jll
if !isdefined(algoimWrapper_jll, :libalgoimwrapper_path)
  error("algoimWrapper_jll not available on this platform")
end

@wrapmodule(libalgoimwrapper)
__init__() = @initcxx

export LevelSetFunction
export ClosureLevelSet
export JuliaFunction2DLevelSet
export JuliaFunction3DLevelSet

to_uvector(x) = to_uvector(x,Val{length(x)}())
to_uvector(x,::Val{2}) = to_2D_uvector(x)
to_uvector(x,::Val{3}) = to_3D_uvector(x)

export to_uvector
export AlgoimUvector

@inline to_const_array(p) = ConstArray(data_array(p),length(p))

export to_const_array

# Low-level methods to fill quadrature data

## Implementation depends on elementary type

using StaticArrays # SVector
using LinearAlgebra: norm

const IN = -1
const OUT = 1
const CUT = 0

## Specialise inline methods for other elementary types

@inline to_array(x::Vector{T}) where {T} = x
@inline to_array(x::SVector{T,N}) where {T,N} = collect(x)

export to_array

# AlgoimCallLevelSetFunction

function julia_function_wrap(x,id::Float32,params::Ptr{Cvoid})
  f = unsafe_pointer_to_objref(params)::Function
  _x = to_const_array(x)
  f(_x)
end

struct AlgoimCallLevelSetFunction{A,B,C,D} <: LevelSetFunction
  φ::A
  ∇φ::B
  cache_φ::C
  cache_∇φ::D
end

function AlgoimCallLevelSetFunction(φ::Function,∇φ::Function)
  AlgoimCallLevelSetFunction{typeof(φ),typeof(∇φ),typeof(nothing),typeof(nothing)}(φ,∇φ,nothing,nothing)
end

## Evaluate differential operators of the level set function

(ls::AlgoimCallLevelSetFunction)(x) = ls.φ(x)

gradient(ls::AlgoimCallLevelSetFunction) = ls.∇φ

function normal(ls::AlgoimCallLevelSetFunction,x,cell_id::Int=1)
  gx = ls.∇φ(x)
  gx/norm(gx)
end

function normal(phi::AlgoimCallLevelSetFunction,x::AbstractVector{<:AbstractVector},cell_id::Int=1)
  map(xi->normal(phi,xi,cell_id),x)
end

export normal
export AlgoimCallLevelSetFunction

## Generic interface

fill_quad_data(phi::AlgoimCallLevelSetFunction,xmin::V,xmax::V,phase::Int,degree::Int,cell_id::Int=1) where {V} =
  fill_quad_data(phi,xmin,xmax,phase,degree,cell_id,Val{length(xmin)}())

function fill_quad_data(phi,xmin,xmax,phase,degree,cell_id,::Val{2})
  ls_value_wrap_2D_c = @safe_cfunction(julia_function_wrap,Float64,(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  ls_grad_wrap_2D_c = @safe_cfunction(julia_function_wrap,ConstCxxRef{AlgoimUvector{Float64,2}},(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  cpp_f = ClosureLevelSet{Int32(2)}(ls_value_wrap_2D_c,phi.φ)
  cpp_g = ClosureLevelSet{Int32(2)}(ls_grad_wrap_2D_c,phi.∇φ)
  jls = JuliaFunction2DLevelSet{Int32(2)}(cpp_f,cpp_g)
  coords, weights = fill_quad_data_in_unit_cube(jls,xmin,xmax,phase,degree,cell_id)
  coords, weights = to_physical_domain!(coords,weights,phi,xmin,xmax,phase,cell_id)
end

function fill_quad_data(phi,xmin,xmax,phase,degree,cell_id,::Val{3})
  ls_value_wrap_3D_c = @safe_cfunction(julia_function_wrap,Float64,(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  ls_grad_wrap_3D_c = @safe_cfunction(julia_function_wrap,ConstCxxRef{AlgoimUvector{Float64,3}},(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  cpp_f = ClosureLevelSet{Int32(3)}(ls_value_wrap_3D_c,phi.φ)
  cpp_g = ClosureLevelSet{Int32(3)}(ls_grad_wrap_3D_c,phi.∇φ)
  jls = JuliaFunction3DLevelSet{Int32(3)}(cpp_f,cpp_g)
  coords, weights = fill_quad_data_in_unit_cube(jls,xmin,xmax,phase,degree,cell_id)
  coords, weights = to_physical_domain!(coords,weights,phi,xmin,xmax,phase,cell_id)
end

function fill_quad_data(julials::LevelSetFunction,phi::AlgoimCallLevelSetFunction,
                        xmin::V,xmax::V,phase::Int,degree::Int,cell_id::Int=1) where {V}
  coords, weights = fill_quad_data_in_unit_cube(julials,xmin,xmax,phase,degree,cell_id)
  coords, weights = to_physical_domain!(coords,weights,phi,xmin,xmax,phase,cell_id)
end

function fill_quad_data_in_unit_cube(phi,xmin::V,xmax::V,phase,degree,cell_id::Int=1) where {V}
  T = eltype(xmin)
  coords = T[]; weights = T[]
  fill_quad_data_cpp(phi,coords,weights,to_array(xmin),to_array(xmax),degree,phase,Float32(cell_id))
  nd = length(xmin); np = length(weights)
  coords = reshape(coords,(nd,np))
  coords = V[ coords[:,i] for i in 1:np ]
  coords, weights
end

function to_physical_domain!(coords::Vector{V},weights::Vector{T},phi::LevelSetFunction,
                             xmin::V,xmax::V,phase::Int,cell_id::Int=1) where {T,V}
  range = xmax - xmin
  coords = map( ci -> xmin .+ ci .* range, coords )
  detJ = prod(range)
  if phase == IN || phase == OUT
    weights = detJ * weights
  elseif phase == CUT 
    if !isempty(weights)
      n = map( ci -> normal(phi,ci,cell_id) .* range, coords )
      detJ_Γ = detJ ./ norm.(n) # = j*sqrt(n⋅inv(c)⋅n)) assuming from-to cuboids
      weights = detJ_Γ .* weights
    end
  else
    error()
  end
  coords, weights
end

export fill_quad_data
export fill_quad_data_in_unit_cube
export to_physical_domain!

fill_cpp_data_raw(phi::AlgoimCallLevelSetFunction,partition::D,xmin::V,xmax::V,degree::Int) where {D,V} =
  fill_cpp_data_raw(phi,partition,xmin,xmax,degree,Val{length(xmin)}())

function fill_cpp_data_raw(phi,partition,xmin,xmax,degree,::Val{2})
  ls_value_wrap_2D_c = @safe_cfunction(julia_function_wrap,Float64,(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  ls_grad_wrap_2D_c = @safe_cfunction(julia_function_wrap,ConstCxxRef{AlgoimUvector{Float64,2}},(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  cpp_f = ClosureLevelSet{Int32(2)}(ls_value_wrap_2D_c,phi.φ)
  cpp_g = ClosureLevelSet{Int32(2)}(ls_grad_wrap_2D_c,phi.∇φ)
  jls = JuliaFunction2DLevelSet{Int32(2)}(cpp_f,cpp_g)
  coords = eltype(xmin)[]
  _fill_cpp_data_degree_dispatch(jls,to_array(partition),to_array(xmin),to_array(xmax),to_array(coords),degree)
  coords
end

function fill_cpp_data_raw(phi,partition,xmin,xmax,degree,::Val{3})
  ls_value_wrap_3D_c = @safe_cfunction(julia_function_wrap,Float64,(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  ls_grad_wrap_3D_c = @safe_cfunction(julia_function_wrap,ConstCxxRef{AlgoimUvector{Float64,3}},(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  cpp_f = ClosureLevelSet{Int32(3)}(ls_value_wrap_3D_c,phi.φ)
  cpp_g = ClosureLevelSet{Int32(3)}(ls_grad_wrap_3D_c,phi.∇φ)
  jls = JuliaFunction3DLevelSet{Int32(3)}(cpp_f,cpp_g)
  coords = eltype(xmin)[]
  _fill_cpp_data_degree_dispatch(jls,to_array(partition),to_array(xmin),to_array(xmax),to_array(coords),degree)
  coords
end

fill_cpp_data(phi::AlgoimCallLevelSetFunction,partition::D,xmin::V,xmax::V,degree::Int=2,trim::Bool=false,limitstol::Float64=1.0e-8) where {D,V} =
  fill_cpp_data(phi,partition,xmin,xmax,degree,trim,limitstol,Val{length(xmin)}())

function trim_to_limits!(coords::Matrix{T},xmin,xmax,limitstol) where {T<:Number}
  map(eachcol(coords)) do cd
    for i in eachindex(cd)
      if ( cd[i] < xmin[i] ) && ( cd[i] > xmin[i] - limitstol )
        cd[i] = xmin[i]
      elseif ( cd[i] > xmax[i] ) && ( cd[i] < xmax[i] + limitstol )
        cd[i] = xmax[i]
      end
    end
  end
end

function fill_cpp_data(phi,partition,xmin,xmax,degree,trim,limitstol,::Val{2})
  coords = fill_cpp_data_raw(phi,partition,xmin,xmax,degree,Val{2}())
  np = (partition[1]+1)*(partition[2]+1)
  coords = reshape(coords,(2,np))
  trim && trim_to_limits!(coords,xmin,xmax,limitstol)
  typeof(xmin)[eachcol(coords)...]
end

function fill_cpp_data(phi,partition,xmin,xmax,degree,trim,limitstol,::Val{3})
  coords = fill_cpp_data_raw(phi,partition,xmin,xmax,degree,Val{3}())
  np = (partition[1]+1)*(partition[2]+1)*(partition[3]+1)
  coords = reshape(coords,(3,np))
  trim && trim_to_limits!(coords,xmin,xmax,limitstol)
  typeof(xmin)[eachcol(coords)...]
end

function _fill_cpp_data_degree_dispatch(phi,partition,xmin,xmax,coords,degree)
  if degree == 2
    fill_cpp_data_taylor_2(phi,partition,xmin,xmax,coords)
  elseif degree == 3
    fill_cpp_data_taylor_3(phi,partition,xmin,xmax,coords)
  elseif degree == 4
    fill_cpp_data_taylor_4(phi,partition,xmin,xmax,coords)
  elseif degree == 5
    fill_cpp_data_taylor_5(phi,partition,xmin,xmax,coords)
  elseif degree == -1
    fill_cpp_data_cubic(phi,partition,xmin,xmax,coords)
  else
    error("Not implemented")
  end
end

export fill_cpp_data
export fill_cpp_data_raw

end

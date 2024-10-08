module Algoim

# Load algoimWrapper module

using CxxWrap

using algoimWrapper_jll
if !isdefined(algoimWrapper_jll, :libalgoimwrapper_path)
  error("algoimWrapper_jll not available on this platform")
end

@wrapmodule(() -> libalgoimwrapper)
__init__() = @initcxx

export LevelSetFunction

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

# Wrapper to implement the strategy presented in the blog post
# https://julialang.org/blog/2013/05/callback/#passing_closures_via_pass-through_pointers
# in order to avoid @cfunction closures which do not work on ARM architectures
# https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/#Closure-cfunctions

function julia_function_wrap(x,id::Float32,params::Ptr{Cvoid})
  f = unsafe_pointer_to_objref(params)::Function
  _x = to_const_array(x)
  f(_x,id)
end

mutable struct CachedLevelSetValue{A,B} <: Function
  f::A
  c::B
end

(φ::CachedLevelSetValue{<:Function})(p,i::Float32) = φ.f(p)

mutable struct CachedLevelSetGradient{A,B} <: Function
  f::A
  c::B
end

(φ::CachedLevelSetGradient{<:Function})(p,i::Float32) = φ.f(p)

export julia_function_wrap
export CachedLevelSetValue
export CachedLevelSetGradient

# AlgoimCallLevelSetFunction

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

function JuliaFunctionLevelSet(phi::AlgoimCallLevelSetFunction,::Val{2})
  ls_v_wrap_c = @safe_cfunction(julia_function_wrap,
    Float64,(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  ls_g_wrap_c = @safe_cfunction(julia_function_wrap,
    ConstCxxRef{AlgoimUvector{Float64,2}},(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  CxxWrap.gcprotect(ls_v_wrap_c); CxxWrap.gcprotect(ls_g_wrap_c);
  _cpp_f = CachedLevelSetValue(phi.φ,phi.cache_φ); CxxWrap.gcprotect(_cpp_f)
  _cpp_g = CachedLevelSetGradient(phi.∇φ,phi.cache_∇φ); CxxWrap.gcprotect(_cpp_g)  
  cpp_f = ClosureLevelSet{Int32(2)}(ls_v_wrap_c,_cpp_f)
  cpp_g = ClosureLevelSet{Int32(2)}(ls_g_wrap_c,_cpp_g)
  JuliaFunction2DLevelSet{Int32(2)}(cpp_f,cpp_g)
end

function JuliaFunctionLevelSet(phi::AlgoimCallLevelSetFunction,::Val{3})
  ls_v_wrap_c = @safe_cfunction(julia_function_wrap,
    Float64,(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  ls_g_wrap_c = @safe_cfunction(julia_function_wrap,
    ConstCxxRef{AlgoimUvector{Float64,3}},(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  CxxWrap.gcprotect(ls_v_wrap_c); CxxWrap.gcprotect(ls_g_wrap_c);
  _cpp_f = CachedLevelSetValue(phi.φ,phi.cache_φ); CxxWrap.gcprotect(_cpp_f)
  _cpp_g = CachedLevelSetGradient(phi.∇φ,phi.cache_∇φ); CxxWrap.gcprotect(_cpp_g)  
  cpp_f = ClosureLevelSet{Int32(3)}(ls_v_wrap_c,_cpp_f)
  cpp_g = ClosureLevelSet{Int32(3)}(ls_g_wrap_c,_cpp_g)
  JuliaFunction3DLevelSet{Int32(3)}(cpp_f,cpp_g)
end

function fill_quad_data(phi::AlgoimCallLevelSetFunction,
                        xmin::V,xmax::V,phase::Int,degree::Int,cell_id::Int=1) where {V}
  jls = JuliaFunctionLevelSet(phi,Val{length(xmin)}())
  fill_quad_data(jls,phi,xmin,xmax,phase,degree,cell_id)
end

function fill_quad_data(jls::LevelSetFunction,phi::AlgoimCallLevelSetFunction,
                        xmin::V,xmax::V,phase::Int,degree::Int,cell_id::Int=1) where {V}
  coords, weights = fill_quad_data_in_unit_cube(jls,xmin,xmax,phase,degree,cell_id)
  coords, weights = to_physical_domain!(coords,weights,phi,xmin,xmax,phase,cell_id)
end

function fill_quad_data_in_unit_cube(phi,xmin::V,xmax::V,phase,degree,cell_id::Int=1) where {V}
  T = eltype(xmin)
  coords = T[]; weights = T[]
  fill_quad_data_cppcall(phi,coords,weights,
                         to_array(xmin),to_array(xmax),
                         degree,phase,Float32(cell_id))
  nd = length(xmin); np = length(weights)
  coords = reshape(coords,(nd,np))
  coords = V[ coords[:,i] for i in 1:np ]
  coords, weights
end

function fill_quad_data(phi1::AlgoimCallLevelSetFunction,phi2::AlgoimCallLevelSetFunction,
                        xmin::V,xmax::V,phase1::Int,phase2::Int,degree::Int,cell_id::Int=1) where {V}
  jls1 = JuliaFunctionLevelSet(phi1,Val{length(xmin)}())
  jls2 = JuliaFunctionLevelSet(phi2,Val{length(xmin)}())
  fill_quad_data(jls1,jls2,phi1,phi2,xmin,xmax,phase1,phase2,degree,cell_id)
end

function fill_quad_data(jls1::LevelSetFunction,jls2::LevelSetFunction,
                        phi1::AlgoimCallLevelSetFunction,phi2::AlgoimCallLevelSetFunction,
                        xmin::V,xmax::V,phase1::Int,phase2::Int,degree::Int,cell_id::Int=1) where {V}
  coords, weights = fill_quad_data_in_unit_cube(jls1,jls2,xmin,xmax,phase1,phase2,degree,cell_id)
  if phase1 == CUT
    coords, weights = to_physical_domain!(coords,weights,phi1,xmin,xmax,phase1,cell_id)
  elseif phase2 == CUT
    coords, weights = to_physical_domain!(coords,weights,phi2,xmin,xmax,phase2,cell_id)
  else # phase1,phase2 ∈ {IN,OUT}; we can choose any of the two
    coords, weights = to_physical_domain!(coords,weights,phi1,xmin,xmax,phase1,cell_id)
  end
end

function fill_quad_data_in_unit_cube(phi1,phi2,xmin::V,xmax::V,phase1,phase2,degree,cell_id::Int=1) where {V}
  T = eltype(xmin)
  coords = T[]; weights = T[]
  fill_quad_data_cppcall(phi1,phi2,coords,weights,
                         to_array(xmin),to_array(xmax),
                         degree,phase1,phase2,Float32(cell_id))
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

export JuliaFunctionLevelSet
export fill_quad_data
export fill_quad_data_in_unit_cube
export to_physical_domain!

fill_cpp_data(phi::AlgoimCallLevelSetFunction,partition::D,xmin::V,xmax::V,
              degree::Int=2,trim::Bool=false,limitstol::Float64=1.0e-8;order::Int=1,
              rmin=zeros(eltype(partition),length(partition)),rmax=partition) where {D,V} =
  fill_cpp_data(phi,partition,xmin,xmax,rmin,rmax,order,degree,trim,limitstol,Val(length(xmin)))

fill_cpp_data(phivals::AbstractVector,partition::D,xmin::V,xmax::V,
              degree::Int=2,trim::Bool=false,limitstol::Float64=1.0e-8;
              rmin=zeros(eltype(partition),length(partition)),rmax=partition) where {D,V} =
  fill_cpp_data(phivals,partition,xmin,xmax,rmin,rmax,degree,trim,limitstol,Val(length(xmin)))

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

function fill_cpp_data(phi,partition,xmin,xmax,rmin,rmax,order,degree,trim,limitstol,::Val{2})
  ls_v_wrap_c = @safe_cfunction(julia_function_wrap,
    Float64,(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  ls_g_wrap_c = @safe_cfunction(julia_function_wrap,
    ConstCxxRef{AlgoimUvector{Float64,2}},(ConstCxxRef{AlgoimUvector{Float64,2}},Float32,Ptr{Cvoid}))
  _cpp_f = CachedLevelSetValue(phi.φ,phi.cache_φ)
  _cpp_g = CachedLevelSetGradient(phi.∇φ,phi.cache_∇φ) 
  cpp_f = ClosureLevelSet{Int32(2)}(ls_v_wrap_c,_cpp_f)
  cpp_g = ClosureLevelSet{Int32(2)}(ls_g_wrap_c,_cpp_g)
  jls = JuliaFunction2DLevelSet{Int32(2)}(cpp_f,cpp_g) # JuliaFunction2DLevelSet
  coords = eltype(xmin)[]
  fill_cpp_data_cppcall(jls,Val(Cint(degree)),to_array(partition),
                        to_array(xmin),to_array(xmax),
                        to_array(rmin),to_array(rmax),
                        to_array(coords),Int32(order))
  np = num_cpps(rmin,rmax,Val(2))
  coords = reshape(coords,(2,np))
  trim && trim_to_limits!(coords,xmin,xmax,limitstol)
  typeof(xmin)[eachcol(coords)...]
end

function fill_cpp_data(phi,partition,xmin,xmax,rmin,rmax,order,degree,trim,limitstol,::Val{3})
  ls_v_wrap_c = @safe_cfunction(julia_function_wrap,
    Float64,(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  ls_g_wrap_c = @safe_cfunction(julia_function_wrap,
    ConstCxxRef{AlgoimUvector{Float64,3}},(ConstCxxRef{AlgoimUvector{Float64,3}},Float32,Ptr{Cvoid}))
  _cpp_f = CachedLevelSetValue(phi.φ,phi.cache_φ)
  _cpp_g = CachedLevelSetGradient(phi.∇φ,phi.cache_∇φ)
  cpp_f = ClosureLevelSet{Int32(3)}(ls_v_wrap_c,_cpp_f)
  cpp_g = ClosureLevelSet{Int32(3)}(ls_g_wrap_c,_cpp_g)
  jls = JuliaFunction3DLevelSet{Int32(3)}(cpp_f,cpp_g) # JuliaFunction3DLevelSet
  coords = eltype(xmin)[]
  fill_cpp_data_cppcall(jls,Val(Cint(degree)),to_array(partition),
                        to_array(xmin),to_array(xmax),
                        to_array(rmin),to_array(rmax),
                        to_array(coords),Int32(order))
  np = num_cpps(rmin,rmax,Val(3))
  coords = reshape(coords,(3,np))
  trim && trim_to_limits!(coords,xmin,xmax,limitstol)
  typeof(xmin)[eachcol(coords)...]
end

function fill_cpp_data(phivals,partition,xmin,xmax,rmin,rmax,degree,trim,limitstol,::Val{D}) where {D}
  coords = eltype(xmin)[]
  fill_cpp_data_cppcall(Val(Cint(D)),Val(Cint(degree)),to_array(partition),
                        to_array(xmin),to_array(xmax),
                        to_array(rmin),to_array(rmax),
                        to_array(phivals),to_array(coords))
  np = num_cpps(rmin,rmax,Val(D))
  coords = reshape(coords,(D,np))
  trim && trim_to_limits!(coords,xmin,xmax,limitstol)
  typeof(xmin)[eachcol(coords)...]
end

num_cpps(rmin,rmax,::Val{2}) = (rmax[1]-rmin[1]+1)*(rmax[2]-rmin[2]+1)
num_cpps(rmin,rmax,::Val{3}) = (rmax[1]-rmin[1]+1)*(rmax[2]-rmin[2]+1)*(rmax[3]-rmin[3]+1)

export fill_cpp_data

end

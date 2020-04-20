import Polynomials.hasneg
import Polynomials.isneg
abstract type AbstractLagrangePolynomial{T<:Number}<:AbstractPolynomial{T} end
coeffs(p::AbstractLagrangePolynomial) = p.x # temp fix, make this work with methods that require coefficients
Base.convert(P::Type{<:AbstractLagrangePolynomial}, p::AbstractLagrangePolynomial) where {T} = P(p.x, p.y, domain(p), p.var)

function showterm(io::IO, ::Type{<:AbstractLagrangePolynomial{T}}, pj::T, var, j, first::Bool, mimetype) where {N, T}
    !first &&  print(io, " ")
    print(io, "$(pj)")
    return true
end

mutable struct LagrangePoly{T<:Number} <: AbstractLagrangePolynomial{T}
    x::Vector{T}
    y::Vector{T}
    weights::Vector
    domain::Interval
    var::Symbol
    function LagrangePoly{T}(x::Vector{T}, y::Vector{T}, domain::Interval, var::Symbol) where {T <: Number}
        return new{T}(x,y,lagrange_bary_weights(x),domain, var)
    end
end
@register LagrangePoly
LagrangePoly(x::AbstractVector{T},y::AbstractVector{T}, domain::Interval, var::SymbolLike = :x) where {T} = LagrangePoly{T}(x,y,domain,Symbol(var))

function lagrange_bary_weights(x)
    numPoints = length(x)
    weights = [prod(1/(x[i] - x[j]) for j in 1:numPoints if i ≠ j) for i in 1:numPoints]
    return weights
end

function lagrange_eval_weights(p::AbstractLagrangePolynomial, xeval, y)
    x = p.x
    w = p.w
    num = mapreduce((x,y,w)-> w*y / (xeval - x),+,x,y,w)
    denom = mapreduce((x,w)-> w / (xeval - x),+,x,w)
    val = num/denom
    return isnan(val) ? y[x.==xeval][1] : val
end

domain(p::AbstractLagrangePolynomial) = p.domain
LagrangePoly(x::Vector{T}, y::Vector{T}, var::Symbol=:x; lower = min(x...), upper = max(x...)) where {T} = LagrangePoly(x,y,Interval(lower, upper),var)
function fit(P::Type{<:AbstractLagrangePolynomial}, x::AbstractVector{T}, y::AbstractVector{T}, var = :x; lower = min(x...), upper = max(x...)) where {T}
    LagrangePoly(x,y,Interval(lower, upper),var)
end

update!(p::AbstractLagrangePolynomial{T}, y::Vector{T}) where {T} = p.y = y

# Legendre-Gauss-Radau

mutable struct LGRPoly{T<:Number} <: AbstractLagrangePolynomial{T}
    x::Vector{T}
    y::Vector{T}
    weights::Vector
    var::Symbol
    function LGRPoly{T}(y::Vector{T}, var::Symbol) where {T <: Number}
        lgrpoints = lgr_points(length(y))
        return new{T}(lgrpoints,y,lagrange_bary_weights(lgrpoints), var)
    end
end
@register LGRPoly
LGRPoly(y::AbstractVector{T}, var::SymbolLike = :x) where {T} = LGRPoly{T}(y, Symbol(var))

function lgr_points(order) # does not include endpoints
    if order>1
        N = 1:order-2
        a = zeros(order-1)
        b = zeros(order-2)
        a[1] = 1 / 3 
        a[2:end] = map((N)-> 1 / (4*N^2 + 8N +3),N)
        b = map((N)-> (N^2 + N)^0.5 / (2N+1),N)
        J = SymTridiagonal(a,b)
        return pushfirst!(eigvals(J),-1)
    end
end

domain(P::Type{<:LGRPoly}) = Interval(-1,1,true, false)
domain(p::LGRPoly) = Interval(-1,1,true, false)

function fit(P::Type{<:AbstractLagrangePolynomial}, y::AbstractVector{T}, var = :x) where {T}
    LGRPoly(y,var)
end
Base.convert(P::Type{<:LGRPoly}, p::LGRPoly) where {T} = P(p.y, p.var)

import Polynomials.hasneg
import Polynomials.isneg

export AbstractLagrangePolynomial
export LagrangePoly
export LGRPoly

abstract type AbstractLagrangePolynomial{T<:Number}<:AbstractPolynomial{T} end
coeffs(p::AbstractLagrangePolynomial) = p.y # temp fix, make this work with methods that require coefficients
Base.convert(P::Type{<:AbstractLagrangePolynomial}, p::AbstractLagrangePolynomial) where {T} = P(p.x, p.y, domain(p), p.var)

function (lp::AbstractLagrangePolynomial{T})(x::S) where {T,S<:Number}
    x ∉ domain(lp) && error("$x outside of domain")
    return lagrange_eval_weights(lp,x,lp.y)
end

function showterm(io::IO, ::Type{<:AbstractLagrangePolynomial{T}}, pj::T, var, j, first::Bool, mimetype) where {N, T}
    !first &&  print(io, " ")
    print(io, "$(pj)")
    return true
end

mutable struct LagrangePoly{T<:Number} <: AbstractLagrangePolynomial{T}
    x::Vector{T}
    y::Vector{T}
    weights::Vector{T}
    domain::Interval
    var::Symbol
    function LagrangePoly{T}(x::Vector{T}, y::Vector{T}, domain::Interval, var::Symbol) where {T <: Number}
        return new{T}(x,y,lagrange_bary_weights(x,T),domain, var)
    end
end
@register LagrangePoly
LagrangePoly(x::AbstractVector{T},y::AbstractVector{T}, domain::Interval, var::SymbolLike = :x) where {T} = LagrangePoly{T}(x,y,domain,Symbol(var))

function lagrange_bary_weights(x::AbstractVector{T}, ::Type{T}) where {T}
    numPoints = length(x)
    w = ones(T,numPoints)
    for j in 2:numPoints
        @inbounds xj = x[j]
        for k in 1:j-1
            @inbounds w[k] *= x[k] - xj
        end
        for k in 1:j-1
            @inbounds w[j] *= xj - x[k]
        end
    end
    return 1 ./ w
end

function lagrange_eval_weights(p::AbstractLagrangePolynomial{T}, xeval, y) where {T}
    x = p.x
    w = p.weights
    num = zero(T)
    denom = zero(T)
    for j = 1:length(x)
        @inbounds wx = w[j] / (xeval - x[j])
        isinf(wx) && return y[j]
        @inbounds num += wx*y[j]
        denom +=wx
    end
    return num/denom
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
    weights::Vector{T}
    var::Symbol
    function LGRPoly{T}(y::Vector{T}, var::Symbol) where {T <: Number}
        lgrpoints = lgr_points(length(y))
        return new{T}(lgrpoints,y,lagrange_bary_weights(lgrpoints, T), var)
    end
end
@register LGRPoly
LGRPoly(y::AbstractVector{T}, var::SymbolLike = :x) where {T} = LGRPoly{T}(y, Symbol(var))

function lgr_points(order) # does not include endpoints
    order<2 && error("order must be greater than one")
    N = 1:order-2
    a = zeros(order-1)
    b = zeros(order-2)
    a[1] = 1 / 3 
    a[2:end] = map((N)-> 1 / (4*N^2 + 8N +3),N)
    b = map((N)-> (N^2 + N)^0.5 / (2N+1),N)
    J = SymTridiagonal(a,b)
    return pushfirst!(eigvals(J),-1)
end

domain(P::Type{<:LGRPoly}) = Interval(-1,1,true, false)
domain(p::LGRPoly) = Interval(-1,1,true, false)

function fit(P::Type{<:AbstractLagrangePolynomial}, y::AbstractVector{T}, var = :x) where {T}
    LGRPoly(y,var)
end
Base.convert(P::Type{<:LGRPoly}, p::LGRPoly) where {T} = P(p.y, p.var)

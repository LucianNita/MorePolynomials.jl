import Polynomials.hasneg
import Polynomials.isneg

export AbstractLagrangePolynomial
export LagrangePoly
export LGRPoly
export derivmatrix
export lgr_points
export lgr_weights

abstract type AbstractLagrangePolynomial{T<:Number}<:AbstractPolynomial{T} end

coeffs(p::AbstractLagrangePolynomial) = p.y # temp fix, make this work with methods that require coefficients

convert(P::Type{<:AbstractLagrangePolynomial}, p::AbstractLagrangePolynomial) where {T} = P(p.x, p.y, domain(p), p.var)

Base.eachindex(p::AbstractLagrangePolynomial) = 1:length(p)

Base.getindex(p::AbstractLagrangePolynomial, i::Int) = (p.x[i],p.y[i])
Base.setindex!(p::AbstractLagrangePolynomial, y::Number, i::Int) = p.y[i] = y # add case for updating x and y
Base.firstindex(p::AbstractLagrangePolynomial) = 1
Base.lastindex(p::AbstractLagrangePolynomial) = length(p.x)

Base.length(p::AbstractLagrangePolynomial) = length(p.x)

function (lp::AbstractLagrangePolynomial{T})(x::S) where {T,S<:Number}
    if x ∉ domain(lp)
        pDomain = domain(lp)
        throw(DomainError(x,"outside of domain $pDomain"))
    end
    return lagrange_eval_weights(lp,x,lp.y)
end

function showterm(io::IO, ::Type{<:AbstractLagrangePolynomial{T}}, pj::Any, var, j, first::Bool, mimetype) where {N, T}
    !first &&  print(io, " ")
    print(io, "$(pj)")
    return true
end

mutable struct LagrangePoly{T<:Number} <: AbstractLagrangePolynomial{T} # consider making this immutable for performance
    x::Vector{T}
    y::Vector{T}
    weights::Vector{T}
    domain::Interval{T}
    var::Symbol
    function LagrangePoly{T}(x::Vector{T}, y::Vector{T}, domain::Interval{T}, var::Symbol) where {T <: Number}
        return new{T}(x,y,lagrange_bary_weights(x),domain, var)
    end
end

@register LagrangePoly

function LagrangePoly(x::Vector{T}, y::Vector{T}, lower, upper, var::Symbol=:x; containlower::Bool=true, containupper::Bool=true) where {T}
    lower = convert(T, lower)
    upper = convert(T, upper)
    return LagrangePoly{T}(x,y,Interval(lower, upper, containlower, containupper),var)
end

function LagrangePoly(x::Vector{T}, y::Vector{T}, var::Symbol=:x; kwargs...) where {T}
    lower = min(x...)
    upper = max(x...)
    return LagrangePoly(x,y,lower,upper,var;kwargs...)
end

function lagrange_bary_weights(x::AbstractVector{T}) where {T}
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


function fit(P::Type{<:AbstractLagrangePolynomial}, x::AbstractVector{T}, y::AbstractVector{T}, var = :x) where {T}
    return LagrangePoly(x,y,var)
end


function derivmatrix(p::AbstractLagrangePolynomial{T}) where {T}
    w = p.weights
    x = p.x
    numPoints = length(x)
    D = zeros(T, numPoints, numPoints)
    for i in 1:numPoints
        diag = zero(T)
        xi = x[i]
        wi = w[i]
        for j = 1:i-1
            @inbounds Dij = w[j] / ( wi * (xi - x[j]) )
            diag -= Dij
            @inbounds D[i,j] = Dij
        end
        for j = i+1:numPoints
            @inbounds Dij = w[j] / ( wi * (xi - x[j]) )
            diag -= Dij
            @inbounds D[i,j] = Dij
        end
        @inbounds D[i,i] = diag
    end
    return D
end

function derivative(p::AbstractLagrangePolynomial{T}, order::Integer = 1) where {T}
    pder = deepcopy(p)
    pder.y = derivmatrix(p) * p.y
    return order == 1 ? pder : derivative(pder, order - 1) # make this better by numerically calculating higher derivatives (you know you can)
end


# Legendre-Gauss-Radau

# todo precompile first few (hudred)? weights (probs doesn't matter)

mutable struct LGRPoly{T<:Number} <: AbstractLagrangePolynomial{T}
    x::Vector{T}
    y::Vector{T}
    weights::Vector{T}
    var::Symbol
    function LGRPoly{T}(y::Vector{T}, var::Symbol) where {T <: Number}
        lgrpoints = lgr_points(length(y))
        return new{T}(lgrpoints,y,lagrange_bary_weights(lgrpoints), var)
    end
end
@register LGRPoly
LGRPoly(y::AbstractVector{T}, var::SymbolLike = :x) where {T} = LGRPoly{T}(y, Symbol(var)) # add function where flipped = true

function lgr_points(numPoints::Int) # find when this becomes slower than fast gaussradau and use that instead
    numPoints<2 && error("number of points must be greater than one")
    N = 1:numPoints-2
    a = zeros(numPoints-1)
    b = zeros(numPoints-2)
    a[1] = 1 / 3 
    a[2:end] = map((N)-> 1 / (4*N^2 + 8N +3),N)
    b = map((N)-> (N^2 + N)^0.5 / (2N+1),N)
    J = SymTridiagonal(a,b)
    return pushfirst!(eigvals(J),-1)
end

domain(P::Type{<:LGRPoly}) = Interval(-1,1)
domain(p::LGRPoly{T}) where {T} = Interval{T}(-1,1)

function fit(P::Type{<:AbstractLagrangePolynomial}, y::AbstractVector{T}, var::SymbolLike = :x) where {T}
    LGRPoly(y,var)
end

Base.convert(P::Type{<:LGRPoly}, p::LGRPoly) where {T} = P(p.y, p.var)

function lgr_weights(p::LGRPoly{T}) where {T} # definitely precompute these
    xlen = length(p)
    return gaussradau(xlen)[2]
end

function lgr_weights(numPoints::Int) # definitely precompute these
    return gaussradau(numPoints)[2]
end

function integrate(p::LGRPoly) 
    return sum(lgr_weights(p).*p.y)
end


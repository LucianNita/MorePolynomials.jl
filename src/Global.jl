import Base.push!
import Base.getindex
import Base.setindex
import Base.firstindex
import Base.lastindex

export domains
export GlobalPoly


mutable struct GlobalPoly{T<:Number} <: AbstractPolynomial{T}
    polys::Vector{AbstractPolynomial{T}}#array of polynomials
    domains::Vector{Interval{T}}
    var::Symbol
    function GlobalPoly{T}(p::AbstractPolynomial{T}, gloDomain::Interval{T}, var::Symbol) where {T<:Number}
        checkends(gloDomain,p)
        checkinfs(gloDomain, domain(p))
        return new{T}([p], [gloDomain], var) # new initializing globalPoly must be created with an initial poly to avoid checks further down the line
    end
end

function GlobalPoly{T}(p::N, lower, upper, var::SymbolLike=:x;  containlower::Bool=containslower(domain(p)), containupper::Bool=containsupper(domain(p))) where {T,N<:AbstractPolynomial}
    lower = convert(T, lower)
    upper = convert(T, upper)
    p = convert(Base.typename(N).wrapper{T},p)
    gloD = Interval(lower, upper, containlower, containupper)
    return GlobalPoly{T}(p, gloD, var)
end
function GlobalPoly{T}(p::AbstractPolynomial, var::SymbolLike=:x;  kwargs...) where {T}
    lower=first(domain(p))
    upper=last(domain(p))
    if T<:Int && (lower==Inf || upper==Inf)
        error("GlobalyPoly of type $T does not support infinite endpoints, consider setting upper and lower bounds or converting to type <:AbstractFloat")
    end
    GlobalPoly{T}(p, lower, upper, var; kwargs...)
end
function GlobalPoly(p::AbstractPolynomial{T}, lower, upper, var::SymbolLike=:x; kwargs...) where {T}
    return GlobalPoly{T}(p, lower, upper, var; kwargs...)
end
function GlobalPoly(p::AbstractPolynomial{T}, var::SymbolLike=:x; kwargs...) where {T}
    return GlobalPoly{T}(p, var; kwargs...)
end

function showterm(io::IO, ::Type{<:GlobalPoly{T}}, pj::AbstractPolynomial{T}, var, j, first::Bool, mimetype) where {T}
    !first &&  print(io, " ")
    print(io, "$(pj)")
    return true
end

length(gloP::GlobalPoly) = length(gloP.polys)

eachindex(p::GlobalPoly) = 1:length(p)
# indexing global polys
getindex(gloP::GlobalPoly, i::Int) = gloP.polys[i]
#setindex!(gloP::GlobalPoly, p::AbstractPolynomial, i) = #todo
firstindex(gloP::GlobalPoly, i::Int) = gloP.polys[begin]
lastindex(gloP::GlobalPoly, i::Int) = gloP.polys[end]

domains(gloP::GlobalPoly) = gloP.domains

function push!(gloP::GlobalPoly, locP::AbstractPolynomial, gloDomain::Interval) 
    gloDomains = domains(gloP)
    checkends(gloDomain,locP)
    #checkinfs(gloDomain, domain(locP))
    for i in 1:length(gloP) 
        if !isempty(intersect(gloDomain, gloDomains[i]))
            conflict = intersect(gloDomain, gloDomains[i])
            error("the domain intersects with an existing polynomial at location [$i]. Intersection is $conflict .")
        end
    end
    push!(gloP.polys, locP)
    push!(gloP.domains, gloDomain)
end
function push!(gloP::GlobalPoly{T}, locP::AbstractPolynomial, gloLower, gloUpper; containlower=containslower(domain(locP)), containupper=containsupper(domain(locP))) where {T}
    if T<:Int && (lower==Inf || upper==Inf)
        error("GlobalyPoly of type $T does not support infinite endpoints, consider setting upper and lower bounds or converting to type <:AbstractFloat")
    end
    gloD = Interval(gloLower, gloUpper, containlower, containupper)
    checkinfs(gloD, domain(locP))
    return push!(gloP, locP, gloD) 
end
function push!(gloP::GlobalPoly, locP::AbstractPolynomial; kwargs...)
    lower=first(domain(locP))
    upper=last(domain(locP))
    return push!(gloP, locP, lower, upper; kwargs...)
end 

function domain(p::GlobalPoly)
    low,up = first(domains(p)[1]), last(domains(p)[end])
    return Interval(first(low), last(up), containslower(low), containsupper(up))
end

function (p::GlobalPoly{T})(x::S) where {T,S<:Number}
    gloDomains = domains(p)
    for i in 1:length(p)
        if x in gloDomains[i]
            return p[i](global2local(domain(p[i]), gloDomains[i], x))
        end
    end
    if x ∉ domain(p)
        pDomain = domain(p)
        throw(DomainError(x,"outside of domain $pDomain"))
    end
    closestMatch = argmin([min(abs.(x .- [first(domains(p)[i]), last(domains(p)[i])])...) for i in 1:length(p)])
    throw(DomainError(x,"not found in global polynomial domains. Closest polynomial has domain $(domains(p)[closestMatch]) at index [$closestMatch]."))
end

function checkinfs(gloD::Interval, locD::Interval)
    lLower = first(locD)
    lUpper = last(locD)
    gLower = first(gloD)
    gUpper = last(gloD)
    if isinf(gLower) || isinf(gUpper) || isinf(lLower) || isinf(lUpper)
        if isinf(gLower) && isinf(gUpper) 
            if xor(isinf(lLower),isinf(lUpper))
                throw(DomainError(gloD,"this global domain is not compatible with local domain $locD"))
            else
                return true
            end
        elseif isinf(lLower)&&isinf(lUpper)
            return true
        elseif !issubset(gloD,locD)
            if !isinf(gLower)&&!isinf(gUpper)
                throw(DomainError(gloD,"this global domain is not a subset of local domain $locD"))
            end
            if issubset(Interval(-Inf,Inf, false, false),Interval(min(lLower,gLower), max(lUpper,gUpper))) 
                throw(DomainError(gloD,"this global domain is not compatible with local domain $locD"))
            end
        end
    end
    return true
end

function checkends(gloD::Interval,locP::AbstractPolynomial)
        !containslower(domain(locP)) && containslower(gloD) && throw(DomainError(gloD,"global lower endpoints incompatible with local polynomial domain $(domain(locP)). Consider either opening the global endpoint or closing the local endpoint."))#check first that local and global bound endpoints are compatible
        !containsupper(domain(locP)) && containsupper(gloD) && throw(DomainError(gloD,"global upper endpoints incompatible with local polynomial $(domain(locP)). Consider either opening the global endpoint or closing the local endpoint."))
end
# following Julia convert() convention where thing being converted to comes first
mapendpoints(from::Interval, to::Interval, x::Number) = first(to) + (last(to) - first(to)) * (x - first(from)) / ( last(from) - first(from) )
function global2local(locD::Interval, gloD::Interval, x::Number)
    lLower = first(locD)
    lUpper = last(locD)
    gLower = first(gloD)
    gUpper = last(gloD)
    if isinf(gLower) || isinf(gUpper) || isinf(lLower) || isinf(lUpper)
        if (!isinf(gLower)&&!isinf(gUpper)) || (isinf(lLower) && isinf(lUpper)) 
            return x
        elseif (!isinf(lLower)&&!isinf(lUpper))
            if (isinf(gLower) && isinf(gUpper))
                return mapendpoints(Interval(-1,1),locD,2*atan(x) / pi)
            elseif isinf(gUpper)
                return mapendpoints(Interval(0,1),locD,2*atan(x-gLower) / pi)
            else
                return mapendpoints(Interval(-1,0),locD,2*atan(x-gUpper) / pi)
            end
        elseif isinf(gUpper)
            return x - gLower + lLower
        else
            return x - gUpper + lUpper
        end

    end
    return lLower + (lUpper - lLower) * (x - gLower) / ( gUpper - gLower )
end

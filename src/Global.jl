import Base.push!
import Base.length
import Base.getindex
import Base.setindex
import Base.firstindex
import Base.lastindex

# some Intervals functions that make my life easier

includeslower(i::Interval) = first(inclusivity(i))
includesupper(i::Interval) = last(inclusivity(i))


mutable struct GlobalPoly{T<:Number} <: AbstractPolynomial{T}
    polys::Vector{AbstractPolynomial{T}}#array of polynomials
    domains::Vector{Interval{T}}
    GlobalPoly{T}(p::AbstractPoly, domain::Interval) where {T<:Number} = new{T}([p], [domain]) # new initializing globalPoly must be created with an initial poly to avoid checks further down the line
end

GlobalPoly(p::AbstractPoly{T}, d::Interval) where {T} = GlobalPoly{T}(p,d)

length(gloP::GlobalPoly) = length(gloP.polys)

# indexing global polys
getindex(gloP::GlobalPoly, i) = gloP.polys[i]
#setindex!(gloP::GlobalPoly, p::AbstractPolynomial, i) = #todo
firstindex(gloP::GlobalPoly) = gloP.polys[begin]
lastindex(gloP::GlobalPoly) = gloP.polys[end]

domains(gloP::GlobalPoly) = gloP.domains

function push!(gloP::GlobalPoly{T}, locP::AbstractPolynomial{T}, gloDomain::Interval) where {T} 
    !includeslower(domain(locP)) && includeslower(gloBound) && throw(DomainError(golDomain,"global lower endpoints incompatible with local polynomial. Consider either opening the global endpoint or closing the local endpoint."))#check first that local and global bound endpoints are compatible
    !includesupper(domain(locP)) && includesupper(gloBound) && throw(DomainError(gloDomain,"global upper endpoints incompatible with local polynomial. Consider either opening the global endpoint or closing the local endpoint."))
    gloDomains = domains(gloP)
    for i in 1:length(gloP)
        if !isempty(intersect(gloDomain, gloDomains[i]))
            conflict = intersect(gloDomain, gloDomains[i])
            error("the domain intersects with an existing polynomial at location [$i]. Intersection is $conflict .")
        end
    end
    push!(gloP.polys, locP)
    push!(gloP.domain, gloBound)
end
push!(gloP::GlobalPoly, locP::AbstractPoly, gloLower, gloUpper) = push!(gloP, locP, Interval(gloLower, gloUpper)) 
push!(gloP::GlobalPoly, locP::AbstractPoly) = push!(gloP, locP, domain(locP))

function domain(p::GlobalPoly)
    low,up = first(domains(p)[1]), last(domains(p)[end])
    return Interval(first(low), last(up), includeslower(low), includesupper(up))
end

function (p::GlobalPoly{T})(x::S) where {T,S<:Number}
    gloDomains = domains(p)
    for i in 1:length(p)
        if x in gloDomains[i]
            return p[i](global_to_local(domain(p[i]), gloDomains[i], x))
        end
    end
    if x ∉ domain(p)
        pDomain = domain(p)
        throw(DomainError(x,"outside of domain $pDomain"))
    end
    closestMatch = argmin([min(abs.(x .- [first(domains(p)[i]), last(domains(p)[i])])...) for i in 1:length(p)])
    throw(DomainError(x,"not found in global polynomial domains. Closest polynomial has domain $(domains(p)[closestMatch]) at index [$closestMatch]."))
end

# following Julia convert() convention where thing being converted to comes first
local_to_global(gloD::Interval, locD::Interval, x::Number) = first(gloD) + (last(gloD) - first(gloD)) * (x - first(locD)) / ( last(locD) - first(locD) )
global_to_local(locD::Interval, gloD::Interval, x::Number) = first(locD) + (last(locD) - first(locD)) * (x - first(gloD)) / ( last(gloD) - first(gloD) )

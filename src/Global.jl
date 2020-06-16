export domains
export PiecewisePoly
export update!
export lengthpoints

mutable struct PiecewisePoly{T<:Number} <: AbstractPolynomial{T}
    polys::Vector{AbstractPolynomial{T}}#array of polynomials
    domains::Vector{Interval{T}}
    var::Symbol
    function PiecewisePoly{T}(subP::AbstractPolynomial{T}, pieDomain::Interval{T}, var::Symbol) where {T<:Number}
        checkends(pieDomain,subP,true)
        checkinfs(pieDomain, domain(subP))
        return new{T}([subP], [pieDomain], var) # new initializing globalPoly must be created with an initial poly to avoid checks further down the line
    end
end

function PiecewisePoly{T}(subP::N, lower, upper; var::SymbolLike=:x,  containlower::Bool=containslower(domain(subP)), containupper::Bool=containsupper(domain(subP))) where {T,N<:AbstractPolynomial}
    lower = convert(T, lower)
    upper = convert(T, upper)
    subP = convert(Base.typename(N).wrapper{T},subP)
    pieDomain = Interval(lower, upper, containlower, containupper)
    return PiecewisePoly{T}(subP, pieDomain, var)
end
function PiecewisePoly{T}(subP::AbstractPolynomial; kwargs...) where {T}
    lower=first(domain(subP))
    upper=last(domain(subP))
    if T<:Int && (lower==Inf || upper==Inf)
        error("PiecewisePoly of type $T does not support infinite endpoints, consider setting upper and lower bounds or converting to type <:AbstractFloat")
    end
    PiecewisePoly{T}(subP, lower, upper; kwargs...)
end
function PiecewisePoly(subP::AbstractPolynomial{T}, lower, upper;kwargs...) where {T}
    return PiecewisePoly{T}(subP, lower, upper;kwargs...)
end
function PiecewisePoly(subP::AbstractPolynomial{T};kwargs...) where {T}
    return PiecewisePoly{T}(subP;kwargs...)
end

function showterm(io::IO, ::Type{<:PiecewisePoly{T}}, pj::Any, var, j, first::Bool, mimetype) where {T}
    !first &&  print(io, " ")
    print(io, "$(pj)")
    return true
end

length(pieP::PiecewisePoly) = length(pieP.polys)

function lengthpoints(pieP::PiecewisePoly)
    lengthp = 0
    pieDomains = domains(pieP)
    lengthp += length(pieP[1])
    for i in 2:length(pieP)
        if iscontinuous(pieP) && last(pieDomains[i-1]) == first(pieDomains[i])
            lengthp -= 1
        end
        @inbounds lengthp += length(pieP[i])
    end
    return lengthp
end

Base.eachindex(pieP::PiecewisePoly) = 1:length(pieP)
# indexing global polys
Base.getindex(pieP::PiecewisePoly, i::Int) = pieP.polys[i]
#Base.setindex!(pieP::PiecewisePoly, p::AbstractPolynomial, i) = #todo
Base.firstindex(pieP::PiecewisePoly) = 1
Base.lastindex(pieP::PiecewisePoly) = length(pieP.polys)

domains(pieP::PiecewisePoly) = pieP.domains

iscontinuous(pieP::PiecewisePoly) = true # update this in future

function push!(pieP::PiecewisePoly, subP::AbstractPolynomial, pieDomain::Interval)
    pieDomains = domains(pieP)
    checkends(pieDomain,subP,iscontinuous(pieP))
    subDomain = checkinfs(pieDomain, domain(subP)) # check this mapping for infinite local domains mapping to non infinite global domians, see internal comment
    if last(pieDomains[end]) <= first(pieDomain)
        if last(pieDomains[end]) == first(pieDomain) && iscontinuous(pieP) && pieP(last(pieDomains[end])) != subP(first(subDomain)) 
            error("new poly \"N\" lower endpoint y value is not coincidence with PiecewisePoly \"G\" (G=$(pieP(last(pieDomains[end]))) vs N=$(subP(first(domain(subP))))). Use a non continuous PiecewisePoly if this is what you want.")
        end
        push!(pieP.polys, subP)
        push!(pieP.domains, pieDomain)
        return pieP
    end
    if first(pieDomains[1]) >= last(pieDomain)
        if first(pieDomains[1]) == last(pieDomain) && iscontinuous(pieP) && pieP(first(pieDomains[1])) != subP(last(subDomain)) 
            error("new poly \"N\" upper endpoint y value is not coincidence with PiecewisePoly \"G\" (G=$(pieP(first(pieDomains[1]))) vs N=$(subP(last(domain(subP))))). Use a non continuous PiecewisePoly if this is what you want.")
        end
        pushfirst!(pieP.polys, subP)
        return pushfirst!(pieP.domains, pieDomain)
    end
    for i in 1:length(pieP)-1 
        if last(pieDomains[i+1]) > first(pieDomain)
            if  last(pieDomains[i]) > first(pieDomain)
                conflict = intersect(pieDomain, pieDomains[i])
                error("the domain intersects with an existing polynomial at location [$(i)]. Intersection is $conflict .")
            end
            if first(pieDomains[i+1]) < last(pieDomain) 
                conflict = intersect(pieDomain, pieDomains[i+1])
                error("the domain intersects with an existing polynomial at location [$(i+1)]. Intersection is $conflict .")
            end
            if last(pieDomains[i]) == first(pieDomain) && iscontinuous(pieP) && pieP(last(pieDomains[i])) != subP(first(subDomain)) 
                error("new poly \"N\" lower endpoint y value is not coincidence with PiecewisePoly \"G\" (G=$(pieP(last(pieDomains[i]))) vs N=$(subP(first(domain(subP))))). Use a non continuous PiecewisePoly if this is what you want.")
            end
            if first(pieDomains[i+1]) == last(pieDomain) && iscontinuous(pieP) && pieP(first(pieDomains[i+1])) != subP(last(subDomain)) 
                error("new poly \"N\" upper endpoint y value is not coincidence with PiecewisePoly \"G\" (G=$(pieP(first(pieDomains[i+1]))) vs N=$(subP(last(domain(subP))))). Use a non continuous PiecewisePoly if this is what you want.")
            end
            insert!(pieP.polys, i+1, subP)
            insert!(pieP.domains, i+1, pieDomain)
            return pieP
        end
    end

    conflict = intersect(pieDomain, pieDomains[end])
    error("new domain endpoints intersect existing polynomial. Intersection is $conflict")
end
function push!(pieP::PiecewisePoly{T}, subP::AbstractPolynomial, gloLower, gloUpper; containlower=containslower(domain(subP)), containupper=containsupper(domain(subP))) where {T}
    if T<:Int && (gloLower==Inf || gloUupper==Inf)
        error("PiecewisePoly of type $T does not support infinite endpoints, consider setting upper and lower bounds or converting to type <:AbstractFloat")
    end
    pieDomain = Interval(gloLower, gloUpper, containlower, containupper)
    return push!(pieP, subP, pieDomain) 
end
function push!(pieP::PiecewisePoly, subP::AbstractPolynomial; kwargs...)
    lower=first(domain(subP))
    upper=last(domain(subP))
    return push!(pieP, subP, lower, upper; kwargs...)
end 

function domain(pieP::PiecewisePoly)
    low,up = domains(pieP)[1], domains(pieP)[end]
    return Interval(first(low), last(up), containslower(low), containsupper(up))
end

function (pieP::PiecewisePoly{T})(x::S) where {T,S<:Number}
    pieDomain = domains(pieP)
    for i in 1:length(pieP)
        if x in pieDomain[i]
            return pieP[i](global2local(domain(pieP[i]), pieDomain[i], x))
        end
    end
    if x ∉ domain(pieP)
        pDomain = domain(pieP)
        throw(DomainError(x,"outside of domain $pDomain"))
    end
    closestMatch = argmin([min(abs.(x .- [first(domains(pieP)[i]), last(domains(pieP)[i])])...) for i in 1:length(pieP)])
    throw(DomainError(x,"not found in PiecewisePoly domains. Closest polynomial has domain $(domains(pieP)[closestMatch]) at index [$closestMatch]."))
end

function update!(pieP::PiecewisePoly{T}, y::AbstractVector{T}) where {T}
    pieDomains = domains(pieP)
    continuouspLengths = [ iscontinuous(pieP) && last(pieDomains[i-1]) == first(pieDomains[i]) ? length(pieP[i]) - 1 : length(pieP[i]) for i in 2:length(pieP)]
    pushfirst!(continuouspLengths,length(pieP[1]))
    sum(continuouspLengths) !=length(y) && throw(DimensionMismatch("tried to assign $(length(y)) elements to $(sum(continuouspLengths)) destinations"))
    Threads.@threads for i in 1:length(pieP)
        if continuouspLengths[i] == length(pieP[i]) # case where endpoints are not coincident
            pieP[i][:] = y[sum(continuouspLengths[1:i-1])+1:sum(continuouspLengths[1:i])]
        else
            pieP[i][:] = y[sum(continuouspLengths[1:i-1]):sum(continuouspLengths[1:i])]
        end
    end
    return pieP
end

function checkinfs(pieDomain::Interval, subDomain::Interval)
    subLower = first(subDomain)
    subUpper = last(subDomain)
    pieLower = first(pieDomain)
    pieUpper = last(pieDomain)
    if isinf(pieLower) || isinf(pieUpper) || isinf(subLower) || isinf(subUpper)
        if isinf(pieLower) && isinf(pieUpper) 
            if xor(isinf(subLower),isinf(subUpper))
                throw(DomainError(pieDomain,"this piecewise domain is not compatible with sub domain $subDomain"))
            else
                return subDomain
            end
        elseif isinf(subLower)&&isinf(subUpper)
            return pieDomain # return global domain so domain mapping doesn't get confused, check if this needs to be done elsewhere too
        elseif !issubset(pieDomain,subDomain)
            if !isinf(pieLower)&&!isinf(pieUpper)
                throw(DomainError(pieDomain,"this piecewise domain is not a subset of sub domain $subDomain"))
            end
            if issubset(Interval(-Inf,Inf, false, false),Interval(min(subLower,pieLower), max(subUpper,pieUpper))) 
                throw(DomainError(pieDomain,"this piecewise domain is not compatible with sub domain $subDomain"))
            end
        end
    end
    return subDomain
end


function checkends(pieDomain::Interval,subP::AbstractPolynomial, iscontinuous::Bool)
    subDomain = domain(subP)
    if iscontinuous
        ( !containslower(subDomain) || !containsupper(subDomain) ) && throw(DomainError(subP),"continuous PiecewisePoly can only contain sub polynomial that have closed endpoints")
    else
        !containslower(subDomain) && containslower(pieDomain) && throw(DomainError(pieDomain,"piecewise lower endpoints incompatible with sub polynomial domain $(subDomain). Consider either opening the piecewise endpoint or closing the sub endpoint."))#check first that local and global bound endpoints are compatible
        !containsupper(subDomain) && containsupper(pieDomain) && throw(DomainError(pieDomain,"piecewise upper endpoints incompatible with sub polynomial domain $(subDomain). Consider either opening the piecewise endpoint or closing the sub endpoint."))
    end
end

# following Julia convert() convention where thing being converted to comes first
mapendpoints(from::Interval, to::Interval, x::Number) = first(to) + (last(to) - first(to)) * (x - first(from)) / ( last(from) - first(from) )
function global2local(subDomain::Interval, pieDomain::Interval, x::Number)
    subLower = first(subDomain)
    subUpper = last(subDomain)
    pieLower = first(pieDomain)
    pieUpper = last(pieDomain)
    if isinf(pieLower) || isinf(pieUpper) || isinf(subLower) || isinf(subUpper)
        if (!isinf(pieLower)&&!isinf(pieUpper)) || (isinf(subLower) && isinf(subUpper)) 
            return x
        elseif (!isinf(subLower)&&!isinf(subUpper))
            if (isinf(pieLower) && isinf(pieUpper))
                return mapendpoints(Interval(-1,1),subDomain,2*atan(x) / pi)
            elseif isinf(pieUpper)
                return mapendpoints(Interval(-1,1),subDomain,(x-1)/(1+x)) # allow user to specify custom function here?
            else
                return mapendpoints(Interval(-1,1),subDomain,(-x-1)/(x-1))
            end
        elseif isinf(pieUpper)
            return x - pieLower + subLower
        else
            return x - pieUpper + subUpper
        end

    end
    return subLower + (subUpper - subLower) * (x - pieLower) / ( pieUpper - pieLower )
end

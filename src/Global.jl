export domains
export GlobalPoly
export update!
export lengthpoints

mutable struct GlobalPoly{T<:Number} <: AbstractPolynomial{T}
    polys::Vector{AbstractPolynomial{T}}#array of polynomials
    domains::Vector{Interval{T}}
    var::Symbol
    function GlobalPoly{T}(p::AbstractPolynomial{T}, gloDomain::Interval{T}, var::Symbol) where {T<:Number}
        checkends(gloDomain,p,true)
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

function showterm(io::IO, ::Type{<:GlobalPoly{T}}, pj::Any, var, j, first::Bool, mimetype) where {T}
    !first &&  print(io, " ")
    print(io, "$(pj)")
    return true
end

length(gloP::GlobalPoly) = length(gloP.polys)

function lengthpoints(gloP::GlobalPoly)
    lengthp = 0
    gDomains = domains(gloP)
    lengthp += length(gloP[1])
    for i in 2:length(gloP)
        if iscontinuous(gloP) && last(gDomains[i-1]) == first(gDomains[i])
            lengthp -= 1
        end
        @inbounds lengthp += length(gloP[i])
    end
    return lengthp
end

Base.eachindex(gloP::GlobalPoly) = 1:length(gloP)
# indexing global polys
Base.getindex(gloP::GlobalPoly, i::Int) = gloP.polys[i]
#Base.setindex!(gloP::GlobalPoly, p::AbstractPolynomial, i) = #todo
Base.firstindex(gloP::GlobalPoly, i::Int) = gloP.polys[begin]
Base.lastindex(gloP::GlobalPoly, i::Int) = gloP.polys[end]

domains(gloP::GlobalPoly) = gloP.domains

iscontinuous(p::GlobalPoly) = true # update this in future

function push!(gloP::GlobalPoly, locP::AbstractPolynomial, gloDomain::Interval) 
    gloDomains = domains(gloP)
    checkends(gloDomain,locP,iscontinuous(gloP))
    checkinfs(gloDomain, domain(locP))
    if last(gloDomains[end]) <= first(gloDomain)
        if last(gloDomains[end]) == first(gloDomain) && iscontinuous(gloP) && gloP(last(gloDomains[end])) != locP(first(domain(locP))) 
            error("new poly \"N\" lower endpoint y value is not coincidence with GlobalPoly \"G\" (G=$(gloP(last(gloDomains[end]))) vs N=$(locP(first(domain(locP))))). Use a non continuous GlobalPoly if this is what you want.")
        end
        push!(gloP.polys, locP)
        return push!(gloP.domains, gloDomain)
    end
    if first(gloDomains[1]) >= last(gloDomain)
        if first(gloDomains[1]) == last(gloDomain) && iscontinuous(gloP) && gloP(first(gloDomains[1])) != locP(last(domain(locP))) 
            error("new poly \"N\" upper endpoint y value is not coincidence with GlobalPoly \"G\" (G=$(gloP(first(gloDomains[1]))) vs N=$(locP(last(domain(locP))))). Use a non continuous GlobalPoly if this is what you want.")
        end
        pushfirst!(gloP.polys, locP)
        return pushfirst!(gloP.domains, gloDomain)
    end
    for i in 1:length(gloP)-1 
        if last(gloDomains[i+1]) > first(gloDomain)
            if  last(gloDomains[i]) > first(gloDomain)
                conflict = intersect(gloDomain, gloDomains[i])
                error("the domain intersects with an existing polynomial at location [$(i)]. Intersection is $conflict .")
            end
            if first(gloDomains[i+1]) < last(gloDomain) 
                conflict = intersect(gloDomain, gloDomains[i+1])
                error("the domain intersects with an existing polynomial at location [$(i+1)]. Intersection is $conflict .")
            end
            if last(gloDomains[i]) == first(gloDomain) && iscontinuous(gloP) && gloP(last(gloDomains[i])) != locP(first(domain(locP))) 
                error("new poly \"N\" lower endpoint y value is not coincidence with GlobalPoly \"G\" (G=$(gloP(last(gloDomains[i]))) vs N=$(locP(first(domain(locP))))). Use a non continuous GlobalPoly if this is what you want.")
            end
            if first(gloDomains[i+1]) == last(gloDomain) && iscontinuous(gloP) && gloP(first(gloDomains[i+1])) != locP(last(domain(locP))) 
                error("new poly \"N\" upper endpoint y value is not coincidence with GlobalPoly \"G\" (G=$(gloP(first(gloDomains[i+1]))) vs N=$(locP(last(domain(locP))))). Use a non continuous GlobalPoly if this is what you want.")
            end
            insert!(gloP.polys, i+1, locP)
            return insert!(gloP.domains, i+1, gloDomain)
        end
    end

    conflict = intersect(gloDomain, gloDomains[end])
    error("new domain endpoints intersect existing polynomial. Intersection is $conflict")
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
    low,up = domains(p)[1], domains(p)[end]
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

function update!(gloP::GlobalPoly{T}, y::AbstractVector{T}) where {T}
    gDomains = domains(gloP)
    continuouspLengths = [ iscontinuous(gloP) && last(gDomains[i-1]) == first(gDomains[i]) ? length(gloP[i]) - 1 : length(gloP[i]) for i in 2:length(gloP)]
    pushfirst!(continuouspLengths,length(gloP[1]))
    sum(continuouspLengths) !=length(y) && throw(DimensionMismatch("tried to assign $(length(y)) elements to $(sum(continuouspLengths)) destinations"))
    Threads.@threads for i in 1:length(gloP)
        if continuouspLengths[i] == length(gloP[i]) # case where endpoints are not coincident
            update!(gloP[i],y[sum(continuouspLengths[1:i-1])+1:sum(continuouspLengths[1:i])])
        else
            update!(gloP[i],y[sum(continuouspLengths[1:i-1]):sum(continuouspLengths[1:i])])
        end
    end
    return gloP
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


function checkends(gloD::Interval,locP::AbstractPolynomial, iscontinuous::Bool)
    lDomain = domain(locP)
    if iscontinuous
        ( !containslower(lDomain) || !containsupper(lDomain) ) && throw(DomainError(locP),"Continuous GlobalPoly can only contain local polynomial that have closed endpoints")
    else
        !containslower(lDomain) && containslower(gloD) && throw(DomainError(gloD,"global lower endpoints incompatible with local polynomial domain $(lDomain). Consider either opening the global endpoint or closing the local endpoint."))#check first that local and global bound endpoints are compatible
        !containsupper(lDomain) && containsupper(gloD) && throw(DomainError(gloD,"global upper endpoints incompatible with local polynomial $(lDomain). Consider either opening the global endpoint or closing the local endpoint."))
    end
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
                return mapendpoints(Interval(-1,1),locD,(x-1)/(1+x)) # allow user to specify custom function here?
            else
                return mapendpoints(Interval(-1,1),locD,(-x-1)/(x-1))
            end
        elseif isinf(gUpper)
            return x - gLower + lLower
        else
            return x - gUpper + lUpper
        end

    end
    return lLower + (lUpper - lLower) * (x - gLower) / ( gUpper - gLower )
end

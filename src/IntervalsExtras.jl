
# some Intervals functions that make my life easier

containsfirst(i::Interval) = first(inclusivity(i))
containslast(i::Interval) = last(inclusivity(i))

function Base.convert(::Type{Interval{T}}, i::Interval{N}) where {T,N}
    ifirst = convert(T, first(i))
    ilast = convert(T, last(i))
    return Interval(ifirst, ilast, containsfirst(i), containslast(i) )
end

Base.promote_rule(::Type{Interval{T}}, ::Type{Interval{S}}) where {T,S} = Interval{promote_type(T,S)}

function Base.intersect(a::AbstractInterval{T}, b::AbstractInterval{S}) where {T,S}
    a,b = promote(a,b)
    return intersect(a,b)
end

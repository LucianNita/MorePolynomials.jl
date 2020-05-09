function legendrePoly(::Type{T},order::Int) where {T}
    order<0 && error("order must be >= 0")
    P = [Polynomial{T}(1), Polynomial{T}([0,1])]
    x = Polynomial{T}([0,1])
    order < 2 && return P[order+1]
    for i = 1:order-1
        Ptmp = ( (2*i + 1 ) * x * P[2] - i * P[1]) / (i+1)
        P[1] = P[2]
        P[2] = Ptmp
    end
    return P[2]
end

legendrePoly(order::T) where {T} = legendrePoly(T,convert(Int,order))

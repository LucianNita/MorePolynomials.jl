macro register(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        #Base.convert(P::Type{<:$poly}, p::$poly) where {T} = P(coeffs(p), p.var) # needs replacing
        Base.promote_rule(::Type{$poly{T}}, ::Type{$poly{S}}) where {T,S} =
            $poly{promote_type(T, S)}
        Base.promote_rule(::Type{$poly{T}}, ::Type{S}) where {T,S<:Number} =
            $poly{promote_type(T, S)}

        function (p::$poly)(x::AbstractVector)
            Base.depwarn(
                "Calling p(x::AbstractVector is deprecated. Use p.(x) instead.",
                Symbol("(p::AbstractPolynomial)"),
            )
            return p.(x)
        end

        #$poly(coeffs::AbstractVector{T}, var::SymbolLike = :x) where {T} =
        #    $poly{T}(coeffs, Symbol(var))
        $poly(n::Number, var = :x) = $poly([n], var)
        $poly{T}(n::S, var = :x) where {T,S<:Number} = $poly(T(n), var)
        $poly{T}(x::AbstractVector{S}, var = :x) where {T,S<:Number} = $poly(T.(x), var)
    end
end

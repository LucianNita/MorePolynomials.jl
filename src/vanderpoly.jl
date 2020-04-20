export vandpoly
using Reduce # for genvandpoly macro 
Reduce.Preload() # required otherwise Reduce has a fit

# The Vandermonde matrix method

macro genvandpoly(index,leny,x,y) # generate poly functions using the Vandermonde matrix method
    leny = :($leny)
    symX = []
    symY = []
    A = Array{Expr}(undef, leny, leny)
    for i = 1:leny
        push!(symX, Meta.parse("x$i"))
        push!(symY, Meta.parse("y$i"))
    end
    for i = 1:leny # Assemble A matrix
        for j = 1:leny
            pwr = leny-j
            var = symX[i]
            A[i,j] = :($var^$pwr)
        end
    end
    coefs = Algebra.:*(Algebra.inv(A),symY) # invert to find coefs of polynomial
    for i = 1:leny
        str = repr(coefs[i])
        for j = 1:leny
            str = replace(str,"y"*string(j)=>"$y"*"["*string(j)*"]")
            str = replace(str,"x"*string(j)=>"$x"*"["*string(j)*"]")
        end
        coefs[i] = Meta.parse(chop(str, head=2, tail=1))
    end

    return coefs[:($index)]
end

function vandpoly(x,y)
    xLength = length(x)
    if xLength == 2
        c1 = @genvandpoly(1, 2,:($x),:($y))
        c2 = @genvandpoly(2, 2,:($x),:($y))
        return Polynomial([c2,c1])
    elseif xLength == 3
        c1 = @genvandpoly(1, 3,:($x),:($y))
        c2 = @genvandpoly(2, 3,:($x),:($y))
        c3 = @genvandpoly(3, 3,:($x),:($y))
        return Polynomial([c3,c2,c1])
    elseif xLength == 4
        c1 = @genvandpoly(1, 4,:($x),:($y))
        c2 = @genvandpoly(2, 4,:($x),:($y))
        c3 = @genvandpoly(3, 4,:($x),:($y))
        c4 = @genvandpoly(4, 4,:($x),:($y))
        return Polynomial([c4,c3,c2,c1])
    elseif xLength == 5
        c1 = @genvandpoly(1, 5,:($x),:($y))
        c2 = @genvandpoly(2, 5,:($x),:($y))
        c3 = @genvandpoly(3, 5,:($x),:($y))
        c4 = @genvandpoly(4, 5,:($x),:($y))
        c5 = @genvandpoly(5, 5,:($x),:($y))
        return Polynomial([c5,c4,c3,c2,c1])
    else 
        Polynomial(inv([xval^j for xval in eachindex(x), j = 0:xlength - 1]) * y) # create poly coefficient by inverting array of linear equations
        println("Warning, function behaves slow when exceeding more than 5 datapoints")
    end
end

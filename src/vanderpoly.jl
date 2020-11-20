export vandpoly #export function vandpoly, make it public
Reduce.Preload() ## required otherwise Reduce has a fit -?

## ##The Vandermonde matrix method #https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/05Interpolation/vandermonde/#:~:text=The%20Vandermonde%20method%20is%20the,Gaussian%20elimination%20or%20PLU%20decomposition.


##Create macro genvandpoly with:
# INPUTS:
#  index - coefficient index we want to return
#  leny - number of points available for leny-1 degree polynomial interpolation
#  x - array with x coordinates of leny points
#  y - array with y coordinates of leny points
# OUTPUTS:
#  coefs[index] - the [$index]'s coefficient

macro genvandpoly(index,leny,x,y) ## generate poly functions using the Vandermonde matrix method
    leny = :($leny) #length of y #expression with number leny
    symX = [] #symmetric X matrix declaration
    symY = [] #symmetric Y matrix declaration
    A = Array{Expr}(undef, leny, leny) #Array of expressions (ie. Vandermonde matrix)
##we push an expression xᵢ,yᵢ into symX,symY
    for i = 1:leny #for all integers i from 1 to leny
        push!(symX, Meta.parse("x$i")) #meta.parse - parses a string to an expression
        push!(symY, Meta.parse("y$i")) #meta.parse - parses a string to an expression
    end
## ## Assemble A matrix #we go over all elements in A ∈ℜⁿˣⁿ matrix #we go over all elements in A and construct them as expressions
    for i = 1:leny  #i is the row number (ie. equation or point)
        for j = 1:leny #j is the column number (ie. power of x)
            pwr = leny-j #power is the total number of columns - index of current column
            var = symX[i] #variable that we need to raise to the power, took from symX
            A[i,j] = :($var^$pwr) #each element of A is an expression var^pwr, with var and pwr being numbers
        end
    end
## ## invert Vandermonde matrix "A" to find coefficients of polynomial
    coefs = Algebra.:*(Algebra.inv(A),symY) #coefs stores polynomial coefficients and is coefs=A\symY
    for i = 1:leny #for all integers i in [1,leny] (ie. all elements in coefs)
        str = repr(coefs[i]) #str is a string representation of element i from coefs
        for j = 1:leny
            str = replace(str,"y"*string(j)=>"$y"*"["*string(j)*"]") #replaces into str, the string yⱼ with the numerical value of y[j]
            str = replace(str,"x"*string(j)=>"$x"*"["*string(j)*"]") #replaces into str, the string xⱼ with the numerical value of x[j]
        end
        coefs[i] = Meta.parse(chop(str, head=2, tail=1)) #remove first 2 and last 1 characters from string str and parse it as an expression
    end
##return the coefficient with index number "index"
    return coefs[:($index)]
end


## Function vandpoly
#  INPUTS:
# x - vector with x coordinates from n points
# y - vector with y coordiantes from n points
#  OUTPUTS:
# Polynomial([cn,cn-1,...,c2,c1]) - Polynomial expression of the form c₁xⁿ⁻¹+c₂xⁿ⁻²+...+cₙ₋₁x+cₙ

function vandpoly(x,y) #used to return a polynomial of degree n-1 out of n datapoints in ℜ²
    xLength = length(x) #xLength is the length of x vector

    if xLength == 2 #if we have only 2 points we fit y=c₁x¹+c₂
        c1 = @genvandpoly(1, 2,:($x),:($y)) #c1 is the coefficient of x¹
        c2 = @genvandpoly(2, 2,:($x),:($y)) #c2 is the coefficient of x⁰
        return Polynomial([c2,c1]) #Return the polynomial c₁x¹+c₂

    elseif xLength == 3 #if we have only 3 points we fit y=c₁x²+c₂x¹+c₃
        c1 = @genvandpoly(1, 3,:($x),:($y)) #c1 is the coefficient of x²
        c2 = @genvandpoly(2, 3,:($x),:($y)) #c2 is the coefficient of x¹
        c3 = @genvandpoly(3, 3,:($x),:($y)) #c3 is the coefficient of x⁰
        return Polynomial([c3,c2,c1]) #Return the polynomial c₁x²+c₂x¹+c₃

    elseif xLength == 4 #if we have only 4 points we fit y=c₁x³+c₂x²+c₃x¹+c₄
        c1 = @genvandpoly(1, 4,:($x),:($y)) #c1 is the coefficient of x³
        c2 = @genvandpoly(2, 4,:($x),:($y)) #c2 is the coefficient of x²
        c3 = @genvandpoly(3, 4,:($x),:($y)) #c3 is the coefficient of x¹
        c4 = @genvandpoly(4, 4,:($x),:($y)) #c4 is the coefficient of x⁰
        return Polynomial([c4,c3,c2,c1]) #Return the polynomial c₁x³+c₂x²+c₃x¹+c₄

    elseif xLength == 5 #if we have only 5 points we fit y=c₁x⁴+c₂x³+c₃x²+c₄x¹+c₅
        c1 = @genvandpoly(1, 5,:($x),:($y)) #c1 is the coefficient of x⁴
        c2 = @genvandpoly(2, 5,:($x),:($y)) #c2 is the coefficient of x³
        c3 = @genvandpoly(3, 5,:($x),:($y)) #c3 is the coefficient of x²
        c4 = @genvandpoly(4, 5,:($x),:($y)) #c4 is the coefficient of x¹
        c5 = @genvandpoly(5, 5,:($x),:($y)) #c5 is the coefficient of x⁰
        return Polynomial([c5,c4,c3,c2,c1]) #Return the polynomial c₁x⁴+c₂x³+c₃x²+c₄x¹+c₅

    else #if we have more than 5 points we simply perform C\y
        Polynomial(inv([xval^j for xval in eachindex(x), j = 0:xlength - 1]) * y) # create poly coefficient by inverting array of linear equations
        println("Warning, function behaves slow when exceeding more than 5 datapoints")
    end
end

export legendrePoly #Just a basic export to make legendrePoly function visible (ie. Public)

## Function that computes the legendre polynomial of nᵗʰ order
#  INPUTS:
# Type T - ?????
# order - Integer defining the desired Legendre polynomial order
#  OUTPUTS:
# Pₙ(x) - Legendre polynomial of order n="order"

function legendrePoly(::Type{T},order::Int) where {T} #compute Legendre polynomial
##Throw an error if we ask for a polynomial of negative order
    order<0 && error("order must be >= 0") #If order is less than zero, we throw an error
##Manually input legendre polynomials of orders 0 and 1
    P = [Polynomial{T}(1), Polynomial{T}([0,1])] #initialize our horizon for pseudo-recursion with known Legendre Polynomials: P=[1,x]
    x = Polynomial{T}([0,1]) #x is just an instance of type Polynomial
##For orders 0 or 1 we return 1 or x respectively
    order<2 && return P[order+1] #if order is 1 return x, if order is 0 return 1
## For order>=2: we are using Bonnet recursion formula: nPₙ(x)=(2n-1)Pₙ₋₁(x)-(n-1)Pₙ₋₂(x)
# Compute Pᵢ₊₁ from Pᵢ and Pᵢ₋₁ in a pseudo-recursive way
    for i = 1:order-1 #i is effectively n-1 (ie. current order-1) because we index from 1 to "order-1" instead of the natural 2 to "order"
        Ptmp = ( (2*i + 1 ) * x * P[2] - i * P[1]) / (i+1) #this computes Pᵢ₊₁ in a pseudo recursive manner
        P[1] = P[2] #here we save Pᵢ for next iteration
        P[2] = Ptmp #here we save Pᵢ₊₁ for next iteration
    end
## Return P of order "order"
    return P[2]
end

legendrePoly(order::T) where {T} = legendrePoly(T,convert(Int,order)) #I suspect here he tried to do some proper recursion, but idk

module MorePolynomials

using Polynomials
import Polynomials.fit
import Polynomials.domain
import Polynomials.coeffs
import Polynomials.showterm
using LinearAlgebra
using Intervals
const SymbolLike = Union{AbstractString,Char,Symbol}
include("Lagrange.jl")

end # module

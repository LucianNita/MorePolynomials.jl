module MorePolynomials

using Polynomials
import Polynomials.fit
import Polynomials.domain
import Polynomials.coeffs
import Polynomials.showterm
using LinearAlgebra
using Intervals
using Reduce # for genvandpoly macro 
const SymbolLike = Union{AbstractString,Char,Symbol}

include("register.jl")
include("Lagrange.jl")
include("vanderpoly.jl")
include("polyhelper.jl")

end # module

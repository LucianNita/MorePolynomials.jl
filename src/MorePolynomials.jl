module MorePolynomials

using Polynomials
import Polynomials.fit
import Polynomials.domain
import Polynomials.coeffs
import Polynomials.showterm
import Polynomials.derivative
import Polynomials.integrate
import Base.convert
import Base.eachindex
import Base.length
using LinearAlgebra
using Intervals
using Reduce # for genvandpoly macro 
using FastGaussQuadrature

const SymbolLike = Union{AbstractString,Char,Symbol}

include("legendre.jl")
include("IntervalsExtras.jl")
include("register.jl")
include("Lagrange.jl")
include("Global.jl")
include("vanderpoly.jl")
include("polyhelper.jl")

end # module

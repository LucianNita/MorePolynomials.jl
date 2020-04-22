# STILL UNDER DEVELOPMENT

![MP logo](https://github.com/sprockmonty/MorePolynomials.jl/blob/master/docs/src/assets/MorePolynomials.png "MorePolynomials logo")

---

# MorePolynomials.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sprockmonty.github.io/MorePolynomials.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sprockmonty.github.io/MorePolynomials.jl/dev)
[![Build Status](https://travis-ci.com/sprockmonty/MorePolynomials.jl.svg?branch=master)](https://travis-ci.com/sprockmonty/MorePolynomials.jl)
[![Codecov](https://codecov.io/gh/sprockmonty/MorePolynomials.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sprockmonty/MorePolynomials.jl)

Fastest Julia Lagrange interpolation or your money back!

## Installation 
```julia
(v1.4) pkg> add https://github.com/sprockmonty/MorePolynomials.jl
julia> using Polynomials, MorePolynomials
```
## Usage
MorePolynomials extends from Polynomials.jl and adds additional polynomials and functionality.

For example, to create a Lagrange Polynomial based on x and y values, use the `fit()` function.

```julia
julia> fit(LagrangePoly, [-1,0,1],[25,0,25])
```
Alternatively, you can use the polyhelper function to guess what kind of polynomial you want. Lower orders will return `Polynomial` generated analytically.
```julia
julia> polyhelper([1,2,3],[1,2,3])
Polynomial(1.0*x)
```

Higher orders will be approximated as Lagrange polynomials of type `LagrangePoly`, and a single array input will return a Legendre–Gauss–Radau based polynomial.
```julia
julia> polyhelper([1,2,3])
LGRPoly(1 2 3)
```

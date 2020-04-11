using Documenter, MorePolynomials

makedocs(;
    modules=[MorePolynomials],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/sprockmonty/MorePolynomials.jl/blob/{commit}{path}#L{line}",
    sitename="MorePolynomials.jl",
    authors="Nathan Davey",
    assets=String[],
)

deploydocs(;
    repo="github.com/sprockmonty/MorePolynomials.jl",
)

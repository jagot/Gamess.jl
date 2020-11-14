using Gamess
using Documenter

makedocs(;
    modules=[Gamess],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/jagot/Gamess.jl/blob/{commit}{path}#L{line}",
    sitename="Gamess.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jagot.github.io/Gamess.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jagot/Gamess.jl",
)

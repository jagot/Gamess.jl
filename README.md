# Gamess

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jagot.github.io/Gamess.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jagot.github.io/Gamess.jl/dev)
[![Build Status](https://github.com/jagot/Gamess.jl/workflows/CI/badge.svg)](https://github.com/jagot/Gamess.jl/actions)
[![Coverage](https://codecov.io/gh/jagot/Gamess.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jagot/Gamess.jl)

Simple library to parse selected pieces of output from
[GAMESS-US](https://www.msg.chem.iastate.edu/gamess/).

Run
```julia
using Gamess
Gamess.load(filename)
```
where `filename` is a file containing the textual output from a GAMESS
run. At the moment, it will read some metadata, the coordinates and
charges of the nuclei in the  molecule, and CIS excited states and
dipole transition moments (if present, supports `CITYP=CIS` and `GUGA`).

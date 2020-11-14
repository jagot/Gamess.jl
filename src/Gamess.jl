module Gamess

import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp

using DataFrames
using LinearAlgebra

function read_until(fun::Function, io::IO, pred::Function)
    while !eof(io)
        l = readline(io)
        pred(l) && break
        fun(l)
    end
end

read_until(io::IO, pred::Function) =
    read_until(l -> nothing, io, pred)

function parse_input_section(section)
    name,tail = split(strip(section), limit=2)
    if name == "\$DATA"
        :DATA => tail
    else
        settings = map(split(tail)) do s
            a,b = split(strip(s), "=")
            Symbol(a) => b
        end
        Symbol(lstrip(name, '$')) => (;settings...)
    end
end

function read_input(io::IO)
    input = Vector{String}()
    read_until(io, l -> occursin("ECHO OF THE FIRST FEW INPUT CARDS", l))
    read_until(io, l -> !occursin("INPUT CARD>", l)) do l
        l = lstrip(l[13:end])
        length(l) > 0 && first(l) == '!' && return
        push!(input, l)
    end
    input = filter(!isempty, map(strip, split(join(input, "\n"), "\$END")))
    (; map(parse_input_section, input)...)
end

function find_title(io::IO)
    read_until(io, l -> occursin("RUN TITLE", l))
    eof(io) && return
    readline(io)
    strip(readline(io))
end

prefloat = re"[-+]?([0-9]+\.[0-9]*|[0-9]*\.[0-9]+)"
float    = prefloat | re.cat(prefloat | re"[-+]?[0-9]+", re"[eE][-+]?[0-9]+")

float.actions[:enter] = [:mark]
float.actions[:exit] = [:float]

function display_machine(machine)
    write("/tmp/gamess_machine.dot", Automa.machine2dot(machine))
    run(`dot -Tsvg -o /tmp/gamess_machine.svg /tmp/gamess_machine.dot`)
end

atom_entry_machine = let
    atom_label = re"[A-Za-z]+"

    charge = float
    x = float
    y = float
    z = float

    atom_label.actions[:enter] = [:mark]
    atom_label.actions[:exit] = [:atom_label]

    atom_entry = re.cat(re" *", atom_label, re" +",
                        charge, re" +",
                        x, re" +",
                        y, re" +",
                        z)

    Automa.compile(atom_entry)
end

# display_machine(atom_entry_machine)

atom_entry_actions = Dict(
    :mark => :(mark = p),
    :atom_label => :(atom_label = String(data[mark:p-1])),
    :float => :(push!(atom_data, parse(Float64, data[mark:p-1])))
)

context = Automa.CodeGenContext()
@eval function read_atom(data::String)
    atom_label = nothing
    atom_data = Vector{Float64}()
    mark = 0

    $(Automa.generate_init_code(context, atom_entry_machine))

    # p_end and p_eof were set to 0 and -1 in the init code,
    # we need to set them to the end of input, i.e. the length of `data`.
    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, atom_entry_machine, atom_entry_actions))

    # We need to make sure that we reached the accept state, else the
    # input did not parse correctly
    iszero(cs) || error("failed to parse on byte ", p)

    (label=atom_label, charge=atom_data[1],
     x=atom_data[2], y=atom_data[3], z=atom_data[4])
end

function read_molecule(io::IO)
    point_group = nothing
    principal_axis_order = 0
    read_until(io, l -> occursin("ATOMIC", l)) do l
        if occursin("POINT GROUP", l)
            point_group=last(rsplit(l, limit=2))
        elseif occursin("PRINCIPAL", l)
            principal_axis_order=parse(Int, last(rsplit(l, limit=2)))
        end
    end
    readline(io)

    labels = Vector{String}()
    charges = Vector{Float64}()
    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()

    read_until(io, isempty) do l
        atom = read_atom(l)
        push!(labels, atom.label)
        push!(charges, atom.charge)
        push!(x, atom.x)
        push!(y, atom.y)
        push!(z, atom.z)
    end
    atoms = DataFrame(Atom=labels, Charge=charges,
                         x=x, y=y, z=z)

    (point_group=point_group,
     principal_axis_order=principal_axis_order,
     atoms=atoms)
end

function load(io::IO)
    input = read_input(io)

    title = find_title(io)
    isnothing(title) && throw(ArgumentError("Could not extract run title from $(io)"))
    molecule = read_molecule(io)

    (title=title,
     input=input,
     molecule=molecule)
end

load(filename::String) = open(load, filename)

end

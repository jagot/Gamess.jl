module Gamess

import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp

using DataFrames
using LinearAlgebra

function read_until(fun::Function, io::IO, pred::Function)
    l = ""
    while !eof(io)
        l = readline(io)
        pred(l) && break
        fun(l)
    end
    l
end

read_until(fun::Function, io::IO, s::String) =
    read_until(fun, io, l -> occursin(s, l))

read_until(io::IO, pred) =
    read_until(l -> nothing, io, pred)

# * Input

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
    read_until(io, "ECHO OF THE FIRST FEW INPUT CARDS")
    read_until(io, l -> !occursin("INPUT CARD>", l)) do l
        l = lstrip(l[13:end])
        length(l) > 0 && first(l) == '!' && return
        push!(input, l)
    end
    input = filter(!isempty, map(strip, split(join(input, "\n"), "\$END")))
    (; map(parse_input_section, input)...)
end

function find_title(io::IO)
    read_until(io, "RUN TITLE")
    eof(io) && return
    readline(io)
    strip(readline(io))
end

# * Molecule

dec      = re"[-+]?[0-9]+"
prefloat = re"[-+]?([0-9]+\.[0-9]*|[0-9]*\.[0-9]+)"
float    = prefloat | re.cat(prefloat | re"[-+]?[0-9]+", re"[eE][-+]?[0-9]+")

dec.actions[:enter] = [:mark]
dec.actions[:exit] = [:dec]
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
    read_until(io, "ATOMIC") do l
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

# * Dipole moments

dipole_machine = let
    dipole = re.cat(re".+", re"= +",
                    float, re" +",
                    float, re" +",
                    float, re" +",
                    float, re" +",
                    re"E.BOHR")

    Automa.compile(dipole)
end

dipole_actions = Dict(
    :mark => :(mark = p),
    :float => :(push!(values, parse(Float64, data[mark:p-1])))
)

context = Automa.CodeGenContext()
@eval function read_dipole(data::String)
    values = Vector{Float64}()
    mark = 0

    $(Automa.generate_init_code(context, dipole_machine))
    p_end = p_eof = lastindex(data)
    $(Automa.generate_exec_code(context, dipole_machine, dipole_actions))
    iszero(cs) || error("failed to parse on byte ", p)

    values
end

# * Einstein coefficients

einstein_machine = let
    einstein = re.cat(re".+", re"EINSTEIN COEFFICIENTS: A= +",
                      float, re" +1/SEC; B= +",
                      float, re" +SEC/G")

    Automa.compile(einstein)
end

einstein_actions = Dict(
    :mark => :(mark = p),
    :float => :(push!(values, parse(Float64, data[mark:p-1])))
)

context = Automa.CodeGenContext()
@eval function read_einstein(data::String)
    values = Vector{Float64}()
    mark = 0

    $(Automa.generate_init_code(context, einstein_machine))
    p_end = p_eof = lastindex(data)
    $(Automa.generate_exec_code(context, einstein_machine, einstein_actions))
    iszero(cs) || error("failed to parse on byte ", p)

    values
end

# * CIS


cis_state_machine = let
    sym = re"[A-Z0-9]+"

    sym.actions[:enter] = [:mark]
    sym.actions[:exit] = [:sym]

    cis_state = re.cat(re" *", re"EXCITED STATE", re" +",
                       dec, re" +",
                       re"ENERGY=", re" +",
                       float, re" +",
                       re"S =", re" +",
                       float, re" +",
                       re"SPACE SYM =", re" +",
                       sym, re" *")

    Automa.compile(cis_state)
end

cis_state_actions = Dict(
    :mark => :(mark = p),
    :sym => :(state = String(data[mark:p-1])),
    :dec => :(),
    :float => :(push!(values, parse(Float64, data[mark:p-1])))
)

context = Automa.CodeGenContext()
@eval function read_cis_state(data::String)
    state = ""
    values = Vector{Float64}()
    mark = 0

    $(Automa.generate_init_code(context, cis_state_machine))

    # p_end and p_eof were set to 0 and -1 in the init code,
    # we need to set them to the end of input, i.e. the length of `data`.
    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, cis_state_machine, cis_state_actions))

    # We need to make sure that we reached the accept state, else the
    # input did not parse correctly
    iszero(cs) || error("failed to parse on byte ", p)

    state, values[1], values[2]
end

function read_cis_dipoles(io::IO)
    left = Vector{Int}()
    right = Vector{Int}()
    g_left = Vector{Int}()
    g_right = Vector{Int}()
    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()
    f = Vector{Float64}()
    A = Vector{Float64}()
    B = Vector{Float64}()

    read_until(io, "CIS TRANSITION DIPOLE MOMENTS AND")
    read_until(io, "-CIS- LAGRANGIAN") do l
        if occursin("TRANSITION FROM THE GROUND STATE TO EXCITED STATE", l)
            if !occursin("OPTICALLY INACTIVE", l)
                push!(left, 0)
                push!(right, parse(Int, l[51:end]))
            end
        elseif occursin("TRANSITION BETWEEN EXCITED STATES", l)
            if !occursin("OPTICALLY INACTIVE", l)
                le,ri = split(l[35:end], "AND")
                push!(left, parse(Int, le))
                push!(right, parse(Int, ri))
            end
        elseif occursin("EXPECTATION VALUE DIPOLE MOMENT FOR EXCITED STATE", l)
            le = parse(Int, l[51:end])
            push!(left, le)
            push!(right, le)
            push!(f, 0)
            push!(A, 0)
            push!(B, 0)
        elseif occursin("STATE MULTIPLICITIES", l)
            g_l,g_r = split(strip(split(l, "=")[2]))
            push!(g_left, parse(Int, g_l))
            push!(g_right, parse(Int, g_r))
        elseif occursin("STATE MULTIPLICITY", l)
            g = parse(Int, split(l, "=")[2])
            push!(g_left, g)
            push!(g_right, g)
        elseif occursin(r"TRANSITION DIPOLE =", l) && occursin(r"E\*BOHR", l)
            values = read_dipole(l)
            push!(x, values[1])
            push!(y, values[2])
            push!(z, values[3])
        elseif occursin(r"STATE DIPOLE =", l) && occursin(r"E\*BOHR", l)
            values = read_dipole(l)
            push!(x, values[1])
            push!(y, values[2])
            push!(z, values[3])
        elseif occursin("OSCILLATOR STRENGTH", l)
            push!(f, parse(Float64, strip(split(l, "=")[2])))
        elseif occursin("EINSTEIN COEFFICIENTS", l)
            values = read_einstein(l)
            push!(A, values[1])
            push!(B, values[2])
        end
    end

    DataFrame(left=left, right=right,
              g_left=g_left, g_right=g_right,
              x=x, y=y, z=z,
              f=f, A=A, B=B)
end

function read_cis(io::IO)
    states = Vector{String}()
    energies = Vector{Float64}()
    Ss = Vector{Float64}()

    read_until(io, "CI-SINGLES EXCITATION ENERGIES") do l
        if occursin("EXCITED STATE", l)
            state, energy, S = read_cis_state(l)
            push!(states, state)
            push!(energies, energy)
            push!(Ss, S)
        end
    end

    length_gauge_dipole = read_cis_dipoles(io)

    (excited_states=DataFrame(State = states, Energy = energies, S=Ss),
     length_gauge_dipole=length_gauge_dipole)
end

# * GUGA

guga_state_machine = let
    guga_state = re.cat(re" *", re"STATE #", re" +",
                        dec, re" +",
                        re"ENERGY =", re" +",
                        float)

    Automa.compile(guga_state)
end

guga_state_actions = Dict(
    :mark => :(mark = p),
    :dec => :(state = parse(Int, data[mark:p-1])),
    :float => :(energy = parse(Float64, data[mark:p-1]))
)

context = Automa.CodeGenContext()
@eval function read_guga_state(data::String)
    state = 0
    energy = 0.0
    mark = 0

    $(Automa.generate_init_code(context, guga_state_machine))

    # p_end and p_eof were set to 0 and -1 in the init code,
    # we need to set them to the end of input, i.e. the length of `data`.
    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, guga_state_machine, guga_state_actions))

    # We need to make sure that we reached the accept state, else the
    # input did not parse correctly
    iszero(cs) || error("failed to parse on byte ", p)

    state, energy
end

function read_guga_dipoles(io::IO, stop)
    left = Vector{Int}()
    right = Vector{Int}()
    g_left = Vector{Int}()
    g_right = Vector{Int}()
    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()
    f = Vector{Float64}()
    A = Vector{Float64}()
    B = Vector{Float64}()

    same = false

    read_until(io, stop) do l
        if occursin("CI STATE NUMBER=", l)
            le,ri = parse.(Int, split(split(l, "=", limit=2)[2])[1:2])
            if le ≠ ri
                push!(left, le)
                push!(right, ri)
                g_l,g_r = parse.(Int, split(rsplit(l, "=", limit=2)[2])[1:2])
                push!(g_left, g_l)
                push!(g_right, g_r)
                same = false
            else
                same = true
            end
        elseif occursin(r"E.BOHR", l)
            if !same
                values = read_dipole(l)
                push!(x, values[1])
                push!(y, values[2])
                push!(z, values[3])
            end
        elseif occursin("OSCILLATOR STRENGTH", l)
            # This is nuts; GUGA output has = as divider in the length
            # form and IS in the velocity form.
            divider = occursin("IS", l) ? "IS" : "="
            same || push!(f, parse(Float64, strip(split(l, divider)[2])))
        elseif occursin("EINSTEIN COEFFICIENTS", l)
            if !same
                values = read_einstein(l)
                push!(A, values[1])
                push!(B, values[2])
            end
        end
    end

    df = DataFrame(left=left, right=right,
                   g_left=g_left, g_right=g_right,
                   x=x, y=y, z=z,
                   f=f)
    # No Einstein coefficients in the velocity form listing
    isempty(A) ? df : hcat(df, DataFrame(A=A, B=B))
end

function read_guga(io::IO)
    read_until(io, "NON-ABELIAN CI WAVEFUNCTION STATE SYMMETRY DRIVER")
    for _ = 1:4
        readline(io)
    end
    read_until(io, l -> !occursin("ORBITAL", l))

    states = Vector{Int}()
    energies = Vector{Float64}()

    read_until(io, "END OF CI-MATRIX DIAGONALIZATION") do l
        if occursin("STATE #", l)
            state, energy = read_guga_state(l)
            push!(states, state)
            push!(energies, energy)
        end
    end

    read_until(io, "LENGTH FORM")
    length_gauge_dipole = read_guga_dipoles(io, "VELOCITY FORM")
    velocity_gauge_dipole = read_guga_dipoles(io, "DONE WITH TRANSITION MOMENT")

    (excited_states=DataFrame(State = states, Energy = energies),
     length_gauge_dipole=length_gauge_dipole,
     velocity_gauge_dipole=velocity_gauge_dipole)
end

# * Energies

function read_energies(io::IO)
    read_until(io, "ENERGY COMPONENTS")
    values = Vector{Pair{Symbol,Float64}}()
    l = read_until(io, "VIRIAL RATIO") do l
        if occursin("=", l)
            label, value = split(l, "=")
            push!(values, Symbol(replace(replace(strip(lowercase(label)), " " => "_"), "-" => "_")) => parse(Float64, value))
        end
    end
    push!(values, :virial_ratio => parse(Float64, split(l, "=")[2]))
    (; values...)
end

# * Interface

function read_ci(io::IO, ci_type)
    if ci_type == "CIS"
        read_cis(io)
    elseif ci_type == "GUGA"
        read_guga(io)
    else
        throw(ArgumentError("Cannot parse CI calculation of type $(ci_type)"))
    end
end

function load(io::IO)
    input = read_input(io)

    title = find_title(io)
    isnothing(title) && throw(ArgumentError("Could not extract run title from $(io)"))
    molecule = read_molecule(io)

    data = (title=title,
            input=input,
            molecule=molecule)

    if lowercase(input.CONTRL.SCFTYP) ≠ "none"
        data = merge(data, (energies=read_energies(io),))
    end

    if :CITYP in keys(input.CONTRL)
        data = merge(data, (ci=read_ci(io, input.CONTRL.CITYP),))
    end

    data
end

load(filename::String) = open(load, filename)

end

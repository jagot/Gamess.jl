using Gamess
using Test

import Gamess: Configuration, substitutions

load_file(name) = Gamess.load(joinpath(dirname(@__FILE__), "data", "$(name).out"))

@testset "Gamess.jl" begin
    @testset "Water CIS" begin
        water = load_file("water")
        atoms = water.molecule.atoms
        @test size(atoms) == (3, 5)
        @test atoms.Atom == ["OXYGEN", "HYDROGEN", "HYDROGEN"]

        ci = water.ci
        states = ci.excited_states
        @test size(states) == (100,3)

        dipoles = ci.length_gauge_dipole
        @test size(dipoles) == (2601,10)
        @test dipoles[1, :z] ≈ 0.958273779741262
        @test dipoles[2599,:z] ≈ 0.044333
        @test dipoles[2599,:A] ≈ -5.4312e-41
        @test dipoles[2599,:B] ≈ 797650.0

        expansions = ci.state_expansions
        @test length(expansions) == 100

        e = first(expansions)
        reference = e.reference
        N = 13
        @test length(reference.α) == length(reference.β) == N
        @test reference == Configuration(vcat(trues(5), falses(N-5)), vcat(trues(5), falses(N-5)))

        @test substitutions(reference, e.configurations[1]) == ([5 => 6], [])
        @test substitutions(reference, e.configurations[2]) == ([5 => 9], [])
        @test substitutions(reference, e.configurations[3]) == ([],[5 => 6])
        @test substitutions(reference, e.configurations[4]) == ([],[5 => 9])
    end

    @testset "Neon CIS" begin
        @testset "Ne I" begin
            neon = load_file("ne-cis")
            atoms = neon.molecule.atoms
            @test size(atoms) == (1, 5)
            @test atoms.Atom == ["NE"]

            ci = neon.ci
            states = ci.excited_states
            @test size(states) == (100,3)

            dipoles = ci.length_gauge_dipole
            @test size(dipoles) == (2603,10)
            @test dipoles[4, :z] ≈ -0.250121

            expansions = ci.state_expansions
            @test length(expansions) == 100

            e = first(expansions)
            reference = e.reference
            N = 121
            @test length(reference.α) == length(reference.β) == N
            @test reference == Configuration(vcat(trues(5), falses(N-5)), vcat(trues(5), falses(N-5)))

            @test substitutions(reference, e.configurations[1]) == ([4 => 10], [])
            @test substitutions(reference, e.configurations[2]) == ([4 => 19], [])
            @test substitutions(reference, e.configurations[13]) == ([],[4 => 10])
            @test substitutions(reference, e.configurations[14]) == ([],[4 => 19])
        end

        @testset "Ne II" begin
            neon = load_file("ne+-guga")
            atoms = neon.molecule.atoms
            @test size(atoms) == (1, 5)
            @test atoms.Atom == ["NE"]

            ci = neon.ci
            states = ci.excited_states
            @test size(states) == (5,2)

            dipoles = ci.length_gauge_dipole
            @test size(dipoles) == (15,10)
            @test dipoles[11, :z] ≈ -0.521326

            dipoles = ci.velocity_gauge_dipole
            @test size(dipoles) == (10,8)
        end
    end
end
